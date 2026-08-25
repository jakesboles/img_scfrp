# Runs DoubletFinder independently on each of the 5 GEM wells (iMG1,
# iMG1_redo, iMG2, iMG3, iMG4) -- doublets form during droplet
# partitioning, which happens once per physical 10x reaction, not once per
# multiplexed sample, so this splits by `pool_dir` (the real GEM well),
# not samples.csv's `batch` column (1-4), which iMG1 and iMG1_redo share
# and would incorrectly pool them back together.
#
# Adapted from als_cns_scrnaseq/r_scripts/04_doubletfinder.R, but that
# script processes one individual sample per SLURM array task (that
# project has one GEM well per sample, so "per sample" and "per GEM well"
# are the same thing there -- not true here, where FRP multiplexing puts
# up to 16 samples in one GEM well). Kept as a single script looping over
# the 5 pools instead of a job array, since 5 units doesn't need SLURM
# array parallelism the way 90 (or even this project's 75 samples) would.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(DoubletFinder)
  library(scCustomize)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp")

plots_dir <- "plots/03_doubletfinder/"
dir.create(plots_dir, showWarnings = F, recursive = T)

data_out_dir <- "data/03_doubletfinder/"
dir.create(paste0(data_out_dir, "bpcells_persample/"), showWarnings = F, recursive = T)
dir.create(paste0(data_out_dir, "metadata_persample/"), showWarnings = F, recursive = T)

message2("Reading in QC-filtered object")

counts <- open_matrix_dir("data/02_qc/bpcells")
meta <- readRDS("data/02_qc/metadata.rds")

pools <- sort(unique(meta$pool_dir))
message(paste0("Found ", length(pools), " pools: ", paste(pools, collapse = ", ")))

# Doublet rate uses each pool's PRE-QC cell count on purpose: doublets form
# during the original partitioning step, so the rate should reflect how
# many cells were originally captured in that GEM well, not how many
# survived 02_qc.R's filtering. 01_obj_assembly's metadata has every cell,
# pre-filter, with `pool_dir` already attached.
preqc_meta <- readRDS("data/01_obj_assembly/metadata.rds")
preqc_counts <- preqc_meta %>%
  count(pool_dir, name = "n_cells_preqc")

# Define function to run DoubletFinder ---------------------------------------

run_doubletfinder <- function(s, doublet_rate, pool){

  # Standard normalization and scaling
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Default Seurat clustering and UMAP
  s <- s %>% RunPCA() %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  s <- RunUMAP(s, dims = 1:10)

  p <- DimPlot(s, reduction = "umap", group.by = "orig.ident")
  ggsave(p, filename = paste0(plots_dir, pool, "_umap.png"),
         units = "in", dpi = 600, height = 5, width = 10)

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))

  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate * nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  # Run DoubletFinder without homotypic adjustment
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK,
                     nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.unadj"

  # Run DoubletFinder with homotypic adjustment
  pANN <- colnames(s@meta.data)[grep("^pANN", colnames(s@meta.data))]
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK,
                     nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.adj"

  # Plot unadjusted vs. adjusted doublets in UMAP coordinates
  p1 <- DimPlot_scCustom(s, reduction = "umap", group.by = "DF.unadj",
                        pt.size = 1, shuffle = TRUE, alpha = 0.6) +
    ggtitle("Unadjusted")
  p2 <- DimPlot_scCustom(s, reduction = "umap", group.by = "DF.adj",
                        pt.size = 1, shuffle = TRUE, alpha = 0.6) +
    ggtitle("Adjusted for Homotypic Proportion")
  p <- p1 + p2 + plot_layout(ncol = 2, nrow = 1, guides = "collect")

  ggsave(p, filename = paste0(plots_dir, pool, "_doublet_umap.png"),
         units = "in", dpi = 600, height = 5, width = 10, bg = "white")

  return(s)
}

# Run DoubletFinder on each pool ----------------------------------------------

for (pool in pools) {
  message2(paste0("Running ", pool))

  # metadata.rds has cell barcodes as actual rownames (02_qc.R's last step
  # before saving is column_to_rownames(var = "cell")), not a "cell"
  # column -- base R `[` subsetting is used here rather than dplyr::filter()
  # to avoid any doubt about whether rownames survive the subset, since a
  # dropped/misaligned barcode here would silently mismatch this metadata
  # against the matrix columns selected below.
  meta_pool <- meta[meta$pool_dir == pool, ]

  # DoubletFinder's synthetic-doublet generation samples cell pairs WITH
  # replacement, which BPCells' `[` operator rejects (its lazy/streaming
  # model isn't built to replay the same column twice in one selection).
  # Materializing to a normal in-memory matrix per pool avoids that --
  # NormalizeData()/ScaleData()/RunPCA() below need the data materialized
  # for a single pool's ~5-16k cells anyway, so this costs nothing extra.
  mat_pool <- as(counts[, rownames(meta_pool)], "dgCMatrix")

  n_preqc <- preqc_counts$n_cells_preqc[preqc_counts$pool_dir == pool]
  doublet_rate <- (n_preqc / 10000) * 0.08

  obj <- CreateSeuratObject(counts = mat_pool, meta.data = meta_pool)
  obj <- run_doubletfinder(obj, doublet_rate, pool)

  message2(paste0("Saving ", pool))

  counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")
  write_matrix_dir(mat = counts_out,
                   dir = paste0(data_out_dir, "bpcells_persample/", pool))

  saveRDS(obj@meta.data,
          file = paste0(data_out_dir, "metadata_persample/", pool, ".rds"))
}
