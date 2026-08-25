# Runs as a SLURM job array (see run_03_doubletfinder.sh), one task per
# (pool_dir, sample) capture -- 75 total: 5 GEM wells x 15 multiplexed
# samples each. Each task independently loads its own capture's raw
# per-sample BPCells matrix from 01_obj_creation.R and the shared QC
# metadata from 02_qc.R, so no task ever materializes the whole-cohort
# object. Adapted from als_cns_scrnaseq/r_scripts/04_doubletfinder.R's
# per-sample pattern.
#
# NOTE on doublet rate: because DoubletFinder runs within one demultiplexed
# sample's own cells rather than across the whole GEM well, the doublet
# rate below is computed from that one sample's own recovered cell count --
# matching the ALS script exactly, as asked. But the real physical
# doublet-formation probability for FRP/probe-multiplexed data depends on
# how many total cells were loaded into the WHOLE GEM well (up to 15
# samples' worth), not the much smaller slice recovered for one probe
# barcode, so this likely underestimates the true rate here. Flagged
# rather than silently deviated from the requested pattern -- worth
# revisiting if it matters for your analysis.

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

# Figure out which capture this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == "") {
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see run_03_doubletfinder.sh), one task per capture, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

message2("Loading manifest")

# cell_quantities.csv has one row per (pool_dir, sample) capture -- the same
# 75 distinct-capture identities 01_obj_creation.R wrote to
# bpcells_persample/. Sorting guarantees the same capture <-> array index
# mapping every run.
manifest <- read.csv("cellbender_scripts/cell_quantities.csv") %>%
  arrange(pool_dir, sample)

if (task_id < 1 | task_id > nrow(manifest)) {
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(manifest), " captures -- check the --array range in ",
              "run_03_doubletfinder.sh against nrow(cell_quantities.csv)."))
}

capture <- manifest[task_id, ]
sample_id <- capture$sample   # distinct per-capture id, e.g. GALC_asyn_1_redo
pool_dir <- capture$pool_dir
# 01_obj_creation.R strips _redo before storing this as the object's
# `sample` metadata (both captures of the same biological sample share it);
# pool_dir is what actually distinguishes them, so both are needed to find
# the right cells in 02_qc.R's metadata below.
sample_unified <- str_remove(sample_id, "_redo$")

message2(paste0("Processing ", sample_id, " (", pool_dir, ", task ", task_id,
                "/", nrow(manifest), ")"))

# Load and filter this capture's matrix ---------------------------------------

message2("Loading raw per-capture matrix")

# Cell barcodes in the per-capture matrices from 01_obj_creation.R are NOT
# sample-prefixed (that prefixing happens later, at the whole-cohort merge
# step), so it's reapplied here to match 02_qc.R's metadata, whose cell
# names are "<sample_id>_<barcode>".
mat <- open_matrix_dir(paste0("data/01_obj_assembly/bpcells_persample/", sample_id))
colnames(mat) <- paste0(sample_id, "_", colnames(mat))

message2("Loading QC metadata for this capture")

meta_all <- readRDS("data/02_qc/metadata.rds")

# Matched on pool_dir + sample_unified (both exact categorical matches)
# rather than string-prefix matching sample_id against rownames -- e.g.
# "GALC_asyn_1" is a literal string prefix of "GALC_asyn_1_redo", so
# prefix matching would silently pull in the wrong capture's cells.
meta_sample <- meta_all[meta_all$pool_dir == pool_dir & meta_all$sample == sample_unified, ]

# On purpose: this is the cell count for this demultiplexed sample BEFORE
# 02_qc.R's discard filter (see the doublet-rate note at the top of this
# file), matching als_cns_scrnaseq/r_scripts/04_doubletfinder.R's
# n_cells_preqc <- nrow(meta_sample) exactly.
n_cells_preqc <- nrow(meta_sample)

meta_sample <- meta_sample[meta_sample$discard == F, ]

# DoubletFinder's synthetic-doublet generation samples cell pairs WITH
# replacement, which BPCells' `[` operator rejects (its lazy/streaming
# model isn't built to replay the same column twice in one selection).
# Materializing to a normal in-memory matrix avoids that -- fine memory-wise
# for one sample's few thousand cells, and
# NormalizeData()/ScaleData()/RunPCA() below need it materialized anyway.
mat <- as(mat[, rownames(meta_sample)], "dgCMatrix")

doublet_rate <- (n_cells_preqc / 10000) * 0.08

obj <- CreateSeuratObject(counts = mat, meta.data = meta_sample)

# Need to remove non-probe genes again -- 01_obj_creation.R's
# bpcells_persample matrices were written before the probe-set gene filter
# was applied to the merged object, so they still carry the full gene set.
probes <- read.csv("/projects/p31535/boles/cellranger_references/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                   skip = 5)
genes <- probes$probe_id %>%
  str_split_i(pattern = "[|]", i = 2) %>%
  unique()

idx <- rownames(obj) %in% genes
obj <- obj[idx, ]

# Define function to run DoubletFinder ---------------------------------------

run_doubletfinder <- function(s, doublet_rate, sample_id){

  # Standard normalization and scaling
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Default Seurat clustering and UMAP
  s <- s %>% RunPCA() %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  s <- RunUMAP(s, dims = 1:10)

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

  # Plot unadjusted vs. adjusted doublets in UMAP coordinates. Filenames use
  # the sample_id passed in (the distinct per-capture id), not
  # s$orig.ident/s$sample -- both are the *unified* name (01_obj_creation.R
  # strips _redo), so iMG1's and iMG1_redo's runs of the same biological
  # sample would otherwise overwrite each other's plots.
  p1 <- DimPlot_scCustom(s, reduction = "umap", group.by = "DF.unadj",
                        pt.size = 1, shuffle = TRUE, alpha = 0.6) +
    ggtitle("Unadjusted")
  p2 <- DimPlot_scCustom(s, reduction = "umap", group.by = "DF.adj",
                        pt.size = 1, shuffle = TRUE, alpha = 0.6) +
    ggtitle("Adjusted for Homotypic Proportion")
  p <- p1 + p2 + plot_layout(ncol = 2, nrow = 1, guides = "collect")

  ggsave(p, filename = paste0(plots_dir, sample_id, "_doublet_umap.png"),
         units = "in", dpi = 600, height = 5, width = 10, bg = "white")

  return(s)
}

message2("Running DoubletFinder")

obj <- run_doubletfinder(obj, doublet_rate, sample_id)

# Save this capture's filtered, DoubletFinder-processed matrix + metadata ---
message2("Saving filtered, DoubletFinder-processed capture")

counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")
write_matrix_dir(mat = counts_out,
                 dir = paste0(data_out_dir, "bpcells_persample/", sample_id))

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata_persample/", sample_id, ".rds"))
