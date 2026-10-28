# Rebuilds the full coculture object from 02_qc.R's filtered BPCells
# counts + metadata, attaches DoubletFinder's per-capture classifications
# (DF.unadj/DF.adj) from 03_doubletfinder.R's metadata_persample/ output,
# then normalizes, scales, and runs PCA. Mirrors
# monocultures/r_scripts/04_norm_pca.R's pattern and parameters exactly
# (NormalizeData/FindVariableFeatures/ScaleData defaults, RunPCA
# npcs = 50) -- these are the same parameters recovered from the lost
# pipeline's own obj@commands (see CLAUDE.md's "Coculture arm, next
# stage" section), so no changes needed there. Doublet removal itself
# isn't done here -- DF.unadj/DF.adj are attached as metadata for a later
# stage (06_qc2.R-equivalent) to act on, same as monoculture.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(BPCells)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp/cocultures")

options(future.globals.maxSize = 64000 * (1024^2))

results_dir <- "results/04_norm_pca/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_out_dir <- "data/04_norm_pca/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

# Read in 02's filtered, whole-cohort object --------------------------------

message2("Reading in filtered object from 02_qc.R")

counts <- open_matrix_dir("data/02_qc/bpcells")
meta <- readRDS("data/02_qc/metadata.rds")
obj <- CreateSeuratObject(counts = counts, meta.data = meta)

# Attach DoubletFinder classifications from 03 -------------------------------

message2("Loading per-capture DoubletFinder metadata")

# Same 30-capture manifest 03_doubletfinder.R's SLURM array uses -- checked
# against what's actually on disk rather than trusting it, so a still-
# running or failed array task is caught here instead of silently
# producing an object with NA DF calls for some cells.
manifest <- read.csv("cellbender_scripts/cell_quantities.csv") %>%
  arrange(pool_dir, sample)

df_dir <- "data/03_doubletfinder/metadata_persample/"
df_files <- list.files(df_dir, pattern = "\\.rds$")
captures_found <- str_remove(df_files, "\\.rds$")

missing_captures <- setdiff(manifest$sample, captures_found)
if (length(missing_captures) > 0) {
  stop(paste0(length(missing_captures), " capture(s) missing from ", df_dir,
              " -- 03_doubletfinder.R may not have finished (or failed) for: ",
              paste(missing_captures, collapse = ", ")))
}

df_list <- list()
for (i in seq_along(captures_found)) {
  message(paste0("Loading ", captures_found[i]))

  df_meta <- readRDS(paste0(df_dir, captures_found[i], ".rds"))

  # Only DF.unadj/DF.adj are needed downstream -- pANN_* is dropped
  # because its column name embeds each capture's own optimal_pK (found
  # independently per capture in 03), so keeping it would explode the
  # combined data frame with mostly-NA per-capture columns on rbind
  # below.
  df_list[[i]] <- df_meta %>%
    dplyr::select(DF.unadj, DF.adj) %>%
    rownames_to_column(var = "cell")
}

df_meta_all <- list_rbind(df_list) %>%
  column_to_rownames(var = "cell")

# Cell identity is the rowname (barcode), already prefixed with the
# distinct sample id by both 01_obj_creation.R (at the whole-cohort
# merge) and 03_doubletfinder.R (re-applied to its own per-capture matrix
# before subsetting) -- same convention on both sides, so a straight
# rowname match is reliable.
missing_cells <- setdiff(colnames(obj), rownames(df_meta_all))
extra_cells <- setdiff(rownames(df_meta_all), colnames(obj))
if (length(missing_cells) > 0 | length(extra_cells) > 0) {
  stop(paste0(length(missing_cells), " cell(s) in the 02_qc.R object have no ",
              "DoubletFinder classification, and ", length(extra_cells),
              " DoubletFinder-classified cell(s) aren't in the 02_qc.R ",
              "object -- 02 and 03 have diverged (rerun whichever is stale ",
              "before continuing)."))
}

obj <- AddMetaData(obj, df_meta_all)

# Normalize, scale, run PCA --------------------------------------------------

message2("Normalizing data")
obj <- NormalizeData(obj)

message2("Finding variable features")
obj <- FindVariableFeatures(obj)

message2("Scaling data")
obj <- ScaleData(obj)

message2("Running PCA")
obj <- RunPCA(obj, npcs = 50)

# Assess the PCA --------------------------------------------------------------

message2("Making elbow plot")

p <- ElbowPlot(obj, ndims = 50)
ggsave(p, filename = paste0(results_dir, "elbow.png"),
       units = "in", dpi = 600, height = 4, width = 5, bg = "white")

message2("Finding PCA elbow")

# Geometric elbow finder, same as monocultures/r_scripts/04_norm_pca.R
# (itself ported from als_cns_scrnaseq/r_scripts/08_pca_evaluation.R):
# draws a line from the first to the last PC's stdev, then picks the PC
# whose point sits farthest (perpendicular distance) from that line -- a
# quantitative stand-in for eyeballing the "corner of the L" off the plot
# above.
find_elbow <- function(stdev){
  n <- length(stdev)

  line_vec <- c(n - 1, stdev[n] - stdev[1])
  line_vec <- line_vec / sqrt(sum(line_vec^2))

  vecs <- cbind((1:n) - 1, stdev - stdev[1])
  proj_len <- vecs %*% line_vec
  proj <- proj_len %*% line_vec
  perp <- vecs - proj
  dist <- sqrt(rowSums(perp^2))

  which.max(dist)
}

elbow_pc <- find_elbow(Stdev(obj[["pca"]]))
message(paste0("Suggested # PCs (max-distance elbow): ", elbow_pc))

message2("Making PC loading plots")

Iterate_PC_Loading_Plots(obj, file_path = results_dir, file_name = "pca_loadings")

# This arm's palettes: condition (10 real genotype-pairs), neuron/
# microglia genotype (5 real genotypes each, microglia additionally has
# "CN" for the neuron-alone wells), batch (3 levels -- independent of
# the monoculture arm's own 4-level batch, see CLAUDE.md).
genotype_levels <- c("KOLF", "GALC", "GBA", "GRN", "LRRK2")
genotype_pal <- c("red", "yellow", "green", "blue", "purple")
microglia_pal <- c("gray50", genotype_pal)
batch_pal <- Dark2_Pal()[1:3]
condition_pal <- DiscretePalette_scCustomize(num_colors = nlevels(obj$condition),
                                             palette = "varibow")

pca_plot <- function(dims){

  # Grouped by condition (genotype-pair, 10 levels), not the 30-level
  # `sample` -- same "condition is the color-by stand-in for real sample
  # groups" convention as monoculture's own PCA plots.
  p1 <- DimPlot_scCustom(obj,
                         reduction = "pca",
                         label = F,
                         dims = dims,
                         group.by = "condition",
                         colors_use = condition_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p2 <- DimPlot_scCustom(obj,
                         reduction = "pca",
                         label = F,
                         dims = dims,
                         group.by = "neuron_genotype",
                         colors_use = genotype_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p3 <- DimPlot_scCustom(obj,
                   reduction = "pca",
                   label = F,
                   dims = dims,
                   group.by = "microglia_genotype",
                   colors_use = microglia_pal,
                   raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p4 <- DimPlot_scCustom(obj,
                   reduction = "pca",
                   label = F,
                   dims = dims,
                   group.by = "batch",
                   colors_use = batch_pal,
                   raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p <- p1 + p2 + p3 + p4 +
    plot_layout(ncol = 2)

  ggsave(p,
         filename = paste0(results_dir, "pca_dims", dims[1], "-", dims[2], ".png"),
         units = "in", dpi = 600,
         height = 8,
         width = 10)

}

pca_plot(c(1,2))
pca_plot(c(1,3))
pca_plot(c(1,4))
pca_plot(c(1,5))
pca_plot(c(2,3))
pca_plot(c(2,4))
pca_plot(c(2,5))
pca_plot(c(3,4))
pca_plot(c(3,5))
pca_plot(c(4,5))

# Save --------------------------------------------------------------------

message2("Saving normalized data as BPCells on-disk matrix")

# Removed first if present, so a rerun doesn't fail on "path already
# exists" against a stale/incomplete directory from a previous attempt.
bpcells_data_dir <- paste0(data_out_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)) {
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

message2("Saving metadata as RDS")

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

message2("Saving PCA reduction as RDS")

saveRDS(obj[["pca"]],
        file = paste0(data_out_dir, "pca.rds"))

message2("Saving variable features as RDS")

saveRDS(VariableFeatures(obj),
        file = paste0(data_out_dir, "variable_features.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   data_mat <- open_matrix_dir(paste0(data_out_dir, "bpcells_data"))
#   meta <- readRDS(paste0(data_out_dir, "metadata.rds"))
#   pca <- readRDS(paste0(data_out_dir, "pca.rds"))
#   obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
#   obj[["RNA"]]$data <- data_mat
#   obj[["pca"]] <- pca
