# Integrates the whole-cohort object across samples using Harmony, on top
# of the PCA computed in 04_norm_pca.R. Mirrors
# monocultures/r_scripts/05_integration_harmony.R's structure and
# parameters (dims = 1:10 throughout, same FindNeighbors/RunUMAP
# settings) -- all of this stage's real design choices were confirmed
# directly with the user rather than guessed, since the lost pipeline's
# own version of this step did things differently (see CLAUDE.md's
# "Coculture arm, next stage" section for the full history):
#   - QC method: monoculture's MAD-threshold approach (02_qc.R), not the
#     rank-based scheme found in the lost pipeline's surviving object.
#   - PC dims: 1:10 throughout, not the lost pipeline's inconsistent
#     10-vs-30 usage.
#   - Harmony integration method: Seurat's IntegrateLayers(method =
#     "HarmonyIntegration"), matching monoculture (the lost pipeline's
#     own invocation wasn't logged either way, so this was a free
#     choice).
#   - Split-by variable: the RNA assay is split by the FULL 30-level
#     genotype-pair x batch-replicate `sample` identity (e.g.
#     "KOLF_KOLF_2") before integrating -- NOT by pool (4 levels, what
#     the lost pipeline actually did) and NOT by the 10-level unified
#     `condition` (monoculture's own analogous choice is to split by its
#     60-level `sample`, not its 5-level `pool_dir`; this arm's 30-level
#     `sample` is the direct equivalent). This additionally separates
#     JSB134-3's two mixed batches into their own layers, something a
#     pool-level split structurally couldn't do here.

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

results_dir <- "results/05_integration/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_out_dir <- "data/05_integration/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

# Read in 04's normalized object + 02's matching raw counts ------------------

message2("Reading in normalized object from 04_norm_pca.R")

data_mat <- open_matrix_dir("data/04_norm_pca/bpcells_data")
meta <- readRDS("data/04_norm_pca/metadata.rds")
pca <- readRDS("data/04_norm_pca/pca.rds")
var_features <- readRDS("data/04_norm_pca/variable_features.rds")

message2("Reading in matching raw counts from 02_qc.R")

# 04_norm_pca.R doesn't re-save raw counts (it only adds metadata
# columns/reductions on top of 02's object, never drops or reorders
# cells), so 02's bpcells directory is still the right source for a real
# counts layer here. Checked and explicitly reordered to 04's cell order
# rather than assumed -- same caution 04 itself uses for its own 02/03
# cell reconciliation.
counts_mat <- open_matrix_dir("data/02_qc/bpcells")

missing_cells <- setdiff(rownames(meta), colnames(counts_mat))
extra_cells <- setdiff(colnames(counts_mat), rownames(meta))
if (length(missing_cells) > 0 | length(extra_cells) > 0) {
  stop(paste0(length(missing_cells), " cell(s) in 04_norm_pca.R's metadata ",
              "have no matching raw counts in data/02_qc/bpcells, and ",
              length(extra_cells), " cell(s) in data/02_qc/bpcells aren't ",
              "in 04's metadata -- 02 and 04 have diverged."))
}
counts_mat <- counts_mat[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["pca"]] <- pca
VariableFeatures(obj) <- var_features

# Plot helper, reused for the pre- and post-integration UMAPs ----------------

genotype_levels <- c("KOLF", "GALC", "GBA", "GRN", "LRRK2")
genotype_pal <- c("red", "yellow", "green", "blue", "purple")
microglia_pal <- c("gray50", genotype_pal)
batch_pal <- Dark2_Pal()[1:3]
condition_pal <- DiscretePalette_scCustomize(num_colors = nlevels(obj$condition),
                                             palette = "varibow")

make_umaps <- function(reduction, type){

  p1 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "condition",
                         colors_use = condition_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p2 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "neuron_genotype",
                         colors_use = genotype_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p3 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "microglia_genotype",
                         colors_use = microglia_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p4 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "batch",
                         colors_use = batch_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p <- p1 + p2 + p3 + p4 +
    plot_layout(ncol = 2)

  ggsave(p,
         filename = paste0(results_dir, type, "_umaps.png"),
         units = "in", dpi = 600,
         height = 8,
         width = 10)
}

# Pre-integration UMAP, for before/after comparison --------------------------

message2("Computing pre-integration UMAP")

obj <- obj %>%
  FindNeighbors(reduction = "pca",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "raw_umap",
          return.model = T)

make_umaps("raw_umap", "raw")

# Integrate with Harmony -------------------------------------------------

message2("Splitting layers by full sample identity for integration")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

message2("Integrating samples using Harmony")

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony",
                       dims = 1:10)

# IntegrateLayers() only adds the new reduction -- it doesn't touch/join
# the assay's existing layers, so counts/data are still split by sample
# here.
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

message2("Computing post-integration UMAP")

obj <- obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

make_umaps("harmony_umap", "harmony")

# Save ----------------------------------------------------------------------

message2("Saving counts matrix as BPCells on-disk matrix")

counts_dir <- paste0(data_out_dir, "bpcells_counts")
if (dir.exists(counts_dir)) {
  unlink(counts_dir, recursive = T)
}

counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")
write_matrix_dir(mat = counts_out, dir = counts_dir)

message2("Saving normalized data matrix as BPCells on-disk matrix")

data_dir <- paste0(data_out_dir, "bpcells_data")
if (dir.exists(data_dir)) {
  unlink(data_dir, recursive = T)
}

write_matrix_dir(mat = obj[["RNA"]]$data, dir = data_dir)

message2("Saving metadata as RDS")

saveRDS(obj@meta.data, file = paste0(data_out_dir, "metadata.rds"))

message2("Saving Harmony reduction as RDS")

saveRDS(obj[["harmony"]], file = paste0(data_out_dir, "harmony.rds"))

message2("Saving Harmony UMAP as RDS")

saveRDS(obj[["harmony_umap"]], file = paste0(data_out_dir, "harmony_umap.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts_mat <- open_matrix_dir(paste0(data_out_dir, "bpcells_counts"))
#   data_mat <- open_matrix_dir(paste0(data_out_dir, "bpcells_data"))
#   meta <- readRDS(paste0(data_out_dir, "metadata.rds"))
#   harmony <- readRDS(paste0(data_out_dir, "harmony.rds"))
#   harmony_umap <- readRDS(paste0(data_out_dir, "harmony_umap.rds"))
#   obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
#   obj[["RNA"]]$data <- data_mat
#   obj[["harmony"]] <- harmony
#   obj[["harmony_umap"]] <- harmony_umap
