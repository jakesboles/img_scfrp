# Clusters the whole-cohort object at two resolutions and finds markers for
# each. Loads 06_qc2.R's BPCells/RDS pieces (not a monolithic object).
#
# Performance fix for FindAllMarkers(): the original draft handed
# obj[["RNA"]]$data straight from open_matrix_dir() (a lazy, on-disk BPCells
# IterableMatrix, stored column-major/cell-major since that's how
# 06_qc2.R wrote it from a standard Seurat `data` layer) directly to
# FindAllMarkers(). Seurat has a documented, known-slow (in some versions,
# warned-about -- satijalab/seurat#10365) code path for exactly this case:
# Wilcoxon marker testing needs per-GENE access across all cells, which is
# the *wrong* direction for a column-major/cell-major matrix (each query has
# to touch every column's compressed data) -- and it runs TWICE here, once
# per resolution.
#
# Considered (and rejected) the "transpose before testing" fix that BPCells'
# own docs recommend for very large datasets that can't fit in memory
# (BPCells implements an efficient disk-backed transpose specifically for
# this): the actual bottleneck for a lazy IterableMatrix is running the test
# *without* ever materializing it, so each per-gene query re-touches disk.
# Once the matrix is fully materialized into an in-memory dgCMatrix -- fine
# at this scale (~40-60k cells), and already this project's established
# pattern for handing BPCells data to non-BPCells-aware tools (see
# wgcna.R/wgcna_stats.R/milo.R) -- a single one-pass read already happens
# regardless of orientation, and presto (the fast backend Seurat's
# FindAllMarkers uses automatically when installed) expects genes-as-rows/
# cells-as-columns anyway, i.e. no transpose relative to Seurat's native
# layout. So: materialize `data` once, up front, outside the resolution
# loop, and skip the transpose. `counts` is left as lazy BPCells --
# clustering/marker-finding never touch it (FindAllMarkers defaults to the
# `data` layer), so materializing it would just cost memory for nothing.
#
# See the chat reply for two further options NOT applied here without
# confirmation, since they change what gets reported/computed:
# restricting FindAllMarkers to VariableFeatures(obj) (~10x fewer genes
# tested, but genuinely excludes non-variable markers), and splitting the
# two resolutions into a 2-task SLURM array for wall-clock parallelism.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(patchwork)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp")

results_dir <- "results/07_clustering/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

data_out_dir <- "data/07_clustering/"
dir.create(data_out_dir,
           showWarnings = F,
           recursive = T)

# FindAllMarkers() uses the much faster presto backend automatically when
# it's installed, and the much slower base-R per-gene/per-cluster loop
# otherwise -- checked explicitly since this is the single biggest lever
# on runtime here, bigger than anything else in this script.
if (!requireNamespace("presto", quietly = TRUE)) {
  warning("Package 'presto' is not installed -- FindAllMarkers() will fall ",
          "back to a much slower implementation. Install with: ",
          "remotes::install_github('immunogenomics/presto')")
}

# Build object from 06 --------------------------------------------------

message2("Reading in integrated object from 06_qc2.R")

counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
meta <- readRDS("data/06_qc2/metadata.rds")
harmony <- readRDS("data/06_qc2/harmony.rds")
harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")
var_features <- readRDS("data/06_qc2/variable_features.rds")

# Materialized once, up front -- see header comment. `counts` is left as
# lazy BPCells since nothing below touches it.
message2("Materializing normalized data matrix")
data_mat <- as(data_mat, "dgCMatrix")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap
VariableFeatures(obj) <- var_features

# Find sNN and clusters ---------------------------------------------------

message2("Building neighbor graph on the harmony embedding")

obj <- obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T)

for (i in c(2, 0.4)){

  colname <- paste0("res", str_replace_all(i, "[.]", "-"))

  message2(paste0("Clustering at resolution = ", i))

  obj <- obj %>%
    FindClusters(method = "igraph",
                 algorithm = 4,
                 resolution = i,
                 cluster.name = colname,
                 graph.name = "RNA_snn")

  message2(paste0("Finding markers for resolution = ", i))

  markers <- FindAllMarkers(obj)

  write.csv(markers,
            file = paste0(results_dir, "cluster_", colname, "_markers.csv"))

  p1 <- DimPlot_scCustom(obj,
                         reduction = "harmony_umap",
                         label = F,
                         group.by = colname)
  ggsave(p1,
         filename = paste0(results_dir, "cluster_", colname, "_umap.png"),
         units = "in", dpi = 600,
         height = 6, width = 8)

  p2 <- dittoBarPlot(obj,
                     var = "genotype",
                     group.by = colname)
  ggsave(p2,
         filename = paste0(results_dir, "clusters_", colname, "_by_genotype.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)

  p3 <- dittoBarPlot(obj,
                     var = "batch",
                     group.by = colname)
  ggsave(p3,
         filename = paste0(results_dir, "clusters_", colname, "_by_batch.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)

  p4 <- dittoBarPlot(obj,
                     var = "treatment",
                     group.by = colname)
  ggsave(p4,
         filename = paste0(results_dir, "clusters_", colname, "_by_treatment.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)

  p5 <- dittoBarPlot(obj,
                     var = "condition",
                     group.by = colname)
  ggsave(p5,
         filename = paste0(results_dir, "clusters_", colname, "_by_condition.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
}

message2("Saving clustered metadata as RDS")

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "clustered_metadata.rds"))

# Downstream scripts should reconstruct the object from 06_qc2.R's pieces
# plus this stage's own metadata (adds res2/res0-4 cluster columns):
#   counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
#   data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
#   meta <- readRDS("data/07_clustering/clustered_metadata.rds")
#   harmony <- readRDS("data/06_qc2/harmony.rds")
#   harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")
#   obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
#   obj[["RNA"]]$data <- data_mat
#   obj[["harmony"]] <- harmony
#   obj[["harmony_umap"]] <- harmony_umap
