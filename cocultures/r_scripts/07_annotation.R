# Broad microglia-vs-neuron cell-type annotation -- the one pipeline step
# this arm needs that monoculture has no counterpart for (monoculture is
# microglia-only; this arm's wells contain both cell types mixed
# together). Clusters the whole-cohort object fresh on 06_qc2.R's harmony
# reduction, plots canonical microglia and neuron marker expression per
# cluster, and labels clusters by eye -- the same by-eye pattern
# 06_qc2.R uses for monoculture's doublet-cluster/non-microglia-cluster
# calls, just applied here to assign a label rather than remove cells.
#
# Confirmed with the user: the lost pipeline did the same cluster +
# manual-marker-labeling, just never persisted the resulting label
# column anywhere the surviving reference object retained (see
# CLAUDE.md's "Coculture arm, next stage" section).
#
# Marker choice: microglia markers (AIF1, CSF1R) match monoculture's own
# 06_qc2.R exactly. Neuron markers (RBFOX3/NeuN, MAP2, TUBB3/Tuj1, SYP,
# STMN2) are canonical cortical-neuron markers chosen as a reasonable
# starting default -- NOT confirmed against what the lost pipeline
# actually used (unrecoverable; a manual labeling step doesn't log to
# obj@commands, see CLAUDE.md). Worth a second look against the actual
# per-cluster expression in results_dir before finalizing the cell_type
# assignment below, and worth telling Claude/noting here if different
# markers turn out to work better.
#
# Meant to be run interactively, same as 02_qc.R/06_qc2.R: inspect the
# marker plots, fill in `cluster_labels` below by eye, then rerun from
# that point down. The script hard-stops if any cluster is left
# unlabeled, rather than silently leaving some cells as NA cell_type.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(patchwork)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp/cocultures")

results_dir <- "results/07_annotation/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_out_dir <- "data/07_annotation/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

# Read in 06's re-integrated object ------------------------------------------

message2("Reading in object from 06_qc2.R")

counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
meta <- readRDS("data/06_qc2/metadata.rds")
harmony <- readRDS("data/06_qc2/harmony.rds")
harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")
var_features <- readRDS("data/06_qc2/variable_features.rds")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap
VariableFeatures(obj) <- var_features

# Cluster fresh on 06's harmony reduction ------------------------------------

message2("Clustering on 06's harmony reduction")

obj <- obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  FindClusters(method = "igraph",
               algorithm = 4,
               resolution = 1,
               cluster.name = "annotation_cluster",
               graph.name = "RNA_snn")

p <- DimPlot_scCustom(obj,
                      reduction = "harmony_umap",
                      group.by = "annotation_cluster",
                      label = T)
ggsave(p, filename = paste0(results_dir, "clusters.png"),
       units = "in", dpi = 600, height = 6, width = 7)

# Plot marker expression per cluster, by eye ---------------------------------

message2("Making marker plots")

microglia_markers <- c("AIF1", "CSF1R")
neuron_markers <- c("RBFOX3", "MAP2", "TUBB3", "SYP", "STMN2")

for (i in c(microglia_markers, neuron_markers, "nCount_RNA", "nFeature_RNA")) {

  p <- VlnPlot_scCustom(obj,
                        group.by = "annotation_cluster",
                        features = i,
                        pt.size = 0) +
    NoLegend()

  if (i %in% c("nCount_RNA", "nFeature_RNA")) { p <- p + scale_y_log10() }

  ggsave(p,
         filename = paste0(results_dir, i, "_vln.png"),
         units = "in", dpi = 600,
         height = 4, width = 10)

  p <- FeaturePlot_scCustom(obj,
                            reduction = "harmony_umap",
                            features = i)
  ggsave(p,
         filename = paste0(results_dir, i, "_umap.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)
}

# Label clusters by eye -------------------------------------------------------

# Fill in after inspecting the marker plots above -- every cluster id in
# `annotation_cluster` must appear as a name here, mapped to "Microglia"
# or "Neuron" (add a third category here, e.g. "Doublet"/"Ambiguous", if
# a cluster doesn't clearly resolve to either -- don't force a bad call
# into one of the two just to satisfy the stop() guard below).
cluster_labels <- c()

cluster_ids <- levels(obj$annotation_cluster)
missing_labels <- setdiff(cluster_ids, names(cluster_labels))
if (length(missing_labels) > 0) {
  stop(paste0("cluster_labels is missing an entry for cluster(s): ",
              paste(missing_labels, collapse = ", "),
              " -- inspect ", results_dir, "<marker>_vln.png/_umap.png ",
              "and fill in cluster_labels before continuing."))
}

obj$cell_type <- factor(cluster_labels[as.character(obj$annotation_cluster)],
                        levels = unique(cluster_labels))

p <- DimPlot_scCustom(obj,
                      reduction = "harmony_umap",
                      group.by = "cell_type")
ggsave(p, filename = paste0(results_dir, "cell_type_umap.png"),
       units = "in", dpi = 600, height = 6, width = 7)

write.csv(obj@meta.data %>% count(annotation_cluster, cell_type),
          file = paste0(results_dir, "cluster_cell_type_mapping.csv"),
          row.names = F)

# Save ------------------------------------------------------------------------

message2("Saving metadata as RDS")

# Only metadata changes at this stage (new annotation_cluster + cell_type
# columns) -- counts/data/harmony/harmony_umap/variable_features are
# unchanged from 06_qc2.R, so downstream scripts reach back there for
# those rather than having them re-saved here.
saveRDS(obj@meta.data, file = paste0(data_out_dir, "metadata.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
#   data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
#   meta <- readRDS("data/07_annotation/metadata.rds")
#   harmony <- readRDS("data/06_qc2/harmony.rds")
#   harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")
#   var_features <- readRDS("data/06_qc2/variable_features.rds")
#   obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
#   obj[["RNA"]]$data <- data_mat
#   obj[["harmony"]] <- harmony
#   obj[["harmony_umap"]] <- harmony_umap
#   VariableFeatures(obj) <- var_features
