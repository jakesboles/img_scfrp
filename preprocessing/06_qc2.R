suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(BPCells)
})

setwd("/projects/b1169/boles/img_scfrp")

plots_dir <- "plots/06_qc2/"
dir.create(plots_dir,
           showWarnings = F,
           recursive = T)
data_dir <- "data/06_qc2/"

counts_mat <- open_matrix_dir("data/05_integration/bpcells_counts")
data_mat <- open_matrix_dir("data/05_integration/bpcells_data")
meta <- readRDS("data/05_integration/metadata.rds")
harmony <- readRDS("data/05_integration/harmony.rds")
harmony_umap <- readRDS("data/05_integration/harmony_umap.rds")
var_features <- readRDS("data/04_norm_pca/variable_features.rds")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap
VariableFeatures(obj) <- var_features

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
  FindClusters(method = "igraph",
               algorithm = 4,
               resolution = 2,
               cluster.name = "cluster",
               graph.name = "RNA_snn")


p1 <- DimPlot_scCustom(obj,
                 reduction = "harmony_umap")

p2 <- DimPlot_scCustom(obj,
                 reduction = "harmony_umap",
                 group.by = "DF.adj")

p3 <- dittoBarPlot(obj,
             var = "DF.adj",
             group.by = "cluster")

design = "
AABB
CCCC
"

p <- p1 + p2 + p3 + 
  plot_layout(design = design)

ggsave(p,
       filename = paste0(plots_dir, "clusters_by_doublets.png"),
       units = "in", dpi = 600,
       height = 10, width = 12)

obj <- obj %>% 
  filter(!(cluster %in% c(31, 35, 36)))