suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(patchwork)
  library(BPCells)
})

setwd("/projects/b1169/boles/img_scfrp")

results_dir <- "results/07_clustering/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

data_out_dir <- "data/07_clustering/"
dir.create(data_out_dir,
           showWarnings = F,
           recursive = T)

# Build object from 06 ----------------------------------------------------

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


# Find sNN and clusters ---------------------------------------------------

obj <- obj %>% 
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T)

for (i in c(2, 0.4)){
  
  colname <- paste0("res", str_replace_all(i, "[.]", "-"))
  
  obj <- obj %>% 
    FindClusters(method = "igraph",
                 algorithm = 4,
                 resolution = i,
                 cluster.name = colname,
                 graph.name = "RNA_snn")
  
  markers <- FindAllMarkers(obj)
  
  write.csv(markers,
            file = paste0(results_dir, "cluster_", colname, "_markers.csv"))
  
  p1 <- DimPlot_scCustom(obj,
                         reduction = "harmony_umap",
                         label = F,
                         group.by = colname)
  ggsave(p1,
         filename = paste0(results_dir, "/cluster_", colname, "_umap.png"),
         units = "in", dpi = 600,
         height = 6, width = 8)
  
  p2 <- dittoBarPlot(obj,
                     var = "genotype",
                     group.by = colname)
  ggsave(p2,
         filename = paste0(results_dir, "/clusters_", colname, "_by_genotype.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
  
  p3 <- dittoBarPlot(obj,
                     var = "batch",
                     group.by = colname)
  ggsave(p3,
         filename = paste0(results_dir, "/clusters_", colname, "_by_batch.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
  
  p4 <- dittoBarPlot(obj,
                     var = "treatment",
                     group.by = colname)
  ggsave(p4,
         filename = paste0(results_dir, "/clusters_", colname, "_by_treatment.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
  
  p5 <- dittoBarPlot(obj,
                     var = "condition",
                     group.by = colname)
  ggsave(p5,
         filename = paste0(results_dir, "/clusters_", colname, "_by_condition.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
}

saveRDS(obj@meta.data,
        file = paste0(data_dir, "clustered_metadata.rds")) 
