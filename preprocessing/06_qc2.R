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

data_out_dir <- "data/06_qc2/"
dir.create(data_out_dir,
           showWarnings = F,
           recursive = T)

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

# Normalize, scale, run PCA --------------------------------------------------

message2("Normalizing data")
obj <- NormalizeData(obj)

message2("Finding variable features")
obj <- FindVariableFeatures(obj)

message2("Scaling data")
obj <- ScaleData(obj)

message2("Running PCA")
obj <- RunPCA(obj, npcs = 50)

# Integrate with Harmony -------------------------------------------------

message2("Splitting layers by sample for integration")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

message2("Integrating samples using Harmony")

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony",
                       dims = 1:10)

# IntegrateLayers() only adds the new reduction -- it doesn't touch/join the
# assay's existing layers, so counts/data are still split by sample here.
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
                         group.by = "genotype",
                         colors_use = genotype_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p3 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "treatment",
                         colors_use = treatment_pal,
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
         filename = paste0(plots_dir, type, "_umaps.png"),
         units = "in", dpi = 600,
         height = 8,
         width = 10)
}

# Save ----------------------------------------------------------------------

message2("Saving counts matrix as BPCells on-disk matrix")

# Removed first if present, so a rerun doesn't fail on "path already exists"
# against a stale/incomplete directory from a previous attempt.
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