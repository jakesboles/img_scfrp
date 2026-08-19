library(Seurat)
library(tidyverse)
library(scCustomize)
library(patchwork)

setwd("/projects/b1169/boles/img_scfrp")

data_out_dir <- "data/05_integration_subsets/"
dir.create(data_out_dir, F, T)

plots_dir <- "plots/05_integration_subsets/"
dir.create(plots_dir, F, T)

obj <- readRDS("data/04_sct_pca/obj.rds")

# Pools 2 and 3 -----------------------------------------------------------

sub23 <- subset(obj,
              subset = batch %in% c("pool2", "pool3"))

# Normalize data

sub23 <- sub23 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# Run PCA and examine 

sub23 <- sub23 %>%
  RunPCA(npcs = 50)

p <- ElbowPlot(sub23,
               ndims = 50)
ggsave(p,
       filename = paste0(plots_dir, "pools2-3_elbow.png"),
       units = "in", dpi = 600,
       height = 4, width = 5,
       bg = "white")

Iterate_PC_Loading_Plots(sub23,
                         file_path = plots_dir,
                         file_name = "pools2-3_pca_loadings")

npcs <- 10

# Integrate with Harmony 

sub23 <- IntegrateLayers(sub23,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony_pca",
                       dims = 1:npcs)

sub23 <- sub23 %>%
  FindNeighbors(reduction = "harmony_pca",
                dims = 1:npcs,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony_pca",
                dims = 1:npcs,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T)

sub23 <- sub23 %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

p1 <- DimPlot_scCustom(sub23,
                       reduction = "harmony_umap",
                       group.by = "sample",
                       colors_use = DiscretePalette_scCustomize(num_colors = 15,
                                                                palette = "varibow"),
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p2 <- DimPlot_scCustom(sub23,
                       reduction = "harmony_umap",
                       group.by = "genotype",
                       colors_use = c("red", "yellow", "green", "blue", "purple"),
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p3 <- DimPlot_scCustom(sub23,
                       reduction = "harmony_umap",
                       group.by = "treatment",
                       colors_use = JCO_Four()[1:3],
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p4 <- DimPlot_scCustom(sub23,
                       reduction = "harmony_umap",
                       group.by = "batch",
                       # colors_use = Dark2_Pal()[1:4],
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 2)

ggsave(p,
       filename = paste0(plots_dir, "pools2-3_umaps.png"),
       units = "in", dpi = 600,
       height = 8,
       width = 10)

saveRDS(sub23,
        paste0(data_out_dir, "pools2-3_harmony_obj.rds"))

# FeaturePlot_scCustom(obj,
#                      reduction = "harmony_umap",
#                      features = c("CD163", "P2RY12", "GPNMB", "C1QA"))
# 
# DimPlot_scCustom(obj,
#                  reduction = "harmony_umap",
#                  group.by = "genotype",
#                  split.by = "batch",
#                  colors_use = c("red", "yellow", "green", "blue", "purple"),
#                  raster = F)

# Pools 1 and 4 -----------------------------------------------------------

sub14 <- subset(obj,
                subset = batch %in% c("pool1", "pool4"))

# Normalize data

sub14 <- sub14 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# Run PCA and examine 

sub14 <- sub14 %>%
  RunPCA(npcs = 50)

p <- ElbowPlot(sub14,
               ndims = 50)
ggsave(p,
       filename = paste0(plots_dir, "pools1-4_elbow.png"),
       units = "in", dpi = 600,
       height = 4, width = 5,
       bg = "white")

Iterate_PC_Loading_Plots(sub14,
                         file_path = plots_dir,
                         file_name = "pools1-4_pca_loadings")

npcs <- 10

# Integrate with Harmony

sub14 <- IntegrateLayers(sub14,
                         method = "HarmonyIntegration",
                         orig.reduction = "pca",
                         new.reduction = "harmony_pca",
                         dims = 1:npcs)

sub14 <- sub14 %>%
  FindNeighbors(reduction = "harmony_pca",
                dims = 1:npcs,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony_pca",
                dims = 1:npcs,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T)

sub14 <- sub14 %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

p1 <- DimPlot_scCustom(sub14,
                       reduction = "harmony_umap",
                       group.by = "sample",
                       colors_use = DiscretePalette_scCustomize(num_colors = 15,
                                                                palette = "varibow"),
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p2 <- DimPlot_scCustom(sub14,
                       reduction = "harmony_umap",
                       group.by = "genotype",
                       colors_use = c("red", "yellow", "green", "blue", "purple"),
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p3 <- DimPlot_scCustom(sub14,
                       reduction = "harmony_umap",
                       group.by = "treatment",
                       colors_use = JCO_Four()[1:3],
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p4 <- DimPlot_scCustom(sub14,
                       reduction = "harmony_umap",
                       group.by = "batch",
                       # colors_use = Dark2_Pal()[1:4],
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p <- p1 + p2 + p3 + p4 + 
  plot_layout(ncol = 2)

ggsave(p,
       filename = paste0(plots_dir, "pools1-4_umaps.png"),
       units = "in", dpi = 600,
       height = 8,
       width = 10)

saveRDS(sub14,
        paste0(data_out_dir, "pools1-4_harmony_obj.rds"))

# Individual pools -------------------------------------------------------------

for (i in c("pool1", "pool2", "pool3", "pool4")){
  
  message(paste0("Running ", i))
  
  sub <- subset(obj,
                  subset = batch == i)
  
  # Normalize data
  
  sub <- sub %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData()
  
  # Run PCA and examine 
  
  sub <- sub %>%
    RunPCA(npcs = 50)
  
  p <- ElbowPlot(sub,
                 ndims = 50)
  ggsave(p,
         filename = paste0(plots_dir, i, "_elbow.png"),
         units = "in", dpi = 600,
         height = 4, width = 5,
         bg = "white")
  
  Iterate_PC_Loading_Plots(sub,
                           file_path = plots_dir,
                           file_name = paste0(i, "_pca_loadings"))
  
  npcs <- 10
  
  sub <- sub %>%
    FindNeighbors(reduction = "pca",
                  dims = 1:npcs,
                  k.param = 15,
                  nn.method = "annoy",
                  annoy.metric = "euclidean",
                  return.neighbor = T) %>%
    FindNeighbors(reduction = "pca",
                  dims = 1:npcs,
                  k.param = 15,
                  nn.method = "annoy",
                  annoy.metric = "euclidean",
                  compute.SNN = T)
  
  sub <- sub %>%
    RunUMAP(umap.method = "uwot",
            # reduction = "cca_pca",
            nn.name = "RNA.nn",
            metric = "euclidean",
            min.dist = 0.5,
            n_neighbors = 15L,
            reduction.name = "umap",
            return.model = T)
  
  p1 <- DimPlot_scCustom(sub,
                         reduction = "umap",
                         group.by = "sample",
                         colors_use = DiscretePalette_scCustomize(num_colors = 15,
                                                                  palette = "varibow"),
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p2 <- DimPlot_scCustom(sub,
                         reduction = "umap",
                         group.by = "genotype",
                         colors_use = c("red", "yellow", "green", "blue", "purple"),
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p3 <- DimPlot_scCustom(sub,
                         reduction = "umap",
                         group.by = "treatment",
                         colors_use = JCO_Four()[1:3],
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p4 <- DimPlot_scCustom(sub,
                         reduction = "umap",
                         group.by = "batch",
                         # colors_use = Dark2_Pal()[1:4],
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p <- p1 + p2 + p3 + p4 + 
    plot_layout(ncol = 2)
  
  ggsave(p,
         filename = paste0(plots_dir, i, "_raw_umaps.png"),
         units = "in", dpi = 600,
         height = 8,
         width = 10)
  
  # Integrate with Harmony
  
  sub <- IntegrateLayers(sub,
                           method = "HarmonyIntegration",
                           orig.reduction = "pca",
                           new.reduction = "harmony_pca",
                           dims = 1:npcs)
  
  sub <- sub %>%
    FindNeighbors(reduction = "harmony_pca",
                  dims = 1:npcs,
                  k.param = 15,
                  nn.method = "annoy",
                  annoy.metric = "euclidean",
                  return.neighbor = T) %>%
    FindNeighbors(reduction = "harmony_pca",
                  dims = 1:npcs,
                  k.param = 15,
                  nn.method = "annoy",
                  annoy.metric = "euclidean",
                  compute.SNN = T)
  
  sub <- sub %>%
    RunUMAP(umap.method = "uwot",
            # reduction = "cca_pca",
            nn.name = "RNA.nn",
            metric = "euclidean",
            min.dist = 0.5,
            n_neighbors = 15L,
            reduction.name = "harmony_umap",
            return.model = T)
  
  p1 <- DimPlot_scCustom(sub,
                         reduction = "harmony_umap",
                         group.by = "sample",
                         colors_use = DiscretePalette_scCustomize(num_colors = 15,
                                                                  palette = "varibow"),
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p2 <- DimPlot_scCustom(sub,
                         reduction = "harmony_umap",
                         group.by = "genotype",
                         colors_use = c("red", "yellow", "green", "blue", "purple"),
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p3 <- DimPlot_scCustom(sub,
                         reduction = "harmony_umap",
                         group.by = "treatment",
                         colors_use = JCO_Four()[1:3],
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p4 <- DimPlot_scCustom(sub,
                         reduction = "harmony_umap",
                         group.by = "batch",
                         # colors_use = Dark2_Pal()[1:4],
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p <- p1 + p2 + p3 + p4 + 
    plot_layout(ncol = 2)
  
  ggsave(p,
         filename = paste0(plots_dir, i, "_harmony_umaps.png"),
         units = "in", dpi = 600,
         height = 8,
         width = 10)
  
  saveRDS(sub,
          paste0(data_out_dir, i, "_harmony_obj.rds"))
}
