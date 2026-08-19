library(scCustomize)
library(Seurat)
library(tidyverse)
library(patchwork)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

data_out_dir <- paste0(proj_dir, "data/05_integration/")
dir.create(data_out_dir, F, T)

plots_dir <- paste0(proj_dir, "plots/05_integration/")
dir.create(plots_dir, F, T)

obj <- readRDS(paste0(proj_dir, "data/04_sct_pca/obj.rds"))

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
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "raw_umap",
          return.model = T)

make_umaps <- function(reduction, type){
  
  p1 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "sample",
                         colors_use = DiscretePalette_scCustomize(num_colors = 15,
                                                                  palette = "varibow"),
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p2 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "genotype",
                         colors_use = c("red", "yellow", "green", "blue", "purple"),
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p3 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "treatment",
                         colors_use = JCO_Four()[1:3],
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p4 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "batch",
                         colors_use = Dark2_Pal()[1:4],
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

make_umaps("raw_umap", "raw")

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony_pca",
                       dims = 1:10)

obj <- obj %>%
  FindNeighbors(reduction = "harmony_pca",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony_pca",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

make_umaps("harmony_umap", "harmony")

saveRDS(obj,
        paste0(data_out_dir, "harmony_obj.rds"))