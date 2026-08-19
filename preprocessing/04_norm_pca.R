library(Seurat)
library(tidyverse)
library(scCustomize)
library(patchwork)

options(future.globals.maxSize = 64000 * (1024^2))

proj_dir <- "/projects/b1169/boles/img_scfrp/"

obj_dir <- paste0(proj_dir, "data/02_qc/")

meta_dir <- paste0(proj_dir, "data/03_doubletfinder/")
meta_files <- list.files(meta_dir)
meta_list <- list()

data_out_dir <- paste0(proj_dir, "data/04_sct_pca/")
dir.create(data_out_dir,
           recursive = T,
           showWarnings = F)

plots_dir <- paste0(proj_dir, "plots/04_sct_pca/")
dir.create(plots_dir,
           showWarnings = F,
           recursive = T)

# Read in object ----------------------------------------------------------

obj <- list.files(obj_dir)
obj <- readRDS(paste0(obj_dir, obj))

# Read in metadata and collapse into one dataframe ------------------------

for (i in seq_along(meta_files)){
  
  message(str_remove_all(meta_files[i], ".rds"))
  
  meta_list[[i]] <- readRDS(paste0(meta_dir, meta_files[i]))
  
  meta_list[[i]] <- meta_list[[i]]@meta.data %>%
    dplyr::select(c(orig.ident:percent_mito, DF.unadj:DF.adj))
}

# meta_list

meta <- list_rbind(meta_list)

setdiff(colnames(obj), rownames(meta))

obj <- AddMetaData(obj, meta)

# Normalize data  ----------------------------------------------------------------

obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

# Run PCA and examine -----------------------------------------------------

obj <- obj %>%
  RunPCA(npcs = 50)

p <- ElbowPlot(obj,
               ndims = 50)

ggsave(p,
       filename = paste0(plots_dir, "elbow.png"),
       units = "in", dpi = 600,
       height = 4, width = 5,
       bg = "white")

Iterate_PC_Loading_Plots(obj,
                         file_path = plots_dir,
                         file_name = "pca_loadings")

pca_plot <- function(dims){
  
  p1 <- DimPlot_scCustom(obj,
                         reduction = "pca",
                         label = F,
                         dims = dims,
                         group.by = "sample",
                         colors_use = DiscretePalette_scCustomize(num_colors = 15,
                                                                  palette = "varibow"),
                         raster = F) +
    # guides(color = guide_legend(position = "inside")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p2 <- DimPlot_scCustom(obj,
                         reduction = "pca",
                         label = F,
                         dims = dims,
                         group.by = "genotype",
                         colors_use = c("red", "yellow", "green", "blue", "purple"),
                         raster = F) +
    # guides(color = guide_legend(position = "inside")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p3 <- DimPlot_scCustom(obj,
                   reduction = "pca",
                   label = F,
                   dims = dims,
                   group.by = "treatment",
                   colors_use = JCO_Four()[1:3],
                   raster = F) +
    # guides(color = guide_legend(position = "inside")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p4 <- DimPlot_scCustom(obj,
                   reduction = "pca",
                   label = F,
                   dims = dims,
                   group.by = "batch",
                   colors_use =Dark2_Pal()[1:4],
                   raster = F) +
    # guides(color = guide_legend(position = "inside")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
  
  p <- p1+ p2 + p3 + p4 + 
    plot_layout(ncol = 2)
  
  ggsave(p,
         filename = paste0(plots_dir, "pca_dims", dims[1], "-", dims[2], ".png"),
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


# pca_plot2 <- function(dims, grouping){
# 
#   
#   # Idents(obj) <- grouping
#   
#   if (length(unique(Idents(obj))) == 3) {
#     pal <- JCO_Four()[1:3]
#   } else {
#     if (length(unique(Idents(obj))) == 5) {
#       pal <- c("red", "yellow", "green", "blue", "purple")
#     } else {
#       if (length(unique(Idents(obj))) == 4) { 
#         pal <- Dark2_Pal()[1:4]
#       } else {
#         if (length(unique(Idents(obj))) == 15) {
#           pal <- DiscretePalette_scCustomize(num_colors = 15,
#                                              palette = "varibow")
#         }
#       }
#     }
#   }
#   
#   DimPlot_scCustom(obj,
#                    reduction = "pca",
#                    label = F,
#                    dims = dims,
#                    colors_use = pal) +
#     # guides(color = guide_legend(position = "inside")) +
#     theme(axis.text = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           legend.justification = "left")
# }
# 
# pca_grid <- function(grouping){
#   
#   # Idents(obj) <- grouping
#   
#   design <- "
#   AA######
#   BBCC####
#   DDEEFF##
#   GGHHIIJJ
#   "
#   
#   pca_plot2(c(1,2)) +
#     pca_plot2(c(1,3)) + pca_plot2(c(2,3)) +
#     pca_plot2(c(1,4)) + pca_plot2(c(2,4)) + pca_plot2(c(3,4)) +
#     pca_plot2(c(1,5)) + pca_plot2(c(2,5)) + pca_plot2(c(3,5)) + pca_plot2(c(4,5)) +
#     plot_layout(design = design,
#                 guides = "collect")
# }

# Idents(obj) <- "batch"
# p <- pca_grid("batch")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_batch.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 12)
# 
# Idents(obj) <- "genotype"
# p <- pca_grid("genotype")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_genotype.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 12)
# 
# Idents(obj) <- "treatment"
# p <- pca_grid("treatment")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_treatment.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 12)
# 
# Idents(obj) <- "sample"
# p <- pca_grid("sample")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_sample.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 12)


# Save object -------------------------------------------------------------

saveRDS(obj,
        file = paste0(data_out_dir, "obj.rds"))
