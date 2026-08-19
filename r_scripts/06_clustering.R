library(Seurat)
library(scCustomize)
library(dittoSeq)
library(reticulate)
library(patchwork)
library(speckle)

setwd("/projects/b1169/boles/img_scfrp")

obj <- readRDS("data/05_integration/harmony_obj.rds")

plots_dir <- "plots/06_clustering"
dir.create(plots_dir,
           F, T)

data_dir <- "data/06_clustering"
dir.create(data_dir,
           F, T)

obj
colnames(obj@meta.data)

# DimPlot_scCustom(obj,
#                  reduction = "harmony_umap",
#                  label = F,
#                  group.by = "sample")

obj <- obj %>%
  FindClusters(method = "igraph",
               algorithm = 4,
               resolution = 2,
               cluster.name = "cluster",
               graph.name = "RNA_snn")

p1 <- DimPlot_scCustom(obj,
                       reduction = "harmony_umap",
                       label = F,
                       group.by = "cluster")
ggsave(p1,
       filename = paste0(plots_dir, "/cluster_umap.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

p2 <- dittoBarPlot(obj,
                   var = "genotype",
                   group.by = "cluster")
ggsave(p2,
       filename = paste0(plots_dir, "/clusters_by_genotype.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

p3 <- dittoBarPlot(obj,
                   var = "batch",
                   group.by = "cluster")
ggsave(p3,
       filename = paste0(plots_dir, "/clusters_by_batch.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

p4 <- dittoBarPlot(obj,
                   var = "treatment",
                   group.by = "cluster")
ggsave(p4,
       filename = paste0(plots_dir, "/clusters_by_treatment.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

# p4 <- DimPlot_scCustom(obj,
#                        reduction = "harmony_umap",
#                        group.by = "batch",
#                        label = F,
#                        colors_use = JCO_Four())
# 
# p1 + p2 + p3 + p4 + 
#   plot_layout(ncol = 2)

# Cluster_Highlight_Plot(obj, cluster_name = "19",
#                        reduction = "harmony_umap")

props <- getTransformedProps(obj$cluster, obj$orig.ident, transform="logit")

genotype <- rep(unique(obj$genotype), each = 3)
treatment <- rep(unique(obj$treatment), each = 5)

model <- model.matrix(~ 0 + genotype + treatment)

propeller.anova(prop.list=props, 
                design=model, 
                coef = c(1,2,3,4,5), 
                robust=TRUE, 
                trend=FALSE, 
                sort=TRUE)

propeller(clusters = obj$cluster,
          sample = obj$orig.ident,
          group = obj$genotype)

saveRDS(obj,
        file = paste0(data_dir, "/clustered_integrated_obj.rds"))
