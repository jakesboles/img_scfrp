library(tidyverse)
library(scCustomize)
library(Seurat)
library(ggalluvial)

plots_dir <- "plots/06_clustering"
dir.create(plots_dir,
           F, T)

data_dir <- "data/06_clustering"
dir.create(data_dir,
           F, T)

tab_dir <- "tab_data/06_clustering"
dir.create(tab_dir,
           F, T)

obj <- readRDS("/projects/b1169/boles/img_scfrp/data/06_clustering/clustered_integrated_obj.rds")

obj <- obj %>%
  FindClusters(method = "igraph",
               algorithm = 4,
               resolution = 0.3,
               cluster.name = "cluster_low_res",
               graph.name = "RNA_snn")

# DimPlot_scCustom(obj,
#                  group.by = "cluster_low_res",
#                  reduction = "harmony_umap")
# 
# df <- obj@meta.data %>% 
#   dplyr::select(c(cluster_low_res, cluster))
# 
# df <- df %>% 
#   group_by(cluster_low_res, cluster) %>% 
#   summarize(n = n())
# 
# pal <- DiscretePalette_scCustomize(num_colors = length(unique(df$cluster)),
#                                    palette = "varibow")
# 
# ggplot(data = df,
#        aes(axis1 = cluster, axis2 = cluster_low_res, y = n)) +
#   scale_x_discrete(limits = c("Resolution = 2", "Resolution = 0.4")) +
#   xlab("Cluster") +
#   geom_alluvium(aes(fill = cluster)) +
#   geom_stratum() +
#   scale_fill_manual(values = pal) +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#   theme_minimal() + 
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank())

obj <- JoinLayers(obj)

Idents(obj) <- "cluster"

markers1 <- FindAllMarkers(obj)

write.csv(markers1, 
          file = paste0(tab_dir, "/markers_res2.csv"))

Idents(obj) <- "cluster_low_res"

markers2 <- FindAllMarkers(obj)

write.csv(markers2,
          file = paste0(tab_dir, "/markers_res0-4.csv"))

saveRDS(obj,
        file = paste0(data_dir, "clustered_integrated_joined_obj.rds"))