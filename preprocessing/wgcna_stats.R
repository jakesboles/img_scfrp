library(Seurat)
library(scCustomize)
library(tidyverse)
library(lme4)
library(emmeans)
library(multcomp)
library(UCell)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

tab_dir <- paste0(proj_dir, "tab_data/wgcna/")
plots_dir <- paste0(proj_dir, "plots/wgcna/")

obj <- readRDS(paste0(proj_dir, "data/wgcna/obj_consensus.rds"))

mods <- read.csv(paste0(proj_dir, "tab_data/wgcna/module_members_consensus.csv"))

gene_sets <- list()

module_names <- unique(mods$module)
module_names <- module_names[!(module_names %in% "grey")]

for (i in module_names){
  gene_sets[[i]] <- mods %>%
    filter(module == i) %>%
    dplyr::pull(gene_name)
}

maxrank <- max(unlist(lapply(gene_sets, length)))

# sub <- subset(seurat_obj,
#               PredictedCellType == cell_file)

obj <- AddModuleScore_UCell(obj,
                            features = gene_sets,
                            maxRank = maxrank)

obj <- SmoothKNN(obj,
                 signature.names = paste0(names(gene_sets), "_UCell"),
                 reduction = "harmony_pca")

scores <- obj@meta.data %>%
  dplyr::select(matches("UCell_kNN|genotype|treatment|batch|orig.ident"))

write.csv(scores,
          file = paste0(proj_dir, "tab_data/wgcna/module_scores_ucell_consensus.csv"))

stats_for_plotting <- list()

cols <- colnames(scores)[str_detect(colnames(scores), "_UCell_kNN")]

colors <- str_remove_all(cols, "_UCell_kNN")

for (j in seq_along(cols)){
  
  df2 <- scores %>%
    dplyr::rename("active_module" = cols[j])
  
  lm <- lmer(active_module ~ genotype * treatment + (1|batch),
             data = df2)
  
  # summary(lm)
  
  # write.csv(summary(lm)$coefficients,
  #           file = paste0(tab_dir, colors[j], "_lm.csv"))
  
  emm <- emmeans(lm, pairwise ~ treatment | genotype)
  
  stats_for_plotting[[j]] <- multcomp::cld(emm, Letters = letters) %>%
    mutate(.group = str_remove_all(.group, " "))
  
  sink(paste0(tab_dir, "lmer_stats_", colors[j], ".txt"))
  cat("Model output")
  cat("\n")
  print(summary(lm))
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat("\n")
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat("Joint tests on the model")
  cat("\n")
  print(joint_tests(lm))
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat("\n")
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat("")
  cat("Post-hoc comparisons & estimated marginal means")
  cat("\n")
  print(emm)
  sink()
  closeAllConnections()
  
  stats_for_plotting[[j]]$module <- colors[j]
  colnames(stats_for_plotting[[j]]) <- c("treatment", "genotype", "emmean", "SE", "df", "lower_CL", "upper_CL", ".group", "module")
}

stats_df <- list_rbind(stats_for_plotting)

p <- stats_df %>%
  ggplot(aes(x = genotype,
             y = emmean)) + 
  geom_crossbar(aes(ymin = lower_CL,
                    ymax = upper_CL,
                    fill = treatment),
                position = position_dodge(width = 1)) + 
  geom_text(aes(label = .group,
                y = upper_CL * 1.03,
                color = treatment),
            position = position_dodge(width = 1)) + 
  # scale_fill_manual(values = c("dodgerblue1", "magenta1", "chartreuse1")) + 
  facet_wrap(. ~ module,
             scales = "free_y",
             ncol = 5) +
  # ggtitle(cell_type) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"))
ggsave(p,
       filename = paste0(plots_dir, "lmer_stats.png"),
       units = "in", dpi = 600,
       height = 7,
       width = 12)

VlnPlot_scCustom(obj,
                 features = c("greenyellow_UCell_kNN"),
                 group.by = "sample",
                 pt.size = 0) + 
  NoLegend()

# modules <- read.csv(paste0(tab_dir, "hmes_consensus.csv"))
# 
# meta <- obj@meta.data
# 
# meta <- meta %>%
#   rownames_to_column(var = "X") %>%
#   left_join(modules, by = "X")
# 
# module_names <- colnames(modules)[-1]
# module_names <- module_names[str_detect(module_names, "grey", negate = T)]

# for (i in "magenta"){
#   df <- meta %>%
#     dplyr::rename("active" = i) 
#   
#   fit <- lmer(active ~ genotype * treatment + (1|orig.ident),
#               data = df)
#   
#   emm <- emmeans(fit,
#                  ~ genotype * treatment)
#   
#   cld <- multcomp::cld(emm, 
#              Letters = letters)
#   
#   cld %>%
#     mutate(.group = str_remove_all(.group, " ")) %>%
#     ggplot(aes(x = genotype,
#                y = emmean)) + 
#     geom_crossbar(aes(ymin = asymp.LCL,
#                       ymax = asymp.UCL,
#                       color = treatment),
#                   position = position_dodge(width = 1))
#   
#   
# }

modules <- modules %>%
  column_to_rownames(var = "X")

obj <- AddMetaData(obj,
                  modules)

VlnPlot_scCustom(obj,
                 features = "turquoise",
                 group.by = "sample")

gene_sets <- list()

module_names <- unique(modules$module)
module_names <- module_names[!(module_names %in% "grey")]

for (i in module_names){
  gene_sets[[i]] <- modules %>%
    filter(module == i) %>%
    pull(gene_name)
}

maxrank <- max(unlist(lapply(gene_sets, length)))


obj <- AddModuleScore_UCell(obj,
                            features = gene_sets,
                            maxRank = maxrank)

obj <- SmoothKNN(obj,
                 signature.names = paste0(names(gene_sets), "_UCell"),
                 reduction = "")

