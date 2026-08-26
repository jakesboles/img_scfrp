# Scores each WGCNA module's expression per cell (UCell, kNN-smoothed over
# the Harmony embedding), then fits genotype x treatment models (random
# intercept for batch) per module and plots estimated marginal means with
# post-hoc letters. Loads the whole-cohort object from
# 05_integration_harmony.R's BPCells/RDS pieces (not a monolithic
# obj_consensus.rds) plus the module gene membership table wgcna.R writes,
# and recomputes its own UCell scores directly from raw counts rather than
# reusing wgcna.R's hdWGCNA-internal module scores.

suppressMessages({
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(lme4)
  library(emmeans)
  library(multcomp)
  library(UCell)
  library(BPCells)
})

setwd("/projects/b1169/boles/img_scfrp")

tab_dir <- "tab_data/wgcna/"
plots_dir <- "plots/wgcna/"

# Read in 05's integrated object ---------------------------------------------

counts_mat <- open_matrix_dir("data/05_integration/bpcells_counts")
meta <- readRDS("data/05_integration/metadata.rds")
harmony <- readRDS("data/05_integration/harmony.rds")

# UCell ranks genes within each cell, which needs efficient random access
# across the whole matrix rather than simple column slicing -- same
# "materialize before handing to a tool that isn't BPCells-aware" caution
# already established for DoubletFinder in 03_doubletfinder.R.
counts_mat <- as(counts_mat, "dgCMatrix")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["harmony"]] <- harmony

mods <- read.csv(paste0(tab_dir, "module_members_consensus.csv"))

gene_sets <- list()

module_names <- unique(mods$module)
module_names <- module_names[!(module_names %in% "grey")]

for (i in module_names){
  gene_sets[[i]] <- mods %>%
    filter(module == i) %>%
    dplyr::pull(gene_name)
}

maxrank <- max(unlist(lapply(gene_sets, length)))

obj <- AddModuleScore_UCell(obj,
                            features = gene_sets,
                            maxRank = maxrank)

obj <- SmoothKNN(obj,
                 signature.names = paste0(names(gene_sets), "_UCell"),
                 reduction = "harmony")

# `sample`, not `orig.ident` -- both hold the identical unified sample id
# throughout this pipeline, but `sample` is this project's documented name.
scores <- obj@meta.data %>%
  dplyr::select(matches("UCell_kNN|genotype|treatment|batch|sample"))

write.csv(scores,
          file = paste0(tab_dir, "module_scores_ucell_consensus.csv"))

stats_for_plotting <- list()

cols <- colnames(scores)[str_detect(colnames(scores), "_UCell_kNN")]

colors <- str_remove_all(cols, "_UCell_kNN")

for (j in seq_along(cols)){

  df2 <- scores %>%
    dplyr::rename("active_module" = cols[j])

  lm <- lmer(active_module ~ genotype * treatment + (1|batch),
             data = df2)

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
  facet_wrap(. ~ module,
             scales = "free_y",
             ncol = 5) +
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

# Quick eyeball plot of one module's smoothed score by sample -- not saved,
# meant to be run interactively; swap in whichever module color you want a
# closer look at.
VlnPlot_scCustom(obj,
                 features = c("greenyellow_UCell_kNN"),
                 group.by = "sample",
                 pt.size = 0) +
  NoLegend()
