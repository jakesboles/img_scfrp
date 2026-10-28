# QC for the assembled coculture object: raw QC plots, per-sample
# MAD-based thresholds, discard flagging, and a filtered object saved
# on-disk with BPCells. Mirrors monocultures/r_scripts/02_qc.R's method
# exactly (same MAD multipliers, same scCustomize-based percent_mito/
# percent_ribo) -- an explicit, user-confirmed choice to use monoculture's
# established QC approach for this arm rather than the rank-based scheme
# this arm's own lost pipeline used (see CLAUDE.md's "Coculture arm, next
# stage" section for that history). Differences from monoculture's
# version are all about the metadata schema, not the QC method itself: no
# treatment axis, a 3-level (not 4-level) batch, and `condition` is the
# genotype-pair itself (neuron_genotype x microglia_genotype, 10 levels)
# rather than a genotype x treatment cross.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(BPCells)
  library(janitor)
  library(ggbeeswarm)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp/cocultures")

results_dir <- "results/02_qc/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_out_dir <- "data/02_qc/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

message2("Reading in object")

counts <- open_matrix_dir("data/01_obj_assembly/bpcells")
meta <- readRDS("data/01_obj_assembly/metadata.rds")
obj <- CreateSeuratObject(counts = counts, meta.data = meta)

# Add QC metrics and grouping variables ------------------------------------

message2("Organizing object for plotting")

genotype_levels <- c("KOLF", "GALC", "GBA", "GRN", "LRRK2")
microglia_levels <- c("CN", genotype_levels)

# The 10 real genotype-pair conditions in this design, in a fixed display
# order (KOLF-neuron conditions first, then the matched homozygous
# pairs) -- most genotype x genotype combinations never occur, so this
# is an explicit list rather than a full cross like monoculture's
# genotype x treatment `condition`.
condition_levels <- c("KOLF_CN", "KOLF_KOLF", "KOLF_GALC", "KOLF_GBA",
                      "KOLF_GRN", "KOLF_LRRK2", "GALC_GALC", "GBA_GBA",
                      "GRN_GRN", "LRRK2_LRRK2")

obj@meta.data <- obj@meta.data %>%
  mutate(neuron_genotype = factor(neuron_genotype, levels = genotype_levels),
         microglia_genotype = factor(microglia_genotype, levels = microglia_levels),
         # This arm's batch is a 3-level genotype-pair replicate number
         # (see CLAUDE.md) -- wholly independent of the monoculture arm's
         # own 4-level batch, despite the shared column name.
         batch = factor(batch, levels = 1:3, labels = paste("Batch", 1:3)),
         condition = factor(paste0(neuron_genotype, "_", microglia_genotype),
                            levels = condition_levels))

condition_pal <- DiscretePalette_scCustomize(num_colors = length(condition_levels), palette = "varibow")
batch_pal <- Dark2_Pal()[1:3]

message2("Adding cell complexity and mito/ribo percentages")

obj <- obj %>%
  Add_Cell_Complexity() %>%
  Add_Mito_Ribo(species = "human")

# Raw QC plots, faceted by batch --------------------------------------------

message2("Making raw QC plots")

juxt_vln <- function(feature, title, xlab, file_suffix){
  df <- obj@meta.data

  p <- df %>%
    ggplot(aes(y = condition, x = {{feature}})) +
    geom_violin(aes(fill = condition)) +
    scale_fill_manual(values = condition_pal) +
    ggtitle(title) +
    scale_y_discrete(limits = rev) +
    xlab(xlab) +
    facet_wrap(. ~ batch, ncol = 3) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_text(hjust = 1))

  if (str_detect(deparse(substitute(feature)), "RNA")) {
    p <- p +
      scale_x_log10()
  }

  ggsave(p, filename = paste0(results_dir, "raw_", file_suffix, ".png"),
         units = "in", dpi = 600, height = 5, width = 8)
}

juxt_vln(nFeature_RNA, "Genes per cell", "# unique genes (x 1000)", "nfeature")
juxt_vln(nCount_RNA, "UMIs per cell", "# UMI", "ncount")
juxt_vln(log10GenesPerUMI, "Cell complexity", "log10(# genes) / log10(# UMIs)", "complexity")
juxt_vln(percent_mito, "Mitochondrial gene content", "% mitochondrial genes", "mito")

# Per-sample median stats, faceted by batch ---------------------------------

message2("Computing raw median stats")

sample_lookup <- obj@meta.data %>%
  distinct(sample, neuron_genotype, microglia_genotype, batch, condition)

median_stats_plot <- function(stats, filename){
  stats %>%
    left_join(sample_lookup, by = "sample") %>%
    pivot_longer(!c(sample, neuron_genotype, microglia_genotype, batch, condition),
                 names_to = "stat", values_to = "value") %>%
    ggplot(aes(y = condition, x = value)) +
    geom_point(aes(fill = batch), shape = 21, size = 3) +
    scale_fill_manual(values = batch_pal) +
    scale_y_discrete(limits = rev) +
    facet_wrap(. ~ stat, ncol = 5, scales = "free_x") +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_blank(),
          legend.title = element_blank())
  ggsave(paste0(results_dir, filename), units = "in", dpi = 600, height = 4, width = 12)
}

stats <- Median_Stats(obj, group.by = "sample")
counts_tbl <- obj@meta.data %>%
  group_by(sample) %>%
  dplyr::summarize(Cell_count = n()) %>%
  adorn_totals(name = "Totals (All Cells)")
stats <- stats %>% left_join(counts_tbl, by = "sample")

write.csv(stats, file = paste0(results_dir, "raw_median_stats.csv"), row.names = F)

# adorn_totals()'s summary row has no sample_lookup match and would plot
# as NA -- drop it rather than relying on a hardcoded row index.
median_stats_plot(stats %>% filter(sample != "Totals (All Cells)"),
                  "raw_median_stats.png")

# Set MAD-based thresholds per sample ---------------------------------------

message2("Setting thresholds")

obj@meta.data <- obj@meta.data %>%
  mutate(log_nFeature = log10(nFeature_RNA),
         log_nCount = log10(nCount_RNA))

samples <- sort(unique(as.character(obj$sample)))

thresh_df <- data.frame(sample = samples,
                        umi_med = rep(NA, length(samples)),
                        umi_mad = rep(NA, length(samples)),
                        feature_med = rep(NA, length(samples)),
                        feature_mad = rep(NA, length(samples)),
                        mito_med = rep(NA, length(samples)),
                        mito_mad = rep(NA, length(samples)))

meta <- obj@meta.data

for (i in seq_along(samples)){
  message(paste0("Getting cutoffs for ", samples[i]))

  df <- meta %>% filter(sample == samples[i])

  thresh_df$umi_med[i] <- median(df$log_nCount)
  thresh_df$umi_mad[i] <- stats::mad(df$log_nCount)
  thresh_df$feature_med[i] <- median(df$log_nFeature)
  thresh_df$feature_mad[i] <- stats::mad(df$log_nFeature)
  thresh_df$mito_med[i] <- median(df$percent_mito)
  thresh_df$mito_mad[i] <- stats::mad(df$percent_mito)
}

# Same thresholds as monocultures/r_scripts/02_qc.R, per the user's
# explicit choice to use this arm's established QC method -- not
# re-derived for coculture either, same inherited-not-validated caveat
# as monoculture's own version carries.
thresh_df <- thresh_df %>%
  mutate(umi_lower = if_else(umi_med - 2 * umi_mad < log10(2000), log10(2000), umi_med - 2 * umi_mad),
         feature_lower = if_else(feature_med - 2 * feature_mad < log10(2000), log10(2000), feature_med - 2 * feature_mad),
         mito_upper = 5,
         umi_upper = umi_med + 2.5 * umi_mad,
         feature_upper = feature_med + 2.5 * feature_mad)

meta <- meta %>%
  rownames_to_column(var = "cell") %>%
  left_join(thresh_df, by = "sample") %>%
  mutate(mito_discard = if_else(percent_mito > mito_upper, T, F),
         umi_discard = if_else(log_nCount < umi_lower | log_nCount > umi_upper, T, F),
         gene_discard = if_else(log_nFeature < feature_lower | log_nFeature > feature_upper, T, F))

# Discard-flagged plots, faceted by batch ------------------------------------

message2("Making discard-flag plots")

discard_plot <- function(feature, discard_flag, title, xlab, file_suffix){
  p <- meta %>%
    ggplot(aes(y = condition, x = {{feature}})) +
    geom_quasirandom(aes(color = {{discard_flag}}), size = 0.2) +
    scale_color_manual(values = c("gray60", "red")) +
    ggtitle(title) +
    scale_y_discrete(limits = rev) +
    xlab(xlab) +
    facet_wrap(. ~ batch, ncol = 3) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_text(hjust = 1))

  if (str_detect(deparse(substitute(feature)), "RNA")) {
    p <- p +
      scale_x_log10()
  }

  ggsave(p, filename = paste0(results_dir, "thresholds_", file_suffix, ".png"),
         units = "in", dpi = 600, height = 5, width = 8)
}

discard_plot(nFeature_RNA, gene_discard, "Discarded due to low/high genes?",
            "# unique genes (x 1000)", "nfeature")
discard_plot(nCount_RNA, umi_discard, "Discarded due to low/high UMIs?",
            "log10(# UMIs)", "ncount")
discard_plot(percent_mito, mito_discard, "Discarded due to high mitochondrial gene content?",
            "% mitochondrial genes", "mito")

# Mark cells for removal and save stats --------------------------------------

message2("Marking cells and saving statistics")

meta <- meta %>%
  mutate(discard = if_else(mito_discard == T | umi_discard == T | gene_discard == T, T, F))

stats <- meta %>%
  group_by(sample, discard) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = "discard", values_from = "n") %>%
  mutate(`TRUE` = if_else(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(retained_percent = (`FALSE` / (`FALSE` + `TRUE`)) * 100) %>%
  dplyr::rename("retained_count" = "FALSE", "discarded_count" = "TRUE")

write.csv(thresh_df, file = paste0(results_dir, "cutoffs.csv"), row.names = F)
write.csv(stats, file = paste0(results_dir, "filtered_counts.csv"), row.names = F)

barcodes <- meta %>%
  select(cell, nCount_RNA, nFeature_RNA, percent_mito,
         mito_discard, umi_discard, gene_discard, discard)
write.csv(barcodes, file = paste0(results_dir, "barcode_qc.csv"), row.names = F)

# Filter object and save as BPCells -----------------------------------------

message2("Filtering object of low quality cells")

obj@meta.data <- obj@meta.data %>%
  rownames_to_column(var = "cell") %>%
  left_join(meta %>% select(cell, mito_discard, umi_discard, gene_discard, discard),
           by = "cell") %>%
  column_to_rownames(var = "cell")

obj_filtered <- subset(obj, discard == F)

message2("Computing filtered median stats")

stats <- Median_Stats(obj_filtered, group.by = "sample")
counts_tbl <- obj_filtered@meta.data %>%
  group_by(sample) %>%
  dplyr::summarize(Cell_count = n()) %>%
  adorn_totals(name = "Totals (All Cells)")
stats <- stats %>% left_join(counts_tbl, by = "sample")

write.csv(stats, file = paste0(results_dir, "filtered_median_stats.csv"), row.names = F)
median_stats_plot(stats %>% filter(sample != "Totals (All Cells)"),
                  "filtered_median_stats.png")

message2("Saving filtered object as BPCells on-disk matrix")

write_matrix_dir(mat = obj_filtered[["RNA"]]$counts,
                 dir = paste0(data_out_dir, "bpcells"))

message2("Saving filtered metadata as RDS")

saveRDS(obj_filtered@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts <- open_matrix_dir("data/02_qc/bpcells")
#   meta <- readRDS("data/02_qc/metadata.rds")
#   obj <- CreateSeuratObject(counts = counts, meta.data = meta)
