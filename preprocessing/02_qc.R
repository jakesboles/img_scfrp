library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(janitor)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

data_in_dir <- paste0(proj_dir, "data/01_obj_creation/")

plots_dir <- paste0(proj_dir, "plots/02_qc/")
dir.create(plots_dir, showWarnings = F,
           recursive = T)

csv_dir <- paste0(proj_dir, "tab_data/02_qc/")
dir.create(csv_dir, showWarnings = F,
           recursive = T)

data_out_dir <- paste0(proj_dir, "data/02_qc/")
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

obj <- readRDS(paste0(data_in_dir, "obj.rds"))


obj@meta.data %>%
  str()

obj@meta.data <- obj@meta.data %>%
  separate(orig.ident,
           into = c("batch", "genotype", "treatment"),
           sep = "_",
           remove = F) %>%
  mutate(sample = paste0(genotype, "_", treatment))

obj@meta.data <- obj@meta.data %>%
  mutate(sample = factor(sample,
                         levels = paste0(rep(c("KOLF", "GALC", "GBA", "GRN", "LRRK2"), each = 3), 
                                         "_", 
                                         rep(c("control", "myelin", "asyn"), times = 5))),
         genotype = factor(genotype,
                           levels = c("KOLF", "GALC", "GBA", "GRN", "LRRK2")),
         treatment = factor(treatment,
                            levels = c("control", "myelin", "asyn")))

treatment_pal <- JCO_Four()[1:3]
genotype_pal <- c("red", "yellow", "green", "blue", "purple")
batch_pal <- Dark2_Pal()[1:4]
sample_pal <- DiscretePalette_scCustomize(num_colors = 15,
                                          palette = "varibow")

obj <- obj %>%
  Add_Cell_Complexity() %>%
  Add_Mito_Ribo(species = "human")

obj@meta.data %>%
  str()

juxt_vln <- function(feature,
                     title,
                     xlab){
  
  df <- obj@meta.data %>%
    mutate(batch = factor(batch,
                          labels = c("Batch 1", "Batch 2", "Batch 3", "Batch 4")),
           nFeature_RNA = nFeature_RNA/1000,
           nCount_RNA = log10(nCount_RNA))
  
  # if (deparse(substitute(feature)) %in% c("nFeature_RNA", "nCount_RNA")){
  #   df[, deparse(substitute(feature))] <- log10(df[, deparse(substitute(feature))])
  # }
  
  p <- df %>%
    ggplot(aes(y = sample,
               x = {{feature}})) + 
    geom_quasirandom(aes(color = sample),
                     size = 0.2) + 
    scale_color_manual(values = sample_pal) + 
    ggtitle(title) +
    scale_y_discrete(limits = rev) +
    # scale_x_continuous(limits = c(0, max),
    #                    expand = c(0, NA)) + 
    xlab(xlab) + 
    facet_wrap(. ~ batch,
               ncol = 4) + 
    theme_bw() + 
    theme(axis.text = element_text(color = "black"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_text(hjust = 1))
  
  # max <- max(df[, deparse(substitute(feature))])
  # 
  # plot <- function(pool){
  #   title <- str_to_title(pool)
  #   
  #   p <- df %>%
  #     filter(batch == pool) %>%
  #     ggplot(aes(y = sample,
  #                x = {{feature}})) + 
  #     geom_quasirandom(aes(color = sample),
  #                      size = 0.2) + 
  #     scale_color_manual(values = sample_pal) + 
  #     ggtitle(title) + 
  #     scale_y_discrete(limits = rev) +
  #     scale_x_continuous(limits = c(0, max),
  #                        expand = c(0, NA)) + 
  #     theme_bw() + 
  #     theme(axis.text = element_text(color = "black"),
  #           legend.position = "none",
  #           plot.title = element_text(hjust = 0.5),
  #           axis.title = element_blank(),
  #           axis.text.y = element_text(hjust = 0.5))
  #         
  #   return(p)
  # }
  # 
  # p1 <- plot("pool1")
  # p2 <- plot("pool2")  
  #   # theme(axis.text.y = element_blank(),
  #   #       axis.ticks.y = element_blank())
  # p3 <- plot("pool3") 
  #   # theme(axis.text.y = element_blank(),
  #   #       axis.ticks.y = element_blank())
  # p4 <- plot("pool4") 
  #   # theme(axis.text.y = element_blank(),
  #   #       axis.ticks.y = element_blank())
  # 
  # p <- p1 + p2 + p3 + p4 + 
  #   plot_layout(ncol = 4,
  #               axes = "collect") + 
  #   plot_annotation(title = deparse(substitute(feature)),
  #                   theme = theme(plot.title = element_text(size = 16, hjust = 0.5),
  #                                 axis.text.y = element_text(hjust = 1)))
  
  ggsave(p,
         filename = paste0(plots_dir, "raw_",
                           str_to_lower(deparse(substitute(feature))), 
                           ".png"),
         units = "in", 
         dpi = 600,
         height = 5,
         width = 10)
}

juxt_vln(nFeature_RNA,
         "Genes per cell",
         "# unique genes (x 1000)")
juxt_vln(nCount_RNA,
         "UMIs per cell",
         "log10(# UMIs)")
juxt_vln(log10GenesPerUMI,
         "Cell complexity",
         "log10(# genes) / log10(# UMIs")
juxt_vln(percent_mito,
         "Mitochondrial gene content",
         "% mitochondrial genes")


stats <- Median_Stats(obj,
                      group_by = "orig.ident")

counts <- obj@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::summarize(Cell_count = n()) %>%
  adorn_totals(name = "Totals (All Cells)")

stats <- stats %>%
  left_join(counts, by = "orig.ident")

stats <- stats[-61, ] %>%
  separate(orig.ident,
           into = c("batch", "genotype", "treatment"),
           sep = "_",
           remove = F) %>%
  mutate(sample = paste0(genotype, "_", treatment))

stats %>%
  mutate(sample = factor(sample,
                         levels = paste0(rep(c("KOLF", "GALC", "GBA", "GRN", "LRRK2"), each = 3), 
                                         "_", 
                                         rep(c("control", "myelin", "asyn"), times = 5)))) %>%
  mutate(batch = factor(batch,
                        labels = c("Batch 1", "Batch 2", "Batch 3", "Batch 4"))) %>%
  pivot_longer(!c(orig.ident, batch, genotype, treatment, sample),
               names_to = "stat", values_to = "value") %>%
  ggplot(aes(y = sample, 
             x = value)) + 
  geom_point(aes(fill = batch),
             shape = 21,
             size = 3) + 
  scale_fill_manual(values = batch_pal) + 
  scale_y_discrete(limits = rev) +
  facet_wrap(. ~ stat,
             ncol = 5,
             scales = "free_x") + 
  theme_bw() + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank())
ggsave(paste0(plots_dir, "raw_median_stats.png"),
       units = "in", dpi = 600,
       height = 4, width = 12)

write.csv(stats, 
          file = paste0(csv_dir, "raw_median_stats.csv"),
          row.names = F)


# Mark cells based on MAD -----------------------------------------------------

obj@meta.data <- obj@meta.data %>%
  mutate(log_nFeature = log10(nFeature_RNA),
         log_nCount = log10(nCount_RNA))

samples <- unique(obj$orig.ident)

thresh_df <- data.frame(orig.ident = samples,
                        umi_med = c(rep(NA, length(samples))),
                        umi_mad = c(rep(NA, length(samples))),
                        feature_med = c(rep(NA, length(samples))),
                        feature_mad = c(rep(NA, length(samples))),
                        mito_med = c(rep(NA, length(samples))),
                        mito_mad = c(rep(NA, length(samples))))

meta <- obj@meta.data

for (i in seq_along(samples)){
  message(paste0("Getting centrality measures for ", samples[i]))
  
  df <- meta %>%
    filter(orig.ident == samples[i])
  
  thresh_df$umi_med[i] <- median(df$log_nCount)
  
  thresh_df$umi_mad[i] <- stats::mad(df$log_nCount)
  
  thresh_df$feature_med[i] <- median(df$log_nFeature)
  
  thresh_df$feature_mad[i] <- stats::mad(df$log_nFeature)
  
  thresh_df$mito_med[i] <- median(df$percent_mito)
  
  thresh_df$mito_mad[i] <- stats::mad(df$percent_mito)
}

# Pretty "strict" cutoffs of 2 x MAD still discards very few cells based 
# on number of counts and unique genes
thresh_df <- thresh_df %>%
  mutate(umi_lower = umi_med - 2*umi_mad,
         feature_lower = feature_med - 2*feature_mad,
         mito_upper = 5,
         umi_upper = umi_med + 2*umi_mad,
         feature_upper = feature_med + 2*feature_mad)

meta <- meta %>%
  rownames_to_column(var = "cell") %>%
  left_join(thresh_df,
            by = "orig.ident") %>% 
  mutate(mito_discard = if_else(percent_mito > mito_upper, T, F),
         umi_discard = if_else(log_nCount < umi_lower | 
                                 log_nCount > umi_upper, T, F),
         gene_discard = if_else(log_nFeature < feature_lower |
                                  log_nFeature > feature_upper, T, F))

discard_plot <- function(feature, 
                         discard_flag,
                         title, 
                         xlab){
  p <- meta %>%
    mutate(batch = factor(batch,
                          labels = c("Batch 1", "Batch 2", "Batch 3", "Batch 4")),
           nFeature_RNA = nFeature_RNA/1000,
           nCount_RNA = log10(nCount_RNA)) %>%
    # dplyr::select(c(sample, batch, percent_mito, log_nFeature, log_nCount,
    #                 mito_discard, umi_discard, gene_discard)) %>%
    # pivot_longer(c(log_nFeature, log_nCount, percent_mito),
    #              names_to = "stat",
    #              values_to = "value") %>%
    ggplot(aes(y = sample, 
               x = {{feature}})) + 
    geom_quasirandom(aes(color = {{discard_flag}}),
                     size = 0.2) +
    scale_color_manual(values = c("gray60", "red")) +
    ggtitle(title) +
    scale_y_discrete(limits = rev) +
    # scale_x_continuous(limits = c(0, max),
    #                    expand = c(0, NA)) + 
    xlab(xlab) + 
    facet_wrap(. ~ batch,
               ncol = 4) + 
    guides(color = guide_legend(override.aes = list(size = 2))) +
    theme_bw() + 
    theme(axis.text = element_text(color = "black"),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_text(hjust = 1))
  
  ggsave(p,
         filename = paste0(plots_dir, "thresholds_",
                           str_to_lower(deparse(substitute(feature))),
                           ".png"),
         units = "in",
         dpi = 600,
         height = 5,
         width = 10)
}

discard_plot(nFeature_RNA,
             gene_discard,
             "Discarded due to low/high genes?",
             "# unique genes (x 1000)")
discard_plot(nCount_RNA,
             umi_discard,
             "Discarded due to low/high UMIs?",
             "log10(# UMIs)")
discard_plot(percent_mito,
             mito_discard,
             "Discarded due to high mitochondrial gene content?",
             "% mitochondrial genes")

meta <- meta %>%
  mutate(discard = if_else(mito_discard == T |
                             umi_discard == T |
                             gene_discard == T,
                           T, F))


stats <- meta %>%
  group_by(orig.ident, discard) %>% 
  dplyr::summarize(n = n()) %>%
  pivot_wider(names_from = "discard",
              values_from = "n") %>%
  mutate(`TRUE` = if_else(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(retained_percent = (`FALSE` / (`FALSE` + `TRUE`)) * 100) %>%
  # separate(orig.ident, 
  #          into = c("id", "tissue"),
  #          sep = "_") %>%
  # mutate(tissue = case_when(tissue == "b" ~ "Motor cortex",
  #                           tissue == "s" ~ "Cervical spinal cord",
  #                           tissue == "m" ~ "Skeletal muscle")) %>%
  dplyr::rename("retained_count" = "FALSE",
                "discarded_count" = "TRUE")

write.csv(thresh_df,
          file = paste0(csv_dir, "cutoffs.csv"),
          row.names = F)

write.csv(stats,
          file = paste0(csv_dir, "filtered_counts.csv"),
          row.names = F)

barcodes <- meta %>%
  dplyr::select(c(cell, nCount_RNA, nFeature_RNA, percent_mito,
                  mito_discard, umi_discard, gene_discard, discard))

write.csv(barcodes,
          file = paste0(csv_dir, "barcode_qc.csv"),
          row.names = F)

# Filtering out cells -----------------------------------------------------

obj@meta.data <- obj@meta.data %>%
  rownames_to_column(var = "cell") %>%
  left_join(meta) %>%
  column_to_rownames(var = "cell")

obj_filtered <- subset(obj, discard == F)

stats <- Median_Stats(obj_filtered,
                      group_by = "orig.ident")

counts <- obj_filtered@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::summarize(Cell_count = n()) %>%
  adorn_totals(name = "Totals (All Cells)")

stats <- stats %>%
  left_join(counts, by = "orig.ident")

stats <- stats[-61, ] %>%
  separate(orig.ident,
           into = c("batch", "genotype", "treatment"),
           sep = "_",
           remove = F) %>%
  mutate(sample = paste0(genotype, "_", treatment))

stats %>%
  mutate(sample = factor(sample,
                         levels = paste0(rep(c("KOLF", "GALC", "GBA", "GRN", "LRRK2"), each = 3), 
                                         "_", 
                                         rep(c("control", "myelin", "asyn"), times = 5)))) %>%
  mutate(batch = factor(batch,
                        labels = c("Batch 1", "Batch 2", "Batch 3", "Batch 4"))) %>%
  pivot_longer(!c(orig.ident, batch, genotype, treatment, sample),
               names_to = "stat", values_to = "value") %>%
  ggplot(aes(y = sample, 
             x = value)) + 
  geom_point(aes(fill = batch),
             shape = 21,
             size = 3) + 
  scale_fill_manual(values = batch_pal) + 
  scale_y_discrete(limits = rev) +
  facet_wrap(. ~ stat,
             ncol = 5,
             scales = "free_x") + 
  theme_bw() + 
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.title = element_blank())
ggsave(paste0(plots_dir, "filtered_median_stats.png"),
       units = "in", dpi = 600,
       height = 4, width = 12)

write.csv(stats, 
          file = paste0(csv_dir, "filtered_median_stats.csv"),
          row.names = F)

saveRDS(obj_filtered,
        file = paste0(data_out_dir, "filtered_obj.rds"))
