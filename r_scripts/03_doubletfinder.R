# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
  library(scCustomize)
  library(ggpubr)
  # library(BPCells)
  library(patchwork)
})

# Script again adapted from Anne's projects

proj_dir <- "/projects/b1169/boles/img_scfrp/"

plots_dir <- paste0(proj_dir, "plots/03_doubletfinder/")
dir.create(plots_dir,
           showWarnings = F, recursive = T)

data_out_dir <- paste0(proj_dir, "data/03_doubletfinder/")
dir.create(data_out_dir,
           showWarnings = F, recursive = T)

obj <- readRDS(paste0(proj_dir, "data/02_qc/filtered_obj.rds"))

# Define function to run DoubletFinder
run_doubletfinder <- function(s, doublet_rate) {
  
  # Standard normalization and scaling 
  s <- s %>% 
    JoinLayers() %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData()
  
  # Default Seurat clustering and TSNE 
  s <- s %>% 
    RunPCA() %>% 
    FindNeighbors(dims = 1:10) %>% 
    FindClusters()
  s <- RunUMAP(s, dims = 1:10) # had to set check_duplicates to FALSE for some samples
  
  p <- DimPlot(s,
               reduction = "umap",
               group.by = "orig.ident")
  ggsave(p,
         filename = paste0(plots_dir, unique(s$batch), "_umap.png"),
         units = "in", dpi = 600,
         height = 5, width = 10)
  
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # Run DoubletFinder without homotypic adjustment 
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.unadj"
  
  # Run DoubletFinder with homotypic adjustment
  pANN <- colnames(s@meta.data)[grep("^pANN", colnames(s@meta.data))]
  print(pANN)
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)
  
  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.adj"
  
  # Plot unadjusted vs. adjusted doublets in TSNE coordinates
  plt <- ggarrange(
    (DimPlot_scCustom(s, 
                      reduction = "umap", 
                      group.by = "DF.unadj", 
                      pt.size = 1, 
                      shuffle = TRUE, 
                      alpha = 0.6) + 
       ggtitle("Unadjusted")),
    (DimPlot_scCustom(s, 
                      reduction = "umap", 
                      group.by = "DF.adj", 
                      pt.size = 1, 
                      shuffle = TRUE, 
                      alpha = 0.6) +
       ggtitle("Adjusted for Homotypic Proportion")),
    ncol = 2, nrow = 1,
    common.legend = T, legend = "right"
  )
  
  ggsave(plt,
         filename = paste0(plots_dir, unique(s$batch), "_doublet_umap.png"),
         units = "in", dpi = 600,
         height = 5, width = 10,
         bg = "white")
  
  return(s)
}

# Split into sample objects ----
obj_list <- SplitObject(obj, split.by = "batch")

# Load pre-filter stats to get counts for each sample
# This way, each sample gets a customized doublet formation rate depending on the # of cells detected 

counts <- read.csv("/projects/b1169/boles/img_scfrp/tab_data/02_qc/raw_median_stats.csv")

counts <- counts %>%
  dplyr::group_by(batch) %>%
  summarize(cell_count = sum(Cell_count))

counts$batch %in% names(obj_list)

idx <- match(counts$batch, names(obj_list))
counts <- counts[idx, ] %>%
  pull(cell_count)

doublet_rates <- (counts / 8000) * 0.004

# Run DoubletFinder on all samples ----
for (s in seq_along(obj_list)) {
  message(paste0(
    paste(rep("~", 40), collapse = ""),
    "Running ", names(obj_list)[s],
    paste(rep("~", 40), collapse = "")
  ))
  
  obj_list[[s]] <- run_doubletfinder(obj_list[[s]], doublet_rates[s])
  saveRDS(obj_list[[s]], 
          paste0(data_out_dir, names(obj_list)[s], ".rds"))
  
}
