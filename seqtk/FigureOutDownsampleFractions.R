suppressMessages({
  library(plyr)
  library(data.table)
  library(tidyverse) 
  library(Seurat)
  library(rlist)
  library(doMC)
  library(UpSetR)
  library(DElegate)
  library(DESeq2)
  library(edgeR)
  library(sparseMatrixStats)
  library(SeuratData)
  library(ComplexHeatmap)
  library(ggrepel)
  library(rsvd)
  library(plotly)
  library(ggthemes)
  library(ggplot2)
  library(cowplot)
  library(ggsci)
  library(ggdark)
  library(viridis)
  library(car)
  library(e1071)
  library(parallel)
  library(aricode)
  library(knitr)
  library(rlist)
  # library(SCP)
  # library(stardust)
  library(DescTools)
  library(ppcor)
  library(fastDummies)
  # library(SPOTlight)
  library(gridExtra)
  library(MAST)
  library(Matrix.utils)
  library(msigdbr)
  library(fgsea)
  library(SCPA)
})

options(future.globals.maxSize=1048576000000)
load("/projects/b1169/projects/sea_ad_hypothalamus/results/preprocessing/qc/out_TW_05-04-2023/helperfunctions.RData")
setwd("/gpfs/projects/b1042/Gate_Lab/boles/img_scfrp/cellranger")

batches <- list.files()

res <- res2 <- list()

for(batch in batches){
  
  samples <- list.files(paste0(batch, "/outs/per_sample_outs"))
  
  for(sample in samples){
    
    dat <- read.csv(paste0(batch, "/outs/per_sample_outs/", sample, "/metrics_summary.csv"))
    
    res[[batch]][[sample]] <- data.frame("Batch" = batch, 
                                         "Cells" = as.numeric(dat$Metric.Value[1]), 
                                         "Reads" = as.numeric(dat$Metric.Value[16]),
                                         "ReadsInCells" = as.numeric(dat$Metric.Value[7]))
    
  }
  
  res[[batch]] <- list.rbind(res[[batch]])

  res2[[batch]] <- data.frame("Batch" = batch, 
                              "Cells" = sum(res[[batch]]$Cells),
                              "Reads" = sum(res[[batch]]$Reads),
                              "ReadsInCells" = sum(res[[batch]]$ReadsInCells))
  
}

res2 <- list.rbind(res2)

res2$RPC <- res2$ReadsInCells/res2$Cells

targetRPC <- min(res2$RPC)

res2$DF <- targetRPC/res2$RPC

write.csv(res2, "DownSampleFractions.csv")