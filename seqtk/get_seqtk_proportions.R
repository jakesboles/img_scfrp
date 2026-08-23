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

setwd("/gpfs/projects/b1169/boles/img_scfrp/")

batches <- list.dirs("cellranger",
                     recursive = F,
                     full.names = F)

numbers <- list()

for(batch in batches){
  
  qc <- read.csv(paste0("cellranger/", batch, "/outs/qc_library_metrics.csv"))
  
  reads <- qc %>% 
    filter(Metric.Name == "Number of reads" & 
             Grouped.By == "Fastq ID") %>% 
    pull(Metric.Value) %>% 
    as.numeric() %>% 
    sum()
  
  rpc <- qc %>%
    filter(Metric.Name == "Mean reads per cell") %>% 
    pull(Metric.Value) %>% 
    as.numeric()
  
  cells <- qc %>% 
    filter(Metric.Name == "Cells") %>% 
    pull(Metric.Value) %>%
    as.numeric()
  
  numbers[[batch]] <- tibble(
    reads = reads,
    rpc = rpc,
    cells = cells,
    batch = batch
  )
  
}

numbers <- list_rbind(numbers)

numbers$target <- min(numbers$rpc)

numbers <- numbers %>% 
  mutate(target_fraction = target / rpc) %>% 
  dplyr::select(c(batch, target_fraction)) %>% 
  filter(target_fraction != 1)

write.table(numbers, "seqtk/downsample_targets.txt",
          row.names = F,
          col.names = F,
          quote = F,
          sep = ",")
