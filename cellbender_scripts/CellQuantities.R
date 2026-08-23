setwd("/projects/b1042/Gate_Lab/thomas/NeilSC/Coculture/")

folders <- list.files()

x <- list()

reslist2 <- list()

library(rlist)


for(folder in folders){
  
  pool_dir <- paste0(folder, "/outs/per_sample_outs")
  
  files_use <- list.files(pool_dir)
  
  files_use <- files_use[!endsWith(files_use, "NotUsed")]
  
  x[[folder]] <- files_use
  
  reslist <- list()
  
  for(donor in x[[folder]]){
    
    dat <- read.csv(paste0(pool_dir, "/", donor, "/metrics_summary.csv"))
    
    insample <- dat$Metric.Value[1]
    
    resdat <- data.frame("Pool" = folder, "Sample" = donor, "Cells" = insample)
    
    reslist[[donor]] <- resdat
    
  }
  
  reslist2[[folder]] <- list.rbind(reslist)
  
}

resdf <- list.rbind(reslist2)

resdf$Cells <- as.numeric(resdf$Cells)

resdf$Droplets <- floor(resdf$Cells*1.75)

write.csv(resdf, "/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/CellQuantities_012126.csv", quote = F)
