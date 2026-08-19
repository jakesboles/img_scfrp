library(tidyverse)
library(fgsea)

setwd("/projects/b1169/boles/img_scfrp")

dir.create("tab_data/deseq2_gsea/")

url <- 'https://github.com/jackbibby1/SCPA/raw/refs/heads/main/gene_sets/h_k_r_go_pid_reg_wik.csv'
pathways <- read.csv(url, 
                     header = F) %>%
  column_to_rownames(var = "V1")

genesets <- list()

for (i in (1:nrow(pathways))){
  genesets[[i]] <- unlist(pathways[i, ]) %>%
    unique()
  
  # get rid of weird empty entry 
  genesets[[i]] <- genesets[[i]][str_detect(genesets[[i]], "[A-Za-z0-9]")]
  
  names(genesets)[i] <- rownames(pathways)[i]
}

files <- list.files("tab_data/deseq2/",
                    pattern = "lfc_shrunk")

for (i in seq_along(files)){
  
  comp <- str_remove_all(files[i],
                         "lfc_shrunk_") %>%
    str_remove_all(".csv")
  
  message(comp)
  message(paste0(i, " out of ", length(files)))
  
  de <- read.csv(paste0("tab_data/deseq2/", files[i]))
  
  de <- de %>%
    arrange(desc(log2FoldChange))
  
  x <- de$log2FoldChange
  names(x) <- de$X
  
  fgseaRes <- fgsea(pathways=genesets, 
                    stats=x)
  
  # For saving as CSV
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  
  write.csv(fgseaRes,
            file = paste0("tab_data/deseq2_gsea/", comp, ".csv"),
            row.names = F)
  
}


