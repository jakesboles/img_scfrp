library(tidyverse)
library(DESeq2)
library(apeglm)
library(IHW)

setwd("/projects/b1169/boles/img_scfrp")

dir.create("tab_data/deseq2/")

obj <- readRDS("data/05_integration/harmony_obj.rds")

bulk <- AggregateExpression(obj,
                            assays = "RNA",
                            group.by = c("orig.ident"),
                            return.seurat = F,
                            slot = "counts")

exp <- bulk$RNA

meta <- obj@meta.data %>%
  dplyr::select(c(orig.ident, batch, genotype, treatment, sample)) %>%
  distinct()

rownames(meta) <- str_replace_all(meta$orig.ident, "_", "-")

rownames(meta)
colnames(exp)

idx <- match(colnames(exp), rownames(meta))
meta <- meta[idx, ]

meta$genotype <- factor(meta$genotype,
                        levels = c("KOLF", "GALC", "GBA", "GRN", "LRRK2"))
meta$treatment <- factor(meta$treatment,
                         levels = c("control", "myelin", "asyn"))

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = meta,
                              design = ~ batch + sample)

keep <- rowSums(counts(dds) >= 10) >= 30

dds <- dds[keep, ]

samples <- levels(meta$sample)

contrast_tbl <- combn(samples, 2)

refs <- unique(contrast_tbl[1, ])

for (i in seq_along(refs)){
  
  message(paste0("REFERENCE: ", refs[i]))
  
  col_idx <- contrast_tbl[1,] %in% refs[i]
  
  contrast_tbl_small <- contrast_tbl[, col_idx] %>%
    as.data.frame()
  
  dds$sample <- relevel(dds$sample, ref = refs[i])
  
  dds <- DESeq(dds,
               fitType = "local")
  
  for (j in 1:ncol(contrast_tbl_small)){
    
    query_group <- contrast_tbl_small[2, j]
    ref_group <- contrast_tbl_small[1, j]
    
    message(paste0(query_group, " vs ", ref_group))
    message(paste0(j, " out of ", ncol(contrast_tbl_small)))
    
    res <- results(dds, 
                   contrast = c("sample", query_group, ref_group),
                   filterFun = ihw,
                   independentFiltering = T)
    
    res_tbl <- as.data.frame(res)
    write.csv(res_tbl, 
              file = paste0("tab_data/deseq2/", query_group, "_vs_", ref_group, ".csv"))
    
    shrunk <- lfcShrink(dds,
                        res = res,
                        coef = paste0("sample_", query_group, "_vs_", ref_group),
                        type = "apeglm")
    
    shrunk_tbl <- as.data.frame(shrunk)
    write.csv(shrunk_tbl,
              file = paste0("tab_data/deseq2/lfc_shrunk_", query_group, "_vs_", ref_group, ".csv"))
  }
  
}
