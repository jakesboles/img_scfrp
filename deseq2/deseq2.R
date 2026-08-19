library(tidyverse)
library(DESeq2)
library(Seurat)
library(SingleCellExperiment)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

tab_dir <- paste0(proj_dir, "tab_data/deseq2/")
dir.create(tab_dir, F, T)

plots_dir <- paste0(proj_dir, "plots/deseq2/")
dir.create(plots_dir, F, T)

obj <- readRDS(paste0(proj_dir, "data/05_integration/harmony_obj.rds"))

obj <- JoinLayers(obj)

bulk <- AggregateExpression(obj,
                            assays = "RNA",
                            group.by = c("orig.ident"),
                            return.seurat = F,
                            slot = "counts")

exp <- bulk$RNA

meta <- s@meta.data %>%
  dplyr::select(c(orig.ident, genotype, treatment, batch)) %>%
  distinct()

rownames(meta) <- meta$id

rownames(meta)
colnames(exp)

idx <- match(colnames(exp), rownames(meta))
meta <- meta[idx, ]

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = meta,
                              design = ~ genotype*treatment + batch)

dds <- DESeq(dds)

# res <- results(dds,
#                contrast = c("Group", "sALS", "Control"))

res_tbl <- as.data.frame(res)

res %>%
  as.data.frame() %>%
  arrange(padj) %>%
  head(30)

res %>%
  as.data.frame() %>%
  mutate(deg = case_when(padj > 0.05 ~ "Not DEG",
                         padj < 0.05 & log2FoldChange > log2(1.5) ~ "Upregulated",
                         padj < 0.05 & log2FoldChange < -log2(1.5) ~ "Downregulated")) %>%
  group_by(deg) %>%
  summarize(num_genes = n())

# VlnPlot(s,
#         features = c("RASSF5", "NR4A3", "TRA2B", "PRDM1"),
#         group.by = "Group")
# 
# plotCounts(dds,
#            gene = "NR4A3",
#            intgroup = "Batch")