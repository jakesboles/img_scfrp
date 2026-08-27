# Pseudobulk DESeq2 across all pairwise genotype x treatment `condition`
# comparisons (105 = 15 choose 2), with `batch` as a covariate. Loads
# 06_qc2.R's BPCells counts + metadata (not a monolithic harmony_obj.rds) --
# only raw counts are needed here (for AggregateExpression(layer = "counts")
# pseudobulk), not the normalized data layer or any reduction.
#
# Contrasts are between `condition` levels, not individual `sample`s: each
# condition has exactly 4 real batch replicates (one sample per batch), so
# design = ~ batch + condition leaves real residual degrees of freedom for
# DESeq2 to estimate dispersion from. The previous draft's design (~ batch +
# sample) was both rank-deficient (`sample` already fully determines
# `batch`, so the two terms are perfectly collinear) and, even without that,
# unfittable in principle -- with exactly one pseudobulk sample per `sample`
# level, a `sample`-level design has zero residual df regardless.
#
# Per-reference-relevel-and-refit loop (rather than one fit + arbitrary
# `results(contrast = ...)` calls) is needed because lfcShrink(type =
# "apeglm") only shrinks a named model coefficient (coef =), not an
# arbitrary contrast -- getting a properly-named, apeglm-shrinkable
# coefficient for every pairwise comparison requires releveling the
# reference and refitting for each one.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(DESeq2)
  library(apeglm)
  library(IHW)
})

setwd("/projects/b1169/boles/img_scfrp")

tab_dir <- "tab_data/deseq2/"
dir.create(tab_dir, showWarnings = F, recursive = T)

# Read in 06's integrated, doublet-cluster-filtered counts + metadata -------

counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
meta <- readRDS("data/06_qc2/metadata.rds")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")

bulk <- AggregateExpression(obj,
                            assays = "RNA",
                            group.by = "sample",
                            return.seurat = F,
                            layer = "counts")

exp <- bulk$RNA

# genotype/treatment/batch/condition are already correctly-leveled factors
# as of 02_qc.R (Title Case treatment labels, explicit genotype/batch/
# condition level order) -- carried through unchanged to 06_qc2.R's
# metadata, so no re-factoring needed here. (The previous draft re-factored
# treatment against lowercase levels -- "control"/"myelin"/"asyn" -- that
# only ever existed in 01_obj_creation.R's raw metadata; against the actual
# Title Case values here, every cell would have silently become NA.)
meta <- obj@meta.data %>%
  dplyr::select(c(sample, batch, genotype, treatment, condition)) %>%
  distinct()

# AggregateExpression() sanitizes "_" to "-" in its output column names
# (Seurat's own convention, to avoid colliding with the "_" it uses to
# join multiple group.by variables) -- matched here so the two can be
# lined up.
rownames(meta) <- str_replace_all(meta$sample, "_", "-")

idx <- match(colnames(exp), rownames(meta))
meta <- meta[idx, ]

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = meta,
                              design = ~ batch + condition)

keep <- rowSums(counts(dds) >= 10) >= 30

dds <- dds[keep, ]

conditions <- levels(meta$condition)

contrast_tbl <- combn(conditions, 2)

refs <- unique(contrast_tbl[1, ])

for (i in seq_along(refs)){

  message(paste0("REFERENCE: ", refs[i]))

  col_idx <- contrast_tbl[1,] %in% refs[i]

  contrast_tbl_small <- contrast_tbl[, col_idx] %>%
    as.data.frame()

  dds$condition <- relevel(dds$condition, ref = refs[i])

  dds <- DESeq(dds,
               fitType = "local")

  for (j in 1:ncol(contrast_tbl_small)){

    query_group <- contrast_tbl_small[2, j]
    ref_group <- contrast_tbl_small[1, j]

    message(paste0(query_group, " vs ", ref_group))
    message(paste0(j, " out of ", ncol(contrast_tbl_small)))

    res <- results(dds,
                   contrast = c("condition", query_group, ref_group),
                   filterFun = ihw,
                   independentFiltering = T)

    res_tbl <- as.data.frame(res)
    write.csv(res_tbl,
              file = paste0(tab_dir, query_group, "_vs_", ref_group, ".csv"))

    shrunk <- lfcShrink(dds,
                        res = res,
                        coef = paste0("condition_", query_group, "_vs_", ref_group),
                        type = "apeglm")

    shrunk_tbl <- as.data.frame(shrunk)
    write.csv(shrunk_tbl,
              file = paste0(tab_dir, "lfc_shrunk_", query_group, "_vs_", ref_group, ".csv"))
  }

}
