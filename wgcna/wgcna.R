# Consensus hdWGCNA network on the whole-cohort integrated object: builds
# metacells, constructs a co-expression network with a separate network per
# `batch` reconciled into one consensus network, and computes module
# eigengenes/scores. Loads 06_qc2.R's BPCells/RDS pieces (counts, normalized
# data, harmony reduction, harmony UMAP, metadata, variable features)
# instead of a single monolithic obj_consensus.rds, and only re-saves what
# THIS stage adds on top of 06's output (the hdWGCNA network/module data,
# and metadata with module scores attached) -- counts/data/harmony/
# harmony_umap are unchanged from 06_qc2.R's cell set, so
# wgcna_stats.R/wgcna_viz.R reach back into data/06_qc2/ directly for those
# rather than having them re-saved here.
#
# hdWGCNA's own misc payload (network, module table, module eigengenes) is
# metacell/module-scale (thousands of metacells x tens of modules), not
# cell x gene scale -- small enough to save as one RDS, same as this
# project already does for other non-matrix pieces (pca.rds, harmony.rds).
# It's the raw/normalized expression matrices specifically that the BPCells
# convention exists to keep off of a single RDS, and this script doesn't
# re-save those.

suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(UCell)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp")

data_out_dir <- "wgcna/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

plots_dir <- "wgcna/"
dir.create(plots_dir, showWarnings = F, recursive = T)

tab_out_dir <- "wgcna/"
dir.create(tab_out_dir, showWarnings = F, recursive = T)

# Read in 06's integrated, doublet-cluster-filtered object -------------------

message2("Reading in integrated object from 06_qc2.R")

counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
meta <- readRDS("data/06_qc2/metadata.rds")
harmony <- readRDS("data/06_qc2/harmony.rds")
harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")

# 06_qc2.R saves its own fresh variable_features.rds (04's copy reflects
# the pre-doublet-cluster-filter cohort).
var_features <- readRDS("data/06_qc2/variable_features.rds")

# hdWGCNA isn't written with BPCells/lazy-matrix awareness in mind (unlike
# Seurat's own NormalizeData()/ScaleData()/RunPCA()/IntegrateLayers(),
# which do support BPCells directly, as 04_norm_pca.R/05_integration_
# harmony.R/06_qc2.R already demonstrate) -- SetupForWGCNA()'s gene_select =
# "fraction" step walks the assay with base-R indexing that a lazy
# IterableMatrix doesn't support, failing with "invalid subscript type
# 'S4'". Same "materialize before handing to a non-BPCells-aware tool"
# caution already established for DoubletFinder (03_doubletfinder.R) and
# UCell (wgcna_stats.R), applied here to the whole assay since hdWGCNA's
# pipeline (gene selection, metacell construction, network/TOM) touches it
# throughout, not just at this one step. Full-cohort scale here is a
# standard in-memory sparse matrix for a single-cell analysis (this is
# hdWGCNA's normal expected input) -- BPCells backing was the atypical
# part, not the memory cost of materializing.
counts_mat <- as(counts_mat, "dgCMatrix")
data_mat <- as(data_mat, "dgCMatrix")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap
VariableFeatures(obj) <- var_features

obj <- ScaleData(obj)

obj$dummy_group <- 1

# Set up WGCNA ------------------------------------------------------------

obj <- SetupForWGCNA(obj,
                     gene_select = "fraction",
                     fraction = 0.01,
                     wgcna_name = "wgcna_consensus")

# Construct meta cells ----------------------------------------------------

# Grouping by `sample` (the unified 60-level biological sample) rather than
# `orig.ident` -- both hold the identical value throughout this pipeline
# (01_obj_creation.R sets orig.ident from the same sample_unified value used
# for `sample`), but `sample` is this project's documented name for it.
# genotype/treatment/batch/condition are redundant with `sample` for
# partitioning purposes (each sample maps to exactly one combination of all
# four), but are included so they land as real columns on the metacell
# object's own metadata, not just implied by its `sample` value.
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("dummy_group", "batch", "genotype", "treatment", "condition", "sample"),
  reduction = 'harmony', # matches the reduction name 05/06_qc2 actually save
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'sample' # set the Idents of the metacell seurat object
)

# Process meta cells as a Seurat dataset ----------------------------------

obj <- NormalizeMetacells(obj)
obj <- ScaleMetacells(obj,
                      features = VariableFeatures(obj))
obj <- RunPCAMetacells(obj,
                       features = VariableFeatures(obj))
obj <- RunHarmonyMetacells(obj,
                           group.by.vars = "sample")
obj <- RunUMAPMetacells(obj,
                        reduction = "harmony",
                        dims = 1:10)

# Same palettes as 02_qc.R/04_norm_pca.R/05_integration_harmony.R's plots.
# First panel groups by `condition` (15 levels), not the 60-level `sample`
# -- the old version of this script grouped by `sample` while only
# supplying a 15-color palette, the same mismatch already fixed in
# 05_integration_harmony.R/06_qc2.R.
genotype_pal <- c("red", "yellow", "green", "blue", "purple")
treatment_pal <- JCO_Four()[1:3]
batch_pal <- Dark2_Pal()[1:4]
condition_pal <- DiscretePalette_scCustomize(num_colors = nlevels(obj$condition),
                                             palette = "varibow")

p1 <- DimPlotMetacells(obj,
                       group.by = "condition",
                       cols = condition_pal,
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p2 <- DimPlotMetacells(obj,
                       group.by = "genotype",
                       cols = genotype_pal,
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p3 <- DimPlotMetacells(obj,
                       group.by = "treatment",
                       cols = treatment_pal,
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p4 <- DimPlotMetacells(obj,
                       group.by = "batch",
                       cols = batch_pal,
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p <- p1 + p2 + p3 + p4 +
  plot_layout(ncol = 2)

ggsave(p,
       filename = paste0(plots_dir, "metacells_umaps.png"),
       units = "in", dpi = 600,
       height = 8,
       width = 10)

# Setup expression matrix -------------------------------------------------

# multi.group.by = "batch" (not pool_dir): consensus networks are built
# per-batch and reconciled, matching the same granularity decision already
# made for Harmony integration in 05 (split by `sample`, not `pool_dir` --
# see that script's header comment). Keeping `batch` here is consistent
# with that choice rather than introducing a finer iMG1/iMG1_redo
# distinction at this stage that upstream integration didn't make either.
obj <- SetMultiExpr(
  obj,
  group_name = "1", # the name of the group of interest in the group.by column
  group.by = 'dummy_group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data', # using normalized data
  multi.group.by = "batch"
)

# Test soft powers --------------------------------------------------------

obj <- TestSoftPowersConsensus(
  obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(obj)

consensus_groups <- unique(obj$batch)
p_list <- lapply(1:length(consensus_groups), function(i){
  cur_group <- consensus_groups[[i]]
  plot_list[[i]][[1]] + ggtitle(paste0('Batch: ', cur_group)) + theme(plot.title=element_text(hjust=0.5))
})

png(paste0(plots_dir, "soft_powers_consensus.png"),
    height = 8, width = 8,
    res = 600,
    units = "in")
wrap_plots(p_list, ncol=2)
dev.off()

power_table <- GetPowerTable(obj)
write.csv(power_table,
          file = paste0(tab_out_dir, "soft_powers_consensus.csv"))

# Construct TOM -----------------------------------------------------------

obj <- ConstructNetwork(
  obj,
  tom_name = 'tom_consensus', # name of the topoligical overlap matrix written to disk
  overwrite_tom = T,
  tom_outdir = data_out_dir,
  consensus = T
)

png(paste0(plots_dir, "dendrogram_consensus.png"),
    height = 8, width = 8,
    res = 600,
    units = "in")
PlotDendrogram(obj, main= "Dendrogram (consensus)")
dev.off()

# Compute eigengenes, module members --------------------------------------

obj <- ModuleEigengenes(
  obj,
  group.by.vars = "sample"
)

hMEs <- GetMEs(obj)
write.csv(hMEs,
          file = paste0(tab_out_dir, "hmes_consensus.csv"))

obj <- ModuleConnectivity(
  obj,
  group.by = 'dummy_group',
  group_name = '1'
)

p <- PlotKMEs(obj, ncol=4, text_size = 4)
png(paste0(plots_dir, "kmes_consensus.png"),
    height = 12, width = 12,
    res = 600,
    units = "in")
print(p)
dev.off()

obj <- ModuleExprScore(
  obj,
  n_genes = 25,
  method='UCell'
)

plot_list <- ModuleFeaturePlot(
  obj,
  features='MEs', # plot the hMEs
  reduction = "harmony_umap",
  order = TRUE,
  ucell = T
  # order so the points with highest hMEs are on top
)

png(paste0(plots_dir, "eigengenes_umap_consensus.png"),
    height = 8, width = 8,
    res = 600,
    units = "in")
print(wrap_plots(plot_list, ncol=4))
dev.off()

plot_list <- ModuleFeaturePlot(
  obj,
  features='scores', # plot the hMEs
  reduction = "harmony_umap",
  order = "shuffle",
  ucell = T
  # order so the points with highest scores are on top
)

png(paste0(plots_dir, "module_scores_umap_consensus.png"),
    height = 8, width = 8,
    units = "in",
    res = 600)
wrap_plots(plot_list, ncol=4)
dev.off()

# get module membership for each gene and write to csv --------------------
# Needed by wgcna_clusterprofiler.R (pathway enrichment per module) and
# wgcna_stats.R (module gene sets for its own UCell scoring).
mods <- obj@misc[["wgcna_consensus"]][["wgcna_modules"]]
write.csv(mods,
          file = paste0(tab_out_dir, "module_members_consensus.csv"))

# Save ----------------------------------------------------------------------

message2("Saving WGCNA network/module data")

# hdWGCNA's own experiment payload (network, modules, eigengenes, metacell
# object) -- metacell/module-scale, not cell x gene-matrix scale, so one
# RDS is appropriate here (see header comment).
saveRDS(obj@misc[["wgcna_consensus"]],
        file = paste0(data_out_dir, "misc_wgcna_consensus.rds"))

message2("Saving metadata (with module scores) as RDS")

# This stage's own metadata, distinct from 06_qc2/metadata.rds --
# ModuleExprScore() above added per-cell module score/eigengene columns on
# top of what 06_qc2 saved.
saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

# Downstream scripts should reconstruct the object from these on-disk
# pieces -- counts/data/harmony/harmony_umap are unchanged from
# 06_qc2.R's cell set, so reach back there directly rather than have this
# stage re-save them:
#   counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
#   data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
#   harmony <- readRDS("data/06_qc2/harmony.rds")
#   harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")
#   meta <- readRDS("wgcna/metadata.rds")
#   obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
#   obj[["RNA"]]$data <- data_mat
#   obj[["harmony"]] <- harmony
#   obj[["harmony_umap"]] <- harmony_umap
#   obj@misc[["wgcna_consensus"]] <- readRDS("wgcna/misc_wgcna_consensus.rds")
#   obj@misc$active_wgcna <- "wgcna_consensus"
