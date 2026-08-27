# Visualizes the consensus WGCNA network built by wgcna.R: module radar
# plot by genotype, module-module correlogram, and per-sample dot/violin
# plots of module eigengene expression. Loads the whole-cohort object from
# 06_qc2.R's BPCells/RDS pieces plus wgcna.R's own metadata.rds and
# misc_wgcna_consensus.rds (the hdWGCNA network/module payload) -- not a
# monolithic data/wgcna/obj.rds, which no longer exists now that wgcna.R
# saves separate pieces.

suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(scCustomize)
  library(patchwork)
  library(BPCells)
})

setwd("/projects/b1169/boles/img_scfrp")

# wgcna.R writes its own outputs directly to wgcna/ (not data/plots/
# tab_data/wgcna/).
plots_dir <- "wgcna/"
dir.create(plots_dir, showWarnings = F, recursive = T)

# Read in 06's integrated, doublet-cluster-filtered object + wgcna.R's
# network data ----------------------------------------------------------

counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
harmony <- readRDS("data/06_qc2/harmony.rds")
harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")

# wgcna.R's own metadata (module scores from ModuleExprScore() attached on
# top of 06_qc2/metadata.rds).
meta <- readRDS("wgcna/metadata.rds")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap

# hdWGCNA's Get*()/Module*() helpers below expect this populated exactly as
# SetupForWGCNA()/wgcna.R itself left it.
obj@misc[["wgcna_consensus"]] <- readRDS("wgcna/misc_wgcna_consensus.rds")
obj@misc$active_wgcna <- "wgcna_consensus"

# Module radar / correlogram --------------------------------------------

p <- ModuleRadarPlot(obj,
                group.by = "genotype",
                base.size = 1,
                ncol = 6)

p

ModuleCorrelogram(obj)

# Per-cell module eigengene dot/violin plots -------------------------------

MEs <- GetMEs(obj,
              harmonized = T)
modules <- GetModules(obj)

mods <- levels(modules$module)
mods <- mods[mods != 'grey']

obj@meta.data <- cbind(obj@meta.data, MEs)
obj

DotPlot(obj,
        features = mods,
        group.by = "sample") +
  scale_color_gradient2()

Clustered_DotPlot(obj,
                  features = mods,
                  group.by = "sample")

Stacked_VlnPlot(obj,
                features = mods[9:14],
                group.by = "sample")
