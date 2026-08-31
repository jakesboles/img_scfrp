# Milo (miloR) differential neighborhood-abundance testing: each genotype
# vs KOLF (the parental line), collapsed across treatment. Adapted from a
# working script from another project -- that script's contrast variable
# ("Group"), reduced dims ("INTEGRATED_PCA"/"INTEGRATED_UMAP2"), and PC
# count (d = 13) are specific to that project's own object and don't apply
# here; this version is rebuilt against this project's own metadata schema
# (genotype/treatment/batch/sample/condition, see CLAUDE.md's metadata
# table) and established PC count (dims = 1:10, matching every other stage
# from 04_norm_pca.R on).
#
# Loads 06_qc2.R's BPCells counts + metadata + harmony/harmony_umap
# reductions (not a monolithic object) -- only counts are needed to
# satisfy CreateSeuratObject()/as.SingleCellExperiment(), since Milo's
# core machinery (neighbor graph, reduced-dim coordinates, cell metadata)
# never touches expression values directly; the normalized `data` layer
# isn't loaded at all. harmony_umap is only used for plotNhoodGraphDA()'s
# 2D layout at the very end -- neighborhoods/distances are computed on the
# real (10-dim) harmony embedding, not the UMAP.
#
# 06_qc2.R computes its own FindNeighbors() internally but doesn't save the
# resulting Graph object, so it's recomputed here on the same harmony
# reduction/dims/k this project uses everywhere else, immediately before
# handing it to miloR::buildFromAdjacency() -- this reuses the same
# neighbor structure the rest of the pipeline is built on, rather than
# letting miloR build an independent kNN graph via its own buildGraph().

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(miloR)
  library(scater)
  library(igraph)
  library(scales)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp")

results_dir <- "results/milo/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_out_dir <- "data/milo/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

# Read in 06's integrated, doublet-cluster-filtered object -------------------

message2("Reading in integrated object from 06_qc2.R")

counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
meta <- readRDS("data/06_qc2/metadata.rds")
harmony <- readRDS("data/06_qc2/harmony.rds")
harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")

# Milo/SingleCellExperiment/scater aren't BPCells-aware -- same
# "materialize before handing to a non-BPCells-aware tool" caution already
# established for DoubletFinder/hdWGCNA/UCell in this project.
counts_mat <- as(counts_mat, "dgCMatrix")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap

message2("Building neighbor graph on the harmony embedding")

obj <- FindNeighbors(obj,
                     reduction = "harmony",
                     dims = 1:10,
                     k.param = 15,
                     compute.SNN = T)

# Convert to SingleCellExperiment / Milo --------------------------------------

message2("Converting to SingleCellExperiment/Milo")

sce <- as.SingleCellExperiment(obj)

# as.SingleCellExperiment() upper-cases reduction names (harmony ->
# HARMONY) -- confirmed here rather than assumed, since this is the first
# time this pipeline has gone through an SCE conversion.
message2("Reduced dims available in the SCE (used as reduced_dims/layout below):")
print(reducedDimNames(sce))

milo <- Milo(sce)

miloR::graph(milo) <- miloR::graph(buildFromAdjacency(obj@graphs$RNA_snn,
                                                       k = 15,
                                                       is.binary = F))

# Build neighborhoods ---------------------------------------------------------

message2("Building neighborhoods")

milo <- makeNhoods(milo,
                   prop = 0.05,
                   k = 15,
                   d = 10, # this project's established PC count -- see 04_norm_pca.R/05_integration_harmony.R/06_qc2.R
                   refined = TRUE,
                   reduced_dims = "HARMONY")

p <- plotNhoodSizeHist(milo)
ggsave(p, filename = paste0(results_dir, "nhood_size_hist.png"),
       units = "in", dpi = 600, height = 5, width = 6)

milo <- countCells(milo,
                   meta.data = as.data.frame(colData(milo)),
                   sample = "sample")

milo <- calcNhoodDistance(milo,
                          d = 10,
                          reduced.dim = "HARMONY")

# Design: genotype vs KOLF, collapsed across treatment ------------------------

design <- data.frame(colData(milo))[, c("sample", "genotype", "batch")]
design <- distinct(design)
rownames(design) <- design$sample

# genotype/batch are already correctly-leveled factors as of 02_qc.R,
# carried through unchanged to 06_qc2.R's metadata -- checked explicitly
# here (rather than re-derived) since the SCE round-trip via colData() is a
# new transformation this pipeline hasn't used before, and re-deriving
# levels with a separately-maintained hardcoded list is exactly the kind of
# drift that has caused bugs elsewhere in this project (see deseq2.R's
# lowercase-treatment bug in CLAUDE.md).
expected_genotype_levels <- c("KOLF", "GALC", "GBA", "GRN", "LRRK2")
if (!is.factor(design$genotype) || !setequal(levels(design$genotype), expected_genotype_levels)) {
  stop("genotype didn't survive the SCE round-trip as the expected factor -- check colData(milo)$genotype.")
}
genotype_levels <- levels(design$genotype)

# batch's factor labels ("Batch 1", ...) contain a space, which would
# otherwise land in the model matrix's column names (e.g. "batchBatch 2")
# -- sanitized to avoid any downstream contrast-string-matching issues,
# same reasoning as deseq2.R's "_"-to-"-" sample-name sanitization.
design$batch <- factor(str_remove_all(as.character(design$batch), " "))

# Test each non-KOLF genotype against KOLF ------------------------------------
#
# design = ~ 0 + genotype + batch: the `0 +` gives genotype a full set of
# per-level columns (needed to write pairwise contrasts against KOLF,
# rather than only k-1 columns relative to an implicit reference); `batch`
# is included as a covariate for the same reason deseq2.R's design
# includes it (~ batch + condition) -- each genotype has 4 real batch
# replicates (one per treatment x batch combination, pooled across
# treatment here), so this controls for batch-driven neighborhood-
# abundance shifts rather than attributing them to genotype.

test_genotypes <- setdiff(genotype_levels, "KOLF")

da_results <- list()

for (g in test_genotypes) {

  message2(paste0(g, " vs KOLF"))

  contrast <- paste0("genotype", g, " - genotypeKOLF")

  res <- testNhoods(milo,
                    design = ~ 0 + genotype + batch,
                    design.df = design,
                    model.contrasts = contrast,
                    fdr.weighting = "graph-overlap",
                    norm.method = "TMM")

  message(paste0(g, " vs KOLF -- neighborhoods at SpatialFDR < 0.05:"))
  print(table(res$SpatialFDR < 0.05))

  write.csv(res,
            file = paste0(results_dir, "da_", g, "_vs_KOLF.csv"),
            row.names = F)

  da_results[[g]] <- res
}

# Neighborhood graph + DA plots ------------------------------------------------

message2("Building neighborhood graph for plotting")

milo <- buildNhoodGraph(milo)

for (g in test_genotypes) {
  p <- plotNhoodGraphDA(milo, da_results[[g]], layout = "HARMONY_UMAP", alpha = 0.1) +
    ggtitle(paste0(g, " vs KOLF"))
  ggsave(p,
         filename = paste0(results_dir, "da_umap_", g, "_vs_KOLF.png"),
         units = "in", dpi = 600, height = 6, width = 7)
}

# Save Milo-specific pieces ----------------------------------------------------

message2("Saving Milo neighborhood/graph data")

# Neighborhood/module-scale (nhoods x samples, nhoods x nhoods), not cell x
# gene-matrix scale -- one RDS is appropriate here, same reasoning as
# wgcna.R's misc_wgcna_consensus.rds. Doesn't include the counts matrix,
# unchanged from data/06_qc2/.
saveRDS(list(nhoods = nhoods(milo),
            nhoodCounts = nhoodCounts(milo),
            nhoodIndex = nhoodIndex(milo),
            nhoodDistances = nhoodDistances(milo),
            graph = miloR::graph(milo),
            nhoodGraph = nhoodGraph(milo)),
       file = paste0(data_out_dir, "milo_nhoods.rds"))

saveRDS(design, file = paste0(data_out_dir, "design.rds"))

# Downstream re-plotting should reconstruct from these on-disk pieces
# rather than rerunning makeNhoods()/countCells()/calcNhoodDistance():
#   milo_pieces <- readRDS("data/milo/milo_nhoods.rds")
#   design <- readRDS("data/milo/design.rds")
#   nhoods(milo) <- milo_pieces$nhoods
#   nhoodCounts(milo) <- milo_pieces$nhoodCounts
#   nhoodIndex(milo) <- milo_pieces$nhoodIndex
#   nhoodDistances(milo) <- milo_pieces$nhoodDistances
#   miloR::graph(milo) <- milo_pieces$graph
#   nhoodGraph(milo) <- milo_pieces$nhoodGraph
# (starting from a freshly-rebuilt `milo <- Milo(sce)` on the same object)
