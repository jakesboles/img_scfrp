# Second-round QC: removes whole clusters heavily enriched for
# DoubletFinder-flagged doublets, then re-derives normalization/PCA/Harmony
# integration on the filtered cell set (repeating 04_norm_pca.R/
# 05_integration_harmony.R's steps). This is deliberately a whole-cluster
# removal, not a per-cell DF.adj == "Doublet" filter -- 03_doubletfinder.R's
# own documented caveat is that per-sample DoubletFinder can't catch
# doublets formed between two DIFFERENT samples sharing a GEM well, so
# clusters that don't resolve into any single expected population (visibly
# doublet-enriched in the plot below) are a complementary way to catch what
# that per-sample granularity missed.
#
# Loads 05_integration_harmony.R's BPCells/RDS pieces (not a monolithic
# object), same as every stage since 04. Saves in the same BPCells pattern
# 05 uses (bpcells_counts, bpcells_data, metadata, harmony, harmony_umap),
# plus its own fresh variable_features.rds -- 04's variable_features.rds
# was selected against the pre-filter cell set, so anything reaching back
# to 04 for this would otherwise be using a stale gene selection basis
# once this stage removes cells. wgcna.R and deseq2.R read from
# data/06_qc2/ (repointed once this stage was adopted), not
# data/05_integration/.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(patchwork)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp")

results_dir <- "results/06_qc2/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

data_out_dir <- "data/06_qc2/"
dir.create(data_out_dir,
           showWarnings = F,
           recursive = T)

# Read in 05's integrated object ---------------------------------------------

message2("Reading in integrated object from 05_integration_harmony.R")

counts_mat <- open_matrix_dir("data/05_integration/bpcells_counts")
data_mat <- open_matrix_dir("data/05_integration/bpcells_data")
meta <- readRDS("data/05_integration/metadata.rds")
harmony <- readRDS("data/05_integration/harmony.rds")
harmony_umap <- readRDS("data/05_integration/harmony_umap.rds")
var_features <- readRDS("data/04_norm_pca/variable_features.rds")

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["harmony_umap"]] <- harmony_umap
VariableFeatures(obj) <- var_features

# Same palettes as 02_qc.R/04_norm_pca.R/05_integration_harmony.R's plots.
genotype_pal <- c("red", "yellow", "green", "blue", "purple")
treatment_pal <- JCO_Four()[1:3]
batch_pal <- Dark2_Pal()[1:4]
condition_pal <- DiscretePalette_scCustomize(num_colors = nlevels(obj$condition),
                                             palette = "varibow")

make_umaps <- function(reduction, type){

  p1 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "condition",
                         colors_use = condition_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p2 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "genotype",
                         colors_use = genotype_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p3 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "treatment",
                         colors_use = treatment_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p4 <- DimPlot_scCustom(obj,
                         reduction = reduction,
                         group.by = "batch",
                         colors_use = batch_pal,
                         raster = F) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")

  p <- p1 + p2 + p3 + p4 +
    plot_layout(ncol = 2)

  ggsave(p,
         filename = paste0(results_dir, type, "_umaps.png"),
         units = "in", dpi = 600,
         height = 8,
         width = 10)
}

# Cluster on the existing (05) integration to find doublet-enriched
# clusters -----------------------------------------------------------------

message2("Clustering on 05's harmony reduction to find doublet-enriched clusters")

obj <- obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  FindClusters(method = "igraph",
               algorithm = 4,
               resolution = 2,
               cluster.name = "doublet_qc_cluster",
               graph.name = "RNA_snn")

p1 <- DimPlot_scCustom(obj,
                       reduction = "harmony_umap")

p2 <- DimPlot_scCustom(obj,
                       reduction = "harmony_umap",
                       group.by = "DF.adj")

p3 <- dittoBarPlot(obj,
                   var = "DF.adj",
                   group.by = "doublet_qc_cluster")

design <- "
AABB
CCCC
"

p <- p1 + p2 + p3 +
  plot_layout(design = design)

ggsave(p,
       filename = paste0(results_dir, "clusters_by_doublets.png"),
       units = "in", dpi = 600,
       height = 10, width = 12)

# Remove doublet-enriched clusters -------------------------------------------

# Chosen by eye from clusters_by_doublets.png above (clusters visibly
# dominated by DF.adj == "Doublet" cells) -- a manual QC judgment call,
# not derivable programmatically, and specific to this exact
# resolution = 2 clustering. Re-inspect and update this list if upstream
# data changes and this script is rerun.
doublet_clusters <- c(31, 35, 36)

message2(paste0("Removing clusters: ", paste(doublet_clusters, collapse = ", ")))

removed_counts <- obj@meta.data %>%
  count(doublet_qc_cluster) %>%
  filter(doublet_qc_cluster %in% doublet_clusters)
write.csv(removed_counts,
          file = paste0(results_dir, "removed_doublet_clusters.csv"),
          row.names = F)

obj <- subset(obj, subset = !(doublet_qc_cluster %in% doublet_clusters))

# Normalize, scale, run PCA --------------------------------------------------

message2("Normalizing data")
obj <- NormalizeData(obj)

message2("Finding variable features")
obj <- FindVariableFeatures(obj)

message2("Scaling data")
obj <- ScaleData(obj)

message2("Running PCA")
obj <- RunPCA(obj, npcs = 50)

# Integrate with Harmony -------------------------------------------------

message2("Splitting layers by sample for integration")

# Same split-by-`sample` choice as 05_integration_harmony.R -- kept
# consistent with that stage's own explicit, user-confirmed decision.
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

message2("Integrating samples using Harmony")

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony",
                       dims = 1:10)

# IntegrateLayers() only adds the new reduction -- it doesn't touch/join the
# assay's existing layers, so counts/data are still split by sample here.
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

message2("Computing post-integration UMAP")

obj <- obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:10,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

make_umaps("harmony_umap", "harmony")

# Save ----------------------------------------------------------------------

message2("Saving counts matrix as BPCells on-disk matrix")

# Removed first if present, so a rerun doesn't fail on "path already exists"
# against a stale/incomplete directory from a previous attempt.
counts_dir <- paste0(data_out_dir, "bpcells_counts")
if (dir.exists(counts_dir)) {
  unlink(counts_dir, recursive = T)
}

counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")
write_matrix_dir(mat = counts_out, dir = counts_dir)

message2("Saving normalized data matrix as BPCells on-disk matrix")

data_dir <- paste0(data_out_dir, "bpcells_data")
if (dir.exists(data_dir)) {
  unlink(data_dir, recursive = T)
}

write_matrix_dir(mat = obj[["RNA"]]$data, dir = data_dir)

message2("Saving metadata as RDS")

saveRDS(obj@meta.data, file = paste0(data_out_dir, "metadata.rds"))

message2("Saving Harmony reduction as RDS")

saveRDS(obj[["harmony"]], file = paste0(data_out_dir, "harmony.rds"))

message2("Saving Harmony UMAP as RDS")

saveRDS(obj[["harmony_umap"]], file = paste0(data_out_dir, "harmony_umap.rds"))

message2("Saving variable features as RDS")

# Freshly selected against the filtered cell set -- 04's variable_features
# would otherwise be the only copy available, but it reflects the
# pre-filter cohort.
saveRDS(VariableFeatures(obj), file = paste0(data_out_dir, "variable_features.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts_mat <- open_matrix_dir("data/06_qc2/bpcells_counts")
#   data_mat <- open_matrix_dir("data/06_qc2/bpcells_data")
#   meta <- readRDS("data/06_qc2/metadata.rds")
#   harmony <- readRDS("data/06_qc2/harmony.rds")
#   harmony_umap <- readRDS("data/06_qc2/harmony_umap.rds")
#   var_features <- readRDS("data/06_qc2/variable_features.rds")
#   obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta, assay = "RNA")
#   obj[["RNA"]]$data <- data_mat
#   obj[["harmony"]] <- harmony
#   obj[["harmony_umap"]] <- harmony_umap
#   VariableFeatures(obj) <- var_features
