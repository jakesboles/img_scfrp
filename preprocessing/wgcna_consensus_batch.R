library(hdWGCNA)
library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(UCell)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

data_out_dir <- paste0(proj_dir, "data/wgcna/")
dir.create(data_out_dir, F, T)

plots_dir <- paste0(proj_dir, "plots/wgcna/")
dir.create(plots_dir, F, T)

tab_out_dir <- paste0(proj_dir, "tab_data/wgcna/")
dir.create(tab_out_dir, F, T)

obj <- readRDS(paste0(proj_dir, "data/05_integration/harmony_obj.rds"))

obj <- JoinLayers(obj)

obj$dummy_group <- 1

# Set up WGCNA ------------------------------------------------------------

obj <- SetupForWGCNA(obj, 
                     gene_select = "fraction", 
                     fraction = 0.01,  
                     wgcna_name = "wgcna_consensus")

# Construct meta cells ----------------------------------------------------

obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("orig.ident", "dummy_group", "batch", "genotype", "treatment", "sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony_pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'orig.ident' # set the Idents of the metacell seurat object
)

# Process meta cells as a Seurat dataset ----------------------------------

obj <- NormalizeMetacells(obj)
obj <- ScaleMetacells(obj,
                      features = VariableFeatures(obj))
obj <- RunPCAMetacells(obj,
                       features = VariableFeatures(obj))
obj <- RunHarmonyMetacells(obj,
                           group.by.vars = "orig.ident")
obj <- RunUMAPMetacells(obj,
                        reduction = "harmony",
                        dims = 1:10)

p1 <- DimPlotMetacells(obj,
                       # reduction = "umap",
                       group.by = "sample",
                       cols = DiscretePalette_scCustomize(num_colors = 15,
                                                          palette = "varibow"),
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p2 <- DimPlotMetacells(obj,
                       # reduction = "umap",
                       group.by = "genotype",
                       cols = c("red", "yellow", "green", "blue", "purple"),
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p3 <- DimPlotMetacells(obj,
                       # reduction = "umap",
                       group.by = "treatment",
                       cols = JCO_Four()[1:3],
                       raster = F) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.justification = "left")

p4 <- DimPlotMetacells(obj,
                       # reduction = "umap",
                       group.by = "batch",
                       cols = Dark2_Pal()[1:4],
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

obj <- SetMultiExpr(
  obj,
  group_name = "1", # the name of the group of interest in the group.by column
  group.by='dummy_group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
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

# assemble with patchwork
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
  tom_outdir = paste0(data_out_dir),
  consensus = T
)

png(paste0(plots_dir, "dendrogram_consensus.png"), 
    height = 8, width = 8,
    res = 600,
    units = "in")
PlotDendrogram(obj, main= "Dendrogram (consensus)")
dev.off()

# TOM <- GetTOM(obj)
# saveRDS(TOM,
#         file = paste0(data_out_dir, "tom.rds"))

# Compute eigengenes, module members --------------------------------------

obj <- ModuleEigengenes(
  obj,
  group.by.vars="orig.ident"
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

# stitch together with patchwork
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

# stitch together with patchwork
png(paste0(plots_dir, "module_scores_umap_consensus.png"), 
    height = 8, width = 8,
    units = "in",
    res = 600)
wrap_plots(plot_list, ncol=4)
dev.off()

# get module membership for each gene and write to csv
mods <- obj@misc[["wgcna_consensus"]][["wgcna_modules"]]
write.csv(mods,
          file = paste0(tab_out_dir, "module_members_consensus.csv"))


# save object -------------------------------------------------------------

saveRDS(obj,
        file = paste0(data_out_dir, "obj_consensus.rds"))

