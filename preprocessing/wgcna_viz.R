library(hdWGCNA)
library(Seurat)
library(scCustomize)
library(patchwork)

proj_dir <- "/projects/b1169/boles/img_scfrp/"

plots_dir <- paste0(proj_dir, "plots/wgcna/")
dir.create(plots_dir, F, T)

obj <- readRDS(paste0(proj_dir, "data/wgcna/obj.rds"))

# obj@misc

# pdf(file = paste0(plots_dir, "module_radars.pdf"),
    # width = 6, height = 4)
p <- ModuleRadarPlot(obj,
                group.by = "genotype",
                base.size = 1,
                ncol = 6)
# dev.off()

p

ModuleCorrelogram(obj)

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
