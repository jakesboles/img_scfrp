# Assembles one Seurat object from all CellBender-corrected count matrices
# for the coculture arm, and saves the counts matrix on-disk with BPCells.
# Mirrors monocultures/r_scripts/01_obj_creation.R's pattern, with two
# simplifications specific to this arm: no samples.csv join (sample-level
# metadata is parsed directly out of each sample_id string below, since
# that's the only place it exists for this arm -- see CLAUDE.md's
# "Coculture arm" section) and no _redo-style capture-vs-unified-sample
# distinction (every (pool_dir, sample) pair in cell_quantities.csv is
# both a distinct capture AND a distinct biological sample here -- unlike
# monoculture's iMG1/iMG1_redo, no two captures in this arm ever share a
# sample identity).

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/img_scfrp/cocultures")

data_out_dir <- "data/01_obj_assembly/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

# Manifest of every real CellBender output: pool_dir,sample,cells,droplets.
# Written by cellbender_scripts/make_cellbender_params.R -- reusing it
# here means this script never needs its own copy of which Cell Ranger/
# CellBender output directory belongs to which sample.
manifest <- read.csv("cellbender_scripts/cell_quantities.csv")

# Sample-level metadata lives entirely in the sample_id string itself for
# this arm: <neuron genotype>_<microglia genotype>[_<batch>] -- batch
# absent means batch 1, and "microglia genotype" == "CN" means a
# neuron-alone well (no microglia added). str_split_fixed(..., 2) is safe
# here because none of the genotype tokens (KOLF/GALC/GBA/GRN/LRRK2) or
# "CN" contain an underscore themselves.
manifest <- manifest %>%
  mutate(batch = str_extract(sample, "_(\\d+)$") %>% str_remove("^_"),
         batch = if_else(is.na(batch), "1", batch),
         base = str_remove(sample, "_\\d+$"))

tokens <- str_split_fixed(manifest$base, "_", 2)
manifest$neuron_genotype <- tokens[, 1]
manifest$microglia_genotype <- tokens[, 2]

# Only keep genes in the FRP probe set -- CellBender's output includes
# other genes despite the library being probe-based. Same reference used
# by the monoculture arm (shared across both arms, see CLAUDE.md).
probes <- read.csv("/projects/p31535/boles/cellranger_references/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                   skip = 5)
genes <- probes$probe_id %>%
  str_split_i(pattern = "[|]", i = 2) %>%
  unique()

# Create Seurat objects from corrected counts -----------------------------
message2("Creating Seurat objects")

# Each sample's matrix is written to its own on-disk BPCells directory
# right after reading, so only one sample's matrix is fully in memory at
# a time instead of holding all of them at once through the merge below.
bpcells_persample_dir <- paste0(data_out_dir, "bpcells_persample/")
dir.create(bpcells_persample_dir, showWarnings = F, recursive = T)

create_object <- function(pool_dir, sample, neuron_genotype, microglia_genotype, batch){
  h5_path <- file.path("cellbender", pool_dir, sample, paste0(sample, "_filtered.h5"))
  counts <- Read_CellBender_h5_Mat(h5_path)

  # dgCMatrix always stores its @x slot as R type "double" regardless of
  # the values it holds, which leaves write_matrix_dir()'s compressed
  # writer with an ambiguously-typed matrix. Converting to an explicit
  # integer type first avoids that.
  counts <- convert_matrix_type(counts, type = "uint32_t")

  bp_dir <- paste0(bpcells_persample_dir, sample)
  write_matrix_dir(mat = counts, dir = bp_dir)
  mat <- open_matrix_dir(dir = bp_dir)

  obj <- CreateSeuratObject(counts = mat, project = sample)
  obj$sample <- sample
  obj$pool_dir <- pool_dir
  obj$neuron_genotype <- neuron_genotype
  obj$microglia_genotype <- microglia_genotype
  obj$batch <- batch
  return(obj)
}

obj_list <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  message(paste0("Creating object for ", row$sample))
  obj_list[[row$sample]] <- create_object(row$pool_dir, row$sample,
                                          row$neuron_genotype, row$microglia_genotype,
                                          row$batch)
}

for (sample in names(obj_list)) {
  message(paste0("Removing non-probe list genes from ", sample))
  idx <- rownames(obj_list[[sample]]) %in% genes
  print(table(idx))
  obj_list[[sample]] <- obj_list[[sample]][idx, ]
}

message2("Merging Seurat objects")

# add.cell.ids uses `sample` (= names(obj_list)) directly -- every
# (pool_dir, sample) pair in this arm is both a distinct capture and a
# distinct biological sample, so there's no monoculture-style redo
# capture to additionally disambiguate here.
obj <- Merge_Seurat_List(obj_list, add.cell.ids = names(obj_list))
obj

message2("Joining layers")

obj <- JoinLayers(obj)

message2("Saving counts matrix as BPCells on-disk matrix")

write_matrix_dir(mat = obj[["RNA"]]$counts,
                 dir = paste0(data_out_dir, "bpcells"))

message2("Saving metadata as RDS")

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts <- open_matrix_dir(paste0(data_out_dir, "bpcells"))
#   meta <- readRDS(paste0(data_out_dir, "metadata.rds"))
#   obj <- CreateSeuratObject(counts = counts, meta.data = meta)
