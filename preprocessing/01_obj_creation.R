# Assembles one Seurat object from all CellBender-corrected count matrices
# for the img_scfrp iMG project, and saves the counts matrix on-disk with
# BPCells (rather than a single in-memory RDS) so downstream scripts can
# load it without reading the full matrix into memory. Mirrors the pattern
# used in als_cns_scrnaseq/r_scripts/01_obj_creation.R.

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

setwd("/projects/b1169/boles/img_scfrp")

data_out_dir <- "data/01_obj_assembly/"
dir.create(data_out_dir, showWarnings = F, recursive = T)

# Manifest of every real CellBender output: pool_dir,sample,cells,droplets.
# Written by cellbender_scripts/make_cellbender_params.R -- reusing it here
# means this script never needs its own copy of which cellranger/CellBender
# output directory belongs to which sample.
manifest <- read.csv("cellbender_scripts/cell_quantities.csv")

# iMG1/iMG1_redo (and no other pool) are two independent GEM wells of the
# same biological samples -- kept as separate Cell Ranger/CellBender runs
# to avoid treating same-barcode reads from unrelated cells in each capture
# as one cell (see cellranger_scripts/README.md). `sample` stays the
# distinct per-capture id (e.g. GALC_asyn_1_redo) so it can still be used
# for on-disk paths and add.cell.ids -- cells from the two captures must
# keep unique prefixes, or exactly the barcode-collision risk that
# motivated splitting them in the first place reappears here at merge
# time. `sample_unified` strips the _redo suffix and is what actually gets
# stored in the object's metadata as `sample`, so both captures' cells
# read as the same sample downstream; which capture a cell actually came
# from is preserved separately via `pool_dir`.
manifest$sample_unified <- str_remove(manifest$sample, "_redo$")

meta <- read.csv("samples.csv") %>%
  select(id, genotype, treatment, batch) %>%
  dplyr::rename(sample = id)

manifest <- manifest %>%
  left_join(meta, by = "sample")

if (any(is.na(manifest$genotype))) {
  stop("No samples.csv metadata found for: ",
       paste(manifest$sample[is.na(manifest$genotype)], collapse = ", "))
}

# Only keep genes in the FRP probe set -- CellBender's output includes
# other genes despite the library being probe-based.
probes <- read.csv("/projects/p31535/boles/cellranger_references/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                   skip = 5)
genes <- probes$probe_id %>%
  str_split_i(pattern = "[|]", i = 2) %>%
  unique()

# Create Seurat objects from corrected counts -----------------------------
message2("Creating Seurat objects")

# Each sample's matrix is written to its own on-disk BPCells directory
# right after reading, so only one sample's matrix is fully in memory at a
# time instead of holding all of them at once through the merge below.
bpcells_persample_dir <- paste0(data_out_dir, "bpcells_persample/")
dir.create(bpcells_persample_dir, showWarnings = F, recursive = T)

create_object <- function(pool_dir, sample, sample_unified, genotype, treatment, date, batch){
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

  obj <- CreateSeuratObject(counts = mat, project = sample_unified)
  obj$sample <- sample_unified
  obj$pool_dir <- pool_dir
  obj$genotype <- genotype
  obj$treatment <- treatment
  obj$date <- date
  obj$batch <- batch
  return(obj)
}

obj_list <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  message(paste0("Creating object for ", row$sample))
  obj_list[[row$sample]] <- create_object(row$pool_dir, row$sample, row$sample_unified,
                                          row$genotype, row$treatment, row$date, row$batch)
}

for (sample in names(obj_list)) {
  message(paste0("Removing non-probe list genes from ", sample))
  idx <- rownames(obj_list[[sample]]) %in% genes
  print(table(idx))
  obj_list[[sample]] <- obj_list[[sample]][idx, ]
}

message2("Merging Seurat objects")

# add.cell.ids uses the distinct per-capture sample id (manifest$sample,
# via names(obj_list)), not sample_unified -- see the comment above on why
# iMG1/iMG1_redo cells must keep distinct barcode prefixes even though
# they'll share the same `sample` metadata value.
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
