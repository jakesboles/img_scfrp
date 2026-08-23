# Reads each pool's per-sample Cell Ranger metrics and writes
# cellbender_scripts/cell_quantities.csv (pool,sample,cells,droplets), the
# target list run_all_cellbender.sh's SLURM array reads from -- one row per
# (pool, sample) CellBender needs to run on.
#
# Base R only: the original version of this script loaded ~40 packages
# (Seurat, DESeq2, MAST, ComplexHeatmap, ...) unrelated to reading a few
# CSVs and computing floor(cells * 1.75), which only slowed startup.

setwd("/projects/b1169/boles/img_scfrp")

# Which cellranger/<dir>/outs/per_sample_outs/ to read per pool. iMG2 is
# the sequencing-depth target (see seqtk/downsample_targets.txt) and was
# never re-run downsampled, so it alone keeps its original output dir name.
pool_dirs <- c(
  iMG1      = "iMG1_downsampled",
  iMG1_redo = "iMG1_redo_downsampled",
  iMG2      = "iMG2",
  iMG3      = "iMG3_downsampled",
  iMG4      = "iMG4_downsampled"
)

rows <- list()

for (pool in names(pool_dirs)) {

  pool_dir <- pool_dirs[[pool]]
  samples_dir <- file.path("cellranger", pool_dir, "outs", "per_sample_outs")
  samples <- list.files(samples_dir)
  samples <- samples[!endsWith(samples, "NotUsed")]

  for (sample in samples) {

    metrics_path <- file.path(samples_dir, sample, "metrics_summary.csv")
    metrics <- read.csv(metrics_path)

    # Filter by Metric.Name rather than assuming a fixed row order (this is
    # what seqtk/get_seqtk_proportions.R already does for the same file).
    cells <- as.numeric(metrics$Metric.Value[metrics$Metric.Name == "Cells"])

    # pool_dir (not just `pool`) is written out so run_all_cellbender.sh's
    # bash array task can build the Cell Ranger --input path and the
    # CellBender output path directly from the CSV, without needing its
    # own copy of the pool -> output-directory mapping above.
    rows[[paste(pool, sample)]] <- data.frame(
      pool_dir = pool_dir,
      sample = sample,
      cells = cells,
      droplets = floor(cells * 1.75)
    )
  }
}

result <- do.call(rbind, rows)

write.csv(result, "cellbender_scripts/cell_quantities.csv",
          row.names = FALSE, quote = FALSE)
