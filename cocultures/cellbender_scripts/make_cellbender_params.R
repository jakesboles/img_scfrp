# Reads each coculture pool's per-sample Cell Ranger metrics and writes
# cellbender_scripts/cell_quantities.csv (pool_dir,sample,cells,droplets),
# the target list run_all_cellbender_cpu.sh's SLURM array reads from --
# one row per (pool, sample) CellBender needs to run on. Coculture-arm
# counterpart of monocultures/cellbender_scripts/make_cellbender_params.R;
# see that file and CLAUDE.md for the full rationale behind each step.
#
# Base R only: no object-building/stats libraries needed to read a few
# CSVs and compute floor(cells * 1.75).

setwd("/projects/b1169/boles/img_scfrp/cocultures")

# Which cellranger/<dir>/outs/per_sample_outs/ to read per pool. Unlike
# monoculture, this arm has no downsampling step (not needed for this
# project, per the user), so each pool's Cell Ranger output directory
# name is identical to its own pool name -- no _downsampled-style
# renaming to track here.
pool_dirs <- c(
  JSB127     = "JSB127",
  `JSB134-1` = "JSB134-1",
  `JSB134-2` = "JSB134-2",
  `JSB134-3` = "JSB134-3"
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

    # Same duplicate-row trap as monoculture's metrics_summary.csv --
    # Category == "Cells" picks out the sample-specific row, not the
    # whole GEM well's aggregate total under the same Metric.Name.
    is_match <- metrics$Metric.Name == "Cells" & metrics$Category == "Cells"
    if (sum(is_match) != 1) {
      stop(sprintf(
        "Expected exactly 1 matching 'Cells' row for %s/%s, found %d. Matching metrics_summary.csv rows:\n%s",
        pool_dir, sample, sum(is_match),
        paste(capture.output(print(metrics[metrics$Metric.Name == "Cells", ])), collapse = "\n")
      ))
    }
    cells <- as.numeric(metrics$Metric.Value[is_match])

    # pool_dir (not just `pool`) is written out so
    # run_all_cellbender_cpu.sh's bash array task can build the Cell
    # Ranger --input path and the CellBender output path directly from
    # the CSV, without needing its own copy of the pool -> output-
    # directory mapping above.
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
