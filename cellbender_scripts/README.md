# CellBender

## Environment: what was wrong, what to do

Thomas's copied-in environment couldn't actually be reproduced from the
files as pushed:

- **`CellbenderEnv2/Cellbender2.yml` has no `dependencies:` list.** The
  original `CellbenderEnv2.sh` created the env at
  `.../CellbenderEnv2/env/Cellbender2` but exported from
  `.../CellbenderEnv2` (the parent folder, one level up) — `conda env
  export` against a path that isn't itself an environment silently produces
  a boilerplate yml with no packages in it. `conda env create -f
  Cellbender2.yml` from this file would build an empty environment, not a
  working CellBender install. See the warning comment at the top of that
  file for the full explanation.
- Even with that fixed, **`pip install cellbender` has no version pin**, so
  re-running the create script today installs whatever's currently on
  PyPI — not necessarily what Thomas's env ended up with. There's no
  captured record (no working yml, no `pip freeze`) of the exact versions
  that were actually used.

**Net effect: there's no way to exactly reproduce Thomas's original
environment from what's in this repo.** Two options, in order of
preference:

1. If you still have access to
   `/projects/b1169/thomas/CellbenderEnv2/env/Cellbender2`, run
   `conda env export --prefix /projects/b1169/thomas/CellbenderEnv2/env/Cellbender2 > real_export.yml`
   against it directly, today, before it's ever cleaned up — that's the
   only way to capture the exact packages he used.
2. Otherwise, use the fixed `CellbenderEnv2/CellbenderEnv2.sh` to build a
   fresh environment for the same purpose. It's not guaranteed
   byte-identical to Thomas's, but it fixes the concrete bugs (broken
   export path, an absolute output path instead of a bare filename) and
   builds under this project's own directory
   (`/projects/b1169/boles/img_scfrp/envs/cellbender`) instead of a
   colleague's, so future access doesn't depend on his directory staying
   around. It keeps `python=3.7.12` and building on a GPU node unchanged
   from Thomas's version (pip's cellbender install can resolve a
   CPU-only PyTorch if no GPU is visible at install time, so removing the
   GPU node wasn't safe to do without knowing that for certain) — it just
   shrinks that allocation from 4 cores/50G/2h to 2 cores/16G/30min, since
   creating the env is a few minutes of single-threaded work, not 2 hours
   of heavy compute.

Either way, after running the fixed script, **check that the resulting
`Cellbender2.yml` actually has a `dependencies:` list** before trusting it
— that's the one thing that silently failed before.

## Running CellBender

`CellbenderBatch.sh` (looped individual `sbatch` calls, one per sample) and
`CellbenderInd.sh` (the per-sample worker) are replaced by a single
`run_all_cellbender.sh`, following the manifest + SLURM array pattern
already used for `seqtk/` and consistent with this repo's other pipeline
steps:

1. `Rscript cellbender_scripts/CellQuantities.R` — reads each pool's
   `cellranger/<pool>/outs/per_sample_outs/<sample>/metrics_summary.csv`
   and writes `cell_quantities.csv` (`pool_dir,sample,cells,droplets`,
   `droplets = floor(cells * 1.75)`). Rewritten from Thomas's version to
   fix two bugs: it wrote `write.csv(..., quote = F)` without
   `row.names = FALSE`, so the CSV's first column was a meaningless row
   index that `CellbenderBatch.sh` then used as its sample lookup key
   (worked, but only by accident — it depended on rlist::list.rbind's row
   naming happening to stay unique); and it read
   `dat$Metric.Value[1]` assuming "Cells" is always the first metric row,
   rather than filtering by `Metric.Name` like
   `seqtk/get_seqtk_proportions.R` already correctly does for the same
   file. Also stripped ~40 unrelated library() calls (Seurat, DESeq2,
   MAST, ...) left over from a different script — this one only needs base
   R.

   **Which pool directories it reads matters**: `iMG2` is the
   sequencing-depth target (see `seqtk/downsample_targets.txt`) and was
   never re-run downsampled, so `CellQuantities.R` reads
   `cellranger/iMG2/...` for it while reading
   `cellranger/<pool>_downsampled/...` for the other four pools. This
   mapping is defined once, at the top of `CellQuantities.R`, and written
   into `cell_quantities.csv`'s `pool_dir` column so the array script never
   needs its own copy of it.

2. `bash cellbender_scripts/run_all_cellbender.sh` — run directly (not via
   `sbatch`) the first time; it computes the array size from
   `cell_quantities.csv` and resubmits itself as
   `sbatch --array=1-N%8 run_all_cellbender.sh`. The `%8` throttle caps how
   many `genomics-gpu`/A100 jobs run concurrently — the GPU partition is a
   much scarcer shared resource than the CPU `genomics` partition used
   elsewhere in this repo, so this is more conservative by default; adjust
   to match what you can actually get from `b1042`/`genomics-gpu`.

   Each array task activates the env from `envs/cellbender` (see above),
   `cd`s into its own `cellbender/<pool_dir>/<sample>/` directory, and runs
   `cellbender remove-background` there. **The `cd` before running matters
   and must stay**: CellBender writes its checkpoint file
   (`ckpt.tar.gz`) relative to the current working directory by default,
   so two jobs sharing a `cwd` can lock each other out of their
   checkpoints — this was in Thomas's original script's comments and is
   preserved here, minus the extra `mkdir` of a same-named directory one
   level down that his version also did (CellBender's `--output` argument
   is a file-path *prefix*, not a directory that needs to pre-exist under
   that exact name, so that third `mkdir` was dead code).

3. Outputs land in `cellbender/<pool_dir>/<sample>/<sample>.h5` (plus
   CellBender's other output files: `_filtered.h5`, `.log`,
   `_metrics.csv`, `_report.html`, checkpoint).

`--expected-cells`/`--total-droplets-included`/`--fpr 0.01`/`--epochs 100`
are unchanged from Thomas's script — no evidence those choices were wrong,
so they weren't second-guessed here.
