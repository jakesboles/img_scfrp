# Cell Ranger cluster-mode setup (coculture arm)

Same cluster-mode mechanics as `monocultures/cellranger_scripts/README.md`
(`cellranger multi --jobmode=slurm.template`, `b1042`/`genomics`,
`test_cluster_mode.sh` before any real run) -- not repeated here in full,
read that file first if you haven't already. This README only covers
what's specific to this arm.

## Why JSB134-1 and JSB134-2 need to be rerun

An earlier copy-paste error while preparing this arm's config files
pointed `JSB134-1_config.csv` and `JSB134-2_config.csv`'s `[libraries]`
section at the same underlying FASTQ selection, so the two pools' Cell
Ranger runs processed identical reads despite having different
`[samples]` barcode/sample-id assignments. The configs have since been
corrected (`JSB134-1` now reads `fastq_id = 134-1`, `JSB134-2` reads
`fastq_id = 134-2`, each pulling the matching files out of the shared
`/projects/b1042/Gate_Lab/boles/coculture_scfrp/fastq` directory by that
`fastq_id`), but the **existing `cellranger/JSB134-1/` and
`cellranger/JSB134-2/` output directories on disk still reflect the old,
buggy run** — they need to be regenerated from the corrected configs.
`JSB127` and `JSB134-3` are not implicated in this bug (their configs
were never duplicated) and don't need a rerun on this basis.

## Before rerunning: clear the stale output

Cell Ranger refuses to run `multi --id "JSB134-1"` into a
`cellranger/JSB134-1/` directory that already exists. Move or delete the
stale directories before submitting (`JSB134-1_cellranger.sh`/
`JSB134-2_cellranger.sh`, or `run_all.sh` for both at once):

```
mv cellranger/JSB134-1 cellranger/JSB134-1_stale_bug  # or: rm -rf cellranger/JSB134-1
mv cellranger/JSB134-2 cellranger/JSB134-2_stale_bug  # or: rm -rf cellranger/JSB134-2
```

**Also check CellBender**: `cellbender_scripts/cell_quantities.csv` and
any CellBender output already produced for `JSB134-1`/`JSB134-2` samples
were generated from the stale Cell Ranger output and are stale too --
regenerate `cell_quantities.csv` (`Rscript
cellbender_scripts/make_cellbender_params.R`) and rerun CellBender for
those samples once this rerun finishes. If CellBender is still actively
running against the old output when you read this, those two pools'
in-flight/completed work will need to be redone either way.

## Running it

1. Validate cluster mode first if you haven't already for this arm:
   `sbatch cellranger_scripts/test_cluster_mode.sh` (see monoculture's
   README for what to check), then `rm -rf cellranger/cluster_mode_test`.
2. Clear the stale `JSB134-1`/`JSB134-2` output as above.
3. `bash cellranger_scripts/run_all.sh` (submits only `JSB134-1` and
   `JSB134-2` — see that script's own comment). `JSB127_cellranger.sh`
   and `JSB134-3_cellranger.sh` exist for completeness/future reuse but
   aren't submitted by `run_all.sh`, since neither is known to need a
   rerun.
