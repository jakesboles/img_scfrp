# Downsampling FASTQs to normalize sequencing depth

`get_seqtk_proportions.R` reads each pool's `cellranger/<pool>/outs/qc_library_metrics.csv`,
computes `target_fraction = min(mean reads per cell across pools) / (pool's
mean reads per cell)`, and writes `downsample_targets.txt` as
`pool,fraction` rows. The pool with the lowest depth (iMG2) gets
`target_fraction == 1` and is filtered out — nothing to downsample there,
so it's correctly absent from the file.

## Do the FASTQs need to be concatenated first?

No. `seqtk sample -s<seed> <file> <fraction>` decides independently, per
read, whether to keep it, based only on the read and the seed — it doesn't
need to know the total read count across a pool up front, so applying the
*same* fraction to every FASTQ file that makes up a pool is statistically
equivalent to concatenating those files first and then subsampling once.
This matters for `iMG3` specifically, which spans two source directories
(`Gate38_11.22.2024` and `FixedRNA_Jake`) — every file from both
directories gets `iMG3`'s one fraction independently, with no concatenation
step, and the result is the same as if they'd been merged first.

(This does rely on `mean reads per cell` in the R script's calculation
already being Cell Ranger's own aggregate metric across every FASTQ that
fed into that pool's `multi` run — which it is, since `multi` reports one
number per pool regardless of how many `[libraries]` rows contributed to
it. So the fraction itself is already computed correctly at the pool
level; per-file application just needs to reuse that one number for every
file in the pool.)

The one thing per-file sampling does *not* preserve automatically is R1/R2
pairing *across* different array tasks — see the seed note in `seqtk.sh`:
R1 and R2 of the same lane must use the same seed (they do — `-s100` is
fixed for every task), which is sufficient because seqtk's per-read keep/drop
decision depends only on the read and the seed, not on anything else in the
file.

## Running it

1. **Build the manifest** (login node, no compute needed):
   `bash seqtk/build_manifest.sh`
   This reads `downsample_targets.txt` and, for each pool, pulls its
   `fastq_id`/`fastqs` rows straight out of
   `cellranger_configs/<pool>_config.csv` — the same file Cell Ranger
   itself uses — so the file list can't drift out of sync with what
   Cell Ranger actually processed. It globs the real files in each
   `fastqs` directory and writes `downsample_manifest.txt` as
   `source_path,fraction,dest_path` rows, one per FASTQ file (R1 and R2
   are separate rows). It prints the row count and the exact `sbatch`
   command to run next.
2. **Submit the array job** using the count `build_manifest.sh` reported:
   `sbatch --array=1-N seqtk/seqtk.sh` (append `%20` or similar, e.g.
   `--array=1-N%20`, to cap how many run concurrently if you don't want
   every file downsampling at once)
3. Downsampled files land under `downsampled_fastq/<original subdir
   name>/<original filename>` — the same subdirectory/filename structure
   as `fastq/`, just under a new top-level directory, so a downsampled
   Cell Ranger config only needs `fastq/` swapped for `downsampled_fastq/`
   in its `[libraries]` paths.
4. Re-run `build_manifest.sh` (and re-submit) any time
   `downsample_targets.txt` or the source FASTQs change — it always
   regenerates the manifest from scratch.
