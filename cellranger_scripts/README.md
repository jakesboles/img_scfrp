# Cell Ranger cluster-mode setup

The 5 pool jobs (`iMG1_cellranger.sh`, `iMG1_redo_cellranger.sh`,
`iMG2_cellranger.sh` ... `iMG4_cellranger.sh`) run `cellranger multi` in [cluster
mode](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode)
via `--jobmode=slurm.template`. Instead of doing all the work inside one
SLURM allocation, the SLURM job you submit only hosts the lightweight Martian
(`mrp`) controller process; `mrp` in turn submits one `sbatch` job per
pipeline stage chunk to the `b1042`/`genomics` account/partition using
`slurm.template`, so a single pool's run is spread across many concurrent
cluster jobs instead of being capped by one node's core count.

## Why iMG1 and iMG1_redo are separate runs

`iMG1` (fastqs under `fastq/Gate38_11.22.2024/`) and `iMG1_redo` (fastqs
under `fastq/Gate40/`) both contain the same 15 batch-1 samples and the same
`BC001`-`BC016` probe barcode assignment, but they came from two physically
distinct 10x Chromium captures — a second cell aliquot loaded onto a new
gel bead pool after the first capture. Because probe-barcoded FRP samples
are demultiplexed by (10x cell barcode, probe barcode) within a single GEM
well, and 10x cell barcodes are drawn from the same fixed whitelist across
*any* two independent gel bead pools, a cell barcode from the `iMG1` capture
and a cell barcode from the `iMG1_redo` capture can collide by chance even
though they came from different physical cells. Running both captures
through one `multi` call (as the original single `iMG1_config.csv` briefly
did, pointing `[libraries]` at both fastq dirs) would let Cell Ranger treat
same-barcode reads from each capture as one cell. Keeping them as two
`multi` runs (`iMG1_config.csv`/`iMG1_cellranger.sh` and
`iMG1_redo_config.csv`/`iMG1_redo_cellranger.sh`, each with a single
`[libraries]` row) avoids that entirely — sample_ids in the redo config are
suffixed `_redo` (e.g. `KOLF_untreated_1_redo`) so downstream analysis can't
mix them up with the originals either.

## Why `b1042`/`genomics`, not `b1169`/`b1169`

We tried both. `b1169`/`b1169` only ever ran ~3 sub-jobs concurrently, and
`sacctmgr` confirmed that account/QOS has a hard cap on running jobs per
user, independent of how busy the cluster is. `b1042`/`genomics`
(Northwestern's shared Genomics Compute Cluster queue) has no such per-user
cap — it's limited only by actual traffic/contention on the genomics nodes
at the time you submit, so real concurrency scales with `--maxjobs` instead
of being capped at ~3 regardless of demand. That's the tradeoff: `genomics`
can queue longer when the shared nodes are busy, but it has real headroom
for parallelism that `b1169` structurally does not.

Every parent job (`iMG1_cellranger.sh`, `iMG1_redo_cellranger.sh`,
`iMG2..4_cellranger.sh`, `test_cluster_mode.sh`) and
`slurm.template` (which governs the sub-jobs Martian dispatches) must point
at the same `b1042`/`genomics` account/partition — they were out of sync
earlier (parent jobs switched but `slurm.template` still said `b1169`),
which would have submitted every actual stage job under the wrong
account. If you ever change this again, grep for `--account`/`--partition`
across `cellranger_scripts/` to make sure all of them move together.

## Before running the real jobs: `test_cluster_mode.sh`

`slurm.template` is unsupported/community-maintained by 10x for SLURM, and a
misconfigured template (or a cluster policy that blocks a compute node from
calling `sbatch` itself) will only show up once `mrp` tries to submit its
first stage job. Finding that out after a multi-day, 350G real job has been
queued for hours is expensive, so validate it first with Cell Ranger's
built-in tiny test dataset:

1. Submit the test: `sbatch cellranger_scripts/test_cluster_mode.sh`
2. Watch sub-jobs appear: `squeue -u $USER` should show a `cluster_mode_test`
   job plus a stream of short-lived child jobs (name prefixed `ID.` by
   Martian) getting submitted and completing — `cellranger testrun` exercises
   ~200 small stage jobs, so you should see real churn within a minute or two
   of the parent job starting.
3. Check the outcome:
   - Success: `cellranger/cluster_mode_test/_log` and the job's SLURM log
     under `logs/` end with `Pipestance completed successfully!`.
   - Failure to submit sub-jobs at all (e.g. immediate error mentioning
     `sbatch` or permission denied): most likely your cluster does not allow
     a compute node to submit further SLURM jobs (nested submission). If so,
     don't run the pool jobs via `sbatch` — instead run `cellranger multi
     --jobmode=...` directly on the login node inside `tmux`/`screen` so the
     parent process lives outside the batch system.
   - Sub-jobs submitted but stuck pending/failing on account or partition
     errors: fix `slurm.template`'s `--account`/`--partition` lines (or
     `--mem`/walltime if your cluster rejects the request) and re-run.
4. Once it reports success, clean up before running the real jobs so the
   `--id` doesn't collide with anything:
   `rm -rf cellranger/cluster_mode_test`
5. Only then submit the real pools, e.g. via `run_all.sh`.

## Tuning

- `--maxjobs 24` / `--jobinterval 100` in each script only cap how many
  stage jobs Martian is *willing* to have pending/running in the queue at
  once and how fast it submits them — they can't force the scheduler to run
  more of your jobs concurrently than the partition's own QOS/fairshare
  allows. Raising `--maxjobs` only helps once the actual account/partition
  limit is higher.
- Each pool's parent SLURM job only requests 4 cores / 16G / 48h (`genomics`'s
  walltime limit is 72h; 48h leaves margin) — this covers `mrp` plus any
  Martian stages marked "local" (which always run inside the parent process
  regardless of jobmode). The real per-stage compute happens in the jobs
  `slurm.template` submits, sized per-stage via
  `__MRO_THREADS__`/`__MRO_MEM_GB__`, capped at 24h walltime each in the
  template.
