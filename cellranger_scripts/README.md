# Cell Ranger cluster-mode setup

The 4 pool jobs (`iMG1_cellranger.sh` ... `iMG4_cellranger.sh`) run
`cellranger multi` in [cluster
mode](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode)
via `--jobmode=slurm.template`. Instead of doing all the work inside one
SLURM allocation, the SLURM job you submit only hosts the lightweight Martian
(`mrp`) controller process; `mrp` in turn submits one `sbatch` job per
pipeline stage chunk to the `b1169`/`b1169` account/partition using
`slurm.template`, so a single pool's run is spread across many concurrent
cluster jobs instead of being capped by one node's core count.

`b1169`/`b1169` only runs ~3 sub-jobs concurrently — almost certainly a
per-account/QOS concurrent-job cap. We tried `b1042`/`genomics`
(Northwestern's shared Genomics Compute Cluster queue) expecting more
concurrency, but it made things worse: still only a few jobs running at a
time, each sitting in the queue much longer before starting. That's
consistent with `genomics` being a busy queue shared across many other
groups' jobs — SLURM's fairshare/backfill scheduling decides when *your*
jobs run based on cluster-wide demand, and no `cellranger`/Martian flag can
override that from the client side. So this is back on `b1169`/`b1169`,
where jobs at least start promptly even though only a few run at once.

If you want to push on this further, two read-only, no-admin-needed checks
are worth running before/instead of another account switch:
- `sacctmgr show qos format=name,maxjobspu,maxsubmitpu -p` — shows whether
  `b1169`'s QOS has an explicit `MaxJobsPU` (max running jobs per user) that
  would explain the ~3 concurrent cap directly.
- `squeue -p genomics | wc -l` (or `sinfo -p genomics`) next time you try
  `genomics`, to see how contested that partition is at the time — if it's
  full of other users' jobs, that's the wait, not a config problem.
If `MaxJobsPU` on `b1169` turns out to be the limit, an HPC admin would need
to raise it — but that's a deliberate ask you make if/when you decide it's
worth it, not something required to use cluster mode at all.

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
- Each pool's parent SLURM job only requests 4 cores / 16G / 96h — this
  covers `mrp` plus any Martian stages marked "local" (which always run
  inside the parent process regardless of jobmode). The real per-stage
  compute happens in the jobs `slurm.template` submits, sized per-stage via
  `__MRO_THREADS__`/`__MRO_MEM_GB__`, capped at 24h walltime each in the
  template.
