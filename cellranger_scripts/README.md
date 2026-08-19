# Cell Ranger cluster-mode setup

The 4 pool jobs (`iMG1_cellranger.sh` ... `iMG4_cellranger.sh`) run
`cellranger multi` in [cluster
mode](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-cluster-mode)
via `--jobmode=slurm.template`. Instead of doing all the work inside one
SLURM allocation, the SLURM job you submit only hosts the lightweight Martian
(`mrp`) controller process; `mrp` in turn submits one `sbatch` job per
pipeline stage chunk to the `b1169` account/partition using
`slurm.template`, so a single pool's run is spread across many concurrent
cluster jobs instead of being capped by one node's core count.

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

- `--maxjobs 24` / `--jobinterval 100` in each script throttle how many
  concurrent stage jobs Martian keeps in the SLURM queue and how fast it
  submits them. Raise `--maxjobs` if `b1169` can support more concurrent jobs
  and you want more parallelism; lower it if you hit per-user pending-job
  limits.
- Each pool's parent SLURM job only requests 4 cores / 16G / 96h — this
  covers `mrp` plus any Martian stages marked "local" (which always run
  inside the parent process regardless of jobmode). The real per-stage
  compute happens in the jobs `slurm.template` submits, sized per-stage via
  `__MRO_THREADS__`/`__MRO_MEM_GB__`, capped at 24h walltime each in the
  template.
