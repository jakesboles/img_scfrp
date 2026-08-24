# CellBender

## Environment

Thomas's originally-pushed environment couldn't be reproduced from the files
as pushed (`Cellbender2.yml` had no `dependencies:` list — `conda env
export` was run against the parent project folder instead of the actual
environment path one level down, which silently produces an empty yml
rather than erroring). Since then, `real_export.yml` was captured directly
from Thomas's still-live environment at
`/projects/b1169/thomas/CellbenderEnv2/env/Cellbender2` — that's the
authoritative record of what's actually installed:
**`cellbender==0.3.0`, `torch==1.13.1`, `nvidia-cuda-runtime-cu11==11.7.99`**
(pip-installed; see the `pip:` block in `real_export.yml` for the rest).

`make_cellbender_env.sh` (fixed from Thomas's original `CellbenderEnv2.sh`)
builds a fresh working environment under this project's own directory
(`/projects/b1169/boles/img_scfrp/envs/cellbender`) rather than a
colleague's, in case that's ever needed again — not guaranteed
byte-identical to Thomas's, but `real_export.yml` is there if exact
reproduction matters.

## Running CellBender

`make_cellbender_params.R` reads each pool's per-sample
`cellranger/<pool_dir>/outs/per_sample_outs/<sample>/metrics_summary.csv`
and writes `cell_quantities.csv` (`pool_dir,sample,cells,droplets`,
`droplets = floor(cells * 1.75)`). `iMG2` is the sequencing-depth target
(see `seqtk/downsample_targets.txt`) and was never re-run downsampled, so
it reads `cellranger/iMG2/...` while the other four pools read
`cellranger/<pool>_downsampled/...` — that mapping lives once at the top of
the script and gets written into the `pool_dir` column so
`run_all_cellbender.sh` never needs its own copy of it.

`run_all_cellbender.sh`: run directly (not via `sbatch`) the first time; it
computes the array size from `cell_quantities.csv` and resubmits itself as
`sbatch --array=1-N%8`. The `%8` throttle caps concurrent A100 jobs — adjust
to match what `b1042`/`genomics-gpu` actually gives you. Each task activates
`envs/cellbender`, `cd`s into its own `cellbender/<pool_dir>/<sample>/`
directory (CellBender writes its checkpoint relative to `cwd` by default, so
two jobs sharing one can lock each other out — this must stay), and runs
`cellbender remove-background` there. Outputs land in
`cellbender/<pool_dir>/<sample>/<sample>.h5` plus CellBender's other output
files.

## Job failures seen 2026-08-24 — root causes and fixes

Three representative failures from one array submission (job 4141845),
all now fixed:

**1. `cell_quantities.csv` had a duplicate row for every single sample**
(150 data rows for 75 real `pool,sample` pairs — confirmed by checking:
every sample in `iMG1_downsampled` had a second row with the identical
value `43625`, regardless of which sample). Cell Ranger's per-sample
`metrics_summary.csv` contains *two* rows with `Metric.Name == "Cells"`:
one under `Category == "Cells"` (the real per-sample count) and one that's
actually the whole GEM well's aggregate total, reported under a different
`Category` but confusingly sharing the same `Metric.Name`. The old filter
(`Metric.Name == "Cells"` alone) matched both, and R's `data.frame()`
silently recycled the resulting two-element vector into two rows instead of
erroring. This is what caused **task 2's `ZeroDivisionError`**: it got
`--expected-cells 43625 --total-droplets-included 76343` for a sample that
actually only had ~3-4k real cells, so CellBender's "empty droplet" region
ended up with zero rows to average over.

  Fixed in `make_cellbender_params.R`: the filter is now
  `Metric.Name == "Cells" & Category == "Cells"`, plus a hard `stop()` with
  a diagnostic dump of the matching rows if it ever finds anything other
  than exactly one match, so a similar mismatch fails loudly next time
  instead of silently writing bad rows again. **The corrupted
  `cell_quantities.csv` has been deleted from the repo** — regenerate it
  with `Rscript cellbender_scripts/make_cellbender_params.R` before
  resubmitting; `run_all_cellbender.sh` will refuse to run without it.

  `run_all_cellbender.sh` also now sanity-checks each row's `cells`/
  `droplets` (positive integers, `droplets > cells`) before launching
  CellBender, and skips that task with a clear error instead of burning
  GPU-minutes on obviously-malformed input. This is a generic format check,
  not a fix for the specific corruption above (a mislabeled-but-well-formed
  pool total would still pass it) — the real fix for that is the `Category`
  filter.

**2. Task 1 (`GALC_asyn_1`) actually succeeded** — full training, checkpoint,
posterior, and all output files completed. The traceback at the very end
(`UnicodeDecodeError: 'ascii' codec can't decode byte 0xee...`, `Unable to
create report`) is CellBender's own HTML-report step failing because it
reads its just-written report file without specifying an encoding, and this
cluster's default locale is plain ASCII rather than UTF-8. Cosmetic only —
every actual output file (`.h5`, `_filtered.h5`, `_metrics.csv`, cell
barcodes, checkpoint) was written before that point. Fixed by exporting
`LC_ALL=en_US.UTF-8`/`LANG=en_US.UTF-8` in `run_all_cellbender.sh`.

**3. Task 3 (`GALC_myelin_1`, sane input: 2426 cells) crashed with
`RuntimeError: CUDA error: unrecognized error code`** partway through
training. This is unrelated to the two issues above — its input was a
single, correctly-sized row. `real_export.yml` confirms `torch==1.13.1`
with a bundled CUDA 11.7 runtime, which does support A100 (Ampere/sm_80)
architecturally, so this isn't a "too old for this GPU" problem. "CUDA
error: unrecognized error code" is a known, deliberately vague error class
in PyTorch/CUDA that most often shows up from GPU resource contention
(e.g. two processes landing on the same physical device), a transient
driver/hardware fault, or occasionally memory corruption carried over from
an earlier failing kernel launch in the same process — not something
diagnosable for certain from a single log. **Not yet fixed here** —
recommended next steps:
- Re-run the same (pool, sample) once the CSV bug above is fixed and see if
  it's reproducible. If it only happens occasionally, that points to
  contention/a transient fault rather than a systematic incompatibility.
- If it recurs consistently, run
  `python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.get_device_capability())"`
  inside the activated `envs/cellbender` environment on a `genomics-gpu`
  node to confirm what CUDA build/compute capability PyTorch actually sees
  at runtime on that specific node.
- If it correlates with how many array tasks are running at once, try
  lowering the `%8` throttle in `run_all_cellbender.sh` to rule out GPU
  contention on shared nodes.
