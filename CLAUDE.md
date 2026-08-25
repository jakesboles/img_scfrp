# img_scfrp — working notes for continuing this pipeline

This file is a memory dump of what's been established across a long session
building this pipeline, so future work doesn't have to re-derive it. Written
by Claude; update it as the project evolves rather than letting it go stale.

## What this project is

iPSC-derived microglia (iMG) scRNA-seq, 10x Chromium **Fixed RNA Profiling**
(FRP) with probe-barcode multiplexing. Design: 5 genotypes (KOLF parental +
GALC/GBA/GRN/LRRK2, PD-associated gene edits) x 3 treatments (untreated,
myelin, asyn) x 4 sequencing batches = up to 60 unique biological samples,
each identified like `GALC_asyn_1` (`<genotype>_<treatment>_<batch-number>`
— note the order; this is *not* `<batch>_<genotype>_<treatment>`, a mixup
that has bitten one draft script already).

Each batch is one GEM well (up to 16 samples multiplexed via probe
barcodes) — **except batch 1**, which was captured *twice* as two
physically independent GEM wells (`iMG1` and `iMG1_redo`: same 15 samples,
different cell aliquot, different gel bead pool — a genuine redo, not a
resequencing). This produces **5 real GEM wells total**
(`iMG1`, `iMG1_redo`, `iMG2`, `iMG3`, `iMG4`) and **75 (pool, sample)
"captures"** overall (60 unique biological samples, 15 of which — batch 1 —
exist in both `iMG1` and `iMG1_redo`).

**This distinction between "GEM well" and "biological sample" is the single
most important recurring fact in this codebase.** `iMG1`/`iMG1_redo` cells
must never be pooled for anything barcode-sensitive (Cell Ranger, CellBender,
DoubletFinder) — two independent captures can produce the same 10x cell
barcode by chance, and pooling them risks silently merging two unrelated
cells into one. At the same time, later analysis stages want them treated as
one biological sample. The convention landed on:
- `samples.csv`: `iMG1`'s 15 samples and their `_redo` duplicates are
  separate rows (e.g. `GALC_asyn_1` and `GALC_asyn_1_redo`), both
  `batch = 1`, differing only in `fastq_dirs`.
- Every barcode-sensitive stage (Cell Ranger, CellBender, DoubletFinder)
  processes `iMG1` and `iMG1_redo` as fully separate runs.
- From `01_obj_creation.R` onward, the merged object's `sample` metadata
  column is *unified* (the `_redo` suffix is stripped, so both captures'
  cells read as the same sample) — but a separate `pool_dir` metadata
  column preserves which physical capture each cell actually came from.
  **`sample` alone is therefore ambiguous between the two captures** for
  anything that needs the distinct identity (e.g. per-capture DoubletFinder
  file naming) — use `pool_dir` + `sample` together, or the raw manifest.

## HPC environment (Northwestern Quest)

- Repo/storage lives at `/projects/b1169/boles/img_scfrp` (this is also
  where `data/`, `plots/`, `tab_data/`, `logs/`, `fastq/`,
  `downsampled_fastq/`, `cellranger/`, `cellbender/`, `envs/` live —
  all gitignored, HPC-only).
- Cell Ranger reference/probe-set files:
  `/projects/p31535/boles/cellranger_references/`.
  Cell Ranger binary: `/projects/p31535/boles/cellranger-10.0.0/cellranger`.
- **SLURM accounts, and why**: `b1042`/`genomics` (CPU) is the reliable
  workhorse used for Cell Ranger, seqtk, CellBender-CPU, and R
  preprocessing — high concurrency has never caused problems there.
  `b1169`/`b1169` was tried first for Cell Ranger cluster mode but has a
  **hard per-user running-job cap** (confirmed via `sacctmgr`), so
  everything moved to `b1042`/`genomics`. `b1042`/`genomics-gpu` (A100) is
  **currently unreliable** — see the CellBender GPU section below; CPU
  fallback is in active use instead.
- R modules: `module load R/4.4.0`, `gcc/11.2.0`,
  `hdf5/1.14.1-2-gcc-12.3.0`.
- Python/CellBender env: `module load python-miniconda3/4.10.3` +
  `mamba`, then `source activate /projects/b1169/boles/img_scfrp/envs/cellbender`.

## Pipeline stages, in order

1. **`cellranger_configs/` + `cellranger_scripts/`** — one `cellranger multi`
   config per pool (5 real pools, plus `_downsampled` variants for 4 of
   them — not `iMG2`, see below). Runs in **cluster mode**
   (`--jobmode=slurm.template`, `cellranger_scripts/slurm.template`),
   dispatching per-stage sub-jobs to `b1042`/`genomics` rather than running
   everything in one node allocation. `test_cluster_mode.sh` validates the
   template cheaply via `cellranger testrun` before committing to a real
   run. See `cellranger_scripts/README.md`.

2. **`seqtk/`** — downsamples FASTQs to normalize sequencing depth across
   the 4 numbered batches (not pools — `iMG1`/`iMG1_redo` get the same
   fraction as batch 1). `iMG2` is the depth *target* (lowest reads/cell),
   so it's absent from `downsample_targets.txt`/never downsampled.
   `build_manifest.sh` derives the actual file list per pool straight from
   that pool's `cellranger_configs/*.csv` (single source of truth for
   which FASTQ dirs belong to which pool) rather than hardcoding it again.
   Key point for anyone re-deriving downsampling fractions: **no
   concatenation needed** — `seqtk sample` decides per-read independently,
   so applying one pool's fraction to each of its constituent FASTQ files
   separately is equivalent to concatenating first.

3. **`cellranger_scripts/run_all_downsampled.sh`** — re-runs Cell Ranger on
   downsampled FASTQs for `iMG1`, `iMG1_redo`, `iMG3`, `iMG4` (again, not
   `iMG2` — it's the target, never downsampled, so its original Cell Ranger
   output is used as-is downstream).

4. **`cellbender_scripts/`** — CellBender ambient-RNA removal, one run per
   sample, driven by `cell_quantities.csv` (`pool_dir,sample,cells,droplets`
   — **the canonical manifest of all 75 real (pool, sample) captures**,
   generated by `make_cellbender_params.R` from each pool's
   `metrics_summary.csv`; several other scripts reuse this file rather than
   re-deriving the sample list). See the GPU saga below — CellBender now
   runs via `run_all_cellbender_cpu.sh` (CPU) rather than
   `run_all_cellbender.sh` (GPU, currently blocked).

5. **`preprocessing/01_obj_creation.R`** — assembles all 75 CellBender
   outputs into one Seurat object, on-disk via **BPCells**
   (`data/01_obj_assembly/bpcells` + `metadata.rds`), driven directly off
   `cell_quantities.csv`. Also writes `data/01_obj_assembly/bpcells_persample/<distinct capture id>`
   for each of the 75 captures *before* merging/gene-filtering (raw genes,
   unprefixed barcodes) — later scripts that need per-capture raw matrices
   reload from here rather than re-deriving them. Attaches
   `genotype`/`treatment`/`date`/`batch` from `samples.csv` and sets the
   unified `sample` + distinguishing `pool_dir` metadata described above.

6. **`preprocessing/02_qc.R`** — reads `01`'s BPCells object, computes
   QC metrics (`scCustomize::Add_Cell_Complexity`/`Add_Mito_Ribo`), raw QC
   plots faceted by batch, per-**sample** (unified, 60-level) MAD-based
   thresholds (2x lower / 2.5x upper on log UMI/genes, mito ceiling 5% —
   inherited from `als_cns_scrnaseq`, not yet validated against this
   dataset's actual distributions), flags+filters low-quality cells, saves
   the **filtered** object as BPCells to `data/02_qc/`. Introduces
   `condition` (genotype x treatment, 15-level factor with explicit level
   order) as a plotting/grouping variable distinct from `sample` — don't
   conflate the two.

7. **`preprocessing/03_doubletfinder.R`** — runs DoubletFinder as a
   **75-task SLURM array**, one per (pool, sample) capture (see
   `run_03_doubletfinder.sh`'s self-resubmitting pattern below). Reloads
   raw per-capture matrices from `01`'s `bpcells_persample/` (not `02`'s
   merged object) and re-applies the probe-set gene filter, since those
   per-capture matrices predate it. Matches into `02`'s QC metadata by
   `pool_dir` + unified `sample` (exact categorical match) rather than
   string-prefix-matching the capture id against cell barcodes — capture
   ids can be literal string prefixes of each other (`GALC_asyn_1` is a
   prefix of `GALC_asyn_1_redo`), so prefix matching would silently pull
   the wrong cells. All output/plot filenames key off the distinct capture
   id, never `orig.ident`/`sample` (which is unified and would collide
   `iMG1`/`iMG1_redo`'s files for the same biological sample).
   **Known caveat, flagged in the file itself**: running DoubletFinder per
   individual demultiplexed sample (rather than per whole GEM well) means
   it can't catch doublets formed between two *different* samples sharing
   the same GEM well, and the doublet-rate formula now reflects one
   sample's own recovered cells rather than the whole pool's loaded cells
   — likely an underestimate for this multiplexed design. This was a
   deliberate, explicit tradeoff the user chose (to match
   `als_cns_scrnaseq`'s literal per-sample pattern) after being warned
   about it — don't "fix" it back to per-pool without asking first.

Stages beyond `03` (`04_norm_pca.R` onward, `deseq2/`, `wgcna/`) exist in
the repo as placeholders/drafts carried over from template repos but
haven't been worked on yet in this session.

**Stale duplicates**: `jobs/`, `deseq2/`, and `wgcna/` at the repo root
contain older copies of files that now also live under `preprocessing/`
(from an earlier reorganization commit). Only the `preprocessing/` copies
have been touched/maintained — the root-level duplicates are likely dead
weight worth cleaning up, but haven't been touched since it wasn't asked
for.

## Conventions established in this codebase

- **Self-resubmitting SLURM arrays**: several scripts (`seqtk/seqtk.sh`,
  `cellbender_scripts/run_all_cellbender*.sh`,
  `preprocessing/run_03_doubletfinder.sh`) follow the same pattern: check
  `$SLURM_ARRAY_TASK_ID`; if unset, compute the task count from a manifest
  CSV/file count and `exec sbatch --array=1-N[%throttle] "$0"`; if set, do
  the actual per-task work. Run these with `bash script.sh` (not `sbatch`)
  the first time.
- **BPCells pattern**: `convert_matrix_type(mat, type = "uint32_t")` before
  `write_matrix_dir()` (dgCMatrix's `@x` is always stored as R `double`
  regardless of content, which write_matrix_dir's compressed writer needs
  disambiguated). Materialize to `dgCMatrix` (`as(mat, "dgCMatrix")`)
  before handing a BPCells matrix to DoubletFinder — its synthetic-doublet
  sampling selects the same column more than once, which BPCells' lazy `[`
  operator rejects.
- **Base R `[` subsetting over `dplyr::filter()`** when a metadata subset's
  rownames (cell barcodes) need to stay reliably aligned with a matrix's
  column selection — used deliberately in `03_doubletfinder.R` to avoid
  any doubt about rowname preservation across dplyr verbs.
- **Manifests over re-derivation**: `cell_quantities.csv` and each pool's
  `cellranger_configs/*.csv` are treated as sources of truth that other
  scripts read from, rather than every script re-deriving the same
  pool→sample or pool→fastq-dir mappings independently.
- **PR-only workflow, no direct pushes to `main`** — every change goes
  through a PR on branch `claude/cellranger-multiconfig-setup-z4gdto`.
  Push a commit, open (or update, if the same branch/PR is still open) a
  PR, never push to `main` directly. The user typically **merges PRs
  quickly**, often before the next request — always `git fetch origin main`
  and `git checkout -B claude/cellranger-multiconfig-setup-z4gdto origin/main`
  at the start of new work to check whether the previous PR merged and
  restart the branch cleanly if so. If a PR is still open, keep committing
  to the same branch instead of restarting (the user will say explicitly
  if they're holding off on merging, e.g. "I will not be merging PR #14").
- **The user edits scripts directly on `main` between requests** — expect
  to `git fetch`/diff against the actual current file contents rather than
  assuming the last-pushed version is still current. Their own draft edits
  are useful context (they often encode real, dataset-specific design not
  present in reference repos) but have contained real bugs worth checking
  carefully rather than adopting wholesale — found so far: `separate()`
  column order mismatches, a hardcoded row-index removal
  (`stats[-61, ]`), a metrics filter that matched more rows than intended,
  wrong SLURM account/partition inherited from a different project.
- **Ask before guessing on consequential/ambiguous decisions**, especially
  ones with real scientific stakes (e.g. DoubletFinder granularity) or that
  contradict something just-established as correct. The user explicitly
  invites this ("ask if you're unsure") and has used
  `AskUserQuestion`-style clarifications productively — don't skip it to
  move faster when a wrong guess would mean redoing real compute (CellBender
  GPU-hours, DoubletFinder runs, etc.).
- **`als_cns_scrnaseq`** (github.com/jakesboles/als_cns_scrnaseq,
  cloned read-only for reference at
  `/home/user/jakesboles/als_cns_scrnaseq`) is the user's other project,
  used repeatedly as a style/pattern reference for R preprocessing
  scripts. Important structural difference to keep in mind every time:
  that project has **one GEM well per sample** (non-multiplexed), so
  "per sample" and "per GEM well" are the same unit there. Here they are
  not (FRP multiplexing, up to 16 samples per GEM well) — patterns from
  that repo need deliberate translation, not literal copying, whenever
  the per-sample-vs-per-pool distinction matters (it mattered for both
  `02_qc.R`'s threshold grouping and `03_doubletfinder.R`'s unit of
  computation).

## Debugging narratives worth remembering

- **CellBender environment**: the originally-pushed `Cellbender2.yml` had
  no `dependencies:` list (exported from the wrong directory — the parent
  project folder instead of the actual env path one level down; `conda env
  export` silently produces an empty yml rather than erroring on this). A
  real export was later captured directly from the live environment as
  `real_export.yml`: `cellbender==0.3.0`, `torch==1.13.1`,
  `nvidia-cuda-runtime-cu11==11.7.99`. `make_cellbender_env.sh` now builds
  a fresh env under this project's own directory rather than a colleague's.
- **`cell_quantities.csv` duplicate-row bug (fixed)**: Cell Ranger's
  per-sample `metrics_summary.csv` has *two* rows with `Metric.Name ==
  "Cells"` — the real per-sample count, and the whole GEM well's aggregate
  total under a different `Category` but the same `Metric.Name`.
  `make_cellbender_params.R`'s filter now requires
  `Metric.Name == "Cells" & Category == "Cells"`, plus a hard `stop()` if
  match count is ever not exactly 1.
- **CellBender GPU saga (unresolved, currently worked around)**: jobs on
  `genomics-gpu` intermittently failed with vague CUDA errors. Root cause,
  established via `sacct`/`scontrol` evidence: `qgpu0102` is the *only* one
  of `genomics-gpu`'s four nodes that account `b1042` has fast, non-buy-in
  access to (the other three — `qgpu0101`, `qgpu0519`, `qgpu0520` — showed
  `Reason=Priority`/`QOS=buyin` with 24h+ estimated start times). Every
  `b1042` GPU job funnels onto that one node, and multiple concurrent
  `--gres=gpu:a100:1` allocations were landing on its single physical GPU
  at once — confirmed by `sacct` showing array tasks starting at the
  identical second, only one surviving each time, even at `%1` throttling
  (ruling out self-contention as the sole cause). **Not yet reported to
  HPC support.** Current workaround: `run_all_cellbender_cpu.sh` runs
  CellBender CPU-only (`--cpu-threads` instead of `--cuda`) on the reliable
  `genomics` CPU partition — slower per sample but has completed
  successfully. Full diagnostic writeup in `cellbender_scripts/README.md`.
- **Other CellBender fixes**: report-generation crashes at the very end of
  otherwise-successful runs with `UnicodeDecodeError` — cosmetic, fixed by
  exporting `LC_ALL=en_US.UTF-8`/`LANG=en_US.UTF-8` (this cluster's default
  locale is ASCII).

## Open items / things worth revisiting

- CellBender GPU oversubscription on `qgpu0102` is unresolved — worth an
  HPC support ticket with the `sacct` evidence in
  `cellbender_scripts/README.md` if GPU speed ever becomes worth pursuing
  again.
- `02_qc.R`'s MAD threshold multipliers are inherited from
  `als_cns_scrnaseq`, not re-derived for this dataset — sanity-check
  against the actual raw QC plots once real data has run through it.
- `03_doubletfinder.R`'s doublet-rate formula likely underestimates the
  true rate for this multiplexed design (see caveat above) — revisit if it
  matters for the analysis.
- Root-level `jobs/`, `deseq2/`, `wgcna/` duplicate content already present
  under `preprocessing/` — worth cleaning up at some point.
- Stages `04` onward are untouched placeholders as of this writing.
