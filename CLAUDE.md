# img_scfrp — working notes for continuing this pipeline

This file is a memory dump of what's been established across a long session
building this pipeline, so future work doesn't have to re-derive it. Written
by Claude; update it as the project evolves rather than letting it go stale.

## Top-level project structure: `monocultures/` vs `cocultures/` (as of Oct 28 2026)

The repo was restructured (commit `61c2414`, "New project structure.
Monoculture and coculture analysis will all happen here") into two
top-level arms sharing one project/repo:

- **`monocultures/`** — everything described in the rest of this file below
  (all pipeline stages, scripts, configs) is this arm. It's a pure
  `git mv` of what used to live at repo root directly into `monocultures/`
  — same `r_scripts/`, `jobs/`, `cellbender_scripts/`, `cellranger_configs/`,
  `cellranger_scripts/`, `seqtk/`, `samples.csv` subtree, just one directory
  deeper. **Every path in the rest of this document that looks repo-root-
  relative (`r_scripts/01_obj_creation.R`, `data/02_qc/`,
  `cellbender_scripts/cell_quantities.csv`, etc.) is actually
  `monocultures/<that path>`** — this file keeps the un-prefixed paths in
  prose below only where they already existed pre-restructure, but assume
  the `monocultures/` prefix throughout unless stated otherwise.
- **`cocultures/`** — new arm, iPSC-derived microglia **co-cultured** with
  iPSC-derived cortical neurons (as opposed to the microglia-alone
  monoculture arm). Currently has only `cocultures/configs/*.csv` (Cell
  Ranger `multi` configs); no pipeline built yet. See the dedicated
  "Coculture arm" section below for what's known/open so far.

**Known issue from this restructure — fixed**: every `.R`/`.sh`/manifest
file under `monocultures/` originally still hardcoded the *pre-restructure*
repo root (`setwd("/projects/b1169/boles/img_scfrp")` in all 14
`r_scripts/*.R` files; `REPO="/projects/b1169/boles/img_scfrp"` in
`cellbender_scripts/*.sh`/`seqtk/*.sh`; absolute `Rscript
/projects/b1169/boles/img_scfrp/r_scripts/...` and `#SBATCH --output
.../logs/...` lines in every `jobs/*.sh`/`cellranger_scripts/*.sh`; FASTQ
paths baked into `cellranger_configs/*.csv` and `seqtk/downsample_manifest.txt`)
instead of the new `/projects/b1169/boles/img_scfrp/monocultures`. Fixed by
a mechanical `/monocultures` path-segment insertion across every affected
file (committed directly to this branch, no PR, per the user's explicit
request for this one change) — verified via a repo-wide grep that no
`.../img_scfrp/<subdir>` reference lacking the `monocultures/` segment
remains under `monocultures/`. One pre-existing, unrelated quirk left as-is:
`seqtk/get_seqtk_proportions.R`'s `setwd("/gpfs/projects/b1169/boles/img_scfrp/monocultures/")`
uses a `/gpfs/`-prefixed automounter path (different from every other
script's plain `/projects/...`) — not a restructure artifact, just an
existing inconsistency, so only the `/monocultures` segment was inserted
into it, the `/gpfs/` prefix was left alone.

**Also noted, low-priority**: a stale, byte-identical duplicate of
`monocultures/samples.csv` still sits at repo root (`/samples.csv`) from
before the restructure moved it — harmless (nothing reads it) but worth
deleting next time this area is touched.

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

(All of the above, and everything through the rest of this file up to the
"Coculture arm" section, describes the **monoculture arm specifically** —
i.e. everything under `monocultures/`.)

## BPCells is the standard data I/O pattern — always, not just for 01-03

**Every pipeline stage, from `01_obj_creation.R` onward, must read and write
its Seurat object as separate on-disk pieces (BPCells counts/data matrix +
metadata as its own RDS), never as one monolithic `saveRDS(obj, ...)`/
`readRDS(...)` of a whole `Seurat` object.** This isn't an incidental choice
made for 01-03 — it's the standing convention for this entire codebase,
including every future stage (`05_integration_*.R` onward) once they're
revisited. Reasons this matters enough to be a hard rule here: the full
merged object has up to ~60k+ cells x the probe-set gene count, and every
stage after `01` only ever needs to touch one on-disk matrix layer (raw
counts, or normalized data) plus a metadata data frame — never both layers,
never the whole `Seurat` object with all its reductions/graphs at once.
Loading everything into memory as one RDS defeats that and risks OOM on
shared HPC nodes for no benefit.

The pattern, consistently applied in `01`/`02`/`04`:
- **Write**: `convert_matrix_type(mat, type = "uint32_t")` (only needed for
  integer count matrices, not normalized `data` layers) →
  `write_matrix_dir(mat = ..., dir = "data/<stage>/bpcells[_something]")`.
  Save whatever else the object carries (metadata, PCA/other reductions,
  variable features, DoubletFinder calls, ...) as **separate** `.rds` files
  in the same `data/<stage>/` directory — never bundled back into one RDS
  object.
- **Read**: `open_matrix_dir("data/<stage>/bpcells...")` +
  `readRDS("data/<stage>/metadata.rds")` (+ any other `.rds` pieces) →
  `CreateSeuratObject(counts = ..., meta.data = ...)`, then reattach
  reductions/extra layers with `obj[["pca"]] <- readRDS(...)` etc. as
  needed. Every script so far ends with a comment block showing exactly
  this reconstruction snippet for the next script to copy — keep doing
  that.
- **Per-capture matrices**: `01_obj_creation.R` and `03_doubletfinder.R`
  additionally write/read individual `bpcells_persample/<capture id>`
  directories so later scripts needing one capture's raw (pre-merge,
  pre-QC) matrix don't have to subset it back out of the whole-cohort
  object.
- If a stage's on-disk directory can be regenerated by a rerun (e.g.
  `04_norm_pca.R`'s `bpcells_data`), `unlink(dir, recursive = T)` it first
  if it already exists, so a rerun doesn't fail on "path already exists"
  against a stale/incomplete previous attempt.
- `data/` is fully gitignored (HPC-only, regenerable from upstream stages)
  — never try to commit anything under it.

If a future draft script (yours, mine, or the user's own edit) reads or
writes a whole object via a single `saveRDS`/`readRDS(...)  # full Seurat
object`, that's a regression from this convention and should be flagged/
fixed, not adopted quietly. This convention is scoped to the monoculture
arm's established scripts (`monocultures/r_scripts/`); apply the same
standard when building out the coculture pipeline unless the user says
otherwise.

## Directory structure: `r_scripts/`, `jobs/`, `data/`, `results/` (under `monocultures/`)

As of the user's own repo-wide reorganization (and now nested one level
deeper under `monocultures/` by the Oct 28 2026 restructure above), the
analysis pipeline (everything formerly under `preprocessing/`, `wgcna/`,
`deseq2/`) follows one flat, consistent layout — **the same split applies
to every stage, without exception**:
- **`monocultures/r_scripts/`** — every `.R` script for the analysis
  pipeline, flat, no subfolders (`01_obj_creation.R` through `06_qc2.R`,
  `wgcna*.R`, `deseq2*.R`).
- **`monocultures/jobs/`** — every SLURM `.sh` job script for the analysis
  pipeline, flat, matching names to their `r_scripts/` counterpart
  (`run_01_obj_creation.sh`, `run_wgcna.sh`, etc.).
- **`monocultures/data/<stage>/`** — unchanged from the BPCells convention
  above: each stage's BPCells matrices + RDS pieces (counts, normalized
  data, reductions, metadata, hdWGCNA's network payload, ...). Gitignored,
  HPC-only, regenerable from upstream stages.
- **`monocultures/results/<stage>/`** — **replaces the old separate
  `plots/<stage>/` + `tab_data/<stage>/` split** — every plot (`.png`) and
  every tabular result (`.csv`, `.txt`) a stage produces now lands in one
  `results/<stage>/` folder. Also gitignored (regenerable by rerunning the
  stage). Each script uses a single `results_dir <- "results/<stage>/"`
  variable for both, rather than separate `plots_dir`/`tab_dir`/`csv_dir`
  variables. `cellranger_scripts/`, `cellbender_scripts/`, and `seqtk/`
  were deliberately **not** touched by this reorg — they're a separate,
  already-working part of the pipeline with their own established
  per-purpose folders, out of scope by explicit user choice.

This means `data/<stage>/` and `results/<stage>/` share the exact same
`<stage>` name for a given script (e.g. `06_qc2.R` reads/writes
`data/06_qc2/` for its BPCells/RDS pieces and `results/06_qc2/` for its
plots/CSVs) — a stage's outputs are never split across more than these two
trees. `wgcna.R`/`deseq2.R` are no longer special cases here: they used to
depart from this pattern (a flat `wgcna/` directory, or `tab_data/deseq2/`)
before this reorg, but now follow `data/wgcna/` + `results/wgcna/` and
`data/deseq2/`-if-needed + `results/deseq2/` like everything else. If you
ever see a script writing to a bare `wgcna/`/`deseq2/` path, a separate
`plots/`/`tab_data/` tree, or reading `.R`/`.sh` scripts from
`preprocessing/`/`wgcna/`/`deseq2/`, that's a leftover from before this
reorg and should be fixed to match. (All of this is, as of Oct 28 2026,
*itself* nested one level deeper under `monocultures/` — see the top-level
structure section above — but the `<stage>`-keyed `data/`+`results/` split
described here is otherwise unchanged.)

## HPC environment (Northwestern Quest)

- Repo/storage lives at `/projects/b1169/boles/img_scfrp`. As of the Oct 28
  2026 restructure, `data/`, `results/`, `logs/`, `fastq/`,
  `downsampled_fastq/`, `cellranger/`, `cellbender/` live one level
  deeper, under the relevant arm — `monocultures/data/`,
  `monocultures/logs/`, etc. for the existing pipeline, and (once built)
  `cocultures/data/`, `cocultures/logs/`, etc. for the new arm — matching
  the `.gitignore`'s `*/data/*`-style patterns (anchored one level deep,
  deliberately, for exactly this split). Every script's `setwd()`/`REPO=`
  now points at `.../img_scfrp/monocultures` — see the top-level structure
  section above for the fix. **`envs/` is the one exception** — it stayed
  at the repo root (`/projects/b1169/boles/img_scfrp/envs/cellbender`, one
  level *above* both `monocultures/` and `cocultures/`), shared by both
  arms rather than duplicated per arm, since it's the same CellBender
  install either way. The `.gitignore`'s own `envs/*` pattern (no `*/`
  prefix, unlike every other line in that block) already reflected this —
  an early pass at fixing monoculture's stale paths incorrectly nested
  `envs/` under `monocultures/` anyway (missed that one inconsistency in
  the `.gitignore`); corrected back to the shared root-level path in every
  script that referenced it (`monocultures/cellbender_scripts/
  run_all_cellbender*.sh`, `make_cellbender_env.sh`, `README.md`, and
  `cocultures/cellbender_scripts/run_all_cellbender_cpu.sh`).
- Cell Ranger reference/probe-set files (shared across both arms, unaffected
  by the restructure): `/projects/p31535/boles/cellranger_references/`.
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
  `mamba`, then `source activate /projects/b1169/boles/img_scfrp/envs/cellbender`
  — note this path is **not** under `monocultures/` or `cocultures/` (see
  the bullet above); both arms' `run_all_cellbender_cpu.sh` activate it by
  this same absolute path, independent of their own `REPO` variable.

## Pipeline stages, in order (monoculture arm — `monocultures/`)

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

5. **`r_scripts/01_obj_creation.R`** — assembles all 75 CellBender
   outputs into one Seurat object, on-disk via **BPCells**
   (`data/01_obj_assembly/bpcells` + `metadata.rds`), driven directly off
   `cell_quantities.csv`. Also writes `data/01_obj_assembly/bpcells_persample/<distinct capture id>`
   for each of the 75 captures *before* merging/gene-filtering (raw genes,
   unprefixed barcodes) — later scripts that need per-capture raw matrices
   reload from here rather than re-deriving them. Attaches
   `genotype`/`treatment`/`date`/`batch` from `samples.csv` and sets the
   unified `sample` + distinguishing `pool_dir` metadata described above.

6. **`r_scripts/02_qc.R`** — reads `01`'s BPCells object, computes
   QC metrics (`scCustomize::Add_Cell_Complexity`/`Add_Mito_Ribo`), raw QC
   plots faceted by batch, per-**sample** (unified, 60-level) MAD-based
   thresholds (2x lower / 2.5x upper on log UMI/genes, mito ceiling 5% —
   inherited from `als_cns_scrnaseq`, not yet validated against this
   dataset's actual distributions), flags+filters low-quality cells, saves
   the **filtered** object as BPCells to `data/02_qc/`. Introduces
   `condition` (genotype x treatment, 15-level factor with explicit level
   order) as a plotting/grouping variable distinct from `sample` — don't
   conflate the two.

7. **`r_scripts/03_doubletfinder.R`** — runs DoubletFinder as a
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

8. **`r_scripts/04_norm_pca.R`** — rebuilds the whole-cohort object from
   `02`'s filtered BPCells counts + metadata, attaches DoubletFinder's
   `DF.unadj`/`DF.adj` per-cell classifications from `03`'s
   `metadata_persample/` (matched by exact rowname/cell-barcode, which
   works because both `01` and `03` prefix barcodes with the distinct
   capture id the same way — no `pool_dir`/`sample` join needed here,
   unlike `03`'s own join back into `02`'s metadata). Hard `stop()`s if any
   of the 75 expected `metadata_persample/*.rds` files are missing (an
   unfinished/failed `03` array) or if cells don't reconcile 1:1 between
   `02`'s object and `03`'s DF calls (02/03 have diverged). Does **not**
   remove flagged doublets — `DF.unadj`/`DF.adj` are attached as metadata
   for a later stage to act on. Then `NormalizeData()`/
   `FindVariableFeatures()`/`ScaleData()`/`RunPCA(npcs = 50)`, an elbow
   plot, and a geometric max-distance elbow finder (ported from
   `als_cns_scrnaseq/r_scripts/08_pca_evaluation.R`) to suggest a PC count
   quantitatively rather than eyeballing the plot. PCA scatter plots
   (dims 1-5, all pairs) colored by `condition`/`genotype`/`treatment`/
   `batch` for visual QC of batch/condition structure in the earliest PCs.
   Output moved from a single `data/04_sct_pca/obj.rds` to separate BPCells
   (`bpcells_data`, the normalized `data` layer) + `metadata.rds` +
   `pca.rds` + `variable_features.rds` pieces, matching the BPCells
   standard above; directory renamed `data/04_sct_pca` → `data/04_norm_pca`
   to match (this does log-normalization via `NormalizeData()`, not
   `SCTransform`).

9. **`r_scripts/05_integration_harmony.R`** — the only integration
   method kept (the user removed `05_integration_cca.R` and
   `05_integration_harmony_2-3.R` as unneeded; `run_05_integration_cca.sh`
   was deleted alongside it as a now-dangling job script). Loads `04`'s
   normalized `data`/PCA/variable-features BPCells/RDS pieces, plus a
   **real raw counts matrix pulled from `02_qc.R`'s `data/02_qc/bpcells`**
   (`04` never re-saves counts of its own — its cell set is identical to
   `02`'s filtered set, so `02`'s bpcells directory is still the source of
   truth for counts at this point) — checked/reordered to `04`'s cell
   order with the same hard-`stop()`-on-mismatch caution `04` uses for its
   own 02/03 reconciliation, rather than assumed. **Splits the RNA assay
   by the unified `sample` column (60 levels), not `pool_dir`**, before
   `IntegrateLayers(method = "HarmonyIntegration")` — an explicit choice
   the user made after being shown the tradeoff via `AskUserQuestion`
   (splitting by `pool_dir` would additionally have corrected for any
   technical difference between `iMG1`'s and `iMG1_redo`'s captures of the
   same sample; splitting by `sample` does not, since both captures share
   that value and land in the same integration layer). **Don't change this
   split-by variable without asking first** — same standing rule as `03`'s
   doublet-rate caveat. Computes a pre-integration UMAP (off `pca`) and a
   post-integration UMAP (off the new `harmony` reduction, `dims = 1:10`
   throughout, matching this project's established PC count) for a
   before/after comparison, both faceted the same way as `04`'s PCA plots
   (`condition`/`genotype`/`treatment`/`batch` — fixed a bug inherited from
   the old placeholder script, which grouped the first panel by the
   60-level `sample` while only supplying a 15-color palette; the
   15-level `condition` is this project's established stand-in for that
   panel). **Writes BPCells counts (`bpcells_counts`) and normalized data
   (`bpcells_data`) matrices, `metadata.rds`, `harmony.rds`, and
   `harmony_umap.rds` to `data/05_integration/`** — both matrices are
   written (not just `data`, as `04` does) specifically so `06_qc2.R`
   could build directly off this one stage's output without reaching
   further back. (The old `06_clustering.R`/`06b_clustering_markers.R`
   placeholders that used to have a dangling-reference problem here were
   removed by the user and replaced by `06_qc2.R` below; `wgcna.R` and
   `deseq2.R` have both since been repointed to `data/06_qc2/` instead of
   reaching back to `05` directly — see their entries below.)

10. **`r_scripts/06_qc2.R`** — a second QC pass, addressing
    `03_doubletfinder.R`'s own documented caveat that per-sample
    DoubletFinder can't catch doublets formed between two *different*
    samples sharing a GEM well, **and**, as of two follow-up commits
    (`c4b4c35`/`2f7d1b7`), also removing non-microglial cell clusters.
    Loads `05`'s BPCells/RDS pieces, clusters on `05`'s own `harmony`
    reduction (`resolution = 2`, same as the old `06_clustering.R`), now
    with `cluster.name = "qc_cluster"` (renamed from `doublet_qc_cluster`
    once this stage started removing more than just doublet clusters —
    **if you see `doublet_qc_cluster` anywhere, that's the pre-`c4b4c35`
    name**), and plots per-cluster `DF.adj` composition
    (`clusters_by_doublets.png`) so doublet-enriched clusters can be
    spotted **by eye** — this is a manual QC judgment call, not
    programmatic, and the resulting cluster numbers are specific to that
    exact clustering. Also plots, per `qc_cluster`, violin + UMAP panels
    for `AIF1`/`CSF1R` (core microglial markers) plus
    `nCount_RNA`/`nFeature_RNA`/`percent_mito`, to spot marker-low clusters
    that are likely contaminating non-microglial cell types (also by eye).
    Removes those *whole clusters* — both the doublet-enriched ones
    (`doublet_clusters <- c(31, 35, 36)`) and the non-microglial ones
    (`non_microglia_clusters <- c(26, 33)`) — together in one `subset()`
    call (`remove <- c(doublet_clusters, non_microglia_clusters)`; note a
    one-commit-lived typo, `non_microglial_clusters`, introduced in
    `c4b4c35` and fixed one commit later in `2f7d1b7` — if a diff/blame
    ever surfaces the broken intermediate version, it's already fixed on
    `main`), logging what was removed to `removed_clusters.csv` (renamed
    from `removed_doublet_clusters.csv` for the same reason as the cluster
    column rename above), then repeats `04`/`05`'s normalize→PCA→Harmony
    steps on the filtered cell set and saves in the same BPCells pattern as
    `05` (`bpcells_counts`, `bpcells_data`, `metadata.rds`, `harmony.rds`,
    `harmony_umap.rds`) plus its own fresh `variable_features.rds` — `04`'s
    copy was selected against the pre-filter cohort, so this stage needs
    its own rather than leaving downstream users of the old one (like
    `wgcna.R`) working off a stale gene selection basis. **Now the
    canonical upstream source for `wgcna/` and `r_scripts/deseq2.R`
    scripts** — both were repointed from `data/05_integration/` to
    `data/06_qc2/` once this stage was finished (see their entries below).
    Meant to be run interactively (no `.sh` job script — though the user
    did add one, `run_06_qc2.sh`, for the batch-runnable back half once the
    doublet/non-microglial clusters are known), matching `02_qc.R`'s
    precedent for a script with an eyeballed QC decision in the middle of
    it. **The specific cluster numbers removed (`31, 35, 36` doublets;
    `26, 33` non-microglial) are, like the original doublet-only version,
    tied to this exact `resolution = 2` clustering on this exact upstream
    object — re-eyeball and update them if `05`'s output ever changes.**

11. **`r_scripts/wgcna.R`, `r_scripts/wgcna_stats.R`, `r_scripts/wgcna_viz.R`,
    `r_scripts/wgcna_clusterprofiler.R`** — consensus hdWGCNA network on the
    whole-cohort object (one network per `batch`, reconciled into one
    consensus network — kept consistent with `05`'s own choice to
    integrate at `sample`, not `pool_dir`, granularity: don't introduce a
    finer `pool_dir`/iMG1-vs-iMG1_redo split here that upstream integration
    didn't make either). `wgcna.R` loads `06_qc2.R`'s BPCells/RDS pieces
    (not a monolithic `obj_consensus.rds`) — including `06`'s own fresh
    `variable_features.rds` (needed for `ScaleMetacells()`/
    `RunPCAMetacells()`; originally reached back to `04_norm_pca.R` for
    this before `06_qc2.R` existed, repointed once `06` started saving its
    own copy). Builds metacells via `MetacellsByGroups()`, grouped by
    `sample` (not `orig.ident` — identical value throughout this pipeline,
    but `sample` is this project's documented name for it) plus
    `genotype`/`treatment`/`batch`/`condition` for convenience (redundant
    with `sample` for partitioning, since each sample implies exactly one
    of each, but carried through so they land as real columns on the
    metacell object's own metadata). Fixed a stale `reduction =
    'harmony_pca'` (05/06 actually save the reduction as `"harmony"`) and
    the same 60-level-`sample`-vs-15-color-palette bug already fixed in
    `05` (metacell UMAP now grouped by `condition`). **hdWGCNA isn't
    BPCells-aware**: unlike Seurat's own `NormalizeData()`/`ScaleData()`/
    `RunPCA()`/`IntegrateLayers()` (which do operate directly on BPCells
    matrices, as `04`/`05`/`06` already demonstrate),
    `SetupForWGCNA(gene_select = "fraction")` crashes with `invalid
    subscript type 'S4'` against a lazy `IterableMatrix` — counts/data are
    materialized to `dgCMatrix` right after loading, before anything
    hdWGCNA-specific runs (same "materialize before handing to a
    non-BPCells-aware tool" caution as `03_doubletfinder.R`, just applied
    to the whole assay here since hdWGCNA's pipeline touches it
    throughout, not one step). **Only saves what this stage adds** —
    `misc_wgcna_consensus.rds` (hdWGCNA's own network/module/eigengene
    payload; metacell/module-scale, not cell x gene-matrix scale, so one
    RDS is appropriate, same reasoning as `pca.rds`/`harmony.rds`
    elsewhere) and `metadata.rds` (with `ModuleExprScore()`'s per-cell
    module score columns attached) — counts/data/harmony/harmony_umap are
    unchanged from `06`'s cell set, so downstream scripts reach back into
    `data/06_qc2/` directly rather than have them re-saved here. **Output
    directory**: `wgcna.R`'s RDS/network pieces (`misc_wgcna_consensus.rds`,
    `metadata.rds`) go to `data/wgcna/`; its plots/CSVs go to
    `results/wgcna/` — this used to be a special-cased flat `wgcna/`
    directory for everything (a departure from the rest of the project's
    `data/<stage>/`+`results/<stage>/` split), but the repo-wide
    `r_scripts/`/`jobs/`/`data/`/`results/` reorg brought it in line with
    every other stage; it's no longer an exception. `wgcna_stats.R`/
    `wgcna_viz.R`/`wgcna_clusterprofiler.R` were updated to match so they
    can still find `wgcna.R`'s output — if this convention changes again,
    all four scripts' path variables need to move together.
    `wgcna_stats.R` (UCell module scoring, kNN-smoothed over `harmony`,
    genotype x treatment `lmer` models with a `batch` random intercept,
    `emmeans`/`multcomp` post-hoc letters) reconstructs its own object from
    `06`'s counts (materialized to `dgCMatrix` before `AddModuleScore_UCell`
    — same reasoning as above) and recomputes UCell scores itself
    from `wgcna.R`'s `module_members_consensus.csv`, rather than reusing
    hdWGCNA's own internal module scores — dropped a chunk of dead code
    past the original script's working section that referenced an
    undefined `modules` variable (a leftover half-finished edit). Fixed
    the same `orig.ident`-vs-`sample` and `harmony_pca`-vs-`harmony`
    naming issues here too. Its trailing quick-eyeball `VlnPlot_scCustom()`
    call (an unsaved, interactive-only sanity check against a hardcoded
    example module name) is now guarded by `if (interactive())` so
    `jobs/run_wgcna_stats.sh` can run it headless via `Rscript` without
    that line erroring out (no display to print to, and the hardcoded
    module color may not exist in a given run's actual module set) after
    everything that actually gets saved has already been written.
    `jobs/run_wgcna_stats.sh`'s resources have since been bumped further
    (now `b1042`/`genomics`, `200G` mem, `24:00:00`) — ambient, user-made
    edit, not something to revert.
    `wgcna_viz.R` (module radar plot,
    module-module correlogram, per-sample module eigengene dot/violin
    plots) reconstructs from `06`'s pieces plus `wgcna.R`'s
    `metadata.rds`/`misc_wgcna_consensus.rds`, explicitly setting
    `obj@misc$active_wgcna <- "wgcna_consensus"` since that's normally set
    automatically by `SetupForWGCNA()` and needs to be replicated by hand
    on a manually-reconstructed object. `wgcna_clusterprofiler.R`
    (`clusterProfiler::enricher()` pathway enrichment per module, against
    combined MSigDB C2 canonical-pathway + C5 GO gene sets) only depends
    on `results/wgcna/module_members_consensus.csv` — no object reconstruction
    needed, otherwise unchanged from its pre-cleanup version.

12. **`r_scripts/deseq2.R`** — pseudobulk DESeq2 across all pairwise
    genotype x treatment `condition` comparisons (105 = 15 choose 2), with
    `batch` as a covariate. Loads `06_qc2.R`'s BPCells counts + metadata
    only (no normalized data/reductions needed for
    `AggregateExpression()` pseudobulk). **Rewrote the
    design/contrast structure, not just the data source** — the prior
    draft's `design = ~ batch + sample` with pairwise `sample`-vs-`sample`
    contrasts (all 1770 pairs) had two compounding problems: `batch` is
    fully determined by `sample` (perfectly collinear, rank-deficient
    design matrix) and, independent of that, a `sample`-level design has
    exactly one pseudobulk column per level (zero residual df for
    dispersion estimation) regardless. Switched to contrasting `condition`
    (each level has 4 real batch replicates, one sample per batch, giving
    real residual df) with `batch` as a covariate — an explicit choice the
    user confirmed via `AskUserQuestion` after the design was flagged as
    unfittable, not something to change back without asking first. Kept
    the per-reference-relevel-and-refit loop structure (needed because
    `lfcShrink(type = "apeglm")` only shrinks a named coefficient, not an
    arbitrary `results(contrast = ...)`), just re-targeted at `condition`
    instead of `sample`. Also fixed: the prior draft re-factored
    `treatment` against lowercase levels (`"control"/"myelin"/"asyn"`)
    that only ever existed in `01`'s raw metadata — by `02_qc.R` onward
    the real values are Title Case (`Control`/`Myelin`/`Asyn`), so every
    cell would have silently become `NA`; deleted that re-factoring
    entirely rather than fixing the levels, since `06_qc2.R`'s metadata
    already carries correctly-leveled `genotype`/`treatment`/`batch`/
    `condition` through from `02` and doesn't need re-deriving. Also
    `samples <- levels(meta$sample)` was `NULL` (`sample` is a plain
    character column, never converted to a factor anywhere in this
    pipeline, so `levels()` on it returns `NULL` and `combn(NULL, 2)`
    would have errored) — moot now that the loop is keyed on `condition`,
    which *is* a real factor from `02_qc.R` onward. **Note on
    `AggregateExpression()`'s `layer`/`slot` argument — unsettled, flagged
    again below**: an earlier pass changed `slot = "counts"` to
    `layer = "counts"` here (matching Seurat v5's current argument name,
    the same class of `slot`-vs-`layer` API drift that broke hdWGCNA — see
    the debugging narrative below), but the Oct 28 2026 restructure commit
    (`61c2414`) removed the `layer = "counts"` argument from this call
    entirely, leaving a bare `AggregateExpression(obj, assays = "RNA",
    group.by = "sample", return.seurat = F)` that falls back to whatever
    `AggregateExpression()`'s own default layer is. **Not yet verified
    whether that default is actually `"counts"`** (vs. `"data"`, which
    would silently pseudobulk normalized/log values instead of raw counts
    — wrong for DESeq2, which expects integer counts) — this script hasn't
    been run since that change. Check this before trusting `deseq2.R`'s
    output; re-add `layer = "counts"` explicitly if the default turns out
    to be anything else. `r_scripts/deseq2_gsea.R` (fgsea per pairwise
    comparison against a fixed external gene-set CSV) reads `deseq2.R`'s
    `lfc_shrunk_*.csv` output from `results/deseq2/` and writes its own
    results to `results/deseq2_gsea/` — repointed for the reorg but
    otherwise untouched/not yet reviewed for correctness the way `deseq2.R`
    itself was. `deseq2_final.R`, an exact pre-fix duplicate of `deseq2.R`
    (same `~ batch + sample` / lowercase-`treatment` bugs), was deleted
    outright as dead weight rather than moved.

13. **`r_scripts/milo.R`** — Milo (`miloR`) differential neighborhood-
    abundance testing, two analyses on the same neighborhoods (built once
    and reused): **(1)** each non-KOLF genotype vs `KOLF`, collapsed across
    treatment (`design = ~ 0 + genotype + batch`, `batch` as a covariate
    for the same reason as `deseq2.R`'s design — each genotype has 4 real
    batch replicates once pooled across treatment); **(2)** genotype x
    treatment interaction — for each non-KOLF genotype and each non-Control
    treatment, whether that genotype's *shift* in neighborhood abundance
    under treatment differs from KOLF's shift under the same treatment
    (`design = ~ 0 + condition + batch`, contrast `(conditionG_T -
    conditionG_Control) - (conditionKOLF_T - conditionKOLF_Control)` per
    genotype x treatment pair — 4 genotypes x 2 treatments = 8 tests). This
    was proposed via `AskUserQuestion`-style writeup in chat before being
    built, since it's a real modeling-choice fork (baseline genotype
    difference vs. genotype-specific treatment response) — the user
    confirmed this framing before it was implemented. First-pass script,
    adapted from a working script from an unrelated other project (its
    contrast variable, reduced-dim names, and PC count were specific to
    that project's own object and don't apply here) — rebuilt against this
    project's own metadata schema and established PC count (`dims = 1:10`).
    Loads `06_qc2.R`'s BPCells counts + `harmony`/`harmony_umap` (no
    normalized `data` layer — Milo's core machinery never touches
    expression values, only the neighbor graph/reduced-dim coordinates/
    metadata). Recomputes `FindNeighbors()` on the `harmony` reduction
    (`06_qc2.R` doesn't save its own `Graph` object) immediately before
    handing it to `buildFromAdjacency()`, reusing the same neighbor
    structure the rest of the pipeline is built on rather than letting
    `miloR::buildGraph()` build an independent one. **`as.SingleCellExperiment()`
    upper-cases reduction names** (`harmony` → `"HARMONY"`) — confirmed via
    a printed `reducedDimNames(sce)` check rather than assumed, since this
    is the first script in this pipeline to go through an SCE conversion;
    if a future `miloR`/Seurat version changes that casing behavior, this
    is the place to look. Explicit `stop()` guard on `genotype`/`condition`
    actually surviving the `colData()`/`data.frame()` round-trip as the
    expected factor, rather than re-deriving levels from a separately-
    maintained hardcoded list (the class of bug that bit `deseq2.R`'s
    `treatment` levels) — `condition`'s expected levels are derived from
    the same `expected_genotype_levels` vector analysis 1 already checked
    against, not a second independent hardcoded list. `batch`'s
    `"Batch 1"`-style labels are sanitized (spaces stripped) before going
    into either design matrix, since `testNhoods()`'s internal
    `edgeR`/`limma` contrast building wants syntactically valid column
    names. Saves only what this stage adds — `data/milo/milo_nhoods.rds`
    (neighborhoods/counts/graph/distances; neighborhood-scale, not cell x
    gene-matrix scale, so one RDS is appropriate), `design.rds`, and
    `design_interaction.rds` — not the counts matrix, unchanged from
    `06_qc2.R`. **Not yet run** — some `miloR` accessor names
    (`nhoodDistances()`, `nhoodGraph()`) and the `HARMONY`/`HARMONY_UMAP`
    capitalization assumption are believed correct but unverified; expect
    to need small fixes on first execution. The user plans to first try
    running just analysis 1 interactively, and only attempt analysis 2 if
    that succeeds. No `.sh` job script yet — not asked for.

14. **`r_scripts/07_clustering.R`** — clusters the whole-cohort object at
    two resolutions (`res2`, `res0-4`) on the `harmony` reduction and runs
    `FindAllMarkers()` for each, plus per-resolution UMAP/composition
    plots (by genotype/batch/treatment/condition). Numbered `07` in the
    pipeline's own file naming (comes right after `06_qc2.R`) but
    documented here as entry 14 since that's the order it was actually
    worked on in this file's chronological memory-dump style — don't read
    the numbering as meaning it runs after `wgcna`/`deseq2`/`milo`; it
    doesn't depend on any of them, only on `06_qc2.R`.
    **Performance fix**: the user's own first draft handed
    `obj[["RNA"]]$data` straight from `open_matrix_dir()` — a lazy,
    on-disk, column-major/cell-major BPCells matrix — directly to
    `FindAllMarkers()`, called twice (once per resolution). This is a
    documented slow (in some Seurat versions, warned-about —
    `satijalab/seurat#10365`) path: Wilcoxon marker testing needs per-gene
    access across all cells, the wrong direction for a matrix stored
    cell-major. The user suggested investigating a transpose-based fix
    before `FindAllMarkers()`; researched via `WebSearch` rather than
    guessed (presto's `wilcoxauc()` — the fast backend Seurat's
    `FindAllMarkers()` uses automatically when installed — expects
    genes-as-rows/cells-as-columns, i.e. no transpose relative to
    Seurat's native layout; BPCells' own docs recommend a transpose
    specifically for keeping a matrix lazy/on-disk at a scale that can't
    fit in memory, which doesn't apply once you fully materialize).
    Landed on: materialize `data` once, up front, outside the resolution
    loop (a single one-pass read happens regardless of orientation once
    you're materializing anyway) — same established pattern as
    `wgcna.R`/`wgcna_stats.R`/`milo.R` for handing BPCells data to a
    non-BPCells-aware tool. `counts` is deliberately left as lazy
    BPCells — clustering/marker-finding never touch it, so materializing
    it would cost memory for nothing. Also added an explicit
    `requireNamespace("presto")` check with a `warning()` if missing,
    since that's a bigger lever on `FindAllMarkers()` speed than anything
    else in the script (presto vs. base-R fallback) — **superseded one
    commit later** (`3d8a36f`, "Tweaked library loading in 07"): the
    `requireNamespace()`/`warning()` check was removed and replaced with a
    plain `library(presto)` alongside the script's other
    `suppressMessages({ library(...) })` imports — presto is now assumed
    installed (which it is, on this project's R environment) rather than
    soft-checked. Also fixed a real
    bug in the user's draft: the final `saveRDS()` referenced an
    undefined `data_dir` (should have been `data_out_dir`) — would have
    errored at the very last line, after all the expensive clustering/
    marker-finding work was already done. **Two further speed options
    raised with the user but not applied without confirmation** (both
    change what gets computed/reported, not just how fast): restricting
    `FindAllMarkers(features = VariableFeatures(obj))` (~10x fewer genes
    tested, but genuinely excludes non-variable markers from the output);
    splitting the two resolutions into a 2-task SLURM array for
    wall-clock parallelism (real gain, adds the complexity of computing/
    caching the neighbor graph once for both tasks to share rather than
    recomputing it twice). `jobs/run_07_clustering.sh` added to run this
    non-interactively; its memory ask has since been bumped from the
    original `128G` to **`300G`** (`8329d9b`, "Upped resource ask due to
    OOM error") — i.e. the materialization fix above did *not* avoid an
    OOM in practice at `128G`/`b1169`/`b1169`, 12h; worth remembering if
    further memory tuning comes up, since `300G` is already a large ask
    for this partition.

**No more stale duplicates as of the `r_scripts/`/`jobs/`/`results/` reorg**:
`preprocessing/`, `wgcna/`, and `deseq2/` (as directories) no longer exist
at all — every `.R` script from all three moved into `r_scripts/`, every
`.sh` job script moved into `jobs/`. This also finally resolved the
long-standing stale-duplicate problem this file used to track here: `jobs/`
used to hold an old, untouched, out-of-date copy of several job scripts
(wrong SLURM accounts, dangling references to deleted `.R` files) alongside
the real maintained ones under `preprocessing/`/`wgcna/` — those old copies
were deleted outright rather than merged, since `jobs/` is now the single
canonical location and there's nothing left to duplicate against. If you
ever see a `.R`/`.sh` file living outside `r_scripts/`/`jobs/` for this
pipeline, or a stray `preprocessing/`/`wgcna/`/`deseq2/` directory, that's
a regression — move it, don't resurrect the old split. (As of Oct 28 2026,
`r_scripts/`/`jobs/`/`data/`/`results/` are themselves nested one level
deeper under `monocultures/` — see the top-level structure section at the
top of this file.)

## Coculture arm (`cocultures/`) — status as of Oct 28 2026

New arm added alongside `monocultures/` in the same restructure. iPSC-
derived microglia **co-cultured** with iPSC-derived cortical neurons
(as opposed to the monoculture arm's microglia-alone design). `cellranger
multi` has already been run on all 4 pools (confirmed via the user's own
`tree -L 4` of the HPC `cellranger/` output dir) and CellBender's
manifest/run scripts now exist (`cocultures/cellbender_scripts/`, see
below) — no R pipeline built yet beyond that.

`cocultures/configs/*.csv` (4 files: `JSB127_config.csv`,
`JSB134-1_config.csv`, `JSB134-2_config.csv`, `JSB134-3_config.csv`), each
a `cellranger multi` config (`[gene-expression]`/`[libraries]`/`[samples]`
sections, same reference/probe-set paths as the monoculture arm), describe
**4 real GEM wells** — `JSB127`, `JSB134-1`, `JSB134-2`, `JSB134-3`. Cell
Ranger's actual output directories (`cellranger/<dir>/outs/`) use these
exact same names — **no `pool_dirs` renaming needed, unlike monoculture's
4-of-5-pools `_downsampled` suffix** (see CellBender section below).
`fastq_id`s are `127-1`/`127-2` for `JSB127` and `134-1`/`134-2`/`134-3` for
the matching `JSB134-*` pool, all reading from one shared FASTQ directory
(`/projects/b1042/Gate_Lab/boles/coculture_scfrp/fastq`) rather than one
directory per pool (unlike the monoculture arm's per-pool FASTQ dirs — not
confirmed whether this is meaningful or just where the files happen to
currently sit, but not relevant to anything downstream). Note
`JSB134-1_config.csv`'s minor formatting quirk (trailing commas/empty
columns throughout, vs. `JSB127`'s cleaner format) — likely harmless for
`cellranger multi`'s parser (it ran successfully), but worth a second look
if that config is ever the one that fails on a rerun.

**Sample-naming convention — confirmed by the user, no longer a guess**:
`sample_id`s follow `<token1>_<token2>[_<N>]`, e.g. `KOLF_CN`, `KOLF_GALC`,
`GALC_GALC_2`, `GBA_GBA_3`.
- `token1` = **cortical-neuron genotype**, `token2` = **microglia
  genotype** (reversed from this file's original guess, which had it
  backwards). The one deviation: `KOLF_CN` has no real `token2` —
  `CN` here literally stands for "**C**ortical **N**eurons [alone]", i.e.
  a neuron-only well with no microglia added (not a second genotype
  label). Every other `sample_id` follows the real two-genotype pattern,
  including the "matched" homozygous pairs where both tokens are the same
  non-KOLF genotype (`GALC_GALC`, `GBA_GBA`, `GRN_GRN`, `LRRK2_LRRK2` —
  genotype-matched neurons + microglia) and the KOLF-neuron/non-KOLF-
  microglia pairs (`KOLF_GALC`, `KOLF_GBA`, `KOLF_GRN`, `KOLF_LRRK2`, plus
  `KOLF_KOLF` as the all-KOLF control). **10 distinct genotype-pair
  "conditions" total**: `KOLF_CN`, `KOLF_KOLF`, `KOLF_GALC`, `KOLF_GBA`,
  `KOLF_GRN`, `KOLF_LRRK2`, `GALC_GALC`, `GBA_GBA`, `GRN_GRN`,
  `LRRK2_LRRK2`.
- The trailing numeric suffix (`_2`, `_3`, absent = batch 1) is the
  **batch/replicate number** — confirmed: each of the 10 conditions has 3
  independent-batch replicates. **This batch numbering is wholly
  independent of the monoculture arm's `batch`** (different experiments,
  different dates/meaning — don't conflate `cocultures`' batch 1/2/3 with
  `monocultures`' `"Batch 1".."Batch 4"`).
- **No treatment axis** — genotype pairing (the two tokens) is the only
  experimental variable for this arm. Unlike monoculture, there is no
  myelin/asyn-style `treatment` column to carry.

**Important structural difference from monoculture, discovered from the
Cell Ranger output tree — note for the eventual `samples.csv`-equivalent**:
in monoculture, each pool/GEM well maps to exactly one batch (the only
exception being `iMG1`/`iMG1_redo`, where the *same* batch was captured
*twice* in two separate wells). In coculture, it's the opposite kind of
exception: **one well, `JSB134-3`, holds samples from *two different*
batches at once** — its `per_sample_outs/` contains `GBA_GBA_2`/
`GRN_GRN_2`/`LRRK2_LRRK2_2` (batch 2) *and* `GBA_GBA_3`/`GRN_GRN_3`/
`LRRK2_LRRK2_3` (batch 3) side by side, apparently because those 3
homozygous-pair conditions' batch-2 and batch-3 replicates got pooled
together into one well rather than following `JSB134-1`'s (batch 2) /
`JSB134-2`'s (batch 3) per-batch pooling that the other 7 conditions use.
Concretely: batch 1 = all 10 conditions in `JSB127`; batch 2 = 7 conditions
in `JSB134-1` (`KOLF_CN_2`, `KOLF_KOLF_2`, `KOLF_GALC_2`, `KOLF_GBA_2`,
`KOLF_GRN_2`, `KOLF_LRRK2_2`, `GALC_GALC_2`) + 3 conditions in `JSB134-3`
(`GBA_GBA_2`, `GRN_GRN_2`, `LRRK2_LRRK2_2`); batch 3 = the same 7-condition
split in `JSB134-2` (`..._3`) + the same 3-condition split in `JSB134-3`
(`..._3`) — 30 total (pool, sample) captures (10 conditions x 3 batches),
no "redo"-style duplicate captures anywhere in this arm. **Practical
consequence**: unlike monoculture, `pool_dir` alone cannot be used to
infer `batch` for coculture (`JSB134-3` cells span 2 batches) — batch must
be parsed directly from the `sample_id`'s own trailing numeric suffix
whenever the eventual coculture `01_obj_creation.R`-equivalent is built.

**CellBender — built, ready to run** (`cocultures/cellbender_scripts/`,
mirroring the monoculture arm's `make_cellbender_params.R` →
`run_all_cellbender_cpu.sh` pipeline, entry 4 above):
- `make_cellbender_params.R`: identical logic to the monoculture version
  (`Metric.Name == "Cells" & Category == "Cells"` row, `droplets =
  floor(cells * 1.75)`), with `pool_dirs` as a plain identity mapping
  (`JSB127 = "JSB127"`, etc. — no `_downsampled` suffix to track, since no
  downsampling step applies to this arm for now, per the user) and
  `setwd()` pointed at `/projects/b1169/boles/img_scfrp/cocultures`.
  Writes `cocultures/cellbender_scripts/cell_quantities.csv` in the same
  `pool_dir,sample,cells,droplets` format as monoculture's.
- `run_all_cellbender_cpu.sh`: same self-resubmitting CPU-only pattern as
  monoculture's (`b1042`/`genomics`, 16 cpus, 64G, 48h) — activates the
  **shared, repo-root-level CellBender conda env**
  (`/projects/b1169/boles/img_scfrp/envs/cellbender`, one level above both
  arms — see the HPC environment section above) by absolute path, rather
  than building a second, redundant copy under `cocultures/envs/`, since
  it's the identical CellBender/torch/CUDA version either way. Its own
  `REPO` (`.../cocultures`) only governs this arm's own `cellranger/`/
  `cellbender/`/`logs/` paths, not the environment.
- Usage, same as monoculture: `Rscript
  cellbender_scripts/make_cellbender_params.R` then `bash
  cellbender_scripts/run_all_cellbender_cpu.sh`, both run from
  `cocultures/`. **Currently running on the HPC** as of the user's latest
  update — not yet confirmed complete/successful.

**Sequencing-depth diagnostic — run, downsampling decided against** (`cocultures/seqtk/
get_seqtk_proportions.R`, mirroring monoculture's
`seqtk/get_seqtk_proportions.R`): the user explicitly does **not** want a
downsampling step for this arm yet, but asked to see *how different*
sequencing depth actually is across the 4 pools before that decision gets
revisited. Same core computation as monoculture's script (per-pool reads,
mean reads/cell, cells from `cellranger/<pool>/outs/qc_library_metrics.csv`;
`target` = the lowest-`rpc` pool; `target_fraction = target/rpc` for the
rest) — trimmed to `library(tidyverse)` alone rather than monoculture's
~40-package import (Seurat, DESeq2, MAST, ComplexHeatmap, ... — unrelated
to this computation, just slows startup; not fixed in the monoculture
script itself since this is a new file, but worth porting that fix back
if `monocultures/seqtk/get_seqtk_proportions.R` is ever revisited). Writes
three things, all into `cocultures/seqtk/`: `depth_by_pool.csv` (the full
reads/rpc/cells/target/target_fraction table, for inspection),
`downsample_targets.txt` (same `pool,target_fraction` format monoculture's
`seqtk.sh` consumes, written in case the no-downsampling decision changes
later — pools already at the target are dropped, same as monoculture), and
`reads_per_cell_by_pool.png` (a bar chart of mean reads/cell per pool with
the target depth marked — the actual "how bad is it" view, which
monoculture's script never needed since it only ever fed the raw
fractions straight to `seqtk`). **Important caveat, noted in the script's
own comments**: this is necessarily **per-pool**, not per-batch, depth —
`qc_library_metrics.csv` is a whole-GEM-well metric, and `JSB134-3` holds
two different batches (`_2`/`_3`) sequenced together in one well, so their
depths can't be separated from this file. The other 3 pools each map to
exactly one batch, so per-pool depth doubles as per-batch depth for those.
**Run — decision made**: the user ran this and found the depth disparity
across coculture pools much less severe than monoculture's (which is why
monoculture needed `seqtk` at all), and has decided to **forgo
downsampling for this arm** on that basis — not a placeholder "for now"
deferral anymore, but a real decision backed by the actual numbers. CPU
CellBender (`cocultures/cellbender_scripts/run_all_cellbender_cpu.sh`) is
therefore running directly against the original (non-downsampled) Cell
Ranger outputs, same as monoculture's `iMG2`.

## Coculture arm, next stage: recovering the lost preprocessing pipeline

The user's own coculture preprocessing scripts (the `01_obj_creation.R`-
through-`0N`-equivalent stages, including a microglia-vs-neuron broad
annotation step that monoculture has no counterpart for) were lost, but a
**finished, fully cleaned/processed Seurat object from that lost pipeline
still exists** on the HPC. Plan: reconstruct the pipeline's actual
parameters from that object rather than from memory, since Seurat logs
every function call that touches an object — method + every argument —
in `obj@commands`, which amounts to a parameter-level transcript of what
was actually run (`NormalizeData`, `FindVariableFeatures`, `RunPCA`,
`IntegrateLayers`, `FindClusters`, etc.). Combined with the object's
`meta.data` schema/value tables and `VariableFeatures`, that's enough to
rebuild the pipeline without transferring the (likely multi-GB) object
itself anywhere — this session has no access to the HPC filesystem and no
practical way to receive a file that size.

`cocultures/r_scripts/inspect_reference_object.R` — a one-off diagnostic,
deliberately **not** a numbered pipeline stage — does this extraction:
basic structure (assays/layers, reductions, graphs, `misc` slot names),
full `meta.data` schema plus a value table/summary per column (this is
where the microglia-vs-neuron annotation column's actual name and
levels will show up), `VariableFeatures`, and the full `obj@commands`
dump. Writes `object_diagnostic.txt`. **Not yet run** — waiting on the
user to point `OBJ_PATH` at the real object and run it on Quest, then
share the output back.

**How the microglia-vs-neuron annotation was done — confirmed by the
user**: clustering + manual marker-based labeling, the same by-eye
pattern `06_qc2.R` already uses for monoculture's non-microglia-cluster
removal, just repurposed here to *label* clusters as microglia/neuron
rather than *remove* them. Practical consequence for reading the
diagnostic output: `obj@commands` will capture the `FindClusters`
call(s) (algorithm, resolution, graph name) but **not** the labeling step
itself, since a manual `obj$cell_type <- ...`-style assignment after
eyeballing marker plots doesn't log to `obj@commands` the way a Seurat
function call does.

**`inspect_reference_object.R` run — results, as of this writing**:
`object_diagnostic.txt` came back for a `146472`-cell x `18129`-gene
object (gene count consistent with the probe-set panel, same as
monoculture). Findings:
- **Reductions/graphs present**: `pca`, `umap` (off `pca`), `harmony`,
  `HarmonyUMAP` (off `harmony`) — same pre/post-integration UMAP pattern
  as monoculture's `05_integration_harmony.R`. `RNA_nn`/`RNA_snn` graphs
  also present (built off `harmony`).
- **`obj@commands` parameters recovered** (directly reusable):
  `NormalizeData` (defaults: `LogNormalize`, `scale.factor = 10000`);
  `FindVariableFeatures` (`selection.method = "vst"`, `nfeatures = 2000`,
  all other args default); `ScaleData` (defaults); `RunPCA` (`npcs = 50`,
  `seed.use = 42`, otherwise default); `FindClusters` (`algorithm = 1`
  [original Louvain, **not** monoculture's `algorithm = 4`/Leiden],
  `resolution = 0.5` [**not** monoculture's `06_qc2.R`'s `resolution =
  2`], `graph.name = "RNA_snn"`, producing 27 clusters, `0`-`26`) — this
  is the clustering the manual microglia/neuron labeling would have been
  done on top of.
- **`obj@commands` gaps, both genuinely informative**: (1) no
  `IntegrateLayers`/`RunHarmony` entry at all — the `harmony` reduction
  exists but whatever produced it isn't logged, so the exact
  invocation (Seurat's `IntegrateLayers(method = HarmonyIntegration)`,
  which monoculture uses, vs. calling the standalone `harmony::RunHarmony()`
  directly, which doesn't log to `obj@commands`) can't be confirmed from
  this object alone. (2) `FindNeighbors(reduction = "pca")` was logged
  with `dims = 1:10`, but the immediately-following `RunUMAP(reduction =
  "pca")` was logged with `dims = 1:30` — and both the harmony-side
  `FindNeighbors`/`RunUMAP` calls used `dims = 1:30` throughout. So the
  lost pipeline's actual PC-count usage was **inconsistent** (10 for one
  call, 30 for three others), unlike monoculture's single established
  `dims = 1:10` used everywhere. **Resolved**: the rebuild standardizes on
  `dims = 1:10` throughout, matching monoculture, not the lost pipeline's
  literal inconsistency.
- **No microglia/neuron label column exists anywhere in this object's
  28-column `meta.data`** — nor any `DF.unadj`/`DF.adj` doublet columns
  (confirmed by the user: DoubletFinder was never run for coculture, the
  same early omission monoculture itself made before `03_doubletfinder.R`
  was added — to be built in from the start this time, not reverse-
  engineered from this object). Given both are missing, this object looks
  like it may be an **intermediate checkpoint** (roughly monoculture's
  `05_integration_harmony.R`-equivalent stage: normalized, scaled, PCA'd,
  Harmony-integrated, clustered) rather than the fully-finished,
  annotated object the user described. **Resolved**: confirmed by the
  user this is the right object to build from — the microglia/neuron
  split was done separately (e.g. in a notebook) and never written back,
  so the rebuild does fresh marker-based labeling on this object's own
  `RNA_snn_res.0.5` clusters rather than looking for a different object.
- **QC/filtering columns are a completely different scheme than
  monoculture's**, and **not reconstructible from `obj@commands`** since
  whatever produced them wasn't a logged Seurat function (most likely
  pre-filtering applied to barcodes before `CreateSeuratObject()`, or a
  bare `meta.data` column assignment): per-metric `<metric>_High`/
  `<metric>_Low` (fixed value cutoffs), `<metric>_RankCutoff_High`/
  `<metric>_RankCutoff_Low` (rank-based cutoffs — cells ranked by that
  metric, a rank threshold applied), and `<metric>_Outlier` (the per-cell
  verdict, as `"Not Outlier"` or `"<High|Low> Outlier (<n>)"`) for each of
  `nCount_RNA`, `nFeature_RNA`, `percent.mt`, plus an overall `Outlier`
  master column — a `rank`/`index` pair of columns also exists, consistent
  with an external rank-based filtering step. All object-wide `Outlier`
  values are `"Not Outlier"` (object is already filtered to survivors),
  but the per-metric `*_Outlier` sub-flags still show real flagged values
  for `percent.mt` specifically (same "keep the sub-flags even though the
  master flag already resolved to kept" pattern as monoculture's `02_qc.R`
  `mito_discard`/`umi_discard`/`gene_discard` + `discard`). **Monoculture's
  established MAD-threshold method (`02_qc.R`) is NOT what produced
  these** — this is a different, rank-involving scheme whose exact rule
  can't be confidently reverse-engineered from column values alone.
  **Resolved**: the rebuild uses monoculture's established MAD-threshold
  method instead — this object's rank-based scheme is abandoned, not
  replicated. Also notable: only `percent.mt` (no
  `percent.mt`-style ribo counterpart) was computed, unlike monoculture's
  `scCustomize::Add_Mito_Ribo()` (which produces both `percent_mito` and
  `percent_ribo`) — the naming (`percent.mt`, dot not underscore) also
  suggests base Seurat's `PercentageFeatureSet()` was used directly rather
  than `scCustomize`.
- **The real genotype-pair + batch-replicate structure isn't decomposed
  into separate columns** — `meta.data$batch` here actually holds the
  **pool/GEM-well identity** (`JSB127`/`JSB134-1`/`JSB134-2`/`JSB134-3`,
  i.e. monoculture's `pool_dir` concept, confusingly reusing the name
  `batch`) — it is **not** the genotype-pair-replicate 1/2/3 batch axis
  discussed earlier in this file. That real batch/replicate number, along
  with the neuron and microglia genotypes, only exists bundled together
  as one string in `sample_ID` (e.g. `"LRRK2_LRRK2_3"`) and
  `sample_barcode` (`sample_ID` + cell barcode concatenated, likely used
  for some join). **The lost pipeline apparently never parsed `sample_ID`
  into separate `neuron_genotype`/`microglia_genotype`/`batch` columns the
  way monoculture's `01_obj_creation.R` parses `genotype`/`treatment`/
  `batch` from `samples.csv`** — doing that parsing is a clear, concrete
  improvement to make in the rebuild (not ambiguous, no need to ask): the
  "Coculture arm" section above already documents exactly how to split
  `sample_ID` into its two genotype tokens and trailing batch number.

## Metadata schema, stage by stage

Every column below lives in `obj@meta.data` (or the standalone
`metadata.rds` that stands in for it). Grouping-variable bugs so far have
all come from conflating two of these that look similar but aren't
(`sample` vs. `condition`, raw vs. relabeled `treatment`, `batch` vs.
`pool_dir`) — check new plotting/grouping code against this table before
trusting it. (This table describes the **monoculture arm** only — the
coculture arm's schema is still undetermined, see above.)

| Column | Introduced in | Levels / type | What it actually distinguishes |
|---|---|---|---|
| `sample` | `01` (from `samples.csv$id`, `_redo` stripped) | character, 60 unique values, e.g. `GALC_asyn_1` | The **unified** biological sample — `iMG1` and `iMG1_redo` cells of the same sample share this value. **Ambiguous between the two captures on its own.** |
| `pool_dir` | `01` (from `cell_quantities.csv`) | character, 5 values: `iMG1`, `iMG1_redo`, `iMG2`, `iMG3`, `iMG4` | The real GEM well / physical capture a cell came from. Combine with `sample` to get back the distinct 75-capture identity (what `03`'s SLURM array and its output filenames key on). |
| `genotype` | `01` (raw from `samples.csv`), releveled as a factor in `02` | factor, 5 levels: `KOLF, GALC, GBA, GRN, LRRK2` (in that order) | Cell line / gene edit. |
| `treatment` | `01` (raw from `samples.csv`: `untreated`/`myelin`/`asyn`), **relabeled** in `02` | factor in `02` onward, levels/labels `Control`/`Myelin`/`Asyn` (raw values are lowercase, e.g. `"untreated"`) | Treatment condition. **`01`'s raw metadata.rds still has the lowercase raw strings** (`untreated`/`myelin`/`asyn`) as a plain character column, not the `02`-onward factor — any code reading `01`'s output directly (rather than `02`'s or later) must match against the raw values, not `Control`/`Myelin`/`Asyn`. |
| `batch` | `01` (raw from `samples.csv`, integer 1-4), releveled as a factor in `02` | factor in `02` onward, levels `"Batch 1".."Batch 4"` (explicit `levels=1:4` set to avoid a silent mislabel if the underlying integers aren't already sorted 1-4) | The sequencing/GEM batch number. **Identical for `iMG1` and `iMG1_redo`** (both `batch = 1`) — `batch` alone cannot distinguish the two captures; use `pool_dir` for that. |
| `condition` | `02` (new; not present in `01`'s raw metadata) | factor, 15 levels = `genotype x treatment`, explicit level order (5 genotype blocks x 3 treatments each) | Grouping/coloring variable for QC and PCA plots — deliberately **distinct from `sample`** (a naming collision in the user's own early draft would have overwritten the real per-sample identity). Each `(condition, batch)` pair maps to exactly one real sample. |
| `date` | `01` (raw from `samples.csv`) | character/date, as given in `samples.csv` | Collection date; carried through but not currently used as a grouping variable in `02`/`04`. |
| `nCount_RNA`, `nFeature_RNA` | Seurat default, at `01`'s `CreateSeuratObject()` | numeric | Raw UMI/gene counts. |
| `log10GenesPerUMI`, `percent_mito`, `percent_ribo` | `02` (`Add_Cell_Complexity()`/`Add_Mito_Ribo()`) | numeric | QC metrics computed after `01`'s object is loaded. |
| `mito_discard`, `umi_discard`, `gene_discard`, `discard` | `02` | logical | Per-cell QC flags from `02`'s per-`sample` MAD thresholds; `02`'s saved object (`data/02_qc/`) is already filtered to `discard == F`, but the flag columns are retained in its metadata for reference. |
| `DF.unadj`, `DF.adj` | `03` (per-capture), attached to the whole-cohort object in `04` | factor/character, `"Singlet"`/`"Doublet"` (DoubletFinder's own labels) | DoubletFinder's doublet call per cell, unadjusted and homotypic-proportion-adjusted. Computed **per individual capture** (75 separate DoubletFinder runs) — see the `03` caveat above about what this granularity can and can't catch. Only present in `03`'s own `metadata_persample/*.rds` files and in `04`'s reconstructed object; **not** in `01`'s or `02`'s metadata. |
| `qc_cluster` | `06` (renamed from `doublet_qc_cluster`, see entry 10 above) | integer, `resolution = 2` Louvain/Leiden cluster id, specific to `06`'s own clustering | Used only as the basis for `06`'s manual doublet-cluster and non-microglial-cluster removal — not a stable/reusable clustering, don't propagate downstream as "the" cluster identity. |

Practical implication for writing new grouping code: if you need "one row
per real sample" (e.g. for a per-sample summary table or threshold), group
by `sample` alone only when you've already restricted to post-`01` unified
data and don't care which capture a cell came from; group by
`pool_dir` + `sample` (or just use the distinct 75-capture manifest,
`cellbender_scripts/cell_quantities.csv`) whenever the two `iMG1` captures
must stay distinguishable. Use `condition` (not `sample`) for anything
genotype x treatment-level (15 groups), and `batch` (not `pool_dir`) for
anything about sequencing batch structure that doesn't care about the
`iMG1`/`iMG1_redo` split.

## Conventions established in this codebase

- **Self-resubmitting SLURM arrays**: several scripts (`seqtk/seqtk.sh`,
  `cellbender_scripts/run_all_cellbender*.sh`,
  `jobs/run_03_doubletfinder.sh`) follow the same pattern: check
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
  wrong SLURM account/partition inherited from a different project, a
  one-commit-lived undefined-variable typo in `06_qc2.R`
  (`non_microglial_clusters` vs. `non_microglia_clusters`, self-fixed one
  commit later), and the repo-wide restructure itself leaving every
  script's `setwd()`/hardcoded repo-root path stale (see the top-level
  structure section at the top of this file — now fixed).
- **Ask before guessing on consequential/ambiguous decisions**, especially
  ones with real scientific stakes (e.g. DoubletFinder granularity) or that
  contradict something just-established as correct. The user explicitly
  invites this ("ask if you're unsure") and has used
  `AskUserQuestion`-style clarifications productively — don't skip it to
  move faster when a wrong guess would mean redoing real compute (CellBender
  GPU-hours, DoubletFinder runs, etc.). This is exactly what happened with
  the coculture arm's sample-naming convention (see that section above) —
  asked rather than guessed at what `CN`/the numeric suffix/a possible
  treatment axis meant, and the user's answers corrected a wrong guess
  (token order was originally assumed backwards) before anything
  downstream was built on it.
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
- **hdWGCNA + BPCells incompatibility**: `wgcna.R`'s
  `SetupForWGCNA(gene_select = "fraction")` crashed with `Error in
  rownames(x)[i] : invalid subscript type 'S4'` the first time it ran
  against a BPCells-backed assay. Root cause: hdWGCNA's internal gene-
  selection code indexes the assay with plain base-R subsetting, which a
  lazy BPCells `IterableMatrix` doesn't support the way a real
  `dgCMatrix`/`matrix` does — unlike Seurat's own `NormalizeData()`/
  `ScaleData()`/`RunPCA()`/`IntegrateLayers()`, which `04_norm_pca.R`/
  `05_integration_harmony.R` already run directly against BPCells matrices
  without issue, hdWGCNA was never written with lazy/on-disk matrices in
  mind. Fixed by materializing both `counts` and `data` to `dgCMatrix`
  right after loading in `wgcna.R`, before any hdWGCNA call — full-cohort
  scale here (~40-60k cells) is a perfectly ordinary in-memory sparse
  matrix for a single-cell analysis (hdWGCNA's normal expected input); the
  BPCells backing was the atypical part, not the memory cost of
  materializing. If a similar `invalid subscript type 'S4'` (or similar
  BPCells-shaped) error shows up in `wgcna_stats.R`/`wgcna_viz.R` or any
  future stage that hands a BPCells matrix to a non-Seurat-native tool,
  this is the same class of fix.
- **hdWGCNA + Seurat v5 `GetAssayData(slot=...)` incompatibility (separate
  from the BPCells issue above)**: after the materialization fix, `wgcna.R`
  hit a second, unrelated crash inside `MetacellsByGroups()` itself —
  `Error: The 'slot' argument of GetAssayData() was deprecated in
  SeuratObject 5.0.0 and is now defunct`. This is a version mismatch
  between the installed `hdWGCNA` (which still calls the old
  `GetAssayData(slot = ...)` API internally) and the installed
  `SeuratObject` (>=5.0.0, where `slot=` is a hard error, not just a
  deprecation warning) — nothing in this repo's own scripts to fix, since
  the offending call lives inside hdWGCNA's source. Real fix: update
  hdWGCNA from its GitHub repo (it's not on CRAN), which has newer commits
  using `layer=` instead of `slot=`:
  `Rscript -e 'remotes::install_github("smorabit/hdWGCNA", upgrade =
  "never", force = TRUE)'` (install `remotes` first if missing). Verify
  with `packageVersion("hdWGCNA")` after. **This is likely not a
  single-call-site problem** — other hdWGCNA functions later in `wgcna.R`
  (`SetMultiExpr`, `ModuleEigengenes`, `ConstructNetwork`, etc.) may call
  `GetAssayData(slot = ...)` the same way, so expect the identical error
  to potentially resurface further into the script until hdWGCNA is fully
  updated; if it recurs even after reinstalling, check `.libPaths()` for a
  stale personal-library copy shadowing the update rather than assuming
  the reinstall didn't include a fix. A monkey-patch shim
  (`assignInNamespace`-ing a `slot`-to-`layer`-translating
  `GetAssayData.Seurat`) was considered as a stopgap but not used — the
  user chose to fix the actual package version instead.

## Open items / things worth revisiting

- ~~Repo-wide stale `setwd()`/`REPO=` paths under `monocultures/` after the
  restructure~~ — fixed (mechanical `/monocultures` path-segment insertion
  across every affected `.R`/`.sh`/manifest file, committed directly to
  this branch without a PR per the user's explicit request). See the
  top-level structure section at the top of this file.
- Stale, byte-identical duplicate of `samples.csv` still at repo root
  (should be deleted — `monocultures/samples.csv` is the real one).
- `r_scripts/deseq2.R`'s `AggregateExpression()` call lost its explicit
  `layer = "counts"` argument in the Oct 28 2026 restructure commit
  (`61c2414`) — not yet verified whether `AggregateExpression()`'s default
  layer is actually `"counts"` (if it's `"data"` instead, this would
  silently pseudobulk normalized values rather than raw counts, which
  DESeq2 needs). Check/re-add before trusting `deseq2.R`'s next run.
- `r_scripts/07_clustering.R` hasn't been run yet with the materialization
  fix *and* the `300G` memory bump together (the OOM that drove the bump
  happened at `128G`, with the fix already in place) — verify
  `FindAllMarkers()` is actually fast and fits in `300G` on the next run.
  The variable-features-restriction and 2-task-array ideas were raised but
  not applied — implement if the user wants either after seeing how this
  run goes.
- `r_scripts/milo.R` hasn't been run yet — first pass, adapted from another
  project's script, now with both the genotype-vs-KOLF and genotype x
  treatment interaction analyses implemented (see the pipeline-stage entry
  above). Verify on first run: the `HARMONY`/`HARMONY_UMAP` reduction-name
  capitalization assumption after `as.SingleCellExperiment()`, and the
  `nhoodDistances()`/`nhoodGraph()` accessor names. The user plans to try
  analysis 1 alone first, interactively, before attempting both.
- `r_scripts/deseq2_gsea.R` has not been reviewed for correctness the way
  `deseq2.R` was — only its paths were repointed for the `results/` reorg.
  Worth the same scrutiny (it reads `deseq2.R`'s per-`condition`
  `lfc_shrunk_*.csv` files now, so it inherits that design, but hasn't
  been run or checked itself).
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
- `06_qc2.R`'s specific removed-cluster lists (`doublet_clusters <- c(31,
  35, 36)`, `non_microglia_clusters <- c(26, 33)`) are tied to one exact
  upstream clustering — re-eyeball and update if `05_integration_harmony.R`'s
  output changes.
- The **coculture arm's** `samples.csv`-equivalent metadata file still
  needs to be built (naming convention is now confirmed — see the
  "Coculture arm" section above — just not yet turned into an actual CSV/
  `01_obj_creation.R`-equivalent). Remember `JSB134-3`'s mixed-batch
  structure when building it: `batch` must come from the `sample_id`
  suffix, not `pool_dir`.
- `cocultures/cellbender_scripts/make_cellbender_params.R` and
  `run_all_cellbender_cpu.sh` are now **running** on the HPC (not yet
  confirmed complete/successful) — verify once finished (in particular,
  that `metrics_summary.csv` existed per-sample under this Cell Ranger
  version's output the same way it does for monoculture, and that
  `sample_raw_feature_bc_matrix.h5` was named identically).
- `inspect_reference_object.R` has been run — see the dedicated section
  above for the full results. **Three open decisions, asked about
  directly rather than guessed — now resolved**: (1) confirmed this
  *is* the right object to build from; it has no microglia/neuron label
  because that split was done separately (e.g. in a notebook) and never
  written back — the rebuild will do fresh marker-based annotation on
  this object's own clusters (`RNA_snn_res.0.5`, 27 clusters), same as
  the lost pipeline would have. (2) **QC will use monoculture's
  established MAD-threshold method** (`02_qc.R`'s approach), not the
  rank-involving scheme found in the surviving object — the surviving
  object's QC scheme is abandoned for this rebuild, not replicated. (3)
  **PC dims: standardize on `1:10` throughout**, matching monoculture's
  single established convention, not the lost pipeline's inconsistent
  10-vs-30 usage.
- **A fourth design fork, surfaced on closer inspection of the object's
  layer names — asked about, now resolved**: the lost pipeline's Harmony
  integration split the RNA assay by **pool** before integrating
  (`counts.JSB127`, `counts.JSB134-1`, `counts.JSB134-2`,
  `counts.JSB134-3` — exactly 4 layers, one per pool, matching the 4
  levels of the object's misleadingly-named `batch` column). This is the
  *same* design fork monoculture's `05_integration_harmony.R` explicitly
  resolved the other way via `AskUserQuestion` — monoculture splits by
  the unified `sample` column (60 levels), not `pool_dir` (5 levels).
  **Resolved for coculture**: split by the **full 30-level genotype-pair
  x batch-replicate sample identity** (e.g. `KOLF_KOLF_2` as its own
  layer) — the finest-grained of the three options considered (pool: 4
  levels; unified genotype-pair condition: 10 levels; full sample: 30
  levels), and notably *not* what the lost pipeline did (pool-level).
  This corrects for technical variation at the same per-sample
  granularity monoculture's own split achieves, while additionally
  separating `JSB134-3`'s two mixed batches into their own layers —
  something a pool-level split structurally couldn't do here, since pool
  and real sample identity don't align 1:1 for 3 of the 10 conditions.
- Decided without asking (low-stakes, consistent with established
  convention, easily revisited): add `percent_ribo` alongside
  `percent_mito` via `scCustomize::Add_Mito_Ribo()` for consistency with
  monoculture, even though the lost pipeline only ever computed
  `percent.mt`; use `IntegrateLayers(method = HarmonyIntegration)` for the
  rebuild's Harmony step, matching monoculture, since the lost pipeline's
  actual invocation wasn't logged either way; parse `sample_ID` into
  separate `neuron_genotype`/`microglia_genotype`/`batch` columns (see the
  dedicated section above), since the lost pipeline apparently never did
  this and it's a clear improvement, not an ambiguous one.
- DoubletFinder will be built into the coculture rebuild from the start
  (confirmed by the user) — the surviving reference object never had it,
  the same omission monoculture itself initially made before
  `03_doubletfinder.R` was added.
- ~~Root-level `jobs/` duplicates several `.sh` job scripts~~ — resolved by
  the `r_scripts/`/`jobs/`/`data/`/`results/` reorg (see the directory-
  structure section near the top and the note at the end of the
  pipeline-stages section above); `jobs/` is now the single canonical `.sh`
  location, nothing left to duplicate.
- `r_scripts/deseq2.R` and `r_scripts/wgcna*.R` now both read from `data/06_qc2/` (see
  the pipeline-stages entries above) — `06_qc2.R`'s output is the
  canonical upstream source for both as of this writing. If `06_qc2.R` is
  ever rerun with a different set of doublet-flagged clusters (or its own
  logic changes), both downstream consumers pick that up automatically
  since neither caches a copy of `06`'s cell set.
- `05_integration_harmony.R` splits by the unified `sample` column (60
  levels) for Harmony integration, not `pool_dir` — an explicit, informed
  choice the user made after being shown the tradeoff (see the `05` pipeline
  entry above). Worth remembering if `iMG1`/`iMG1_redo` batch structure
  ever looks under-corrected in the integrated UMAP.
- **Standing instruction from the user (2026 session)**: no need to
  proactively watch/monitor PRs for CI status or reviews going forward —
  the user communicates directly in chat instead. Don't spend turns
  polling PR state unless asked.
</content>
