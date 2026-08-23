#!/bin/bash
# Builds seqtk/downsample_manifest.txt: one row per FASTQ file that needs
# downsampling, as `source_path,fraction,dest_path`.
#
# Run this on a login node before submitting seqtk.sh (it only globs files,
# no compute needed), and re-run it any time downsample_targets.txt or the
# source FASTQs change.
#
# Which files belong to a pool is derived straight from that pool's
# cellranger_configs/<pool>_config.csv [libraries] section (fastq_id +
# fastqs dir) rather than hardcoded here again, so there is exactly one
# place that records which files make up a pool. See seqtk/README.md for
# why per-file fractions don't require concatenating FASTQs first.

set -euo pipefail

REPO="/projects/b1169/boles/img_scfrp"
MANIFEST="$REPO/seqtk/downsample_manifest.txt"
TARGETS="$REPO/seqtk/downsample_targets.txt"

> "$MANIFEST"

while IFS=',' read -r pool fraction; do
    [ -z "$pool" ] && continue
    config="$REPO/cellranger_configs/${pool}_config.csv"
    if [ ! -f "$config" ]; then
        echo "WARNING: no config for pool '$pool' at $config, skipping" >&2
        continue
    fi

    # Pull fastq_id,fastqs pairs out of the [libraries] table (there may be
    # more than one row per pool, e.g. iMG3 spans two source directories)
    awk -F',' '
        /^\[libraries\]/ { inlib=1; next }
        /^\[/            { inlib=0 }
        inlib && NF>=2 && $1 != "fastq_id" { print $1","$2 }
    ' "$config" | while IFS=',' read -r fastq_id fastqs_dir; do
        for read_num in 1 2; do
            for src in "$fastqs_dir/${fastq_id}"_S*_L*_R${read_num}_001.fastq.gz; do
                [ -e "$src" ] || continue
                dest_subdir=$(basename "$fastqs_dir")
                dest="$REPO/downsampled_fastq/${dest_subdir}/$(basename "$src")"
                echo "${src},${fraction},${dest}" >> "$MANIFEST"
            done
        done
    done
done < "$TARGETS"

n=$(wc -l < "$MANIFEST")
echo "Wrote $n files to $MANIFEST"
echo "Submit with: sbatch --array=1-${n} seqtk/seqtk.sh"
