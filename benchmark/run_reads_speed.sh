#!/bin/bash
# CNEX reads-mode speed re-run (head2head).
#
# Controlled timing of CNEX validate+assemble at -t 4 on the 0.25x-8x wgsim
# reads, gzipped end-to-end (same conditions as benchmark.sh head2head runs).
#
# CACHE=cold (default): evict R1/R2 from the page cache with
#   posix_fadvise(DONTNEED) right before each run (no root needed) and record
#   the fincore residency (should be 0B = cold). This is page-cache cold, not
#   device cold (RAID controller cache may still hold data).
# CACHE=warm: pre-read R1/R2 into the page cache (cat) before each run.
#
# Wall times are I/O- and contention-noisy on this FS; core-hours (user+sys) is
# the headline metric (see benchmark/reads_speed.py).
#
# Usage:
#   bash benchmark/run_reads_speed.sh                 # all depths, cold, 1 rep
#   CACHE=warm REPS=3 bash benchmark/run_reads_speed.sh
#   CACHE=cold REPS=1 bash benchmark/run_reads_speed.sh 5x 8x
#
# Env: BASE (~/temp/head2head), CNEX (repo root), CACHE (cold|warm),
#      REPS (default 1), THREADS (default 4)
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
BASE=${BASE:-$HOME/temp/head2head}
CNEX=${CNEX:-$HERE/..}
THREADS=${THREADS:-4}
REPS=${REPS:-1}
CACHE=${CACHE:-cold}
DEPTHS=${*:-"0.25x 0.5x 1x 2x 3x 4x 5x 6x 7x 8x"}
OUT="$BASE/stats/speed_${CACHE}"

case "$CACHE" in cold|warm) ;; *) echo "CACHE must be cold|warm" >&2; exit 1;; esac

mkdir -p "$OUT"
cd "$BASE"

evict_cache() {
    python3 -c 'import os, sys
for p in sys.argv[1:]:
    fd = os.open(p, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)' "$@"
}

for d in $DEPTHS; do
    R1="reads/$d/reads_R1.fq.gz"
    R2="reads/$d/reads_R2.fq.gz"
    [ -f "$R1" ] || { echo "ERROR: $R1 not found" >&2; exit 1; }
    WORK="cnex/speed_${CACHE}_$d"
    rm -rf "$WORK"
    mkdir -p "$WORK"
    mkdir -p "$BASE/cnex/$d"

    for r in $(seq 1 "$REPS"); do
        if [ "$CACHE" = cold ]; then
            evict_cache "$R1" "$R2"
            echo "[$(date +%H:%M:%S)] $d rep $r/$REPS: cache evicted"
            fincore --bytes "$R1" "$R2" | tee -a "$OUT/cache_state_${d}.txt"
        else
            echo "[$(date +%H:%M:%S)] $d rep $r/$REPS: warming cache"
            cat "$R1" "$R2" > /dev/null
        fi

        echo "[$(date +%H:%M:%S)] $d rep $r: validate (-t $THREADS)"
        /usr/bin/time -v -o "$OUT/cnex_validate_${d}_r${r}.txt" \
            "$CNEX/validate" "$R1" "$R2" \
            --mers data/mers_table.tsv --type fastq -t "$THREADS" \
            --pigz "$BASE/bin/pigz" \
            --rescue-panels data/Homo_sapiens.fa --rescue-tau 0.90 \
            --t3 "1.5,5,25" --vote-frac 0.05 \
            -o "$WORK/validate_out" > "$OUT/cnex_validate_${d}_r${r}.log" 2>&1

        echo "[$(date +%H:%M:%S)] $d rep $r: assemble"
        /usr/bin/time -v -o "$OUT/cnex_assemble_${d}_r${r}.txt" \
            "$CNEX/assemble" "$WORK/validate_out" \
            --mers data/mers_table.tsv --panels data/Homo_sapiens.fa --trim \
            -o "$WORK/contigs.fa" > "$OUT/cnex_assemble_${d}_r${r}.log" 2>&1

        cp "$WORK/contigs.fa" "$BASE/cnex/$d/contigs.fa"
        md5sum "$WORK/contigs.fa" >> "$OUT/contigs_md5_${d}.txt"
        n=$(awk '/^>/{n++} END{print n+0}' "$WORK/contigs.fa")
        echo "[$(date +%H:%M:%S)] $d rep $r done: contigs=$n"
    done

    rm -rf "$WORK/validate_out"
done

echo "[$(date +%H:%M:%S)] ALL DONE"
