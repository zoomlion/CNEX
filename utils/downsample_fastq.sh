#!/bin/bash
# Parallel downsampling of paired-end FASTQ to multiple depths.
# Usage: bash utils/downsample_fastq.sh <R1.fq.gz> <R2.fq.gz> <out_dir>
set -e

SEQTK=/home/zhengjiangmin/Software/seqtk/seqtk
PIGZ=/home/zhengjiangmin/Software/mercury/src/pigz/pigz
R1=$1; R2=$2; OUT=$3
mkdir -p "$OUT"

echo "=== Parallel downsampling ===" | tee -a "$OUT/sampling.log"

for spec in "0.25x 2500000" "0.5x 5000000" "1x 10000000" "2x 20000000" "3x 30000000" "5x 51662516"; do
    name=$(echo $spec | cut -d' ' -f1)
    pairs=$(echo $spec | cut -d' ' -f2)
    dir="$OUT/$name"; mkdir -p "$dir"
    (
        echo "[$(date +%H:%M:%S)] Sampling $name ($pairs pairs)..."
        $SEQTK sample -s42 "$R1" $pairs 2>/dev/null | $PIGZ -p 2 > "$dir/reads_R1.fq.gz"
        $SEQTK sample -s42 "$R2" $pairs 2>/dev/null | $PIGZ -p 2 > "$dir/reads_R2.fq.gz"
        echo "[$(date +%H:%M:%S)] $name done"
    ) >> "$OUT/sampling.log" 2>&1 &
done

wait
echo "=== All done ===" | tee -a "$OUT/sampling.log"

for d in "$OUT"/*/; do
    n=$(echo $(zcat "${d}reads_R1.fq.gz" | wc -l) / 4 | bc 2>/dev/null || echo "?")
    echo "  $(basename $d): $n reads"
done
