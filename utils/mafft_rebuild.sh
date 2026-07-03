#!/bin/bash
# MAFFT re-alignment pipeline for CNEX mers table rebuild
set -e
HERE=$(cd "$(dirname "$0")"/.. && pwd)
N_THREADS=10
K=13

echo "=== Phase 1: Extract unaligned sequences ==="
python3 "$HERE/utils/mafft_extract.py" extract

echo "=== Phase 2: MAFFT parallel alignment ==="
T0=$(date +%s)
ls /tmp/mafft_blocks/*.fa | parallel -j $N_THREADS --bar \
  "mafft --auto {} 2>/dev/null | awk '/^>/{if(s!=\"\")print s; print; s=\"\"; next} {s=s\"\"\$0} END{if(s!=\"\")print s}' > /tmp/mafft_aligned/{/}"
echo "MAFFT done in $(( ($(date +%s) - T0) / 60 ))m"

echo "=== Phase 3: Merge aligned blocks ==="
python3 "$HERE/utils/mafft_extract.py" merge

echo "=== Phase 4: Rebuild mers table ==="
mv "$HERE/benchmark/dbs/cnex/mers_table.tsv" "$HERE/benchmark/dbs/cnex/mers_table.tsv.bak2" 2>/dev/null || true
"$HERE/bin/mertable" "$HERE/benchmark/data/blocks_10k_mafft.fa" \
  -k $K -c 4 --min-entropy 1.2 \
  -o "$HERE/benchmark/dbs/cnex/mers_table.tsv"
echo "Done."
