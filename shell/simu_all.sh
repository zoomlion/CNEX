#!/bin/bash
# simu_all.sh - Run 5X read simulation for all NEW gnathostome species
# Uses xargs -P for parallel execution
# Usage: ./shell/simu_all.sh [-j N]
set -euo pipefail

MAX_JOBS=3
while getopts "j:" opt; do
    case $opt in
        j) MAX_JOBS=$OPTARG ;;
    esac
done

MERCURY_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$MERCURY_DIR"

echo "=== CNEX 5X Read Simulation ==="
echo "Parallel jobs: $MAX_JOBS"
echo "Start: $(date)"
echo ""

# Build command list for xargs
# Format: each line is: sp genome_fa n_pairs
CMDS_FILE=$(mktemp /tmp/simu_cmds_XXXX.txt)
trap "rm -f $CMDS_FILE" EXIT

# Read the Makefile and find new species
# Parse GENOME_* and N_* from species.mk
declare -A GENOME_MAP
declare -A N_MAP

eval $(grep -E '^GENOME_' shell/species.mk | sed 's/ :=/=/;s/^/GENOME_MAP[/;s/=/]=/')
eval $(grep -E '^N_' shell/species.mk | sed 's/ :=/=/;s/^/N_MAP[/;s/=/]=/')

NEW_SPECIES=$(grep '^NEW_SPECIES' shell/species.mk | sed 's/^NEW_SPECIES := //')

TOTAL=0
for sp in $NEW_SPECIES; do
    genome="${GENOME_MAP[$sp]:-}"
    n="${N_MAP[$sp]:-}"
    if [ -z "$genome" ] || [ -z "$n" ]; then
        echo "WARNING: $sp missing genome or N, skipping" >&2
        continue
    fi
    # Skip if already exists
    if [ -f "reads/${sp}_r1.fq.gz" ] && [ "$(stat -c%s "reads/${sp}_r1.fq.gz")" -gt 1000000 ]; then
        echo "[skip] $sp — already exists"
        continue
    fi
    echo "$sp	$genome	$n" >> "$CMDS_FILE"
    TOTAL=$((TOTAL + 1))
done

echo ""
echo "Species to simulate: $TOTAL"
if [ "$TOTAL" -eq 0 ]; then
    echo "All done!"
    exit 0
fi
echo ""

# Run with xargs -P for parallelism
# xargs reads from CMDS_FILE and runs simu_one.sh for each line
< "$CMDS_FILE" xargs -P "$MAX_JOBS" -I{} bash -c '
    IFS=$'\''\t'\'' read -r sp genome n <<< "$1"
    echo "[$(date +%H:%M)] Starting $sp ($n reads)"
    cd "'"$MERCURY_DIR"'" && bash shell/simu_one.sh "$sp" "$genome" "$n"
' _ {} \;

echo ""
echo "=== Simulation complete ==="
echo "Finish: $(date)"
echo "Read files: $(ls reads/*_r1.fq.gz 2>/dev/null | wc -l)"
echo "Reads size: $(du -sh reads/ | cut -f1)"
echo "Disk free: $(df -h / | tail -1 | awk '{print $4}')"
