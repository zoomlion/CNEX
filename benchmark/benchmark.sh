#!/bin/bash
# CNEX Benchmark — unified multi-species pipeline script
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
PATH="$HOME/Software/bin:$PATH"
CNEX_BIN="$HERE/../bin"

SPECIES=(Homo_sapiens Mus_musculus Gallus_gallus Scyliorhinus_canicula Callorhinchus_milii)
QUERYLIB_FRACS=(0.10 0.20 0.40 0.60 0.80)

usage() { echo "Usage: $0 {all|extract|dbs|search|m8extract|sensitivity|querylib|querylib-extract|querylib-search|querylib-analyze|querylib-plot|analyze|clean} [species...]" >&2; exit 1; }
[ $# -lt 1 ] && usage
CMD=$1; shift
TARGETS=("$@")
[ ${#TARGETS[@]} -eq 0 ] && TARGETS=("${SPECIES[@]}")

case $CMD in
# ─── Step 1: Extract per-target queries/truth/blocks ───
extract)
    for T in "${TARGETS[@]}"; do
        echo "=== Extract (target=$T, seed=42, 10k blocks) ==="
        python3 "$HERE/01_extract.py" --target "$T"
    done
    ;;

# ─── Step 2: Build databases per species ───
dbs)
    for T in "${TARGETS[@]}"; do
        GENOME="$HERE/data/genomes/$T/$T.fa"
        mkdir -p "$HERE/dbs/$T/mmseqs" "$HERE/dbs/$T/blast" "$HERE/dbs/$T/cnex" "$HERE/stats/$T"

        echo "=== [$T] MMseqs2 createdb ==="
        [ -f "$HERE/dbs/$T/mmseqs/$T" ] || \
        /usr/bin/time -v -o "$HERE/stats/$T/time_mmseqs_createdb.txt" \
            mmseqs createdb "$GENOME" "$HERE/dbs/$T/mmseqs/$T"

        echo "=== [$T] makeblastdb ==="
        [ -f "$HERE/dbs/$T/blast/$T.nhr" ] || \
        /usr/bin/time -v -o "$HERE/stats/$T/time_makeblastdb.txt" \
            makeblastdb -in "$GENOME" -dbtype nucl -out "$HERE/dbs/$T/blast/$T"

        echo "=== [$T] CNEX mertable ==="
        [ -f "$HERE/dbs/$T/cnex/mers_table.tsv" ] || \
        /usr/bin/time -v -o "$HERE/stats/$T/time_cnex_mertable.txt" \
            "$CNEX_BIN/mertable" "$HERE/data/blocks/$T.fa" -k 13 -c 4 -o "$HERE/dbs/$T/cnex/mers_table.tsv"
    done
    ;;

# ─── Step 3: Run searches per species ───
search)
    for T in "${TARGETS[@]}"; do
        GENOME="$HERE/data/genomes/$T/$T.fa"
        mkdir -p "$HERE/results/$T/mmseqs" "$HERE/results/$T/blastn" "$HERE/results/$T/cnex" "$HERE/stats/$T"

        echo "=== [$T] Preprocess queries ==="
        python3 -c "
from pathlib import Path
inp=Path('$HERE/data/queries/$T.fa'); out=Path('$HERE/data/queries/${T}_clean.fa')
lines=inp.read_text().split('\n'); c=[]
for l in lines:
    if l.startswith('>'): c.append(l)
    elif l:
        s=l.replace('-','').upper()
        if s: c.append(s)
out.write_text('\n'.join(c)+'\n')
print(f'  $T clean queries: {(len(c)-1)//2} sequences')
"

        echo "=== [$T] CNEX validate ==="
        /usr/bin/time -v -o "$HERE/stats/$T/time_cnex_validate.txt" \
            "$CNEX_BIN/validate" "$GENOME" --mers "$HERE/dbs/$T/cnex/mers_table.tsv" \
            --type genome -t 4 \
            --rescue-panels "$HERE/data/blocks/$T.fa" --rescue-tau 0.90 --t3 "1.5,5,25" \
            -o "$HERE/results/$T/cnex/validate_out"

        echo "=== [$T] CNEX assemble ==="
        /usr/bin/time -v -o "$HERE/stats/$T/time_cnex_assemble.txt" \
            "$CNEX_BIN/assemble" "$HERE/results/$T/cnex/validate_out" \
            --mers "$HERE/dbs/$T/cnex/mers_table.tsv" \
            --panels "$HERE/data/blocks/$T.fa" --trim \
            -o "$HERE/results/$T/cnex/contigs.fa"

        echo "=== [$T] MMseqs2 search (s=5.7) ==="
        rm -rf "$HERE/results/$T/mmseqs/tmp"
        /usr/bin/time -v -o "$HERE/stats/$T/time_mmseqs_search.txt" \
            mmseqs easy-search "$HERE/data/queries/${T}_clean.fa" "$HERE/dbs/$T/mmseqs/$T" \
            "$HERE/results/$T/mmseqs/out_s5.7_noid.m8" "$HERE/results/$T/mmseqs/tmp" \
            --min-seq-id 0.7 -s 5.7 --threads 16 --search-type 3
        rm -rf "$HERE/results/$T/mmseqs/tmp"

        echo "=== [$T] BLASTN search (w15) ==="
        /usr/bin/time -v -o "$HERE/stats/$T/time_blastn_w15.txt" \
            blastn -db "$HERE/dbs/$T/blast/$T" -query "$HERE/data/queries/${T}_clean.fa" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -perc_identity 70 -word_size 15 -num_threads 16 -out "$HERE/results/$T/blastn/out_w15.txt"
    done
    ;;

# ─── Step 4: Extract best hits per block (per species) ───
m8extract)
    python3 "$HERE/02_m8extract.py" "${TARGETS[@]}"
    ;;

# ─── Step 5: MMseqs2 sensitivity sweep (human only) ───
sensitivity)
    T=Homo_sapiens
    echo "=== MMseqs2 sensitivity sweep (human only) ==="
    mkdir -p "$HERE/stats/$T"
    LOG="$HERE/stats/$T/sensitivity_summary.txt"; > "$LOG"
    for S in 1.0 4.0 7.5; do
        [ -f "$HERE/results/$T/mmseqs/out_s${S}.m8" ] && echo "  s=${S} already done, skip" && continue
        echo "  s=${S} ..."
        rm -rf "$HERE/results/$T/mmseqs/tmp_s${S}"
        /usr/bin/time -v -o "$HERE/stats/$T/time_mmseqs_s${S}.txt" \
            mmseqs easy-search "$HERE/data/queries/${T}_clean.fa" "$HERE/dbs/$T/mmseqs/$T" \
            "$HERE/results/$T/mmseqs/out_s${S}.m8" "$HERE/results/$T/mmseqs/tmp_s${S}" \
            --min-seq-id 0.7 -s ${S} --threads 16 --search-type 3
        wall=$(grep "Elapsed" "$HERE/stats/$T/time_mmseqs_s${S}.txt" | awk '{print $NF}')
        hits=$(wc -l < "$HERE/results/$T/mmseqs/out_s${S}.m8")
        echo "s=${S}: wall=${wall}, hits=${hits}"
        echo "s=${S}: wall=${wall}, hits=${hits}" >> "$LOG"
        rm -rf "$HERE/results/$T/mmseqs/tmp_s${S}"
    done
    # append s=5.7 (already run)
    if [ -f "$HERE/results/$T/mmseqs/out_s5.7_noid.m8" ]; then
        wall=$(grep "Elapsed" "$HERE/stats/$T/time_mmseqs_search.txt" | awk '{print $NF}')
        hits=$(wc -l < "$HERE/results/$T/mmseqs/out_s5.7_noid.m8")
        echo "s=5.7: wall=${wall}, hits=${hits}"
        echo "s=5.7: wall=${wall}, hits=${hits}" >> "$LOG"
    fi
    echo "Done. See $LOG"
    ;;

# ─── Step 6: CNEX query-lib subsampling sweep (block-level, nested) ───
querylib-extract)
    for T in "${TARGETS[@]}"; do
        for F in "${QUERYLIB_FRACS[@]}"; do
            echo "=== Extract ($T, block-frac=$F) ==="
            python3 "$HERE/01_extract.py" --target "$T" --block-frac "$F"
        done
    done
    ;;

querylib-search)
    for T in "${TARGETS[@]}"; do
        GENOME="$HERE/data/genomes/$T/$T.fa"
        for F in "${QUERYLIB_FRACS[@]}"; do
            TAG="bf${F}"
            QDIR="$HERE/data/querylib/$T/$TAG"
            MERS="$HERE/dbs/$T/querylib/$TAG/mers_table.tsv"
            OUT="$HERE/results/querylib/$T/$TAG/cnex"
            mkdir -p "$HERE/stats/querylib/$T/$TAG" "$(dirname "$MERS")"
            [ -d "$OUT" ] && rm -rf "$OUT"
            mkdir -p "$OUT"

            echo "=== [$T $TAG] CNEX mertable ==="
            /usr/bin/time -v -o "$HERE/stats/querylib/$T/$TAG/time_cnex_mertable.txt" \
                "$CNEX_BIN/mertable" "$QDIR/blocks.fa" -k 13 -c 4 -o "$MERS"

            echo "=== [$T $TAG] CNEX validate ==="
            /usr/bin/time -v -o "$HERE/stats/querylib/$T/$TAG/time_cnex_validate.txt" \
                "$CNEX_BIN/validate" "$GENOME" --mers "$MERS" \
                --type genome -t 16 -o "$OUT/validate_out"

            echo "=== [$T $TAG] CNEX assemble ==="
            /usr/bin/time -v -o "$HERE/stats/querylib/$T/$TAG/time_cnex_assemble.txt" \
                "$CNEX_BIN/assemble" "$OUT/validate_out" \
                --mers "$MERS" --trim -o "$OUT/contigs.fa"
        done
    done
    ;;

querylib-analyze)
    python3 "$HERE/05_querylib_analyze.py" "${TARGETS[@]}"
    ;;

querylib-plot)
    python3 "$HERE/06_querylib_plot.py" "${TARGETS[@]}"
    ;;

querylib)
    $0 querylib-extract "${TARGETS[@]}"
    $0 querylib-search "${TARGETS[@]}"
    $0 querylib-analyze "${TARGETS[@]}"
    $0 querylib-plot "${TARGETS[@]}"
    ;;

# ─── Step 7: Analyze (multi-species) ───
analyze)
    python3 "$HERE/03_analyze.py" "${TARGETS[@]}"
    ;;

# ─── Full pipeline ───
all)
    $0 extract "${TARGETS[@]}"
    $0 dbs "${TARGETS[@]}"
    $0 search "${TARGETS[@]}"
    $0 sensitivity
    $0 m8extract "${TARGETS[@]}"
    $0 analyze "${TARGETS[@]}"
    ;;

clean)
    rm -rf "$HERE/data/queries" "$HERE/data/truth" "$HERE/data/blocks" \
           "$HERE/dbs" "$HERE/results" "$HERE/stats" "$HERE/figures" \
           "$HERE/data/querylib"
    echo "Cleaned (backup_human_* kept)."
    ;;

*) usage;;
esac
