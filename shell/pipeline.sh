#!/bin/bash
# pipeline.sh - Full CNEX pipeline: 01 → 02 → 03 → 04
# Usage: bash shell/pipeline.sh
set -e

THREADS=8
WGSIM=~/Software/wgsim/wgsim
PIGZ=src/pigz/pigz
ASSEMBLER=src/02.validate
MERS_TABLE=mers_table.tsv
MSA=most-cons-cne.fa
RESULTS=results

# ─── 10 species ───
SPECIES=(Homo Mus Danio Loxodonta Gallus Xenopus Latimeria Callorhinchus Polypterus Rhincodon)

# ─── Genome paths & 3X PE pair counts ───
declare -A GENOME N_PAIRS
GENOME[Homo]="$HOME/Source/Genome/Homo_sapiens/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
N_PAIRS[Homo]=30997508
GENOME[Mus]="$HOME/Source/Genome/Mus_musculus/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
N_PAIRS[Mus]=27282225
GENOME[Danio]="$HOME/Source/Genome/Danio_rerio/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
N_PAIRS[Danio]=13734714
GENOME[Loxodonta]="$HOME/Source/Genome/Loxodonta_africana/Loxodonta_africana.loxAfr3.dna_sm.toplevel.fa"
N_PAIRS[Loxodonta]=31967609
GENOME[Gallus]="$HOME/Source/Genome/Gallus_gallus/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa"
N_PAIRS[Gallus]=5326827
GENOME[Xenopus]="$HOME/Source/Genome/Xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.dna_sm.toplevel.fa"
N_PAIRS[Xenopus]=7201992
GENOME[Latimeria]="$HOME/Source/Genome/Latimeria_chalumnae/Latimeria_chalumnae.LatCha1.dna_sm.toplevel.fa"
N_PAIRS[Latimeria]=14302959
GENOME[Callorhinchus]="$HOME/Source/Genome/Callorhinchus_milii/Callorhinchus_milii.Callorhinchus_milii-6.1.3.dna_sm.toplevel.fa"
N_PAIRS[Callorhinchus]=4872492
GENOME[Polypterus]="$HOME/Source/Genome/Polypterus_senegalus/Polypterus_senegalus.ASM1683550v1_genomic_GCF.fa"
N_PAIRS[Polypterus]=18369400
GENOME[Rhincodon]="$HOME/Source/Genome/Rhincodon_typus/Rhincodon_typus.ASM164234v2_genomic_GCF.fa"
N_PAIRS[Rhincodon]=14657997

echo "============================================"
echo "  CNEX Pipeline: 01 → simu → 02 → 03 → 04"
echo "============================================"

# ─── Step 1: 01.confi_mer ───
echo ""
echo "[Step 1] Building confident k-mer table ..."
if [ ! -f "$MERS_TABLE" ]; then
    python3 src/01.confi_mer.py "$MSA" -k 13 -c 4 -o "$MERS_TABLE"
    echo "  mers_table: $(wc -l < $MERS_TABLE) lines"
else
    echo "  $MERS_TABLE already exists, skipping"
fi

# ─── Step 2: Simulate reads (3X, gzipped) ───
echo ""
echo "[Step 2] Generating 3X simulated reads ..."
mkdir -p reads
for sp in "${SPECIES[@]}"; do
    if [ -f "reads/${sp}_r1.fq.gz" ] && [ $(stat -c%s "reads/${sp}_r1.fq.gz" 2>/dev/null || echo 0) -gt 1000000 ]; then
        echo "  [$sp] already exists, skipping"
        continue
    fi
    N=${N_PAIRS[$sp]}
    FA=${GENOME[$sp]}
    echo "  [$sp] ${N} PE pairs from $FA ..."
    mkfifo /tmp/fifo1_$$ /tmp/fifo2_$$
    $PIGZ < /tmp/fifo1_$$ > "reads/${sp}_r1.fq.gz" &
    $PIGZ < /tmp/fifo2_$$ > "reads/${sp}_r2.fq.gz" &
    $WGSIM -N $N -1 150 -2 150 -e 0.001 -r 0.001 -R 0.0001 "$FA" /tmp/fifo1_$$ /tmp/fifo2_$$
    wait
    rm -f /tmp/fifo1_$$ /tmp/fifo2_$$
    echo "  [$sp] done"
done
echo "  All reads generated"
ls -lhS reads/*.fq.gz

# ─── Step 3: 02.validate ───
echo ""
echo "[Step 3] Validating reads/genome (02.validate) ..."
rm -rf "$RESULTS"
mkdir -p "$RESULTS"
for sp in "${SPECIES[@]}"; do
    mkdir -p "$RESULTS/$sp"
done

# Run in batches of 5 to avoid overload
run_02() {
    local sp=$1
    $ASSEMBLER "reads/${sp}_r1.fq.gz" "reads/${sp}_r2.fq.gz" \
        --mers "$MERS_TABLE" -t $((THREADS/2)) --output_dir "$RESULTS/$sp" \
        --pigz src/pigz/pigz \
        --min-c 5 --vote-frac 0.1 --vote-ratio 3.0 \
        > "$RESULTS/$sp/02.log" 2>&1
    echo "  [$sp] done ($(wc -l < $RESULTS/$sp/Assemble.0.reads) lines)"
}

echo "  Batch 1/2 ..."
for sp in Homo Mus Danio Loxodonta Gallus; do
    run_02 "$sp" &
done; wait

echo "  Batch 2/2 ..."
for sp in Xenopus Latimeria Callorhinchus Polypterus Rhincodon; do
    run_02 "$sp" &
done; wait

# ─── Step 4: 03.debruijn ───
echo ""
echo "[Step 4] Assembling CNE contigs (03.debruijn) ..."
for sp in "${SPECIES[@]}"; do
    python3 src/03.debruijn.py "$RESULTS/$sp/" --mers "$MERS_TABLE" \
        -o "$RESULTS/$sp/assembled.fasta" \
        > "$RESULTS/$sp/03.log" 2>&1 &
done; wait

for sp in "${SPECIES[@]}"; do
    n=$(grep -c "^>" "$RESULTS/$sp/assembled.fasta" 2>/dev/null || echo 0)
    echo "  [$sp] $n contigs"
done

# ─── Step 5: 04.alignfree_phylo ───
echo ""
echo "[Step 5] Building phylogenetic tree (04) ..."

FASTAS=""
for sp in "${SPECIES[@]}"; do
    FASTAS="$FASTAS $RESULTS/$sp/assembled.fasta"
done

LABELS="--labels"
for sp in "${SPECIES[@]}"; do
    LABELS="$LABELS $sp"
done

python3 src/04.alignfree_phylo.py $FASTAS $LABELS \
    -t $THREADS -d mash -o "$RESULTS/tree10.nwk" --dist "$RESULTS/dist10.csv"

echo ""
echo "============================================"
echo "  Pipeline complete!"
echo "  Tree: $RESULTS/tree10.nwk"
echo "  Dist: $RESULTS/dist10.csv"
echo "============================================"
