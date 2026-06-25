#!/bin/bash
# simu.sh - Generate 3X simulated PE reads for 10 species
set -e
WGSIM=~/Software/wgsim/wgsim
PIGZ=src/pigz/pigz
GENOME_DIR=~/Source/Genome
mkdir -p reads

gen() {
    local SP="$1" N="$2" FA="$3"
    local out="reads/${SP}_r1.fq.gz"
    if [ -f "$out" ] && [ $(stat -c%s "$out" 2>/dev/null || echo 0) -gt 1000000 ]; then
        echo "[$SP] already exists, skipping"
        return 0
    fi
    echo "[$SP] ${N} PE pairs (3X) ..."
    rm -f /tmp/fifo_r1_$$ /tmp/fifo_r2_$$
    mkfifo /tmp/fifo_r1_$$ /tmp/fifo_r2_$$
    $PIGZ < /tmp/fifo_r1_$$ > "reads/${SP}_r1.fq.gz" &
    local pid1=$!
    $PIGZ < /tmp/fifo_r2_$$ > "reads/${SP}_r2.fq.gz" &
    local pid2=$!
    $WGSIM -N "$N" -1 150 -2 150 -e 0.001 -r 0.001 -R 0.0001 \
        "$FA" /tmp/fifo_r1_$$ /tmp/fifo_r2_$$
    wait $pid1 $pid2
    rm -f /tmp/fifo_r1_$$ /tmp/fifo_r2_$$
    local sz=$(du -sh "reads/${SP}_r1.fq.gz" | cut -f1)
    echo "[$SP] done: $sz"
}

G="$GENOME_DIR"

gen Homo         30997508 "$G/Homo_sapiens/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gen Mus          27282225 "$G/Mus_musculus/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa"
gen Danio        13734714 "$G/Danio_rerio/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa"
gen Loxodonta    31967609 "$G/Loxodonta_africana/Loxodonta_africana.loxAfr3.dna_sm.toplevel.fa"
gen Gallus        5326827 "$G/Gallus_gallus/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa"
gen Xenopus       7201992 "$G/Xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.dna_sm.toplevel.fa"
gen Latimeria    14302959 "$G/Latimeria_chalumnae/Latimeria_chalumnae.LatCha1.dna_sm.toplevel.fa"
gen Callorhinchus 4872492 "$G/Callorhinchus_milii/Callorhinchus_milii.Callorhinchus_milii-6.1.3.dna_sm.toplevel.fa"
gen Polypterus   18369400 "$G/Polypterus_senegalus/Polypterus_senegalus.ASM1683550v1_genomic_GCF.fa"
gen Rhincodon    14657997 "$G/Rhincodon_typus/Rhincodon_typus.ASM164234v2_genomic_GCF.fa"

echo "=== All done ==="
ls -lhS reads/*.fq.gz
