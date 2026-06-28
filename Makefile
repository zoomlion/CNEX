# CNEX Pipeline Makefile
#   make reads  - generate 3X simulated reads
#   make 01     - build confident k-mer table
#   make 02     - validate reads or genome (parallel: make -j 5 02)
#   make 03     - de Bruijn assembly (parallel: make -j 10 03)
#   make 04     - build phylogenetic tree
#   make all    - full pipeline
#   make clean  - remove results (keeps reads)

THREADS        ?= 8
THREADS_PER_02 ?= 8
MIN_C          ?= 10
VOTE_FRAC      ?= 0.1
VOTE_RATIO     ?= 5.0
WGSIM          := $(HOME)/Software/wgsim/wgsim
PIGZ           := src/pigz/pigz
ASSEMBLER      := src/02.validate
MERS_TABLE     := mers_table.tsv
RESULTS        := results

SPECIES := Homo Mus Danio Loxodonta Gallus Xenopus Latimeria Callorhinchus Polypterus Rhincodon

# Genome paths & 3X PE pair counts
DANIO_GENOME         := $(HOME)/Source/Genome/Danio_rerio/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa
DANIO_N              := 22891190
HOMO_GENOME          := $(HOME)/Source/Genome/Homo_sapiens/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
HOMO_N               := 51662513
MUS_GENOME           := $(HOME)/Source/Genome/Mus_musculus/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
MUS_N                := 45470375
LOXODONTA_GENOME     := $(HOME)/Source/Genome/Loxodonta_africana/Loxodonta_africana.loxAfr3.dna_sm.toplevel.fa
LOXODONTA_N          := 53279348
GALLUS_GENOME        := $(HOME)/Source/Genome/Gallus_gallus/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa
GALLUS_N             := 8878045
XENOPUS_GENOME       := $(HOME)/Source/Genome/Xenopus_tropicalis/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.dna_sm.toplevel.fa
XENOPUS_N            := 12003320
LATIMERIA_GENOME     := $(HOME)/Source/Genome/Latimeria_chalumnae/Latimeria_chalumnae.LatCha1.dna_sm.toplevel.fa
LATIMERIA_N          := 23838265
CALLORHINCHUS_GENOME := $(HOME)/Source/Genome/Callorhinchus_milii/Callorhinchus_milii.Callorhinchus_milii-6.1.3.dna_sm.toplevel.fa
CALLORHINCHUS_N      := 8120820
POLYPTERUS_GENOME    := $(HOME)/Source/Genome/Polypterus_senegalus/Polypterus_senegalus.ASM1683550v1_genomic_GCF.fa
POLYPTERUS_N         := 30615667
RHINCODON_GENOME     := $(HOME)/Source/Genome/Rhincodon_typus/Rhincodon_typus.ASM164234v2_genomic_GCF.fa
RHINCODON_N          := 24429995

# ─── Phony targets ───
.PHONY: all reads 00 01 02 03 04 clean

all: 04

reads: $(foreach sp,$(SPECIES),reads/$(sp)_r1.fq.gz)

00: most-cons-cne.filtered.fa

01: $(MERS_TABLE)

02: $(foreach sp,$(SPECIES),$(RESULTS)/$(sp)/.02_done)

03: $(foreach sp,$(SPECIES),$(RESULTS)/$(sp)/assembled.fasta)

04: $(RESULTS)/tree10.nwk

# ─── Step 00: MSA block filter ───
most-cons-cne.filtered.fa: most-cons-cne.fa
	@echo "[00] Filtering MSA blocks ..."
	python3 src/00.filter_msa.py most-cons-cne.fa --min-species 8 --min-avg-len 50 -o $@

# ─── Step 01: confident k-mer table ───
$(MERS_TABLE): most-cons-cne.filtered.fa
	@echo "[01] Building confident k-mer table ..."
	python3 src/01.confi_mer.py most-cons-cne.filtered.fa -k 13 -c 4 -o $@
	@echo "[01] done ($(shell wc -l < $@) lines)"

# ─── Simulated reads generation ───

reads/Danio_r1.fq.gz: $(DANIO_GENOME)
	@./shell/wgsim_gen.sh Danio "$(DANIO_GENOME)" $(DANIO_N) "$(WGSIM)" "$(PIGZ)"

reads/Homo_r1.fq.gz: $(HOMO_GENOME)
	@./shell/wgsim_gen.sh Homo "$(HOMO_GENOME)" $(HOMO_N) "$(WGSIM)" "$(PIGZ)"

reads/Mus_r1.fq.gz: $(MUS_GENOME)
	@./shell/wgsim_gen.sh Mus "$(MUS_GENOME)" $(MUS_N) "$(WGSIM)" "$(PIGZ)"

reads/Loxodonta_r1.fq.gz: $(LOXODONTA_GENOME)
	@./shell/wgsim_gen.sh Loxodonta "$(LOXODONTA_GENOME)" $(LOXODONTA_N) "$(WGSIM)" "$(PIGZ)"

reads/Gallus_r1.fq.gz: $(GALLUS_GENOME)
	@./shell/wgsim_gen.sh Gallus "$(GALLUS_GENOME)" $(GALLUS_N) "$(WGSIM)" "$(PIGZ)"

reads/Xenopus_r1.fq.gz: $(XENOPUS_GENOME)
	@./shell/wgsim_gen.sh Xenopus "$(XENOPUS_GENOME)" $(XENOPUS_N) "$(WGSIM)" "$(PIGZ)"

reads/Latimeria_r1.fq.gz: $(LATIMERIA_GENOME)
	@./shell/wgsim_gen.sh Latimeria "$(LATIMERIA_GENOME)" $(LATIMERIA_N) "$(WGSIM)" "$(PIGZ)"

reads/Callorhinchus_r1.fq.gz: $(CALLORHINCHUS_GENOME)
	@./shell/wgsim_gen.sh Callorhinchus "$(CALLORHINCHUS_GENOME)" $(CALLORHINCHUS_N) "$(WGSIM)" "$(PIGZ)"

reads/Polypterus_r1.fq.gz: $(POLYPTERUS_GENOME)
	@./shell/wgsim_gen.sh Polypterus "$(POLYPTERUS_GENOME)" $(POLYPTERUS_N) "$(WGSIM)" "$(PIGZ)"

reads/Rhincodon_r1.fq.gz: $(RHINCODON_GENOME)
	@./shell/wgsim_gen.sh Rhincodon "$(RHINCODON_GENOME)" $(RHINCODON_N) "$(WGSIM)" "$(PIGZ)"

# ─── Step 02: read validation ───
$(RESULTS)/%/.02_done: reads/%_r1.fq.gz reads/%_r2.fq.gz $(MERS_TABLE)
	@echo "[02] $* ..."
	mkdir -p $(RESULTS)/$*
	$(ASSEMBLER) reads/$*_r1.fq.gz reads/$*_r2.fq.gz \
		--mers $(MERS_TABLE) -t $(THREADS_PER_02) \
		--output_dir $(RESULTS)/$* --pigz src/pigz/pigz \
		--min-c $(MIN_C) --vote-frac $(VOTE_FRAC) --vote-ratio $(VOTE_RATIO) \
		> $(RESULTS)/$*/02.log 2>&1
	touch $@
	@echo "[02] $* done"

# ─── Step 03: de Bruijn assembly ───
$(RESULTS)/%/assembled.fasta: $(RESULTS)/%/.02_done
	@echo "[03] $* ..."
	python3 src/03.debruijn.py $(RESULTS)/$*/ --mers $(MERS_TABLE) \
		-o $@ --trim > $(RESULTS)/$*/03.log 2>&1
	@echo "[03] $* done ($(shell grep -c '^>' $@ 2>/dev/null || echo 0) contigs)"

# ─── Step 04: phylogenetic tree ───
$(RESULTS)/tree10.nwk: $(foreach sp,$(SPECIES),$(RESULTS)/$(sp)/assembled.fasta)
	@echo "[04] Building phylogenetic tree ..."
	python3 src/04.alignfree_phylo.py $(foreach sp,$(SPECIES),$(RESULTS)/$(sp)/assembled.fasta) \
		-t $(THREADS) -d containment -m per-element \
		--fusion weighted \
		-o $@ --dist $(RESULTS)/dist10.csv
	@echo "[04] done"
	@echo "  Tree: $@"
	@echo "  Dist: $(RESULTS)/dist10.csv"

# ─── Clean ───
clean:
	rm -rf $(RESULTS)
	@echo "Cleaned $(RESULTS)/"
