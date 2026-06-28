# CNEX pipeline Makefile
# Pattern rules — one rule covers all species
# Usage:
#   make -f shell/pipeline.mk 02 -j 2
#   make -f shell/pipeline.mk 03 -j 8
#   make -f shell/pipeline.mk 04

ASSEMBLER  ?= src/02.validate
MERS_TABLE ?= mers_table.tsv
THREADS    ?= 8
MIN_C      ?= 10
VOTE_FRAC  ?= 0.1
VOTE_RATIO ?= 5.0
RESULTS    ?= results

# Auto-discover species with valid reads
VALID_SPECIES := $(foreach f,$(patsubst reads/%_r1.fq.gz,%,$(wildcard reads/*_r1.fq.gz)),$(if $(findstring _,$f),$f))

# ─── Pattern rules ───

results/%/.02_done: reads/%_r1.fq.gz reads/%_r2.fq.gz $(MERS_TABLE)
	@mkdir -p $(dir $@)
	$(ASSEMBLER) $< $(subst _r1,_r2,$<) \
		--mers $(MERS_TABLE) -t $(THREADS) \
		--output_dir $(dir $@) --pigz src/pigz/pigz \
		--min-c $(MIN_C) --vote-frac $(VOTE_FRAC) --vote-ratio $(VOTE_RATIO) \
		> $(dir $@)02.log 2>&1
	touch $@

results/%/assembled.fasta: results/%/.02_done
	python3 src/03.debruijn.py $(dir $@) --mers $(MERS_TABLE) \
		-o $@ --trim > $(dir $@)03.log 2>&1

# ─── Phony targets ───

.PHONY: 02 03 04

02: $(foreach sp,$(VALID_SPECIES),$(RESULTS)/$(sp)/.02_done)

03: $(foreach sp,$(VALID_SPECIES),$(RESULTS)/$(sp)/assembled.fasta)

04: 03
	python3 src/04.alignfree_phylo.py $(RESULTS)/ -o $(RESULTS)/tree_all.nwk -b 100
