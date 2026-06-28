#!/usr/bin/env python3
"""
discover_genomes.py - Auto-discover gnathostome (jawed vertebrate) species
in ~/Source/Genome/, skip non-gnathostomes and already-processed species,
detect genome FASTA files, compute 5X coverage read counts,
and output a Makefile fragment for the CNEX pipeline.

Output: shell/species.mk (Makefile fragment)
"""

import os, sys, glob, math, re, subprocess

GENOME_DIR = os.path.expanduser("~/Source/Genome")
MERCURY_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
READS_DIR = os.path.join(MERCURY_DIR, "reads")
RESULTS_DIR = os.path.join(MERCURY_DIR, "results")
OUTPUT_MK = os.path.join(MERCURY_DIR, "shell", "species.mk")

# Non-gnathostomes (invertebrates + cyclostomes) — exclude
EXCLUDE = {
    "Acropora_millepora",
    "Amphimedon_queenslandica",
    "Asterias_rubens",
    "Branchiostomia_floridae",
    "Caenorhabditis_elegans",
    "Ciona_intestinalis",
    "Ciona_savignyi",
    "Crassostrea_gigas",
    "Drosophila_melanogaster",
    "Eptatretus_burgeri",
    "Hydra_vulgaris",
    "Octopus_bimaculoides",
    "Octopus_sinensis",
    "Petromyzon_marinus",
    "Petromyzon_marinus_new",
    "Rhopilema_esculentum",
    "Protopterus_annectens",  # 38G genome — too large
    "temp",
}

# Already processed 10 species — skip simulation but include in species list
ALREADY_DONE = {
    "Homo_sapiens": "Homo_sapiens",
    "Mus_musculus": "Mus_musculus",
    "Danio_rerio": "Danio_rerio",
    "Loxodonta_africana": "Loxodonta_africana",
    "Gallus_gallus": "Gallus_gallus",
    "Xenopus_tropicalis": "Xenopus_tropicalis",
    "Latimeria_chalumnae": "Latimeria_chalumnae",
    "Callorhinchus_milii": "Callorhinchus_milii",
    "Polypterus_senegalus": "Polypterus_senegalus",
    "Rhincodon_typus": "Rhincodon_typus",
}

# Manual short-name overrides for special cases
SHORT_NAME_OVERRIDE = {
    "Goodes_thornscrub_tortoise": "Gopherus_evgoodei",
}


def get_genome_size(fasta_path):
    """Get genome size from .fai file (fast) or file size (approx)."""
    fai = fasta_path + ".fai"
    if os.path.exists(fai):
        total = 0
        with open(fai) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    total += int(parts[1])
        if total > 0:
            return total
    # Fallback: file size is a close approximation for plain-text FASTA
    return os.path.getsize(fasta_path)


def find_genome_fasta(species_dir):
    """Find the main genome FASTA file in a species directory.
    
    Priority:
    1. *.dna_sm.toplevel.fa  (Ensembl)
    2. *genomic_GCF.fa       (NCBI RefSeq)
    3. *genomic_GCA.fa       (NCBI GenBank)
    4. *.fa / *.fasta / *.fna (other)
    """
    candidates = []
    for fname in os.listdir(species_dir):
        fpath = os.path.join(species_dir, fname)
        if not os.path.isfile(fpath):
            continue
        if fname.endswith(".fai") or fname.endswith(".gff") or fname.endswith(".gff3"):
            continue
        if "proteome" in fname.lower() or ".cds" in fname.lower() or ".pep" in fname.lower():
            continue
        if fname.endswith(".fa") or fname.endswith(".fasta") or fname.endswith(".fna"):
            candidates.append(fpath)
    
    if not candidates:
        return None
    
    # Score-based selection
    def score(path):
        f = os.path.basename(path)
        s = 0
        if "dna_sm" in f or "dna" in f:
            s += 100  # Ensembl-style: whole genome
        if "toplevel" in f or "primary_assembly" in f:
            s += 50
        if "genomic" in f:
            s += 30
        if "GCF" in f or "GCA" in f:
            s += 20
        if "chromosome" in f.lower() or "chromosome" in f.lower():
            s -= 10
        s += os.path.getsize(path) / 1e9  # prefer larger files
        return s
    
    best = max(candidates, key=score)
    return best


def get_short_name(species_name):
    """Get a short binomial name (already binomial)."""
    if species_name in SHORT_NAME_OVERRIDE:
        return SHORT_NAME_OVERRIDE[species_name]
    return species_name


def main():
    os.makedirs(os.path.join(MERCURY_DIR, "shell"), exist_ok=True)
    
    if not os.path.isdir(GENOME_DIR):
        print(f"Error: genome directory {GENOME_DIR} not found", file=sys.stderr)
        sys.exit(1)
    
    all_species = []
    new_species = []
    done_species_names = set(ALREADY_DONE.values())
    
    for entry in sorted(os.listdir(GENOME_DIR)):
        species_dir = os.path.join(GENOME_DIR, entry)
        if not os.path.isdir(species_dir):
            continue
        if entry in EXCLUDE:
            continue
        
        short = get_short_name(entry)
        is_new = short not in done_species_names
        
        genome_fa = find_genome_fasta(species_dir)
        if genome_fa is None:
            print(f"  SKIP {entry}: no genome FASTA found", file=sys.stderr)
            continue
        
        genome_size = get_genome_size(genome_fa)
        # 5X PE150: N = genome_size * 5 / 300
        n_pairs = genome_size * 5 // 300
        
        info = {
            "species_dir": entry,
            "short_name": short,
            "genome_fa": genome_fa,
            "genome_size": genome_size,
            "n_pairs": n_pairs,
            "is_new": is_new,
        }
        all_species.append(info)
        if is_new:
            new_species.append(info)
    
    # Print summary
    print(f"\n{'=' * 60}")
    print(f"Discovered {len(all_species)} gnathostome species total")
    print(f"  Already processed: {len(all_species) - len(new_species)}")
    print(f"  New (need simulation): {len(new_species)}")
    print(f"{'=' * 60}\n")
    
    print(f"{'Species':30s} {'Genome size':>12s} {'5X N':>15s} {'Status':>10s}")
    print(f"{'-' * 70}")
    for info in all_species:
        status = "DONE" if not info["is_new"] else "NEW"
        size_str = f"{info['genome_size'] / 1e9:.2f}G"
        n_str = f"{info['n_pairs']:,}"
        print(f"{info['short_name']:30s} {size_str:>12s} {n_str:>15s} {status:>10s}")
    
    # Generate Makefile fragment
    lines = []
    lines.append("# Auto-generated by discover_genomes.py")
    lines.append(f"# {len(all_species)} gnathostome species total")
    lines.append(f"# {len(new_species)} new species needing simulation")
    lines.append("")
    
    # ALL_SPECIES list
    all_names = [info["short_name"] for info in all_species]
    lines.append(f"ALL_SPECIES := {' '.join(all_names)}")
    lines.append("")
    
    # NEW_SPECIES list (need simulation)
    new_names = [info["short_name"] for info in new_species]
    lines.append(f"NEW_SPECIES := {' '.join(new_names)}")
    lines.append("")
    
    # DONE_SPECIES list (already processed)
    done_names = [info["short_name"] for info in all_species if not info["is_new"]]
    lines.append(f"DONE_SPECIES := {' '.join(done_names)}")
    lines.append("")
    
    # Genome file paths
    lines.append("# Genome file paths")
    for info in all_species:
        lines.append(f"GENOME_{info['short_name']} := {info['genome_fa']}")
    lines.append("")
    
    # 5X N values
    lines.append("# 5X PE150 read pair counts")
    for info in all_species:
        lines.append(f"N_{info['short_name']} := {info['n_pairs']}")
    lines.append("")
    
    # Read file prerequisites
    lines.append("# Read file prerequisites")
    for info in new_species:
        lines.append(f"reads/{info['short_name']}_r1.fq.gz: $(GENOME_{info['short_name']})")
        lines.append(f"\t@./shell/simu_one.sh {info['short_name']} \"$(GENOME_{info['short_name']})\" $(N_{info['short_name']})")
    lines.append("")
    
    # phony read target for new species
    lines.append(f".PHONY: reads_new")
    lines.append(f"reads_new: $(foreach sp,$(NEW_SPECIES),reads/$(sp)_r1.fq.gz)")
    lines.append("")
    
    # phony read target for all species
    lines.append(f".PHONY: reads_all")
    lines.append(f"reads_all: $(foreach sp,$(ALL_SPECIES),reads/$(sp)_r1.fq.gz)")
    lines.append("")
    
    # 02 assembly targets
    lines.append("# 02 assembly targets (all species)")
    for info in all_species:
        sn = info["short_name"]
        lines.append(f"$(RESULTS)/{sn}/.02_done: reads/{sn}_r1.fq.gz reads/{sn}_r2.fq.gz $(MERS_TABLE)")
        lines.append(f"\t@mkdir -p $(RESULTS)/{sn}")
        lines.append(f"\t$(ASSEMBLER) reads/{sn}_r1.fq.gz reads/{sn}_r2.fq.gz \\")
        lines.append(f"\t\t--mers $(MERS_TABLE) -t $(THREADS_PER_02) \\")
        lines.append(f"\t\t--output_dir $(RESULTS)/{sn} --pigz src/pigz/pigz \\")
        lines.append(f"\t\t--min-c $(MIN_C) --vote-frac $(VOTE_FRAC) --vote-ratio $(VOTE_RATIO) \\")
        lines.append(f"\t\t> $(RESULTS)/{sn}/02.log 2>&1")
        lines.append(f"\ttouch $@")
    lines.append("")
    
    lines.append(f".PHONY: 02_all")
    lines.append(f"02_all: $(foreach sp,$(ALL_SPECIES),$(RESULTS)/$(sp)/.02_done)")
    lines.append("")
    
    # 03 assembly targets
    lines.append("# 03 assembly targets (all species)")
    for info in all_species:
        sn = info["short_name"]
        lines.append(f"$(RESULTS)/{sn}/assembled.fasta: $(RESULTS)/{sn}/.02_done")
        lines.append(f"\tpython3 src/03.debruijn.py $(RESULTS)/{sn}/ --mers $(MERS_TABLE) \\")
        lines.append(f"\t\t-o $@ --trim > $(RESULTS)/{sn}/03.log 2>&1")
    lines.append("")
    
    lines.append(f".PHONY: 03_all")
    lines.append(f"03_all: $(foreach sp,$(ALL_SPECIES),$(RESULTS)/$(sp)/assembled.fasta)")
    lines.append("")
    
    # 04 target
    lines.append("# 04 phylogenetic tree (all species)")
    lines.append(f"04_all: $(foreach sp,$(ALL_SPECIES),$(RESULTS)/$(sp)/assembled.fasta)")
    lines.append(f"\tpython3 src/04.alignfree_phylo.py $(RESULTS)/ results/ -o $(RESULTS)/tree_all.nwk")
    lines.append("")
    
    content = "\n".join(lines)
    with open(OUTPUT_MK, "w") as f:
        f.write(content)
    
    print(f"\nMakefile fragment written to: {OUTPUT_MK}")
    print("\nTo simulate new species:")
    print("  make -f shell/species.mk reads_new -j 3")
    print("\nTo run full pipeline on all species:")
    print("  make -f shell/species.mk 02_all -j 3")
    print("  make -f shell/species.mk 03_all -j 3")
    print("  make -f shell/species.mk 04_all")


if __name__ == "__main__":
    main()
