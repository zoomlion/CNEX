#!/usr/bin/env python3
"""Classify CNE elements by genomic type (intergenic/intron) using GFF coordinates.

Input:  blocks_10k.fa (element FASTA blocks with species coordinates)
        gff_dir/       (GFF3 files named by species)

Output: element_tags.tsv (ele_id \\t type \\t chr \\t start \\t end)

Usage:
  python3 utils/classify_elements.py \\
      --msa benchmark/data/blocks_10k.fa \\
      --gff-dir benchmark/data/gff \\
      -o element_tags.tsv
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter


GFF_EXTENSIONS = [".gff3", ".gff3.gz", ".gff", ".gff.gz"]
HEADER_RE = re.compile(r">([^.]+)\.(\S+)([+-]):(\d+)-(\d+)")

ELEMENT_IDS_SAMPLED = None  # set from CLI --element-ids or auto-detect


def parse_args():
    p = argparse.ArgumentParser(description="Classify CNE elements by genomic type")
    p.add_argument("--msa", required=True, help="blocks_10k.fa path")
    p.add_argument("--gff-dir", default=None, help="Directory of GFF3 files (*.gff3, *.gff, .gz)")
    p.add_argument("--ref-species", default=None, help="Override reference species")
    p.add_argument("-o", "--output", default="element_tags.tsv")
    return p.parse_args()


def parse_blocks(fa_path):
    """Read blocks_10k.fa, yield (block_index, [{species, chr, strand, start, end}])."""
    with open(fa_path, errors='replace') as f:
        content = f.read()
    raw_blocks = [b.strip() for b in content.strip().split('###') if b.strip()]
    for bi, block in enumerate(raw_blocks):
        records = []
        for line in block.split('\n'):
            line = line.strip()
            if not line:
                continue
            m = HEADER_RE.match(line)
            if m:
                records.append({
                    "species": m.group(1),
                    "chr": m.group(2),
                    "strand": m.group(3),
                    "start": int(m.group(4)),
                    "end": int(m.group(5)),
                })
        yield bi, records


def find_gff(gff_dir, species_name):
    """Find GFF file for a species in gff_dir."""
    if not gff_dir or not os.path.isdir(gff_dir):
        return None
    for ext in GFF_EXTENSIONS:
        path = os.path.join(gff_dir, species_name + ext)
        if os.path.isfile(path):
            return path
    return None


def classify_with_bedtools(ele_coords, gff_path):
    """Classify elements as intron/intergenic using bedtools intersect.

    bedtools -wao output:
        BED(4) | GFF(9+) | overlap → GFF feature at field index 6 (0-based).
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        bed_path = f.name
        for eid, ch, st, en in ele_coords:
            if st > en:
                st, en = en, st
            f.write(f"{ch}\t{st}\t{en}\t{eid}\n")

    cmd = ["bedtools", "intersect", "-a", bed_path, "-b", gff_path, "-wao"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except Exception as e:
        os.unlink(bed_path)
        print(f"  Warning: bedtools failed ({e}), all elements -> intergenic")
        return {eid: "intergenic" for eid, _, _, _ in ele_coords}

    result = {eid: "intergenic" for eid, _, _, _ in ele_coords}
    for line in r.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 11:
            continue
        try:
            eid = int(parts[3])
        except ValueError:
            continue
        feature = parts[6]  # GFF column 3 (0-indexed: 6 after 4 BED cols)
        if feature in ("CDS", "exon"):
            result[eid] = "exon"
        elif feature in ("gene", "mRNA", "transcript") and result.get(eid) != "exon":
            result[eid] = "intron"

    os.unlink(bed_path)
    return result


def main():
    args = parse_args()

    if not os.path.isfile(args.msa):
        sys.exit(f"MSA file not found: {args.msa}")

    # Phase 1: count CNEs per species
    print("Counting CNEs per species...")
    species_counter = Counter()
    block_coords = {}  # element_id → [(species, chr, start, end)]
    for bi, records in parse_blocks(args.msa):
        seen = set()
        for rec in records:
            sp = rec["species"]
            if sp not in seen:
                species_counter[sp] += 1
                seen.add(sp)
        if records:
            # use first record's coordinate as representative (for counting only)
            block_coords[bi] = records

    # Phase 2: determine reference species
    ref_species = args.ref_species
    if not ref_species:
        # order by CNE count descending
        sorted_species = [sp for sp, _ in species_counter.most_common()]
        if not sorted_species:
            sys.exit("No species found in MSA")
        ref_species = sorted_species[0]
        print(f"  Most abundant: {ref_species} ({species_counter[ref_species]} CNEs)")

    # Find GFF for ref species
    gff_path = find_gff(args.gff_dir, ref_species) if args.gff_dir else None
    if gff_path:
        print(f"  GFF: {gff_path}")
    else:
        print(f"  No GFF found for {ref_species}, trying alternatives..." if args.gff_dir else "  No --gff-dir, all -> intergenic")
        if args.gff_dir:
            sorted_species = [sp for sp, _ in species_counter.most_common()]
            for sp in sorted_species:
                if sp == ref_species:
                    continue
                gff_path = find_gff(args.gff_dir, sp)
                if gff_path:
                    ref_species = sp
                    print(f"  Using GFF for {ref_species} instead")
                    break
            if not gff_path:
                print("  No GFF found for any species, all -> intergenic")

    # Phase 3: extract reference coordinates for each element
    print("Extracting reference coordinates...")
    ele_coords = []  # [(ele_id, chr, start, end)]
    for bi, records in block_coords.items():
        # Find the ref species record
        ref_rec = None
        for rec in records:
            if rec["species"] == ref_species:
                ref_rec = rec
                break
        if not ref_rec:
            # fallback: use any species
            ref_rec = records[0]
        ele_coords.append((bi, ref_rec["chr"], ref_rec["start"], ref_rec["end"]))

    # Phase 4: classify using GFF
    if gff_path:
        print("Classifying with bedtools intersect...")
        types = classify_with_bedtools(ele_coords, gff_path)
    else:
        print("All elements -> intergenic")
        types = {eid: "intergenic" for eid, _, _, _ in ele_coords}

    # Phase 5: output (skip exon elements)
    print(f"Writing {args.output}...")
    kept = 0
    with open(args.output, 'w') as f:
        f.write("ele_id\ttype\tchr\tstart\tend\n")
        for eid, ch, st, en in sorted(ele_coords):
            t = types.get(eid, "intergenic")
            if t == "exon":
                continue
            f.write(f"{eid}\t{t}\t{ch}\t{st}\t{en}\n")
            kept += 1

    # Summary (exon excluded)
    from collections import Counter
    cnt = Counter()
    for t in types.values():
        if t != "exon":
            cnt[t] += 1
    print(f"  {kept} elements (exon removed)")
    for k, v in cnt.most_common():
        print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
