#!/usr/bin/env python3
"""
Step 05: Per-element alignment (FAMSA) + phylogeny.

Methods:
  --method concat (default): align → concat_msa → IQ-TREE 3
  --method astral (fallback): align → FastTree → ASTRAL

Input:  results/fasta/<bunch_id>.fasta  (from step 04)
Output:
  concat: supermatrix.fa + partitions.txt + supermatrix.fa.treefile
  astral: species_tree.nwk + gene_trees.nwk

Usage:
  python3 src/05.element_phylo.py --max-elements 0 --method concat
"""

import argparse
import os
import subprocess
import sys
import glob
import time
from collections import Counter, OrderedDict
import multiprocessing as mp
try:
    mp.set_start_method('spawn')
except RuntimeError:
    pass
Pool = mp.get_context('spawn').Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import config as C
from concat_msa import trim_alignment_by_occupancy, read_fasta


def parse_args():
    p = argparse.ArgumentParser(description="Element phylogeny: FAMSA + concat/IQ-TREE or ASTRAL")
    p.add_argument("--elements-dir", default="results/fasta",
                   help="Directory with per-element FASTAs (default: results/fasta)")
    p.add_argument("--method", default="concat", choices=["concat", "astral"],
                   help="Phylogenetic method (default: concat)")
    p.add_argument("--famsa", default=C.FAMSA, help="Path to famsa binary")
    p.add_argument("--fasttree", default=C.FastTree, help="Path to FastTree binary")
    p.add_argument("--iqtree3", default=C.IQTREE3, help="Path to IQ-TREE 3 binary (--method concat)")
    p.add_argument("--iqtree-threads", type=int, default=C.IQTREE_THREADS,
                   help="Threads for IQ-TREE 3 (default: config value)")
    p.add_argument("--astral-bin", default=C.ASTRAL,
                   help="Path to ASTRAL IV / ASTER binary (--method astral)")
    p.add_argument("--astral-dir", default="",
                   help="Path to ASTRAL III directory (used if --astral-bin not found)")
    p.add_argument("--astral-jar", default="",
                   help="Path to ASTRAL III jar (overrides --astral-dir)")
    p.add_argument("--min-cne-per-species", type=int, default=C.MIN_CNE_PER_SPECIES,
                   help="Min CNE count per species to retain (default: config value)")
    p.add_argument("--min-occupancy", type=float, default=0.3,
                   help="Min species occupancy per element (default: 0.3)")
    p.add_argument("--min-site-occupancy", type=float, default=0.5,
                   help="Min site occupancy per column for trimming (default: 0.5)")
    p.add_argument("-t", "--threads", type=int, default=4,
                   help="Threads for IQ-TREE / ASTRAL (default: 4)")
    p.add_argument("--parallel", type=int, default=1,
                   help="Parallel elements for FAMSA alignment (default: 1)")
    p.add_argument("--max-elements", type=int, default=100,
                   help="Max elements to process (0=all)")
    p.add_argument("--resume", action="store_true", default=True,
                   help="Skip existing .aln files (default: on)")
    p.add_argument("--no-resume", action="store_false", dest="resume",
                   help="Force reprocess all elements")
    return p.parse_args()


def find_fasta_files(elements_dir, max_elements):
    files = sorted(glob.glob(os.path.join(elements_dir, "*.fasta")),
                   key=lambda x: int(os.path.basename(x).replace(".fasta", "")))
    if max_elements > 0:
        files = files[:max_elements]
    return files


def read_fasta_to_dict(path):
    seqs = OrderedDict()
    hdr = None
    buf = []
    with open(path, errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('>'):
                if hdr is not None:
                    seqs[hdr] = ''.join(buf)
                hdr = s[1:]
                buf = []
            else:
                buf.append(s)
        if hdr is not None:
            seqs[hdr] = ''.join(buf)
    return seqs


def write_fasta(seqs, path):
    with open(path, 'w') as f:
        for hdr, seq in seqs.items():
            f.write(f'>{hdr}\n{seq}\n')


def filter_species_by_cne_count(fasta_files, min_cne):
    """Compute species-to-keep based on CNE count across elements."""
    counter = Counter()
    for fp in fasta_files:
        seqs = read_fasta_to_dict(fp)
        counter.update(seqs.keys())
    keep = {sp for sp, cnt in counter.items() if cnt >= min_cne}
    removed = len(counter) - len(keep)
    if removed:
        print(f"  Species filter: {len(keep)} kept, {removed} removed (< {min_cne} CNEs)")
    return keep


def _flatten_fasta(path):
    with open(path, errors='replace') as f:
        lines = f.readlines()
    out = []
    buf = []
    for line in lines:
        if line.startswith(">"):
            if buf:
                out.append("".join(buf) + "\n")
                buf = []
            out.append(line)
        else:
            buf.append(line.strip())
    if buf:
        out.append("".join(buf) + "\n")
    with open(path, "w") as f:
        f.writelines(out)


def _align_one(args):
    fasta_path, famsa_bin, resume, aln_dir, keep_sp = args
    base = os.path.basename(fasta_path).replace(".fasta", "")
    aln_path = os.path.join(aln_dir, base + ".aln")
    if resume and os.path.isfile(aln_path):
        return True, fasta_path, None
    seqs = read_fasta_to_dict(fasta_path)
    if keep_sp is not None:
        seqs = {h: s for h, s in seqs.items() if h in keep_sp}
    if len(seqs) < 2:
        return False, fasta_path, "too_few_species"
    tmp_fa = fasta_path.replace('.fasta', '_filtered.fasta')
    write_fasta(seqs, tmp_fa)
    cmd = [famsa_bin, "-t", "1", tmp_fa, aln_path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return False, fasta_path, "align"
    _flatten_fasta(aln_path)
    os.remove(tmp_fa)
    return True, fasta_path, None


def run_concat_msa(aln_dir, min_occ, min_site_occ, out_dir, iqtree3, iqtree_threads):
    concat_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "utils", "concat_msa.py")
    supermatrix = os.path.join(out_dir, "supermatrix.fa")
    partitions = os.path.join(out_dir, "partitions.txt")
    cmd = ["python3", concat_script, "-i", aln_dir, "--suffix", ".aln",
           "--min-occupancy", str(min_occ),
           "--min-site-occupancy", str(min_site_occ),
           "-o", supermatrix, "-p", partitions]
    if subprocess.run(cmd).returncode != 0:
        print("  concat_msa.py failed!")
        return False
    iq_cmd = [iqtree3, "-s", supermatrix, "-p", partitions,
              "-m", "GTR+G4", "-bb", "1000", "-nt", str(iqtree_threads)]
    subprocess.run(iq_cmd)
    return True


def run_astral_method(fasta_files, aln_dir, out_dir, fasttree_bin,
                      astral_bin, astral_jar, threads, parallel, resume, keep_sp,
                      min_site_occ=0.0):
    nwk_dir = os.path.join(out_dir, "nwk")
    astral_outdir = os.path.join(out_dir, "astral")
    os.makedirs(nwk_dir, exist_ok=True)
    os.makedirs(astral_outdir, exist_ok=True)

    # Detect ASTRAL mode
    use_iv = os.path.isfile(astral_bin) if astral_bin else False
    jar_path = astral_jar
    if not use_iv and not jar_path:
        print("  No ASTRAL binary found, skipping ASTRAL step.")
        return

    # FastTree per element
    print(f"\n--- Gene trees (FastTree) ---")
    ok_count = fail_count = 0
    total = len(fasta_files)
    for i, fp in enumerate(fasta_files):
        base = os.path.basename(fp).replace(".fasta", "")
        nwk_path = os.path.join(nwk_dir, base + ".nwk")
        aln_path = os.path.join(aln_dir, base + ".aln")
        if resume and os.path.isfile(nwk_path):
            ok_count += 1
            continue
        if not os.path.isfile(aln_path):
            fail_count += 1
            continue
        # Trim alignment by site occupancy
        if min_site_occ > 0:
            sq = read_fasta(aln_path)
            trimmed = trim_alignment_by_occupancy(sq, min_site_occ)
            aln_trimmed = aln_path + ".trimmed"
            with open(aln_trimmed, 'w') as tf:
                for h, s in trimmed.items():
                    tf.write(f'>{h}\n{s}\n')
            fasta_input = aln_trimmed
        else:
            fasta_input = aln_path
        # Detect nucleotide type
        with open(fasta_input) as f:
            seq = "".join(line.strip() for line in f if not line.startswith(">"))[:100]
        is_nt = all(c in "ACGTacgtNn-" for c in seq) if seq else True
        cmd = [fasttree_bin]
        if is_nt:
            cmd.append("-nt")
        cmd.extend(["-gtr", "-nosupport"])
        with open(fasta_input) as inp, open(nwk_path, "w") as out:
            r = subprocess.run(cmd, stdin=inp, stdout=out, stderr=subprocess.PIPE, text=True)
        if r.returncode == 0:
            ok_count += 1
        else:
            fail_count += 1
    print(f"  FastTree: {ok_count} OK, {fail_count} FAIL")

    # Collect gene trees
    all_trees = []
    for fp in fasta_files:
        base = os.path.basename(fp).replace(".fasta", "")
        nwk_path = os.path.join(nwk_dir, base + ".nwk")
        if os.path.isfile(nwk_path):
            with open(nwk_path) as f:
                t = f.read().strip()
                if t:
                    all_trees.append(t)
    gene_trees_path = os.path.join(astral_outdir, "gene_trees.nwk")
    with open(gene_trees_path, 'w') as f:
        for t in all_trees:
            f.write(t + '\n')
    print(f"  Collected {len(all_trees)} gene trees")

    # ASTRAL
    if len(all_trees) < 2:
        print("  Too few gene trees for ASTRAL (need >= 2)")
        return
    species_tree_path = os.path.join(astral_outdir, "species_tree.nwk")
    if use_iv:
        cmd = [astral_bin, "-i", gene_trees_path, "-o", species_tree_path, "-t", str(threads)]
    else:
        cmd = ["java", "-Xmx8g", "-jar", jar_path, "-i", gene_trees_path,
               "-o", species_tree_path, "--extraLevel", "0"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode == 0:
        with open(species_tree_path) as f:
            print(f"\nSpecies tree:\n{f.read().strip()}")
    else:
        print(f"  ASTRAL failed: {r.stderr[:300]}")


def main():
    args = parse_args()
    famsa = os.path.expanduser(args.famsa)

    if not os.path.isfile(famsa):
        sys.exit(f"FAMSA not found: {famsa}")

    print(f"Method:   {args.method}")
    print(f"FAMSA:    {famsa}")
    print(f"Min CNE/species: {args.min_cne_per_species}")
    print(f"Max elements:    {'all' if args.max_elements == 0 else args.max_elements}")
    print()

    base_dir = os.path.dirname(os.path.abspath(args.elements_dir))
    aln_dir = os.path.join(base_dir, "aln")
    os.makedirs(aln_dir, exist_ok=True)

    fasta_files = find_fasta_files(args.elements_dir, args.max_elements)
    print(f"Found {len(fasta_files)} element FASTAs")
    if not fasta_files:
        sys.exit("No element FASTAs found. Run step 04 first.")

    keep_sp = filter_species_by_cne_count(fasta_files, args.min_cne_per_species)

    # FAMSA alignment
    print(f"\n--- Aligning elements (FAMSA) ---")
    t0 = time.time()
    work = [(f, famsa, args.resume, aln_dir, keep_sp) for f in fasta_files]
    ok_count = fail_count = 0
    if args.parallel > 1:
        with Pool(args.parallel) as p:
            for ok, _, _ in p.imap_unordered(_align_one, work):
                if ok: ok_count += 1
                else: fail_count += 1
    else:
        for w in work:
            ok, _, _ = _align_one(w)
            if ok: ok_count += 1
            else: fail_count += 1
    print(f"  Done in {time.time() - t0:.1f}s ({ok_count} OK, {fail_count} FAIL)")

    # Phylogeny
    if args.method == "concat":
        iqtree3 = os.path.expanduser(args.iqtree3)
        if not os.path.isfile(iqtree3):
            sys.exit(f"IQ-TREE 3 not found: {iqtree3}")
        print(f"\n--- Concatenation + IQ-TREE 3 ---")
        ok = run_concat_msa(aln_dir, args.min_occupancy, args.min_site_occupancy,
                            base_dir, iqtree3, args.iqtree_threads)
        if ok:
            print(f"  Species tree: {os.path.join(base_dir, 'supermatrix.fa.treefile')}")
    else:
        fasttree = os.path.expanduser(args.fasttree)
        astral_bin = os.path.expanduser(args.astral_bin)
        run_astral_method(fasta_files, aln_dir, base_dir, fasttree,
                          astral_bin, args.astral_jar, args.threads,
                          args.parallel, args.resume, keep_sp,
                          args.min_site_occupancy)

    print("\nDone!")


if __name__ == "__main__":
    main()
