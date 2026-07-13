#!/usr/bin/env python3
"""
Step 05: Per-element alignment (FAMSA) + gene tree (FastTree) + species tree (ASTRAL).

Input:
  results/elements/<bunch_id>.fasta  -- per-element multi-species FASTA (from step 05)

Output:
  results/elements/<bunch_id>.aln    -- aligned FASTA
  results/elements/<bunch_id>.nwk    -- gene tree
  results/species_tree.nwk           -- ASTRAL species tree
  results/gene_trees.nwk             -- all gene trees concatenated

Usage:
  python3 src/06.element_phylo.py --num-threads 8 --parallel 8 --max-elements 0
"""

import argparse
import os
import subprocess
import sys
import glob
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as C


def parse_args():
    p = argparse.ArgumentParser(description="Element phylogeny: famsa + fasttree + astral")
    p.add_argument("--elements-dir", default="results/fasta",
                   help="Directory with per-element FASTAs (default: results/fasta)")
    p.add_argument("--famsa", default=C.FAMSA,
                   help="Path to famsa binary")
    p.add_argument("--fasttree", default=C.FastTree,
                   help="Path to FastTree binary")
    p.add_argument("--astral-bin", default=C.ASTRAL,
                   help="Path to ASTRAL IV / ASTER native binary")
    p.add_argument("--astral-dir", default="",
                   help="Path to ASTRAL III directory (used only if --astral-bin not found)")
    p.add_argument("--astral-jar", default="",
                   help="Path to ASTRAL III jar (overrides --astral-dir)")
    p.add_argument("-t", "--threads", type=int, default=4,
                   help="Threads for ASTRAL IV / ASTER (default: 4)")
    p.add_argument("--parallel", type=int, default=1,
                    help="Number of parallel element processes (famsa always uses 1 thread each)")
    p.add_argument("--max-elements", type=int, default=100,
                    help="Max elements to process (0=all)")
    p.add_argument("--max-gap", type=float, default=0.5,
                    help="Max gap fraction per column in alignment trimming (default: 0.5)")
    p.add_argument("--resume", action="store_true", default=True,
                   help="Skip elements with existing .nwk files (default: on)")
    p.add_argument("--no-resume", action="store_false", dest="resume",
                   help="Force reprocess all elements")
    return p.parse_args()


def find_fasta_files(elements_dir, max_elements):
    files = sorted(glob.glob(os.path.join(elements_dir, "*.fasta")),
                   key=lambda x: int(os.path.basename(x).replace(".fasta", "")))
    if max_elements > 0:
        files = files[:max_elements]
    return files


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


def trim_alignment(aln_path, max_gap=0.5):
    """Remove alignment columns where gap fraction > max_gap."""
    with open(aln_path) as f:
        lines = f.readlines()
    hdrs, seqs = [], []
    for line in lines:
        if line.startswith('>'):
            hdrs.append(line.strip())
        else:
            seqs.append(list(line.strip()))
    if not seqs or not seqs[0]:
        return aln_path
    ncol = len(seqs[0])
    keep = [j for j in range(ncol)
            if sum(seqs[i][j] == '-' for i in range(len(seqs))) / len(seqs) <= max_gap]
    out_path = aln_path.replace('.aln', '.trimmed.aln')
    with open(out_path, 'w') as f:
        for i, h in enumerate(hdrs):
            f.write(f'{h}\n{"".join(seqs[i][j] for j in keep)}\n')
    return out_path


def _process_one(args_tuple):
    fasta_path, famsa_bin, fasttree_bin, resume, aln_dir, nwk_dir, max_gap = args_tuple

    base = os.path.basename(fasta_path).replace(".fasta", "")
    aln_path = os.path.join(aln_dir, base + ".aln")
    nwk_path = os.path.join(nwk_dir, base + ".nwk")

    # Align (FAMSA, always 1 thread per element process)
    if resume and os.path.isfile(aln_path):
        aln_ok = True
    else:
        cmd = [famsa_bin, "-t", "1", fasta_path, aln_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return False, fasta_path, "align"
        _flatten_fasta(aln_path)
        aln_ok = True

    # Alignment trimming (zero deps)
    tree_in = trim_alignment(aln_path, max_gap)

    # Gene tree (FastTree)
    if resume and os.path.isfile(nwk_path):
        tree_ok = True
    else:
        if not os.path.isfile(tree_in):
            return False, fasta_path, "no_aln"
        with open(tree_in) as f:
            seq = ""
            for line in f:
                if not line.startswith(">"):
                    seq += line.strip()
                if len(seq) > 100:
                    break
        is_nucleotide = all(c in "ACGTacgtNn-" for c in seq) if seq else True
        cmd = [fasttree_bin]
        if is_nucleotide:
            cmd.append("-nt")
        cmd.extend(["-gtr", "-nosupport"])
        with open(aln_path) as inp, open(nwk_path, "w") as out:
            result = subprocess.run(cmd, stdin=inp, stdout=out,
                                    stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            return False, fasta_path, "tree"
        tree_ok = True

    return True, fasta_path, None


def run_astral_iv(bin_path, gene_trees_path, output_path, threads=4):
    cmd = [bin_path, "-i", gene_trees_path, "-o", output_path, "-t", str(threads)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ASTRAL IV error: {result.stderr[:500]}", file=sys.stderr)
        return False
    return True


def run_astral_iii(jar_path, gene_trees_path, output_path):
    cmd = ["java", "-Xmx8g", "-jar", jar_path, "-i", gene_trees_path,
           "-o", output_path, "--extraLevel", "0"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ASTRAL III error: {result.stderr[:500]}", file=sys.stderr)
        return False
    return True


def main():
    args = parse_args()

    famsa = os.path.expanduser(args.famsa)
    fasttree = os.path.expanduser(args.fasttree)
    astral_bin = os.path.expanduser(args.astral_bin)
    astral_dir = os.path.expanduser(args.astral_dir)
    astral_jar = os.path.expanduser(args.astral_jar) if args.astral_jar else ""

    if not os.path.isfile(famsa):
        sys.exit(f"FAMSA not found: {famsa}")
    if not os.path.isfile(fasttree):
        sys.exit(f"FastTree not found: {fasttree}")

    # Detect ASTRAL mode: IV (native) > III (jar)
    use_astral_iv = False
    if os.path.isfile(astral_bin):
        use_astral_iv = True
    elif not astral_jar:
        candidates = [os.path.join(astral_dir, "astral.5.7.1.jar")]
        for f in glob.glob(os.path.join(astral_dir, "*.jar")):
            candidates.insert(0, f)
        for c in candidates:
            if os.path.isfile(c):
                astral_jar = c
                break
        if not astral_jar:
            sys.exit(f"No ASTRAL binary found.\n"
                     f"  Try: --astral-bin {astral_bin} (IV/native)\n"
                     f"  Or:  --astral-jar <path> (III/Java)")

    print(f"FAMSA:    {famsa}")
    print(f"FastTree: {fasttree}")
    if use_astral_iv:
        print(f"ASTRAL:   {astral_bin} (IV/native, {args.threads} threads)")
    else:
        print(f"ASTRAL:   {astral_jar} (III/Java)")
    print(f"Parallel:      {args.parallel} elements simultaneously")
    print(f"Resume:   {'on' if args.resume else 'off'}")
    print()

    # Setup output subdirectories based on elements_dir parent
    base_dir = os.path.dirname(os.path.abspath(args.elements_dir))
    aln_dir = os.path.join(base_dir, "aln")
    nwk_dir = os.path.join(base_dir, "nwk")
    astral_outdir = os.path.join(base_dir, "astral")
    for d in (aln_dir, nwk_dir, astral_outdir):
        os.makedirs(d, exist_ok=True)

    fasta_files = find_fasta_files(args.elements_dir, args.max_elements)
    print(f"Found {len(fasta_files)} element FASTAs")
    if not fasta_files:
        sys.exit("No element FASTAs found. Run step 05 first.")

    # Process elements (parallel)
    total = len(fasta_files)
    work = [(f, famsa, fasttree, args.resume, aln_dir, nwk_dir, args.max_gap) for f in fasta_files]

    t0 = time.time()
    ok_count = 0
    fail_count = 0
    done_count = 0

    if args.parallel > 1:
        with Pool(args.parallel) as p:
            for ok, fpath, stage in p.imap_unordered(_process_one, work):
                done_count += 1
                if ok:
                    ok_count += 1
                else:
                    fail_count += 1
                if done_count % 500 == 0 or done_count == total:
                    elapsed = time.time() - t0
                    rate = done_count / elapsed if elapsed > 0 else 0
                    eta = (total - done_count) / rate if rate > 0 else 0
                    print(f"  {done_count}/{total}  OK={ok_count}  "
                          f"FAIL={fail_count}  {rate:.1f}/s  ETA={eta:.0f}s")
    else:
        for w in work:
            ok, fpath, stage = _process_one(w)
            done_count += 1
            if ok:
                ok_count += 1
            else:
                fail_count += 1
            if done_count % 100 == 0 or done_count == total:
                elapsed = time.time() - t0
                rate = done_count / elapsed if elapsed > 0 else 0
                eta = (total - done_count) / rate if rate > 0 else 0
                print(f"  {done_count}/{total}  OK={ok_count}  "
                      f"FAIL={fail_count}  {rate:.1f}/s  ETA={eta:.0f}s")

    t1 = time.time()
    print(f"\n  Done in {t1 - t0:.1f}s  ({ok_count} OK, {fail_count} FAIL)")

    # Collect gene trees
    print(f"\n--- Collecting gene trees ---")
    all_trees = []
    for fasta_path in fasta_files:
        base = os.path.basename(fasta_path).replace(".fasta", "")
        nwk_path = os.path.join(nwk_dir, base + ".nwk")
        if os.path.isfile(nwk_path):
            with open(nwk_path) as f:
                tree = f.read().strip()
                if tree:
                    all_trees.append(tree)

    gene_trees_path = os.path.join(astral_outdir, "gene_trees.nwk")
    with open(gene_trees_path, "w") as f:
        for t in all_trees:
            f.write(t + "\n")
    print(f"Collected {len(all_trees)} gene trees -> {gene_trees_path}")

    # ASTRAL
    version = "IV" if use_astral_iv else "III"
    print(f"\n--- Running ASTRAL {version} ---")
    species_tree_path = os.path.join(astral_outdir, "species_tree.nwk")
    if len(all_trees) < 2:
        print("Too few gene trees for ASTRAL (need >= 2)")
    else:
        t0 = time.time()
        if use_astral_iv:
            ok = run_astral_iv(astral_bin, gene_trees_path, species_tree_path, args.threads)
        else:
            ok = run_astral_iii(astral_jar, gene_trees_path, species_tree_path)
        t1 = time.time()
        if ok:
            with open(species_tree_path) as f:
                tree = f.read().strip()
            print(f"  ASTRAL {version} done in {t1 - t0:.1f}s")
            print(f"\nSpecies tree:\n{tree}")
        else:
            print(f"ASTRAL {version} failed!")

    print("\nDone!")


if __name__ == "__main__":
    main()
