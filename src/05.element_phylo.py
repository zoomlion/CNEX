#!/usr/bin/env python3
"""
Step 05: Per-element alignment (FAMSA) + phylogeny.

Methods:
  --method concat (default): align → concat_msa → script → IQ-TREE 3
  --method astral           : align → FastTree → script → ASTRAL
  --method both             : run concat + astral in a single command

Tags + thresholds: loops over tag × threshold combos, generates run scripts.
  --submit (default): auto-execute scripts after generation.
  --dry-run: generate scripts only (for manual cluster submission).

Usage:
  python3 src/05.element_phylo.py --max-elements 0 --method concat
"""

import argparse
import os
import shutil
import subprocess
import sys
import glob
import time
from collections import Counter, OrderedDict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "utils"))
import config as C
from concat_msa import trim_alignment_by_occupancy, read_fasta


# ─── helpers ───────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="Element phylogeny: FAMSA + concat/IQ-TREE or ASTRAL")
    p.add_argument("--elements-dir", default="results/fasta",
                   help="Directory with per-element FASTAs")
    p.add_argument("--method", default=C.DEFAULT_METHOD, choices=["concat", "astral", "both"])
    p.add_argument("--famsa", default=C.FAMSA)
    p.add_argument("--fasttree", default=C.FastTree)
    p.add_argument("--iqtree3", default=C.IQTREE3)
    p.add_argument("--astral-bin", default=C.ASTRAL)
    p.add_argument("--astral-dir", default="")
    p.add_argument("--astral-jar", default="")
    p.add_argument("--element-tags", default=C.ELEMENT_TAGS_FILE,
                   help="TSV (ele_id\\ttag); empty = no tags")
    p.add_argument("--concat-length-quantiles", type=str, default=C.CONCAT_LENGTH_QUANTILES)
    p.add_argument("--block-gap", type=str, default=C.ASTRAL_BLOCK_GAPS,
                   help="Comma-separated kb thresholds for astral block clustering (default: config)")
    p.add_argument("--min-cne-per-species", type=int, default=C.MIN_CNE_PER_SPECIES)
    p.add_argument("--min-occupancy", type=float, default=0.3)
    p.add_argument("--min-site-occupancy", type=float, default=0.5)
    p.add_argument("-t", "--threads", type=int, default=C.THREADS,
                   help="Worker count: FAMSA/FastTree parallelism, IQ-TREE/ASTRAL threads")
    p.add_argument("--max-elements", type=int, default=100)
    p.add_argument("--resume", action="store_true", default=True)
    p.add_argument("--no-resume", action="store_false", dest="resume")
    p.add_argument("--submit", dest="dry_run", action="store_false",
                   default=C.DRY_RUN, help="Auto-execute run scripts (default)")
    p.add_argument("--dry-run", dest="dry_run", action="store_true",
                   default=C.DRY_RUN, help="Only generate scripts, no execution")
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


def read_element_tags(path):
    """Return (tag_dict, coord_dict).

    tag_dict:  {type_name: set(ele_id)}  — for tag-based filtering
    coord_dict: {ele_id: (type, chr, start, end)} — for gap clustering

    Format: ele_id \\t type \\t chr \\t start \\t end
    """
    tags = {}
    coords = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('ele_id'):
                continue
            parts = line.split('\t')
            try:
                ele = int(parts[0])
            except (ValueError, IndexError):
                continue
            if len(parts) >= 2:
                tag = parts[1].strip()
                tags.setdefault(tag, set()).add(ele)
            if len(parts) >= 5:
                try:
                    coords[ele] = (parts[1], parts[2], int(parts[3]), int(parts[4]))
                except ValueError:
                    pass
    return tags, coords


def filter_species_by_cne_count(fasta_files, min_cne):
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
    fasta_path, famsa_bin, resume, aln_dir, keep_sp, min_site_occ = args
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
    # local trim: remove gap-rich columns
    if min_site_occ > 0:
        sq = read_fasta(aln_path)
        trimmed = trim_alignment_by_occupancy(sq, min_site_occ)
        trim_path = aln_path.replace(".aln", ".trimmed.aln")
        with open(trim_path, 'w') as f:
            for h, s in trimmed.items():
                f.write(f'>{h}\n{s}\n')
    return True, fasta_path, None


def _fasttree_one(args):
    fp, aln_dir, out_dir, fasttree_bin = args[:4]
    base_name = os.path.basename(fp).replace(".fasta", "")
    nwk_path = os.path.join(out_dir, "nwk", base_name + ".nwk")
    fasta_input = os.path.join(aln_dir, base_name + ".trimmed.aln")
    if not os.path.isfile(fasta_input):
        fasta_input = os.path.join(aln_dir, base_name + ".aln")
        if not os.path.isfile(fasta_input):
            return False, base_name
    with open(fasta_input) as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">"))[:100]
    is_nt = all(c in "ACGTacgtNn-" for c in seq) if seq else True
    cmd = [fasttree_bin]
    if is_nt:
        cmd.append("-nt")
    cmd.extend(["-gtr", "-nosupport"])
    with open(fasta_input) as inp, open(nwk_path, "w") as out:
        r = subprocess.run(cmd, stdin=inp, stdout=out, stderr=subprocess.PIPE, text=True)
    return r.returncode == 0, base_name


def _fasttree_cluster(args):
    fasttree_bin, fasta_path, nwk_path = args
    if not os.path.isfile(fasta_path):
        return False, fasta_path
    cmd = [fasttree_bin, "-nt", "-gtr", "-nosupport"]
    with open(fasta_path) as inp, open(nwk_path, "w") as out:
        r = subprocess.run(cmd, stdin=inp, stdout=out, stderr=subprocess.PIPE, text=True)
    return r.returncode == 0, nwk_path


def compute_aligned_lengths(aln_dir, fasta_files):
    """Yield average non-gap length per element (aligned sequences)."""
    for fp in fasta_files:
        base = os.path.basename(fp).replace(".fasta", "")
        ap = os.path.join(aln_dir, base + ".aln")
        if not os.path.isfile(ap):
            continue
        sq = read_fasta(ap)
        if not sq:
            continue
        lens = [sum(1 for c in s if c not in '-Nn?') for s in sq.values()]
        yield sum(lens) / len(lens)


def quantile_cutoffs(values, quantiles):
    """Compute cutoffs for given quantile percentages (0-100)."""
    sv = sorted(values)
    n = len(sv)
    cuts = {}
    for q in quantiles:
        idx = max(0, min(n - 1, int(n * q / 100)))
        cuts[q] = sv[idx]
    return cuts


def bin_cluster(ele_ids, coords, bin_size):
    """Assign element IDs to fixed genomic bins of bin_size bp.

    Bins with fewer than 2 elements are discarded.
    """
    bins = {}
    for eid in ele_ids:
        if eid not in coords:
            continue
        _, ch, st, _ = coords[eid]
        bin_start = (st // bin_size) * bin_size
        bins.setdefault(f"{ch}_{bin_start}", []).append(eid)
    return {k: v for k, v in bins.items() if len(v) >= 2}


def concat_block_alignments(member_ids, aln_dir):
    """Concatenate trimmed alignments of multiple elements into one supermatrix.

    Pads missing species with gaps so all sequences have equal length.
    """
    block_alignments = []
    all_sp = set()
    for eid in member_ids:
        ap = os.path.join(aln_dir, f"{eid}.trimmed.aln")
        if not os.path.isfile(ap):
            ap = os.path.join(aln_dir, f"{eid}.aln")
        if not os.path.isfile(ap):
            continue
        sq = read_fasta(ap)
        block_alignments.append(sq)
        all_sp.update(sq.keys())

    merged = {sp: [] for sp in all_sp}
    for sq in block_alignments:
        alen = len(next(iter(sq.values()))) if sq else 0
        for sp in all_sp:
            seq = sq.get(sp, "")
            merged[sp].append(seq if seq else "-" * alen)

    return {sp: "".join(seqs) for sp, seqs in merged.items()}


def write_script(path, cmds, description=""):
    """Write a run.sh script, return path."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# " + description + "\n")
        f.write('cd "$(dirname "$0")"\n')
        f.write('\n'.join(cmds) + '\n')
    os.chmod(path, 0o755)
    return path


def execute_script(path, submit):
    """Run script if submit mode, else print info."""
    if submit:
        print(f"  Running: {path}")
        subprocess.run(["bash", path])
    else:
        print(f"  Script:  {path}")


# ─── concat pipeline ───────────────────────────────────────

def run_concat_subset(aln_dir, keep_fastas, min_occ, out_dir,
                      iqtree3, iqtree_threads, submit, tag_name, thr_label,
                      partition=False):
    """Filter alignments, run concat_msa, write + execute run script."""
    if not keep_fastas:
        return
    aln_subdir = os.path.join(out_dir, "aln_filtered")
    os.makedirs(aln_subdir, exist_ok=True)

    # symlink filtered alignments
    for fp in keep_fastas:
        base_name = os.path.basename(fp).replace(".fasta", "")
        src = os.path.join(aln_dir, base_name + ".trimmed.aln")
        if not os.path.isfile(src):
            src = os.path.join(aln_dir, base_name + ".aln")
        dst = os.path.join(aln_subdir, base_name + ".aln")
        if os.path.isfile(src) and not os.path.isfile(dst):
            os.symlink(os.path.abspath(src), dst)

    # concat_msa.py
    concat_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "utils", "concat_msa.py")
    supermatrix = os.path.join(out_dir, "supermatrix.fa")
    partitions = os.path.join(out_dir, "partitions.txt")
    cmd = ["python3", concat_script, "-i", aln_subdir, "--suffix", ".aln",
           "--min-occupancy", str(min_occ), "-o", supermatrix]
    if partition:
        cmd += ["-p", partitions]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(supermatrix):
        print(f"  concat_msa.py failed for {tag_name}/{thr_label}")
        return
    # print key lines from concat_msa.py output
    for line in r.stdout.strip().split('\n'):
        if any(k in line for k in ('Total species', 'Occupancy', 'Supermatrix', 'Partitions')):
            print(f"  {line.strip()}")

    # write run.sh
    script_path = os.path.join(out_dir, "run.sh")
    n_threads = iqtree_threads
    iq_cmd = f"iqtree3 -s {os.path.basename(supermatrix)}"
    if partition:
        iq_cmd += f" -p {os.path.basename(partitions)}"
    iq_cmd += f" -m GTR+F+R4 -bb 1000 -nt {n_threads}"
    cmds = [iq_cmd]
    write_script(script_path, cmds,
                 f"IQ-TREE 3: {tag_name} / {thr_label} ({len(keep_fastas)} elements)")
    execute_script(script_path, submit)
    return script_path


# ─── astral pipeline ───────────────────────────────────────

def run_astral_subset(aln_dir, keep_fastas, out_dir, fasttree_bin,
                       astral_bin, astral_jar, astral_threads, min_site_occ,
                       submit, tag_name, thr_label, threads=1,
                       tag_coords=None, block_gap=0):
    """ASTRAL pipeline: block-gap cluster → concat → FastTree → ASTRAL.

    If block_gap > 0 and tag_coords is provided, elements are clustered by
    genomic proximity before concatenation and tree inference.
    """
    if not keep_fastas:
        return

    # detect ASTRAL mode
    use_iv = os.path.isfile(astral_bin) if astral_bin else False
    jar_path = astral_jar

    nwk_dir = os.path.join(out_dir, "nwk")
    os.makedirs(nwk_dir, exist_ok=True)

    # Determine bins: fixed-size bins or per-element
    if block_gap > 0 and tag_coords:
        ele_ids = sorted(int(os.path.basename(fp).replace(".fasta", "")) for fp in keep_fastas)
        bin_items = bin_cluster(ele_ids, tag_coords, block_gap)
        cluster_items = list(bin_items.items())
        print(f"  Bins ({block_gap}bp): {len(cluster_items)} blocks ({tag_name}/{thr_label})")
    else:
        # per-element clustering (each element is its own "cluster")
        cluster_items = [(int(os.path.basename(fp).replace(".fasta", "")), [fp])
                         for fp in keep_fastas]

    # Prepare FastTree work items: (fasta_path, nwk_path) for each cluster
    work = []
    for block_id, members in cluster_items:
        if block_gap > 0 and tag_coords:
            concat = concat_block_alignments(members, aln_dir)
            if len(concat) < 2:
                continue
            base_name = f"{tag_name}_{block_id}"
            fasta_path = os.path.join(out_dir, "nwk", base_name + ".fa")
            with open(fasta_path, 'w') as f:
                for sp, s in sorted(concat.items()):
                    f.write(f">{sp}\n{s}\n")
        else:
            fp = (list(members)[0]) if isinstance(members, (set, list)) else members
            base_name = os.path.basename(fp).replace(".fasta", "") if isinstance(fp, str) else str(members)
            fasta_path = os.path.join(aln_dir, base_name + ".trimmed.aln")
            if not os.path.isfile(fasta_path):
                fasta_path = os.path.join(aln_dir, base_name + ".aln")
                if not os.path.isfile(fasta_path):
                    continue

        nwk_path = os.path.join(nwk_dir, base_name + ".nwk")
        if os.path.isfile(nwk_path):
            continue
        work.append((fasttree_bin, fasta_path, nwk_path))

    # Run FastTree in parallel
    ok_count = fail_count = 0
    if threads > 1:
        with Pool(threads) as p:
            for ok, _ in p.imap_unordered(_fasttree_cluster, work):
                if ok: ok_count += 1
                else: fail_count += 1
    else:
        for w in work:
            ok, _ = _fasttree_cluster(w)
            if ok: ok_count += 1
            else: fail_count += 1
    print(f"  FastTree: {ok_count} OK, {fail_count} FAIL ({tag_name}/{thr_label})")

    # collect gene trees
    all_trees = []
    if block_gap > 0 and tag_coords:
        for block_id, _ in cluster_items:
            base_name = f"{tag_name}_{block_id}"
            nwk_path = os.path.join(nwk_dir, base_name + ".nwk")
            if os.path.isfile(nwk_path):
                with open(nwk_path) as f:
                    t = f.read().strip()
                    if t:
                        all_trees.append(t)
    else:
        for fp in keep_fastas:
            base_name = os.path.basename(fp).replace(".fasta", "")
            nwk_path = os.path.join(nwk_dir, base_name + ".nwk")
            if os.path.isfile(nwk_path):
                with open(nwk_path) as f:
                    t = f.read().strip()
                    if t:
                        all_trees.append(t)
    gene_trees_path = os.path.join(out_dir, "gene_trees.nwk")
    with open(gene_trees_path, 'w') as f:
        for t in all_trees:
            f.write(t + '\n')
    print(f"  Collected {len(all_trees)} gene trees ({tag_name}/{thr_label})")

    if len(all_trees) < 2:
        print("  Too few trees for ASTRAL")
        return

    # write ASTRAL run script
    species_tree_path = os.path.join(out_dir, "species_tree.nwk")
    script_path = os.path.join(out_dir, "run.sh")
    if use_iv:
        cmd = f"{astral_bin} -i {os.path.basename(gene_trees_path)} -o {os.path.basename(species_tree_path)} -t {astral_threads}"
    else:
        cmd = f"java -Xmx8g -jar {jar_path} -i {os.path.basename(gene_trees_path)} -o {os.path.basename(species_tree_path)} --extraLevel 0"
    cmds = [cmd]
    write_script(script_path, cmds,
                 f"ASTRAL: {tag_name} / {thr_label} ({len(all_trees)} gene trees)")
    execute_script(script_path, submit)
    return script_path


# ─── main ──────────────────────────────────────────────────

def main():
    args = parse_args()
    famsa = os.path.expanduser(args.famsa)
    if not shutil.which(famsa):
        sys.exit(f"FAMSA not found: {famsa}")

    print(f"Method:   {args.method}")
    print(f"Dry-run:  {args.dry_run}")
    print()
    submit = not args.dry_run  # --dry-run inverts submit

    base_dir = os.path.dirname(os.path.abspath(args.elements_dir))
    aln_dir = os.path.join(base_dir, "aln")
    os.makedirs(aln_dir, exist_ok=True)

    fasta_files = find_fasta_files(args.elements_dir, args.max_elements)
    print(f"Found {len(fasta_files)} element FASTAs")
    if not fasta_files:
        sys.exit("No element FASTAs found. Run step 04 first.")

    keep_sp = filter_species_by_cne_count(fasta_files, args.min_cne_per_species)

    # ─── FAMSA alignment ─────────────────────────────────
    print(f"\n--- Aligning elements (FAMSA) ---")
    t0 = time.time()
    work = [(f, famsa, args.resume, aln_dir, keep_sp, args.min_site_occupancy) for f in fasta_files]
    ok_count = fail_count = 0
    if args.threads > 1:
        with Pool(args.threads) as p:
            for ok, _, _ in p.imap_unordered(_align_one, work):
                if ok: ok_count += 1
                else: fail_count += 1
    else:
        for w in work:
            ok, _, _ = _align_one(w)
            if ok: ok_count += 1
            else: fail_count += 1
    print(f"  Done in {time.time() - t0:.1f}s ({ok_count} OK, {fail_count} FAIL)")

    # ─── Tags ────────────────────────────────────────────
    if args.element_tags and os.path.isfile(args.element_tags):
        tag_dict, tag_coords = read_element_tags(args.element_tags)
        print(f"Tags: {len(tag_dict)} groups: {', '.join(sorted(tag_dict))}")
    elif args.element_tags and not os.path.isfile(args.element_tags):
        print(f"  Warning: --element-tags ({args.element_tags}) not found.")
        print(f"  Copy the example:  cp element_tags.example.tsv element_tags.tsv")
        print(f"  Or leave ELEMENT_TAGS_FILE empty in config.py for no tags.")
        tag_dict, tag_coords = {}, {}
    else:
        tag_dict, tag_coords = {}, {}
    tag_dict["all"] = None  # always include full set

    # ─── Thresholds ──────────────────────────────────────
    concat_quantiles = [int(x) for x in args.concat_length_quantiles.split(",") if x.strip()] \
                       if args.concat_length_quantiles.strip() else []
    astral_block_gaps = [int(x) for x in args.block_gap.split(",") if x.strip()] \
                        if args.block_gap.strip() else []

    concat_levels = [None] + (concat_quantiles if concat_quantiles else [])
    astral_levels = [None] + (astral_block_gaps if astral_block_gaps else [])

    # ─── Method-dispatch ─────────────────────────────────
    methods = ["concat", "astral"] if args.method == "both" else [args.method]

    for m in methods:
        all_scripts = []

        if m == "concat":
            iqtree3 = os.path.expanduser(args.iqtree3)
            if not shutil.which(iqtree3):
                sys.exit(f"IQ-TREE 3 not found: {iqtree3}")
            base_out = os.path.join(base_dir, "iqtree")
            levels = concat_levels
            quantiles = concat_quantiles
            method_tag_dict = {"all": None}  # concat: no tag splitting
        else:
            fasttree = os.path.expanduser(args.fasttree)
            if not shutil.which(fasttree):
                sys.exit(f"FastTree not found: {fasttree}")
            astral_bin = os.path.expanduser(args.astral_bin)
            base_out = os.path.join(base_dir, "astral")
            levels = astral_levels
            method_tag_dict = {k: v for k, v in tag_dict.items() if k != "all"}
            if not method_tag_dict:
                method_tag_dict = {"all": None}  # fallback when no tags

        for tag_name, tag_ids in sorted(method_tag_dict.items()):
            tag_files = fasta_files
            if tag_ids is not None:
                tag_files = [f for f in fasta_files
                             if int(os.path.basename(f).replace(".fasta", "").split('.')[0]) in tag_ids]
            if len(tag_files) < 3:
                print(f"  Skip tag '{tag_name}': only {len(tag_files)} elements")
                continue

            if m == "concat":
                tag_lens = list(compute_aligned_lengths(aln_dir, tag_files))
                tag_cuts = quantile_cutoffs(tag_lens, concat_quantiles) if concat_quantiles else {}
            else:
                tag_cuts = {}

            for level in levels:
                if m == "concat":
                    thr_label = f"quantile_{level}" if level is not None else "all"
                else:
                    thr_label = f"block_gap_{level}" if level is not None else "all"
                print(f"\n=== {m.upper()}: {tag_name} / {thr_label} ===")
                out_dir = os.path.join(base_out, tag_name, thr_label)
                os.makedirs(out_dir, exist_ok=True)

                if level is not None and m == "concat":
                    cutoff = tag_cuts.get(level, 0)
                    keep = []
                    for fp in tag_files:
                        b = os.path.basename(fp).replace(".fasta", "")
                        ap = os.path.join(aln_dir, b + ".aln")
                        if os.path.isfile(ap):
                            sq = read_fasta(ap)
                            avg = sum(sum(1 for c in s if c not in '-Nn?') for s in sq.values()) / len(sq)
                            if avg >= cutoff:
                                keep.append(fp)
                    if len(keep) < 3:
                        print(f"  Skip {tag_name}/{thr_label}: {len(keep)} elements")
                        continue
                    print(f"  {tag_name}/{thr_label}: {len(keep)}/{len(tag_files)} elements (≥{cutoff:.0f}bp, P{level})")
                else:
                    keep = tag_files
                    print(f"  {tag_name}/{thr_label}: {len(keep)} elements")

                if m == "concat":
                    sp = run_concat_subset(aln_dir, keep, args.min_occupancy,
                                           out_dir, iqtree3, args.threads, submit,
                                           tag_name, thr_label,
                                           C.PARTITION)
                else:
                    block_gap = level * 1000 if level else 0
                    sp = run_astral_subset(aln_dir, keep, out_dir, fasttree,
                                           astral_bin, args.astral_jar,
                                           args.threads, args.min_site_occupancy,
                                           submit, tag_name, thr_label,
                                           args.threads, tag_coords, block_gap)
                if sp:
                    all_scripts.append(sp)

        # write run_all.sh
        all_scr_path = os.path.join(base_out, "run_all.sh")
        cmds = [f"bash {os.path.relpath(s, base_out)}" for s in all_scripts]
        write_script(all_scr_path, cmds, f"Run all {m.upper()} scripts")

    print("\nDone!")


if __name__ == "__main__":
    main()
