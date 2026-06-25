#!/usr/bin/env python3
"""
04.alignfree_phylo.py - Align-free phylogenetic tree construction
from assembled CNE sequences using k-mer presence/absence.

Features:
  1. Canonical k-mers (reverse-complement aware)
  2. Mash distance & Containment distance (linear evolutionary metrics)
  3. Shannon entropy filtering (k=13 LUT + dynamic fallback)
  4. C++ accelerated bitmap operations via hip.popcount
"""

import argparse
import os
import sys
import time
from collections import defaultdict

import numpy as np
from multiprocessing import Pool

BASE2BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BITS_PER_BASE = 2
DEFAULT_MIN_ENTROPY = 1.4

import math as _math

_ENTROPY_LUT = {}
for _a in range(14):
    for _c in range(14 - _a):
        for _g in range(14 - _a - _c):
            _t = 13 - _a - _c - _g
            _key = tuple(sorted([_a, _c, _g, _t]))
            if _key not in _ENTROPY_LUT:
                _s = 0.0
                for _cnt in _key:
                    if _cnt > 0:
                        _p = _cnt / 13.0
                        _s -= _p * _math.log2(_p)
                _ENTROPY_LUT[_key] = _s

_SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))

try:
    sys.path.insert(0, _SRC_DIR)
    from hip import popcount as _hpc
    _CPP_AVAILABLE = True
except (ImportError, ModuleNotFoundError):
    sys.path.insert(0, os.path.join(_SRC_DIR, '..', 'src'))
    try:
        from hip import popcount as _hpc
        _CPP_AVAILABLE = True
    except (ImportError, ModuleNotFoundError):
        _hpc = None
        _CPP_AVAILABLE = False

_kmer_size = 13


def _entropy(a, c, g, t, k):
    if k == 13:
        return _ENTROPY_LUT.get(tuple(sorted([a, c, g, t])), 0.0)
    s = 0.0
    for cnt in (a, c, g, t):
        if cnt > 0:
            p = cnt / float(k)
            s -= p * _math.log2(p)
    return s


def encode_kmer(seq, k, min_entropy=None):
    code = 0
    a = c = g = t = 0
    for base in seq:
        b = BASE2BITS.get(base)
        if b is None:
            return None
        code = (code << BITS_PER_BASE) | b
        if min_entropy is not None:
            if base == 'A':
                a += 1
            elif base == 'C':
                c += 1
            elif base == 'G':
                g += 1
            else:
                t += 1
    if min_entropy is not None:
        if _entropy(a, c, g, t, k) < min_entropy:
            return None

    rc = 0
    tmp = code
    for _ in range(k):
        rc = (rc << BITS_PER_BASE) | ((tmp & 3) ^ 3)
        tmp >>= BITS_PER_BASE
    return code if code < rc else rc


def extract_data(fasta_path, k, with_contigs=False, min_entropy=None):
    global_codes = set()
    contig_codes_list = [] if with_contigs else None
    max_code = 4 ** k

    def process_seq(seq):
        if len(seq) < k:
            return
        seq = seq.upper()
        local_codes = []
        for i in range(len(seq) - k + 1):
            code = encode_kmer(seq[i:i + k], k, min_entropy)
            if code is not None:
                local_codes.append(code)
        if local_codes:
            if with_contigs:
                contig_codes_list.append(local_codes)
            global_codes.update(local_codes)

    with open(fasta_path) as f:
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_parts:
                    process_seq(''.join(seq_parts))
                    seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            process_seq(''.join(seq_parts))

    if _CPP_AVAILABLE:
        bitmap = _hpc.build_bitmap(list(global_codes), max_code)
    else:
        bitmap = _build_bitmap_py(global_codes, max_code)

    if with_contigs:
        return bitmap, len(global_codes), contig_codes_list
    return bitmap, len(global_codes)


def extract_freq_data(fasta_path, k, with_contigs=False, min_entropy=None):
    global_freq = {}
    contig_freqs = [] if with_contigs else None

    def process_seq(seq):
        if len(seq) < k:
            return
        seq = seq.upper()
        local_freq = {}
        for i in range(len(seq) - k + 1):
            code = encode_kmer(seq[i:i + k], k, min_entropy)
            if code is not None:
                local_freq[code] = local_freq.get(code, 0) + 1
        if local_freq:
            if with_contigs:
                contig_freqs.append(local_freq)
            for code, cnt in local_freq.items():
                global_freq[code] = global_freq.get(code, 0) + cnt

    with open(fasta_path) as f:
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seq_parts:
                    process_seq(''.join(seq_parts))
                    seq_parts = []
            else:
                seq_parts.append(line)
        if seq_parts:
            process_seq(''.join(seq_parts))

    n_kmers = len(global_freq)
    if with_contigs:
        return global_freq, n_kmers, contig_freqs
    return global_freq, n_kmers


def _build_bitmap_py(codes, max_code):
    n_bytes = (max_code + 7) // 8
    ba = bytearray(n_bytes)
    for code in codes:
        ba[code >> 3] |= (1 << (code & 7))
    return int.from_bytes(ba, 'little')


def _bit_count(x):
    """Safe popcount for Python < 3.10."""
    if hasattr(x, 'bit_count'):
        return x.bit_count()
    return int.bit_count(x)


def jaccard_distance(a, b):
    if _CPP_AVAILABLE:
        return _hpc.jaccard_distance(a, b)
    inter = _bit_count(a & b)
    union = _bit_count(a | b)
    if union == 0:
        return 1.0
    return 1.0 - inter / union


def bray_curtis_distance(a, b):
    if _CPP_AVAILABLE:
        return _hpc.bray_curtis_distance(a, b)
    cnt_a = _bit_count(a)
    cnt_b = _bit_count(b)
    inter = _bit_count(a & b)
    s = cnt_a + cnt_b
    if s == 0:
        return 1.0
    return 1.0 - 2.0 * inter / s


def mash_distance(a, b):
    if _CPP_AVAILABLE:
        return _hpc.mash_distance(a, b, _kmer_size)
    inter = _bit_count(a & b)
    cnt_a = _bit_count(a)
    cnt_b = _bit_count(b)
    union = cnt_a + cnt_b - inter
    if union == 0:
        return 1.0
    J = inter / union
    val = 2.0 * J / (1.0 + J)
    if val <= 0.0:
        return 1.0
    return min(1.0, -_math.log(val) / _kmer_size)


def containment_distance(a, b):
    if _CPP_AVAILABLE:
        return _hpc.containment_distance(a, b, _kmer_size)
    cnt_a = _bit_count(a)
    cnt_b = _bit_count(b)
    inter = _bit_count(a & b)
    mn = cnt_a if cnt_a < cnt_b else cnt_b
    if mn == 0 or inter == 0:
        return 1.0
    C = inter / mn
    return min(1.0, -_math.log(C) / _kmer_size)


def cosine_distance(a, b):
    if len(a) > len(b):
        a, b = b, a
    dot = 0
    for k, v in a.items():
        dot += v * b.get(k, 0)
    na = _math.sqrt(sum(v * v for v in a.values()))
    nb = _math.sqrt(sum(v * v for v in b.values()))
    if na == 0 or nb == 0:
        return 1.0
    return 1.0 - dot / (na * nb)


_DISTANCE_FUNCTIONS = {
    'jaccard': jaccard_distance,
    'bray-curtis': bray_curtis_distance,
    'mash': mash_distance,
    'containment': containment_distance,
    'cosine': cosine_distance,
}

_dist_func = mash_distance


def _set_dist_func(name):
    global _dist_func
    _dist_func = _DISTANCE_FUNCTIONS[name]


def _compute_row(args):
    i, bitmaps = args
    n = len(bitmaps)
    row = np.zeros(n, dtype=np.float64)
    for j in range(i + 1, n):
        row[j] = _dist_func(bitmaps[i], bitmaps[j])
    return i, row


def build_distance_matrix(bitmaps, n_jobs=1):
    n = len(bitmaps)
    dist = np.zeros((n, n), dtype=np.float64)

    if n_jobs <= 1 or n <= 3:
        for i in range(n):
            for j in range(i + 1, n):
                d = _dist_func(bitmaps[i], bitmaps[j])
                dist[i][j] = d
                dist[j][i] = d
    else:
        tasks = [(i, bitmaps) for i in range(n - 1)]
        with Pool(min(n_jobs, n - 1)) as pool:
            for i, row in pool.imap_unordered(_compute_row, tasks):
                for j in range(i + 1, n):
                    dist[i][j] = row[j]
                    dist[j][i] = row[j]

    return dist


def _get_leaves(node, children):
    if node not in children:
        return [node]
    result = []
    for child_name, _ in children[node]:
        result.extend(_get_leaves(child_name, children))
    return result


def neighbor_joining(dist, labels):
    n = len(labels)
    if n == 1:
        return f"{labels[0]};", {}, {}, labels[0]

    nodes = list(labels)
    D = dist.copy()
    children = {}
    next_id = 0

    while len(nodes) > 2:
        m = len(nodes)
        r = np.sum(D, axis=1) / (m - 2)

        q_min = float('inf')
        best_i = best_j = 0
        best_d = 0.0
        eps = 1e-14
        for i in range(m):
            for j in range(i + 1, m):
                q = D[i][j] - r[i] - r[j]
                if q < q_min - eps:
                    q_min = q
                    best_i, best_j = i, j
                    best_d = D[i][j]
                elif abs(q - q_min) <= eps and D[i][j] < best_d:
                    best_i, best_j = i, j
                    best_d = D[i][j]

        d_i = max(0.0, (D[best_i][best_j] + r[best_i] - r[best_j]) / 2)
        d_j = max(0.0, (D[best_i][best_j] + r[best_j] - r[best_i]) / 2)

        new_name = f"node{next_id}"
        next_id += 1
        children[new_name] = [(nodes[best_i], d_i), (nodes[best_j], d_j)]

        keep = [k for k in range(m) if k != best_i and k != best_j]
        new_m = len(keep) + 1
        new_D = np.zeros((new_m, new_m), dtype=np.float64)
        new_nodes = [new_name] + [nodes[k] for k in keep]

        for ki, k in enumerate(keep):
            d_k = max(0.0, (D[k][best_i] + D[k][best_j] - D[best_i][best_j]) / 2)
            new_D[0][ki + 1] = d_k
            new_D[ki + 1][0] = d_k
            for kj, l in enumerate(keep):
                if l > k:
                    new_D[ki + 1][kj + 1] = D[k][l]
                    new_D[kj + 1][ki + 1] = D[k][l]

        nodes = new_nodes
        D = new_D

    if len(nodes) == 2:
        d = D[0][1]
        root = f"node{next_id}"
        children[root] = [(nodes[0], d / 2), (nodes[1], d / 2)]
    else:
        root = nodes[0]

    clades = {}
    for node in children:
        leaves = frozenset(_get_leaves(node, children))
        if len(leaves) > 1 and leaves != frozenset(labels):
            clades[node] = leaves

    return _to_newick(root, children) + ";", clades, children, root


def _to_newick(node, children):
    if node not in children:
        return node
    parts = []
    for child_name, branch_len in children[node]:
        child_str = _to_newick(child_name, children)
        parts.append(f"{child_str}:{branch_len:.6f}")
    return f"({', '.join(parts)})"


def _to_newick_with_support(node, children, support):
    if node not in children:
        return node
    parts = []
    for child_name, branch_len in children[node]:
        child_str = _to_newick_with_support(child_name, children, support)
        parts.append(f"{child_str}:{branch_len:.6f}")
    result = f"({', '.join(parts)})"
    if node in support:
        result += f"{support[node]:.0f}"
    return result


def _bitmap_from_contig_codes(contig_codes_list, indices, max_code):
    if _CPP_AVAILABLE:
        return _hpc.build_bitmap_from_contigs(contig_codes_list, indices, max_code)
    merged = set()
    for idx in indices:
        merged.update(contig_codes_list[idx])
    return _build_bitmap_py(merged, max_code)


def bootstrap_phylo(all_contig_codes, labels, k, n_jobs, n_bootstrap):
    max_code = 4 ** k

    ref_bitmaps = []
    for contig_codes_list in all_contig_codes:
        n = len(contig_codes_list)
        if n == 0:
            if _CPP_AVAILABLE:
                ref_bitmaps.append(_hpc.build_bitmap([], max_code))
            else:
                ref_bitmaps.append(_build_bitmap_py(set(), max_code))
        else:
            ref_bitmaps.append(_bitmap_from_contig_codes(
                contig_codes_list, list(range(n)), max_code))

    ref_dist = build_distance_matrix(ref_bitmaps, n_jobs)
    _, ref_clades, ref_children, ref_root = neighbor_joining(ref_dist, labels)

    clade_to_node = {}
    for node, clade in ref_clades.items():
        clade_to_node[clade] = node

    if not clade_to_node:
        print("Warning: no internal clades (need >= 4 species)")
        ref_newick = _to_newick(ref_root, ref_children) + ";"
        return ref_newick, {}, ref_children, ref_root

    clade_counts = defaultdict(int)
    for rep in range(n_bootstrap):
        if (rep + 1) % 100 == 0:
            print(f"\r  replicate {rep + 1}/{n_bootstrap}...", end="", flush=True)

        boot_bitmaps = []
        for contig_codes_list in all_contig_codes:
            n_contigs = len(contig_codes_list)
            if n_contigs == 0:
                if _CPP_AVAILABLE:
                    boot_bitmaps.append(_hpc.build_bitmap([], max_code))
                else:
                    boot_bitmaps.append(_build_bitmap_py(set(), max_code))
            else:
                indices = np.random.randint(0, n_contigs, n_contigs)
                boot_bitmaps.append(_bitmap_from_contig_codes(
                    contig_codes_list, indices, max_code))

        boot_dist = build_distance_matrix(boot_bitmaps, max(1, n_jobs // 2))
        _, boot_clades, _, _ = neighbor_joining(boot_dist, labels)

        for clade in set(boot_clades.values()):
            if clade in clade_to_node:
                clade_counts[clade] += 1

    if n_bootstrap >= 100:
        print(f"\r  {n_bootstrap}/{n_bootstrap} replicates done          ")

    support = {}
    for clade, node in clade_to_node.items():
        support[node] = clade_counts[clade] / n_bootstrap * 100

    newick_supp = _to_newick_with_support(ref_root, ref_children, support) + ";"
    return newick_supp, support, ref_children, ref_root


def bootstrap_phylo_cosine(all_contig_freqs, labels, n_jobs, n_bootstrap):
    ref_freqs = []
    for contig_freq_list in all_contig_freqs:
        merged = {}
        for cf in contig_freq_list:
            for code, cnt in cf.items():
                merged[code] = merged.get(code, 0) + cnt
        ref_freqs.append(merged)

    ref_dist = build_distance_matrix(ref_freqs, n_jobs)
    _, ref_clades, ref_children, ref_root = neighbor_joining(ref_dist, labels)

    clade_to_node = {}
    for node, clade in ref_clades.items():
        clade_to_node[clade] = node

    if not clade_to_node:
        print("Warning: no internal clades (need >= 4 species)")
        ref_newick = _to_newick(ref_root, ref_children) + ";"
        return ref_newick, {}, ref_children, ref_root

    clade_counts = defaultdict(int)
    for rep in range(n_bootstrap):
        if (rep + 1) % 100 == 0:
            print(f"\r  replicate {rep + 1}/{n_bootstrap}...", end="", flush=True)

        boot_freqs = []
        for contig_freq_list in all_contig_freqs:
            n_contigs = len(contig_freq_list)
            if n_contigs == 0:
                boot_freqs.append({})
            else:
                indices = np.random.randint(0, n_contigs, n_contigs)
                merged = {}
                for idx in indices:
                    for code, cnt in contig_freq_list[idx].items():
                        merged[code] = merged.get(code, 0) + cnt
                boot_freqs.append(merged)

        boot_dist = build_distance_matrix(boot_freqs, max(1, n_jobs // 2))
        _, boot_clades, _, _ = neighbor_joining(boot_dist, labels)

        for clade in set(boot_clades.values()):
            if clade in clade_to_node:
                clade_counts[clade] += 1

    if n_bootstrap >= 100:
        print(f"\r  {n_bootstrap}/{n_bootstrap} replicates done          ")

    support = {}
    for clade, node in clade_to_node.items():
        support[node] = clade_counts[clade] / n_bootstrap * 100

    newick_supp = _to_newick_with_support(ref_root, ref_children, support) + ";"
    return newick_supp, support, ref_children, ref_root


def plot_heatmap(dist, labels, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_size = max(6, len(labels) * 0.7)
    fig, ax = plt.subplots(figsize=(fig_size + 2, fig_size))
    im = ax.imshow(dist, cmap='YlOrRd', aspect='auto')

    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_title('Evolutionary Distance Matrix', fontsize=11)

    cbar = fig.colorbar(im, ax=ax, shrink=0.85)
    cbar.set_label('Estimated Evolutionary Distance', fontsize=9)

    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{dist[i][j]:.4f}", ha='center', va='center',
                    fontsize=7, color='white' if dist[i][j] > 0.5 else 'black')

    plt.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main():
    global _kmer_size
    parser = argparse.ArgumentParser(
        description="Align-free phylogenetic tree from CNE sequences "
                    "(canonical k-mers, Mash distance, C++ accelerated)"
    )
    parser.add_argument("inputs", nargs="+",
                        help="Per-species assembled CNE fasta files (at least 2)")
    parser.add_argument("-k", "--kmer", type=int, default=13,
                        help="K-mer size (default: 13)")
    parser.add_argument("-o", "--output", default="tree.nwk",
                        help="Output Newick tree file (default: tree.nwk)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Thread count (default: 1)")
    parser.add_argument("--labels", nargs="+",
                        help="Species labels (default: basename of input files)")
    parser.add_argument("--dist", help="Save distance matrix as CSV")
    parser.add_argument("--plot", help="Save distance heatmap as PNG")
    parser.add_argument("-b", "--bootstrap", type=int, default=0,
                        help="Bootstrap replicates (default: 0 = no bootstrap)")
    parser.add_argument("-d", "--distance",
                        choices=["jaccard", "bray-curtis", "mash", "containment", "cosine"],
                        default="mash",
                        help="Distance metric (default: mash)")
    parser.add_argument("--min-entropy", type=float, default=DEFAULT_MIN_ENTROPY,
                        help=f"Minimum Shannon entropy for k-mers "
                             f"(default: {DEFAULT_MIN_ENTROPY}, set 0 to disable)")
    args = parser.parse_args()

    _kmer_size = args.kmer

    if args.min_entropy <= 0:
        min_entropy = None
    else:
        min_entropy = args.min_entropy

    _set_dist_func(args.distance)

    if len(args.inputs) < 2:
        print("Error: at least 2 species required", file=sys.stderr)
        sys.exit(1)
    if args.kmer < 4 or args.kmer > 20:
        print("Error: k-mer size should be between 4 and 20", file=sys.stderr)
        sys.exit(1)

    max_code = 4 ** args.kmer

    labels = args.labels if args.labels else [
        os.path.splitext(os.path.basename(p))[0] for p in args.inputs
    ]
    if len(labels) != len(args.inputs):
        print(f"Error: {len(labels)} labels for {len(args.inputs)} inputs", file=sys.stderr)
        sys.exit(1)

    print(f"K-mer size        : {args.kmer}bp")
    print(f"K-mer space       : 4^{args.kmer} = {max_code:,} "
          f"({args.kmer * BITS_PER_BASE} bits, {max_code / 8 / 1024 / 1024:.1f} MB)")
    print(f"Backend           : {'C++ (hip.popcount)' if _CPP_AVAILABLE else 'Python (slow)'}")
    print(f"Distance metric   : {args.distance}")
    print(f"Canonical k-mers  : yes")
    if min_entropy is not None:
        print(f"Min entropy       : {min_entropy}")
    else:
        print(f"Min entropy       : off")
    print(f"Species count     : {len(args.inputs)}")
    if args.bootstrap:
        print(f"Bootstrap reps    : {args.bootstrap}")
    print()

    need_contigs = args.bootstrap > 0
    is_cosine = (args.distance == 'cosine')
    data_objects = []
    all_contig_data = [] if need_contigs else None

    for path, label in zip(args.inputs, labels):
        t0 = time.time()
        if is_cosine:
            if need_contigs:
                freq, n, contig_freqs = extract_freq_data(
                    path, args.kmer, with_contigs=True, min_entropy=min_entropy)
                all_contig_data.append(contig_freqs)
            else:
                freq, n = extract_freq_data(path, args.kmer, min_entropy=min_entropy)
            data_objects.append(freq)
        else:
            if need_contigs:
                bm, n, contig_codes = extract_data(
                    path, args.kmer, with_contigs=True, min_entropy=min_entropy)
                all_contig_data.append(contig_codes)
            else:
                bm, n = extract_data(path, args.kmer, min_entropy=min_entropy)
            data_objects.append(bm)
        elapsed = time.time() - t0
        fill_pct = n / max_code * 100
        extra = ""
        if need_contigs:
            extra = f"  {len(all_contig_data[-1]):>6,} contigs"
        print(f"  {label:>12s}  {n:>10,} k-mers ({fill_pct:.2f}% filled)  {elapsed:.1f}s{extra}")

    print(f"\nComputing pairwise {args.distance} distances "
          f"({len(data_objects)}x{len(data_objects)}, "
          f"{args.threads} threads)...", end=" ", flush=True)
    t0 = time.time()
    dist = build_distance_matrix(data_objects, args.threads)
    print(f"{time.time() - t0:.1f}s")

    print("\nDistance matrix:")
    col_w = max(12, max(len(l) for l in labels) + 2)
    print(" " * col_w + "".join(f"{l:>{col_w}}" for l in labels))
    for i, l in enumerate(labels):
        row = "".join(f"{dist[i][j]:{col_w}.6f}" for j in range(len(labels)))
        print(f"{l:>{col_w}}" + row)

    if args.dist:
        np.savetxt(args.dist, dist, delimiter=",", header=",".join(labels), comments="")
        print(f"\nDistance matrix saved to {args.dist}")

    if args.plot:
        try:
            plot_heatmap(dist, labels, args.plot)
            print(f"Heatmap saved to {args.plot}")
        except Exception as e:
            print(f"Warning: plot failed: {e}", file=sys.stderr)

    if args.bootstrap:
        print(f"\nBootstrapping ({args.bootstrap} replicates)...")
        t0 = time.time()
        if is_cosine:
            newick, support, children, root = bootstrap_phylo_cosine(
                all_contig_data, labels, args.threads, args.bootstrap)
        else:
            newick, support, children, root = bootstrap_phylo(
                all_contig_data, labels, args.kmer, args.threads, args.bootstrap)
        elapsed = time.time() - t0
        print(f"  done in {elapsed:.1f}s")

        with open(args.output, "w") as f:
            f.write(newick + "\n")

        print(f"\nBootstrap support values:")
        for node, value in sorted(support.items(), key=lambda x: -x[1]):
            clade_leaves = sorted(_get_leaves(node, children))
            bar = "#" * max(1, int(value / 5))
            print(f"  {value:5.0f}% {bar:20s} ({', '.join(clade_leaves)})")

        print(f"\nTree with bootstrap (Newick):\n{newick}")
        print(f"\nSaved to: {args.output}")
    else:
        print("\nBuilding Neighbor-Joining tree...", end=" ", flush=True)
        t0 = time.time()
        newick = neighbor_joining(dist, labels)[0]
        print(f"{time.time() - t0:.1f}s")

        with open(args.output, "w") as f:
            f.write(newick + "\n")

        print(f"\nTree (Newick):\n{newick}")
        print(f"\nSaved to: {args.output}")


if __name__ == "__main__":
    main()
