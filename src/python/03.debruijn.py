#!/usr/bin/env python3
import os, sys, argparse, random, re, string
from collections import defaultdict
from hip import debruijn
from hip.validator import MerQueryManager, validate_read, dna_encoder


def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]


def debruijn_assembler(reads, k, min_count=2, mer_query=None, ele_id=None):
    graph = debruijn.de_bruijn_graph(reads, k, min_count)
    if not graph:
        return ""
    node_scores = {}
    if mer_query is not None and ele_id is not None:
        ns = _build_node_scores(graph, mer_query, int(ele_id))
        if ns is not None:
            node_scores = ns
    return debruijn.assemble_sequence(graph, node_scores)


def _build_node_scores(graph, mer_query, ele_id):
    scores = {}
    for node in graph:
        score = 0
        for i in range(len(node) - 13 + 1):
            try:
                code = dna_encoder(node[i:i+13])
            except Exception:
                continue
            result = mer_query.query_mer(code)
            if result is not None and result[0] == ele_id:
                score += 1
        if score > 0:
            scores[node] = score
    return scores if scores else None


def score_contig(seq, mer_query, ele_id):
    sc = 0
    for i in range(len(seq) - 13 + 1):
        try:
            code = dna_encoder(seq[i:i+13])
        except Exception:
            continue
        result = mer_query.query_mer(code)
        if result is not None and result[0] == ele_id:
            sc += 1
    return sc


GENOME_ID = re.compile(r'^(\S+):(\d+)-(\d+)$')


def reads_generator(f, max_loci_gap=5000):
    for first_line in f:
        parts = first_line.strip().split("\t")
        if len(parts) != 4:
            continue
        if GENOME_ID.match(parts[0]):
            yield from _genome_generator(f, first_line, max_loci_gap)
        else:
            yield from _reads_generator(f, first_line)
        return


def _reads_generator(f, first_line):
    prev = None
    buf = {}
    for line in [first_line] + list(f):
        p = line.strip().split("\t")
        if len(p) != 4: continue
        rid, strand, eid, seq = p
        seq = seq.upper()
        strand = int(strand)
        if prev is None:
            prev = eid
        if eid != prev:
            yield prev, prev, buf
            buf = {}
            prev = eid
        buf[rid] = (strand, seq)
    if prev is not None:
        yield prev, prev, buf


def _genome_generator(f, first_line, max_loci_gap=5000):
    buf = []
    cur_eid = None

    def flush():
        if not buf:
            return
        eid = buf[0]['ele_id']
        buf.sort(key=lambda r: (r['chr'], r['start']))
        groups = defaultdict(list)
        for r in buf:
            groups[(r['chr'], r['strand'])].append(r)
        for (chr_name, strand), grp in groups.items():
            cluster = []
            cstart = grp[0]['start']
            for r in grp:
                if cluster and (r['start'] - cluster[-1]['end']) > max_loci_gap:
                    cend = cluster[-1]['end']
                    s = '+' if cluster[0]['strand'] == 1 else '-'
                    did = f"{eid}.{chr_name}:{cstart}-{cend}({s})"
                    tmp = {x['rid']: (x['strand'], x['seq']) for x in cluster}
                    yield did, int(eid), tmp
                    cluster = []
                    cstart = r['start']
                cluster.append(r)
            if cluster:
                cend = cluster[-1]['end']
                s = '+' if cluster[0]['strand'] == 1 else '-'
                did = f"{eid}.{chr_name}:{cstart}-{cend}({s})"
                tmp = {x['rid']: (x['strand'], x['seq']) for x in cluster}
                yield did, int(eid), tmp

    for line in [first_line] + list(f):
        p = line.strip().split("\t")
        if len(p) != 4: continue
        rid, strand, eid, seq = p
        m = GENOME_ID.match(rid)
        if not m: continue
        r = {
            'rid': rid, 'ele_id': int(eid),
            'strand': int(strand), 'seq': seq.upper(),
            'chr': m.group(1),
            'start': int(m.group(2)),
            'end': int(m.group(3)),
        }
        if cur_eid is None:
            cur_eid = r['ele_id']
        if r['ele_id'] != cur_eid:
            yield from flush()
            buf = []
            cur_eid = r['ele_id']
        buf.append(r)
    yield from flush()


def assemble_reads(temp_reads, k=35, min_count=2, max_reads=200, mer_query=None, ele_id=None):
    reads = []
    n = 0
    for rid, (strand, seq) in temp_reads.items():
        local_seq = seq if strand == 1 else reverse_complement(seq)
        reads.append(local_seq)
        n += 1
        if n >= max_reads: break
    if not reads:
        return ""
    return debruijn_assembler(reads, k, min_count, mer_query, ele_id)


def filter_contig(seq, mer_query, mer_size, ele_id, min_c=3):
    if len(seq) < mer_size + min_c:
        return True
    confi_id, _ = validate_read(seq, mer_query, mer_size, min_c)
    return confi_id == ele_id


def trim_contig(seq, mer_query, ele_id, min_density=0.3):
    if len(seq) < 26:
        return seq
    n = len(seq) - 12
    first_hit, last_hit = -1, -1
    for i in range(n):
        try:
            code = dna_encoder(seq[i:i + 13])
        except Exception:
            continue
        result = mer_query.query_mer(code)
        if result is not None and result[0] == ele_id:
            if first_hit < 0: first_hit = i
            last_hit = i
    if first_hit < 0 or last_hit < 0:
        return seq
    margin = 6
    start = max(0, first_hit - margin)
    end = min(n, last_hit + margin)
    trimmed = seq[start:end + 12]
    return seq if len(trimmed) < 20 else trimmed


def main():
    parser = argparse.ArgumentParser(description="De Bruijn Graph Assembly")
    parser.add_argument("inputs_dir", help="directory holding Assemble.*.reads files")
    parser.add_argument("-k", "--kmer", type=int, default=35, help="k-mer length")
    parser.add_argument("-o", "--output", default="assembled.fasta", help="output fasta")
    parser.add_argument("--mers", required=True, help="mers table for assembly validation")
    parser.add_argument("--min_c", type=int, default=3, help="min k-mer matches for validation")
    parser.add_argument("--min_count", type=int, default=2,
        help="min k-mer occurrence to keep in de Bruijn graph")
    parser.add_argument("--max_reads", type=int, default=200,
                        help="max reads per element/locus for assembly")
    parser.add_argument("--max-loci-gap", type=int, default=5000,
                        help="max gap (bp) between reads on same chr/strand for one locus")
    parser.add_argument("--trim", action="store_true", default=False,
                        help="Trim contigs to high-density confident k-mer regions")
    parser.add_argument("--min-density", type=float, default=0.3,
                        help="minimum confident k-mer density for trim")
    args = parser.parse_args()

    if not os.path.exists(args.inputs_dir):
        raise FileNotFoundError(f"Input directory not found: {args.inputs_dir}")

    print("Loading mers table ...", end=" ", flush=True)
    mer_query = MerQueryManager()
    mer_query.load_from_file(args.mers)
    mer_size = mer_query.get_mer_size()
    print(f"done ({mer_query.size()} mers, k={mer_size})")

    tag = ''.join(random.choices(string.ascii_letters, k=8))
    input_f = f"{tag}.reads"
    os.system(f"cat {args.inputs_dir}/*.reads | sort -k3n -k1 > {input_f}")
    if not os.path.exists(input_f):
        raise FileNotFoundError(f"Sorted reads file not found: {input_f}")

    best = {}
    results = []
    from_genome = None

    with open(input_f) as f:
        for display_id, ele_id, temp_reads in reads_generator(f, args.max_loci_gap):
            if from_genome is None:
                from_genome = isinstance(display_id, str) and '.' in display_id

            seq = assemble_reads(temp_reads, args.kmer, args.min_count, args.max_reads,
                                 mer_query, ele_id)
            if not seq:
                continue
            if args.trim:
                seq = trim_contig(seq, mer_query, int(ele_id), args.min_density)
                if not seq or len(seq) < 20:
                    continue
            if not filter_contig(seq, mer_query, mer_size, int(ele_id), args.min_c):
                continue

            if from_genome:
                sc = score_contig(seq, mer_query, int(ele_id))
                if ele_id not in best or sc > best[ele_id][2]:
                    best[ele_id] = (display_id, seq, sc)
            else:
                results.append((display_id, seq))
            print(f"\rAssembled {len(results) + len(best)} elements", end="", flush=True)
    print()

    if from_genome:
        for ele_id in sorted(best, key=lambda x: int(x)):
            did, seq, _ = best[ele_id]
            results.append((did, seq))

    with open(args.output, "w") as f:
        for seq_id, seq in results:
            f.write(f">{seq_id}\n{seq}\n")

    os.remove(input_f)
    print(f"Done. {len(results)} sequences written to {args.output}")


if __name__ == "__main__":
    main()
