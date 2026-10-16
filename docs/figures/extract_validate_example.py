#!/usr/bin/env python3
"""Extract a real accepted window from a CNEX validate run and independently
recompute the validator's vote / chain / gate numbers (data source for the
validate principle figure).

Read-only with respect to the pipeline: reads `Assemble.*.reads` (accepted
windows) and the mers table only.

The recomputation mirrors src/hip/validator_core.hpp exactly:
  k = 13, window = 150 (nspan = 138), max_diff = 2, vote_frac = 0.1,
  gate rows T1 = (3.0, 7, 25), T2 = (2.0, 10, 30), OR-combined;
  tally order = match order (position ascending, mers-table row order),
  best  = first element reaching the max count, second = next best.

Usage:
  python3 extract_validate_example.py --list            # rank candidates
  python3 extract_validate_example.py --pick <header>   # dump JSON
"""

import argparse
import glob
import json
import random
import sys

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
        "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}


def encode(seq, i, k):
    code = 0
    for j in range(i, i + k):
        b = BASES.get(seq[j])
        if b is None:
            return None
        code = (code << 2) | b
    return code


def rc(seq):
    return "".join(COMP.get(c, "N") for c in reversed(seq))


def match_list(seq, table, k):
    """(window position x, element id, element locus) in validator order."""
    out = []
    for x in range(len(seq) - k + 1):
        code = encode(seq, x, k)
        if code is None:
            continue
        rows = table.get(code)
        if rows:
            out.extend((x, eid, loci) for eid, loci in rows)
    return out


def tally(matches):
    counts = {}
    seen = []
    for _, eid, _ in matches:
        if eid in counts:
            counts[eid] += 1
        else:
            counts[eid] = 1
            seen.append(eid)
    best_id, best_cnt, second = -1, 0, 0
    for eid in seen:
        c = counts[eid]
        if c > best_cnt:
            second, best_cnt, best_id = best_cnt, c, eid
        elif c > second:
            second = c
    return best_id, best_cnt, second, len(matches)


def chain(matches, eid, max_diff, k):
    first = last = prev_x = prev_l = None
    passing = m = 0
    for x, e, l in matches:
        if e != eid:
            continue
        if first is None:
            first = x
        elif abs((x - prev_x) - (l - prev_l)) <= max_diff:
            passing += 1
        prev_x, prev_l, last, m = x, l, x, m + 1
    if m == 0:
        return 0, 0
    return passing, last - first + k


def evaluate(seq, table, k, tiers, max_diff, vote_frac):
    res = {}
    for name in ("fwd", "rc"):
        strand_seq = seq if name == "fwd" else rc(seq)
        matches = match_list(strand_seq, table, k)
        best_id, best_cnt, second, n = tally(matches)
        passing = span = 0
        passed = []
        if best_id >= 0 and best_cnt > n * vote_frac:
            passing, span = chain(matches, best_id, max_diff, k)
            for rv, mc, ms in tiers:
                if best_cnt >= int(rv * second) and passing >= mc and span >= ms:
                    passed.append([rv, mc, ms])
        res[name] = dict(best_id=best_id, best_cnt=best_cnt, second=second,
                         n=n, passing=passing, span=span, passed=passed,
                         matches=matches)
    if not res["fwd"]["passed"] and not res["rc"]["passed"]:
        res["emit_strand"] = None
    elif not res["rc"]["passed"]:
        res["emit_strand"] = "fwd"
    elif not res["fwd"]["passed"]:
        res["emit_strand"] = "rc"
    else:
        res["emit_strand"] = ("fwd" if res["fwd"]["n"] >= res["rc"]["n"]
                              else "rc")
    return res


def collect_candidates(paths, per_file, rng):
    cands = []
    for path in paths:
        with open(path) as fh:
            lines = fh.read().splitlines()
        lines = [ln for ln in lines if ln]
        take = min(per_file, len(lines))
        for ln in rng.sample(lines, take):
            header, strand, eid, seq = ln.split("\t")
            cands.append((header, int(strand), int(eid), seq))
    return cands


def stream_table(path, wanted_codes):
    table = {}
    with open(path) as fh:
        for ln in fh:
            mer, eid, loci, _cnt = ln.rstrip("\n").split("\t")
            code = 0
            ok = True
            for ch in mer:
                b = BASES.get(ch)
                if b is None:
                    ok = False
                    break
                code = (code << 2) | b
            if not ok or code not in wanted_codes:
                continue
            table.setdefault(code, []).append((int(eid), int(loci)))
    return table


def score(ev, cand_strand, tiers):
    if ev["emit_strand"] is None:
        return None
    if cand_strand != (1 if ev["emit_strand"] == "fwd" else -1):
        return None
    me = ev[ev["emit_strand"]]
    if me["passed"] != [list(tiers[1])]:   # T2-only (OR showcase)
        return None
    if not (12 <= me["best_cnt"] <= 45 and me["second"] >= 2):
        return None
    if me["passing"] < 10 or me["span"] < 60:
        return None
    ratio = me["best_cnt"] / me["second"]
    pts = [m for m in me["matches"] if m[1] == me["best_id"]]
    monotone = all(pts[i][2] >= pts[i - 1][2] for i in range(1, len(pts)))
    clean = me["passing"] / (len(pts) - 1) if len(pts) > 1 else 0.0
    return (monotone and clean >= 0.75, clean, len(pts), me["passing"],
            -abs(ratio - 2.3))


def load_bunch_rows(blocks_path, bunch_id):
    """Species rows of one MSA bunch, mirroring 01.mertable.cpp parsing
    (split on '###', drop 'ref*')."""
    with open(blocks_path) as fh:
        bunches = fh.read().split("###\n")
    if bunch_id >= len(bunches):
        return None
    rows, sid, s = [], None, []
    for ln in bunches[bunch_id].split("\n"):
        if ln.startswith(">"):
            if sid is not None:
                rows.append((sid, "".join(s)))
            sid = ln[1:].strip()
            s = []
        elif ln:
            s.append(ln.replace(" ", "").upper())
    if sid is not None:
        rows.append((sid, "".join(s)))
    return [(h, q) for h, q in rows if not h.startswith("ref")]


def load_bunch_consensus(blocks_path, bunch_id):
    """Column-wise majority (non-gap) over all species rows of one bunch."""
    rows = load_bunch_rows(blocks_path, bunch_id)
    if not rows:
        return None
    length = max(len(q) for _, q in rows)
    consensus = []
    for i in range(length):
        counts = {}
        for _, q in rows:
            if i < len(q) and q[i] != "-":
                counts[q[i]] = counts.get(q[i], 0) + 1
        consensus.append(max(counts, key=counts.get) if counts else "-")
    return "".join(consensus)


def match_species_count(rows, mer, locus, k):
    """How many aligned species rows carry `mer` starting at MSA column `locus`."""
    hits = []
    for header, q in rows:
        pure = q.replace("-", "")
        bias_val, pure_idx = 0, None
        for i, ch in enumerate(q):
            if ch == "-":
                bias_val += 1
            else:
                if i == locus:
                    pure_idx = i - bias_val
                    break
        if pure_idx is not None and pure[pure_idx:pure_idx + k] == mer:
            hits.append(header.split(".")[0])
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads-glob", default=(
        "benchmark/results/Callorhinchus_milii/cnex/c6/"
        "validate_out/Assemble.*.reads"))
    ap.add_argument("--mers", default=(
        "benchmark/dbs/Callorhinchus_milii/cnex/mers_table_a3.tsv"))
    ap.add_argument("--per-file", type=int, default=25)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--k", type=int, default=13)
    ap.add_argument("--max-diff", type=int, default=2)
    ap.add_argument("--vote-frac", type=float, default=0.1)
    ap.add_argument("--tiers", default="3.0,7,25 2.0,10,30")
    ap.add_argument("--pick", default=None)
    ap.add_argument("--diag", action="store_true")
    ap.add_argument("--blocks", default=(
        "benchmark/data/blocks/Callorhinchus_milii.fa"))
    ap.add_argument("--out", default="docs/figures/validate_example.json")
    args = ap.parse_args()

    tiers = [tuple(float(v) if i == 0 else int(v)
                   for i, v in enumerate(spec.split(",")))
             for spec in args.tiers.split()]

    rng = random.Random(args.seed)
    paths = sorted(glob.glob(args.reads_glob))
    if not paths:
        sys.exit("no reads files match " + args.reads_glob)

    if args.pick is not None:
        # pick an exact emitted record (no T2-only ranking required)
        pick_line = None
        for p in paths:
            with open(p) as fh:
                for ln in fh:
                    if ln.startswith(args.pick + "\t"):
                        pick_line = ln.rstrip("\n")
                        break
            if pick_line is not None:
                break
        if pick_line is None:
            sys.exit("record not found in emitted files: " + args.pick)
        header, s, e, q = pick_line.split("\t")
        cands = [(header, int(s), int(e), q)]
    else:
        cands = collect_candidates(paths, args.per_file, rng)
    print(f"candidates: {len(cands)} from {len(paths)} files", file=sys.stderr)

    wanted = set()
    for _, _, _, seq in cands:
        for s in (seq, rc(seq)):
            for x in range(len(s) - args.k + 1):
                code = encode(s, x, args.k)
                if code is not None:
                    wanted.add(code)
    print(f"distinct k-mer codes wanted: {len(wanted)}", file=sys.stderr)
    table = stream_table(args.mers, wanted)
    print(f"table rows kept: {len(table)} mers", file=sys.stderr)

    if args.pick is None:
        ranked = []
        mismatches = 0
        for header, strand, eid, seq in cands:
            ev = evaluate(seq, table, args.k, tiers, args.max_diff,
                          args.vote_frac)
            sc = score(ev, strand, tiers)
            rep = (ev["emit_strand"], ev[ev["emit_strand"]]["best_id"]
                   if ev["emit_strand"] else None)
            if ev["emit_strand"] and rep[1] != eid or \
               (ev["emit_strand"] == "fwd") != (strand == 1):
                mismatches += 1
                continue
            if args.diag:
                me_ = ev[ev["emit_strand"]] if ev["emit_strand"] else None
                if me_:
                    print(f"DIAG {header}\tstrand={strand}\tele={eid}"
                          f"\tc1={me_['best_cnt']}\tc2={me_['second']}"
                          f"\tn={me_['n']}\tpassing={me_['passing']}"
                          f"\tspan={me_['span']}\tpassed={me_['passed']}")
            if sc:
                ranked.append((sc, header, strand, eid, seq, ev))
        ranked.sort(key=lambda r: r[0], reverse=True)
        print(f"replication mismatches dropped: {mismatches}", file=sys.stderr)
        print(f"T2-only candidates: {len(ranked)}", file=sys.stderr)

        for sc, header, strand, eid, seq, ev in ranked[:20]:
            m = ev[ev["emit_strand"]]
            print(f"{header}\tstrand={strand}\tele={eid}"
                  f"\tc1={m['best_cnt']}\tc2={m['second']}\tn={m['n']}"
                  f"\tpassing={m['passing']}\tspan={m['span']}")
        return

    # exact pick: evaluate and verify against the emitted record
    header, strand, eid, seq = cands[0]
    ev = evaluate(seq, table, args.k, tiers, args.max_diff, args.vote_frac)
    key = ev["emit_strand"]
    if key is None or ev[key]["best_id"] != eid or \
       ((key == "fwd") != (strand == 1)):
        sys.exit("replication mismatch for pick " + args.pick)
    m = ev[key]
    best = m["best_id"]
    points = [[x, l] for x, e, l in m["matches"] if e == best]

    # Alignment view (figure): read letters vs element consensus, with the
    # constant offset of the winning chain and the columns that break mers.
    cons = load_bunch_consensus(args.blocks, best)
    offsets = [l - x for x, l in points]
    offset = max(set(offsets), key=offsets.count)
    read_full = seq if ev["emit_strand"] == "fwd" else rc(seq)
    x0 = max(0, min(x for x, _ in points) - 6)
    x1 = min(len(read_full) - 1, max(x for x, _ in points) + 6)
    read_letters = read_full[x0:x1 + 1]
    consensus_letters = (cons[x0 + offset:x1 + offset + 1]
                         if cons else "?" * (x1 - x0 + 1))
    consensus_letters = consensus_letters.ljust(x1 - x0 + 1, "-")[:x1 - x0 + 1]
    mismatches = [x for x in range(x0, x1 + 1)
                  if cons and read_full[x] != cons[x + offset]]

    # Provenance: which aligned species rows actually carry the chain k-mers.
    rows = load_bunch_rows(args.blocks, best) or []
    species_counts = [len(match_species_count(rows, read_full[x:x + args.k],
                                              l, args.k)) for x, l in points]

    # Lookup example for the mechanism panel: the matched position with the
    # most table assignments (shows one read position -> several candidates).
    lookup_example = None
    for x, l in points:
        code = encode(read_full, x, args.k)
        assigns = table.get(code) if code is not None else None
        if assigns and (lookup_example is None
                        or len(assigns) > len(lookup_example[2])):
            lookup_example = (x, read_full[x:x + args.k],
                              [(e, lo) for e, lo in assigns])

    data = dict(
        source=dict(reads_file=args.reads_glob, mers=args.mers,
                    blocks=args.blocks),
        k=args.k, window_size=len(seq), max_diff=args.max_diff,
        vote_frac=args.vote_frac, tiers=[list(t) for t in tiers],
        window=dict(header=header, strand=strand, element=eid, seq=seq),
        emitted=dict(strand=1 if ev["emit_strand"] == "fwd" else -1,
                     element=best),
        fwd={kk: vv for kk, vv in ev["fwd"].items() if kk != "matches"},
        rc={kk: vv for kk, vv in ev["rc"].items() if kk != "matches"},
        chain_points=points,
        matches=[list(mm) for mm in m["matches"]],
        alignment=dict(offset=offset, x_range=[x0, x1],
                       read_letters=read_letters,
                       consensus_letters=consensus_letters,
                       mismatch_x=mismatches,
                       species_rows=len(rows),
                       match_species_counts=species_counts,
                       min_species_count=min(species_counts) if species_counts else 0),
        lookup_example=(dict(x=lookup_example[0], mer=lookup_example[1],
                             assignments=lookup_example[2])
                        if lookup_example else None),
    )
    with open(args.out, "w") as fh:
        json.dump(data, fh, indent=1)
    print(f"wrote {args.out}", file=sys.stderr)
    assert data["emitted"]["element"] == eid
    assert data["emitted"]["strand"] == strand


if __name__ == "__main__":
    main()
