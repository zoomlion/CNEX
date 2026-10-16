#!/usr/bin/env python3
"""Convert a UCSC MAF whole-genome alignment into CNEX's internal block-MSA format.

Each MAF alignment block becomes a FASTA block whose headers carry the genomic
coordinates, matching what the pipeline (01.mertable / rescue panels) and the
benchmark tooling expect:

    s <src> <start> <size> <strand> ...   ->   >{src}{strand}:{start+1}-{start+size}

Blocks are emitted in file order and separated by a line "###", so the output
is byte-compatible with the "###"-delimited block MSA used downstream. Only
conversion happens here; block filtering (00.filter_msa.py) is a separate step.

Usage:
    python3 utils/maf_convert.py INPUT.maf [-o OUTPUT.fa]

INPUT may be gzip-compressed (by extension or magic bytes). OUTPUT defaults to
most-cons-cne.fa.
"""
import argparse
import gzip
import os
import sys
import tempfile


def open_text(path):
    """Open a possibly gzip-compressed file as text, streaming."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if path.endswith(".gz") or magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def convert(inp, out):
    """Stream MAF blocks from `inp` to `out`. Returns (n_blocks, n_records)."""
    n_blocks = 0
    n_records = 0
    saw_maf = False
    cur = []  # list of (header, seq) for the current block

    def flush():
        nonlocal n_blocks
        if not cur:
            return
        for header, seq in cur:
            out.write(">" + header + "\n")
            out.write(seq + "\n")
        out.write("###\n")
        n_blocks += 1
        cur.clear()

    for raw in inp:
        line = raw.rstrip("\r\n")
        if not line:
            continue
        first = line[0]

        if line.startswith("a") and (len(line) == 1 or line[1] == " " or line[1] == "\t"):
            saw_maf = True
            flush()
        elif line.startswith("s ") or line.startswith("s\t"):
            parts = line.split(None, 6)
            if len(parts) < 7:
                continue  # malformed s line
            saw_maf = True
            src, start0, size, strand = parts[1], parts[2], parts[3], parts[4]
            try:
                start0i = int(start0)
                sizei = int(size)
            except ValueError:
                continue
            if sizei <= 0:
                continue
            header = "{}{}:{}-{}".format(src, strand, start0i + 1, start0i + sizei)
            cur.append((header, parts[6]))
            n_records += 1
        elif line.startswith("##maf"):
            saw_maf = True
        elif first == "#" or line.startswith("track"):
            continue
        elif first in ("i", "e", "q"):
            continue
        # anything else is ignored

    flush()

    if not saw_maf or n_blocks == 0:
        raise ValueError("input does not look like a MAF (no 'a'/'s' records found)")
    return n_blocks, n_records


def main():
    ap = argparse.ArgumentParser(
        description="Convert a UCSC MAF alignment to CNEX's ### block-MSA format."
    )
    ap.add_argument("input", help="input MAF file (plain or .gz)")
    ap.add_argument("-o", "--output", default="most-cons-cne.fa",
                    help="output block-MSA file (###-delimited; "
                         "default: most-cons-cne.fa)")
    args = ap.parse_args()

    # Write to a temp file in the output directory, then atomically replace the
    # target on success, so a bad input never truncates/clobbers an existing file.
    outdir = os.path.dirname(os.path.abspath(args.output)) or "."
    try:
        inp = open_text(args.input)
    except OSError as exc:
        sys.exit("Error: {}".format(exc))
    with inp:
        fd, tmp = tempfile.mkstemp(prefix=".maf_convert.", dir=outdir, text=True)
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        try:
            with os.fdopen(fd, "w") as out:
                n_blocks, n_records = convert(inp, out)
        except ValueError as exc:
            os.unlink(tmp)
            sys.exit("Error: {}".format(exc))
        except BaseException:
            os.unlink(tmp)
            raise

    os.replace(tmp, args.output)
    print("blocks={} records={} -> {}".format(n_blocks, n_records, args.output),
          file=sys.stderr)


if __name__ == "__main__":
    main()
