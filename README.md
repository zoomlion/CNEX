# CNEX — Conserved Non-coding Element Discovery

CNEX is a pipeline for rapid identification and phylogenetic analysis of conserved non-coding elements (CNEs) from whole-genome sequences. It supports both sequencing reads and direct genome FASTA input with an ultra-fast sliding-window mode (6X coverage, ~15 seconds for human genome).

## Dependencies

| Tool | Version | Purpose | Source |
|------|---------|---------|--------|
| g++ | 9+ (C++17) | Compile validate and assemble binaries | gcc.gnu.org |
| Python | 3.8+ | Pipeline scripts | python.org |
| pigz | 2.6+ | Parallel gzip decompression | [madler/pigz](https://github.com/madler/pigz) |
| robin-map | 1.3+ | C++ hash map for k-mer lookup | [Tessil/robin-map](https://github.com/Tessil/robin-map) |
| wgsim | — | Read simulation (optional) | [lh3/wgsim](https://github.com/lh3/wgsim) |

pigz and robin-map source are bundled under `src/pigz/` and `src/robin-map/`.

## Quick Start

```bash
# Build and install
make && make install

# Validate reads against the confident k-mer table
cnex validate reads.fq.gz --mers mers_table.tsv -t 8 --output-dir out

# Or validate a genome directly (sliding window, ~15s for human)
cnex validate genome.fa --mers mers_table.tsv --type genome -t 8 --output-dir out

# Assemble validated reads into contigs
cnex assemble out/ --mers mers_table.tsv --trim

# Or run both steps in one command
cnex pipeline genome.fa --mers mers_table.tsv --type genome -t 8 --output-dir out
```

## Pipeline Steps

```
00. filter_msa.py     — Filter multiple sequence alignment blocks
01. confi_mer.py      — Build confident k-mer table
02. validate          — C++ k-mer validator (reads or genome FASTA)
03. assemble          — C++ de Bruijn assembler (per-locus best pick)
04. map_contig_to_element.py — Group contigs by CNE element
05. element_phylo.py  — Per-element gene tree reconstruction
```

Steps 00, 01, 04, 05 are Python scripts in `src/`. Steps 02 (validate) and 03 (assemble) are C++ binaries built via `make`.

## CLI Reference

### `cnex validate`

Validate reads or genome against a confident k-mer table.

```bash
cnex validate <files...> --mers <table> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mers` | — | Confident k-mer table (required) |
| `-t, --threads` | 8 | Worker threads |
| `--output-dir` | `out` | Output directory |
| `--type` | auto | Input type: `genome` or `fastq` |
| `--min-c` | 10 | Minimum colinear k-mer matches |
| `--vote-frac` | 0.1 | Minimum vote fraction |
| `--vote-ratio` | 5.0 | Minimum vote ratio (1st/2nd) |
| `--window-size` | 150 | Sliding window size (genome mode) |
| `--step-size` | 25 | Sliding window step (genome mode, ~6X) |

### `cnex assemble`

De Bruijn assembly from validated reads.

```bash
cnex assemble <input_dir> --mers <table> [options]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mers` | — | Confident k-mer table (required) |
| `-o, --output` | `assembled.fasta` | Output FASTA |
| `-k, --kmer` | 35 | K-mer length |
| `--min-c` | 3 | Minimum k-mer matches for validation |
| `--trim` | on | Trim contigs to confident region |
| `--max-loci-gap` | 5000 | Max gap (bp) for locus clustering |

### `cnex pipeline`

Run validate + assemble in sequence.

```bash
cnex pipeline <files...> --mers <table> [options]
```

Accepts all `validate` options. Assembly output goes to `<output-dir>/assembled.fasta`.

## Input Formats

CNEX supports both plain text and gzip-compressed files:

- **FASTQ** (`.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`) — sequencing reads
- **FASTA** (`.fa`, `.fa.gz`, `.fasta`) — genome sequences (auto-detected)

When using `--type genome`, the input is treated as a genome FASTA and processed with a sliding window.

## Output Format (validate)

Tab-separated, 4 columns:
```
<seq_id>  <strand>  <ele_id>  <sequence>
```

In genome mode, `seq_id` encodes genomic coordinates: `chr:start-end`.

## Output Format (assemble)

FASTA with one entry per element. In genome mode, the header includes the best locus:
```
>ele_id
...sequence...
```

## License

MIT — see [LICENSE](LICENSE).
