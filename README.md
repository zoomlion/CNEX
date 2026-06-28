# CNEX — Conserved Non-coding Element Discovery

CNEX is a pipeline for fast identification of conserved non-coding elements (CNEs) from whole-genome sequences. It supports both sequencing reads and direct genome FASTA input with an ultra-fast sliding-window mode (~15 seconds for the human genome at 6X coverage).

## Dependencies

### Build & Runtime

| Tool | Version | Purpose | Source |
|------|---------|---------|--------|
| g++ | 9+ (C++17) | Compile C++ binaries | gcc.gnu.org |
| Python | 3.8+ | CLI and pipeline scripts | python.org |
| pigz | 2.6+ | Parallel gzip decompression (bundled in `src/pigz/`) | [madler/pigz](https://github.com/madler/pigz) |
| robin-map | — | C++ hash map (bundled in `src/robin-map/`) | [Tessil/robin-map](https://github.com/Tessil/robin-map) |

### Phylogenetic Tools (optional — needed for steps 04-05)

| Tool | Purpose | Default Path | Source |
|------|---------|-------------|--------|
| FAMSA | Multiple sequence alignment | `~/Software/famsa` | [github.com/refresh-bio/FAMSA](https://github.com/refresh-bio/FAMSA) |
| FastTree | Gene tree inference | `~/Software/fasttree-2.2.0/FastTree` | [microbesonline.org/fasttree](http://www.microbesonline.org/fasttree/) |
| ASTRAL | Species tree inference | `~/Software/ASTRAL-5.7.1/astral_exe/astral.5.7.1.jar` | [github.com/smirarab/ASTRAL](https://github.com/smirarab/ASTRAL) |

To check your setup:
```bash
cnex setup
```

## Quick Start

### Build

```bash
make && make install
```

### Single Genome Pipeline

```bash
# From a confident k-mer table (recommended — runs 02+03)
cnex pipeline genome.fa --mers mers_table.tsv --type genome -t 8 -o results/

# From raw MSA (runs 00+01 then 02+03)
cnex pipeline genome.fa --cne most-cons-cne.fa -t 8 -o results/
```

### Multi-Genome Pipeline

```bash
# All genomes in a directory
cnex pipeline genomes/ --mers mers_table.tsv -j 4 -t 8 -o results/

# Glob pattern
cnex pipeline "genomes/*.fa" --mers mers_table.tsv -j 4 -t 8 -o results/

# Each genome gets its own subdirectory under results/
# results/sp1/contigs.fa  results/sp2/contigs.fa  results/sp3/contigs.fa
```

### Step-by-Step

```bash
# Validate reads or genome
cnex validate genome.fa --mers mers_table.tsv --type genome -t 8 -o out/

# Assemble validated reads into contigs
cnex assemble out/ --mers mers_table.tsv --trim
```

## CLI Reference

### `cnex pipeline`

Run the full CNE discovery pipeline on one or more genomes.

```
usage: cnex pipeline [-h] [--cne CNE] [--mers MERS] [-t THREADS] [-j JOBS]
                     [-o OUTPUT_DIR] [--no-trim] [--type {genome,fastq}]
                     [--window-size WINDOW_SIZE] [--step-size STEP_SIZE]
                     [--min-c MIN_C] [--vote-frac VOTE_FRAC]
                     [--vote-ratio VOTE_RATIO]
                     files [files ...]
```

Input can be a directory, individual files, or glob patterns.

| Option | Default | Description |
|--------|---------|-------------|
| `--cne` | — | Start from MSA file (runs 00+01) |
| `--mers` | — | Start from existing k-mer table |
| `-t, --threads` | 8 | Worker threads per genome |
| `-j, --jobs` | 4 | Parallel genomes |
| `-o, --output-dir` | `out` | Output root directory |
| `--type` | auto | Input type: `genome` or `fastq` |
| `--no-trim` | — | Disable contig trimming (on by default) |
| `--min-c` | 10 | Minimum colinear k-mer matches |
| `--vote-frac` | 0.1 | Minimum vote fraction |
| `--vote-ratio` | 5.0 | Minimum vote ratio |

### `cnex validate`

Validate reads or genome against a confident k-mer table.

```
usage: cnex validate [-h] --mers MERS [-t THREADS] [-o OUTPUT_DIR]
                     [--type {genome,fastq}] [--min-c MIN_C]
                     [--vote-frac VOTE_FRAC] [--vote-ratio VOTE_RATIO]
                     [--window-size WINDOW_SIZE] [--step-size STEP_SIZE]
                     files [files ...]
```

### `cnex assemble`

De Bruijn assembly from validated reads.

```
usage: cnex assemble [-h] --mers MERS [-o OUTPUT] [-k KMER] [--min-c MIN_C]
                     [--min-count MIN_COUNT] [--max-reads MAX_READS]
                     [--max-loci-gap MAX_LOCI_GAP] [--no-trim]
                     input_dir
```

### `cnex setup`

Check that all binaries and optional phylogenetic tools are available.

## Input Formats

Supports plain text and gzip-compressed files:

- **FASTQ** (`.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`) — sequencing reads
- **FASTA** (`.fa`, `.fa.gz`, `.fasta`) — genome sequences (auto-detected)

## Output

### Validate Output (`Assemble.*.reads`)

Tab-separated, 4 columns:

```
<seq_id>  <strand>  <ele_id>  <sequence>
```

In genome mode, `seq_id` encodes genomic coordinates: `chr:start-end`.

### Assemble Output (`contigs.fa`)

FASTA with one entry per element. In genome mode, the header includes the best locus.

```
>ele_id
...sequence...
```

## Pipeline Steps

| Step | Script/Binary | Description |
|------|--------------|-------------|
| 00 | `00.filter_msa.py` | Filter MSA blocks (invoked by `pipeline --cne`) |
| 01 | `01.confi_mer.py` | Build confident k-mer table (invoked by `pipeline --cne`) |
| 02 | `validate` (C++) | k-mer validation (reads or genome FASTA sliding window) |
| 03 | `assemble` (C++) | De Bruijn assembly (per-locus best pick for genome mode) |
| 04 | `04.map_contig_to_element.py` | Group contigs by CNE element |
| 05 | `05.element_phylo.py` | Per-element gene trees (needs FAMSA, FastTree, ASTRAL) |

## License

MIT — see [LICENSE](LICENSE).
