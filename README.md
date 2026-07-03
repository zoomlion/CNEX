# CNEX — Conserved Non-coding Element Discovery

CNEX is a pipeline for fast identification of conserved non-coding elements (CNEs) from whole-genome sequences. It supports both sequencing reads and direct genome FASTA input with an ultra-fast sliding-window mode (~3 minutes for the human genome at 6X coverage, 16 threads).

## Quick Start

```bash
# Clone — --recurse-submodules is required for robin-map
git clone --recurse-submodules git@github.com:zoomlion/CNEX.git
cd CNEX

# If you already cloned without submodules, run this:
# git submodule update --init

# Build all C++ binaries in one step
make && make install

# Check your setup
cnex setup

# Run pipeline on a genome
cnex pipeline genome.fa --mers mers_table.tsv --type genome -t 8 -o results/
```

## Examples

### Single Genome

```bash
# From a confident k-mer table (runs 02+03)
cnex pipeline genome.fa --mers mers_table.tsv --type genome -t 8 -o results/

# From raw MSA (runs 00+01 then 02+03)
cnex pipeline genome.fa --cne most-cons-cne.fa -t 8 -o results/
```

### Multi-Genome

```bash
# Directory of genomes — parallel with -j
cnex pipeline genomes/ --mers mers_table.tsv -j 4 -t 8 -o results/

# Glob pattern
cnex pipeline "genomes/*.fa" --mers mers_table.tsv -j 4 -t 8 -o results/

# Each genome gets its own subdirectory:
# results/Homo_sapiens/contigs.fa  results/Mus_musculus/contigs.fa
```

### Step-by-Step

```bash
# Validate reads or genome
cnex validate genome.fa --mers mers_table.tsv --type genome -t 8 -o out/

# Assemble validated reads into contigs
cnex assemble out/ --mers mers_table.tsv --trim
```

## Dependencies

### Build & Runtime

| Tool | Version | Purpose | Source |
|------|---------|---------|--------|
| g++ | 9+ (C++17) | Compile C++ binaries (mertable, validate, assemble) | gcc.gnu.org |
| Python | 3.8+ | CLI and pipeline helper scripts | python.org |
| pigz | 2.6+ | Parallel gzip decompression (bundled in `src/pigz/`) | [madler/pigz](https://github.com/madler/pigz) |
| robin-map | — | C++ hash map (git submodule, bundled in `src/robin-map/`) | [Tessil/robin-map](https://github.com/Tessil/robin-map) |

### Phylogenetic Tools (optional — needed for steps 04-05)

| Tool | Purpose | Default Path | Source |
|------|---------|-------------|--------|
| FAMSA | Multiple sequence alignment | `famsa` (via PATH or `--famsa`) | [refresh-bio/FAMSA](https://github.com/refresh-bio/FAMSA) |
| FastTree | Gene tree inference | `FastTree` (via PATH or `--fasttree`) | [microbesonline.org/fasttree](http://www.microbesonline.org/fasttree/) |
| ASTRAL IV (ASTER) | Species tree inference (**recommended**) | `astral` (via PATH or `--astral-bin`) | [chaoszhang/ASTER](https://github.com/chaoszhang/ASTER) |
| ASTRAL III | Species tree inference (fallback) | (not set, use `--astral-jar`) | [smirarab/ASTRAL](https://github.com/smirarab/ASTRAL) |

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
| `--min-c` | 7 | Minimum colinear k-mer matches |
| `--vote-frac` | 0.1 | Minimum vote fraction |
| `--vote-ratio` | 3.0 | Minimum vote ratio |

### `cnex validate`

Validate reads or genome against a confident k-mer table.

```
usage: cnex validate [-h] --mers MERS [-t THREADS] [-o OUTPUT_DIR]
                     [--type {genome,fastq}] [--min-c MIN_C]
                     [--vote-frac VOTE_FRAC] [--vote-ratio VOTE_RATIO]
                     [--window-size WINDOW_SIZE] [--step-size STEP_SIZE]
                     files [files ...]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mers` | — | Confident k-mer table (required) |
| `-t, --threads` | 8 | Worker threads |
| `-o, --output-dir` | `out` | Output directory |
| `--type` | auto | Input type: `genome` or `fastq` |
| `--min-c` | 10 | Minimum colinear k-mer matches |
| `--vote-frac` | 0.1 | Minimum vote fraction |
| `--vote-ratio` | 5.0 | Minimum vote ratio (1st/2nd) |
| `--window-size` | 150 | Sliding window size (genome mode) |
| `--step-size` | 25 | Sliding window step (genome mode, ~6X) |

### `cnex assemble`

De Bruijn assembly from validated reads.

```
usage: cnex assemble [-h] --mers MERS [-o OUTPUT] [-k KMER] [--min-c MIN_C]
                     [--min-count MIN_COUNT] [--max-reads MAX_READS]
                     [--max-loci-gap MAX_LOCI_GAP] [--no-trim]
                     input_dir
```

| Option | Default | Description |
|--------|---------|-------------|
| `--mers` | — | Confident k-mer table (required) |
| `-o, --output` | `assembled.fasta` | Output FASTA |
| `-k, --kmer` | 35 | K-mer length |
| `--min-c` | 3 | Minimum k-mer matches for validation |
| `--trim` | on | Trim contigs to confident region |
| `--max-loci-gap` | 5000 | Max gap (bp) for locus clustering |

### `cnex setup`

Check that all binaries and optional phylogenetic tools are available.

```
cnex setup
```

## Pipeline Steps

| Step | Script/Binary | Description |
|------|--------------|-------------|
| 00 | `00.filter_msa.py` | Filter MSA blocks (invoked by `pipeline --cne`) |
| 01 | `mertable` (C++) | Build confident k-mer table (invoked by `pipeline --cne`) |
| 02 | `validate` (C++) | k-mer validation (reads or genome FASTA sliding window) |
| 03 | `assemble` (C++) | De Bruijn assembly (per-locus best pick for genome mode) |
| 04 | `04.map_contig_to_element.py` | Group contigs by CNE element |
| 05 | `05.element_phylo.py` | Per-element gene trees (needs FAMSA, FastTree, ASTRAL) |

## License

MIT — see [LICENSE](LICENSE).
