# CNEX — Conserved Non-coding Element Discovery

Fast identification of conserved non-coding elements (CNEs) from whole-genome sequences. Supports both sequencing reads and direct genome FASTA input. Key features:

- **Ultra-fast sliding-window genome mode** (~3 min for human genome at 6X, 16 threads)
- **k-mer guided assembly** — De Bruijn graph with element-specific k-mer prior
- **SNP/INDEL detection** via De Bruijn bubble scanning
- **Phylogeny pipeline** — concat/IQ-TREE 3 or FastTree + ASTRAL, with configurable threshold and tag filtering

## Install

```bash
git clone --recurse-submodules git@github.com:zoomlion/CNEX.git
cd CNEX
make && make install
cnex setup
```

## Configuration

Copy the example and edit:

```bash
cp config.example.py config.py
```

Key settings in `config.py`:

| Setting | Default | Description |
|---------|---------|-------------|
| `ALIGNMENT_JOBS` | 20 | Parallel FAMSA / FastTree processes |
| `IQTREE_THREADS` | 20 | Threads for IQ-TREE 3 |
| `ASTRAL_THREADS` | 20 | Threads for ASTRAL |
| `MIN_CNE_PER_SPECIES` | 100 | Minimum CNEs per species to retain |
| `CONCAT_LENGTH_QUANTILES` | `25,50,75` | Length quantile thresholds for concat method |
| `ASTRAL_LENGTH_QUANTILES` | `25,50,75` | Length quantile thresholds for astral method |
| `ELEMENT_TAGS_FILE` | `""` | Path to TSV for per-tag filtering |
| `DRY_RUN` | `False` | Only generate scripts, skip execution |

## Quick Start

### 1. Validate reads/genome against a k-mer table

```bash
# Genome mode (sliding windows)
cnex validate genome.fa --mers mers_table.tsv --type genome -t 16 -o out/

# Reads mode (FASTQ)
cnex validate reads.fq.gz --mers mers_table.tsv --type fastq -t 16 -o out/
```

### 2. Assemble contigs

```bash
# Basic assembly (--trim default on, --repeat-ratio 1.1 default)
cnex assemble out/ --mers mers_table.tsv -o contigs.fa

# With SNP/INDEL detection
cnex assemble out/ --mers mers_table.tsv --snp -o contigs.fa
```

Output: `contigs.fa`, `variants.tsv`, `snp_elements.gfa`

### 3. Phylogeny

```bash
# Default: concat + IQ-TREE 3
python3 src/04.map_contig_to_element.py
python3 src/05.element_phylo.py

# Or ASTRAL method
python3 src/05.element_phylo.py --method astral

# Both methods
python3 src/05.element_phylo.py --method both

# Dry-run: only generate scripts, no execution
python3 src/05.element_phylo.py --dry-run
```

## Phylogeny Pipeline

### Methods

| Method | Command | Pipeline |
|--------|---------|----------|
| **concat** (default) | `--method concat` | FAMSA → concat_msa → IQ-TREE 3 |
| **astral** | `--method astral` | FAMSA → FastTree → ASTRAL |
| **both** | `--method both` | Both methods in one run |

### Thresholds and Tags

Elements are filtered by aligned length quantiles (P25/P50/P75). Each combination generates:
- A subdirectory with `run.sh`
- A top-level `run_all.sh` to execute all

```
results/
├── iqtree/
│   ├── all/run.sh
│   ├── all_quantile_25/run.sh
│   ├── all_quantile_50/run.sh
│   └── run_all.sh
├── astral/
│   ├── all/run.sh
│   └── run_all.sh
├── aln/       # alignments (shared)
└── fasta/     # element FASTAs (from step 04)
```

Optional tag filtering via TSV:

```bash
cp element_tags.example.tsv element_tags.tsv
# edit element_tags.tsv with your element IDs
python3 src/05.element_phylo.py --method both
```

### Run Modes

| Mode | Behavior |
|------|----------|
| **submit** (default) | Generate run.sh → execute immediately |
| **--dry-run** | Generate run.sh only → `bash run_all.sh` manually |

## SNP / INDEL Detection

```bash
cnex assemble out/ --mers mers_table.tsv --snp -o contigs.fa
```

Output files:

| File | Description |
|------|-------------|
| `contigs.fa` | Assembled contigs |
| `variants.tsv` | SNP/INDEL candidates (ele_id, pos, ref, alt, coverage, frequency) |
| `snp_elements.gfa` | GFA 1.0 graph with bubble paths |

## Dependencies

### Build & Runtime

| Tool | Version | Source |
|------|---------|--------|
| g++ | 9+ (C++17) | gcc.gnu.org |
| Python | 3.8+ | python.org |
| pigz | 2.6+ | bundled in `src/pigz/` |
| robin-map | — | bundled in `src/robin-map/` |

### Phylogenetic Tools

| Tool | Method | Default path |
|------|--------|-------------|
| FAMSA | Both | `famsa` (PATH) |
| FastTree | **astral** only | `FastTree` (PATH) |
| ASTRAL (ASTER) | **astral** only | `astral` (PATH or `--astral-bin`) |
| IQ-TREE 3 | **concat** only | `iqtree3` (PATH or `--iqtree3`) |

## Benchmark

Latest results (human genome, 7772 ground truth CNEs):

| Method | Recall | Precision |
|--------|--------|-----------|
| **CNEX** | **83.1%** | **98.3%** |
| MMseqs2 | 88.2% | 98.4% |
| BLASTN | 89.1% | 98.6% |

## Pipeline Steps

| Step | Description |
|------|-------------|
| 00 | `00.filter_msa.py` — Filter MSA blocks |
| 01 | `mertable` — Build k-mer table |
| 02 | `validate` — k-mer validation |
| 03 | `assemble` — De Bruijn assembly |
| 04 | `04.map_contig_to_element.py` — Group by element |
| 05 | `05.element_phylo.py` — Phylogeny (IQ-TREE / ASTRAL) |

## License

MIT
