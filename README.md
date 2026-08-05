# CNEX — Conserved Non-coding Element Discovery

Fast identification of conserved non-coding elements (CNEs) from whole-genome sequences. Supports both sequencing reads and direct genome FASTA input. Key features:

- **Ultra-fast sliding-window genome mode** (~3 min for human genome at 6X, 16 threads)
- **k-mer guided assembly** — De Bruijn graph with element-specific k-mer prior
- **SNP/INDEL detection** via De Bruijn bubble scanning
- **Phylogeny pipeline** — wASTRAL/ASTRAL (default) with block-gap clustering, or concat/IQ-TREE 3

## Install

```bash
git clone --recurse-submodules git@github.com:zoomlion/CNEX.git
cd CNEX
make && make install
cnex setup
```

## Configuration

```bash
cp config.example.py config.py
```

Key settings in `config.py`:

| Setting | Default | Description |
|---------|---------|-------------|
| `THREADS` | 20 | Worker count: MAFFT/FastTree parallelism, IQ-TREE / ASTRAL threads |
| `MIN_CNE_PER_SPECIES` | 100 | Minimum CNEs per species to retain |
| `DEFAULT_METHOD` | `astral` | `astral`, `concat`, or `both` |
| `ASTRAL_BLOCK_GAPS` | `1000,2000` | kb thresholds for astral block clustering |
| `ELEMENT_TAGS_FILE` | `""` | Path to element_tags.tsv for per-type filtering |
| `PARTITION` | `False` | Output partition file for IQ-TREE |
| `DRY_RUN` | `True` | Generate scripts only; use `--submit` to execute |

## Quick Start

```bash
# 1. Validate + Assemble
cnex validate genome.fa --mers mers_table.tsv --type genome -t 16 -o out/
cnex assemble out/ --mers mers_table.tsv -o contigs.fa

# 2. Element FASTA → Phylogeny
python3 src/04.map_contig_to_element.py
python3 src/05.element_phylo.py              # dry-run (scripts only)
python3 src/05.element_phylo.py --submit     # execute
```

## Phylogeny Pipeline

### Methods

| Method | Flag | Pipeline |
|--------|------|----------|
| **astral** (default) | `--method astral` | MAFFT → block-gap cluster → concat → FastTree → wASTRAL/ASTRAL |
| **concat** | `--method concat` | MAFFT → concat_msa → IQ-TREE 3 |
| **both** | `--method both` | both in one pass |

### ASTRAL: Block-Binning

Elements are grouped into fixed genomic bins (default 1000kb / 2000kb). Each bin is concatenated into a super-locus before FastTree inference. Requires `element_tags.tsv` for coordinates; without it, each element builds an independent gene tree.

```
element_tags.tsv: ele_id→(type, chr, start, end)
    ↓ bin_cluster(bin_size=1000kb)
{chr_0: [ele_1, ...], chr_1000000: [ele_2, ...], ...}
    ↓ concat_block_alignments + FastTree
block_{type}_0.nwk  block_{type}_1.nwk  ...
    ↓ ASTRAL
species_tree.nwk
```

Output structure:

```
results/astral/{type}/
├── block_1000kb/run.sh
├── block_2000kb/run.sh
└── run_all.sh
```

### Concat: Per-Tag Supermatrix

IQ-TREE builds one supermatrix per tag (all/intergenic/intron), without quantile filtering:

```
results/iqtree/{type}/
├── all/run.sh
└── run_all.sh
```

### Element Tags (Optional)

Generate per-element classification (intron vs intergenic) from GFF:

```bash
# Requires GFF3 files named by species in a directory
python3 utils/classify_elements.py \
    --msa blocks_10k.fa \
    --gff-dir gff/ \
    -o element_tags.tsv
```

Then set `ELEMENT_TAGS_FILE = "element_tags.tsv"` in `config.py`, or pass `--element-tags element_tags.tsv` on the CLI.

With tags, each type (`all`, `intergenic`, `intron`) gets its own subdirectory and independent analysis.

### Run Modes

| Mode | Behavior |
|------|----------|
| **dry-run** (default) | Generate `run.sh` + `run_all.sh` only; execute via `bash run_all.sh` |
| **`--submit`** | Generate and execute immediately |

## Output

| File | Description |
|------|-------------|
| `results/aln/*.trimmed.aln` | Per-element trimmed alignments (trimal) |
| `results/astral/{type}/{bin}/species_tree.nwk` | ASTRAL species tree |
| `results/iqtree/{type}/all/supermatrix.fa.treefile` | IQ-TREE tree |
| `variants.tsv` | SNP/INDEL candidates (from `--snp`) |
| `snp_elements.gfa` | GFA 1.0 graph with bubble paths |

## Dependencies

| Tool | Method | Default path |
|------|--------|-------------|
| MAFFT | Both | `mafft` (default aligner, G-INS-i: `--globalpair --maxiterate 1000`) |
| FastTree | **astral** | `FastTree` (`-gtr -gamma -spr 4 -mlacc 2 -slownni -boot 1000`) |
| wASTRAL | **astral** | `wastral` (preferred; falls back to ASTRAL if absent) |
| ASTRAL (ASTER) | **astral** | `astral` |
| IQ-TREE 3 | **concat** | `iqtree3` |
| trimal | Both | `trimal` (alignment trimming) |
| bedtools | classification | `bedtools` (for `classify_elements.py`) |

## Benchmark

| Method | Recall | Precision |
|--------|--------|-----------|
| **CNEX** | **83.1%** | **98.3%** |
| MMseqs2 | 88.2% | 98.4% |
| BLASTN | 89.1% | 98.6% |

## License

MIT
