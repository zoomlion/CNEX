# CNEX — Conserved Non-coding Element Discovery

Fast identification of conserved non-coding elements (CNEs) from whole-genome sequences. Supports both sequencing reads and direct genome FASTA input. Key features:

- **Ultra-fast sliding-window genome mode** (~33 s validation for human at 6X/16 threads; full validate+assemble pipeline ~3 min)
- **k-mer guided assembly** — De Bruijn graph with element-specific k-mer prior
- **SNP/INDEL detection** via De Bruijn bubble scanning
- **Phylogeny pipeline** — concat with IQ-TREE 3 / RAXML-NG

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

`config.py` is a local file (git-ignored); `config.example.py` holds the tracked
defaults.

Key settings in `config.py`:

| Setting | Default | Description |
|---------|---------|-------------|
| `THREADS` | 20 | Worker count: MAFFT parallelism, IQ-TREE / RAXML-NG threads |
| `MIN_CNE_PER_SPECIES` | 100 | Minimum CNEs per species to retain |
| `DEFAULT_METHOD` | `concat` | `concat`, `astral`, or `both` |
| `MAFFT_MODE` | `auto` | MAFFT mode: `ginsi` (high precision) / `auto` / `fftn2` (fastest) |
| `IQTREE_MODEL` | `MFP` | IQ-TREE model when `--species-file` is used (small set, model search) |
| `IQTREE_MODEL_FULL` | `GTR+F+R4` | IQ-TREE model for the full species run (large set, FreeRate) |
| `RAXML_MODEL` | `GTR+R4` | RAXML-NG single model (unpartitioned) |
| `RAXML_BS` | 200 | RAXML-NG bootstrap replicates |
| `ELEMENT_TAGS_FILE` | `""` | Path to element_tags.tsv for per-type filtering |
| `PARTITION` | `True` | Output partition file for IQ-TREE |
| `DRY_RUN` | `True` | Generate scripts only; use `--submit` to execute |

## Quick Start

```bash
# 1. Multi-species: validate + assemble
#    --cne builds the k-mer table from a reference CNE MSA; alternatively
#    pass an existing table with --mers mers_table.tsv
cnex pipeline genomes/ --cne blocks_10k.fa -t 16 -j 8 -o results --trim

# 2. Group per-species contigs into per-element multi-species FASTAs
#    (--max-elements 0 = all elements; the default processes only the first 100)
python3 src/04.map_contig_to_element.py --max-elements 0 --parallel 8

# 3. Align + build trees
python3 src/05.element_phylo.py --max-elements 0            # dry-run: write run.sh scripts
python3 src/05.element_phylo.py --max-elements 0 --submit   # execute
```

For a single genome instead of a whole set:

```bash
cnex validate genome.fa --mers mers_table.tsv -t 16 -o out/
cnex assemble out/ --mers mers_table.tsv -o contigs.fa
```

For paired-end reads, the validated reads are written as `Assemble.*.reads`
and can be assembled directly:

```bash
cnex validate reads_R1.fq.gz reads_R2.fq.gz \
  --mers mers_table.tsv --type fastq -t 4 \
  --rescue-panels blocks.fa --rescue-tau 0.90 --t3 "1.5,5,25" \
  -o out/
cnex assemble out/ --mers mers_table.tsv -o contigs.fa
```

`blocks.fa` is the leave-one-out block MSA used for the T3 panel rescue;
`--rescue-tau` is required when `--rescue-panels` is given.

### Key CLI parameters

`cnex` passes only explicit overrides to the C++ binaries; the defaults below
are the binary defaults.

| Parameter | Default | Meaning |
|---|---:|---|
| `validate -t` | 4 | threads |
| `validate --type` | auto | `genome` or `fastq` |
| `validate --vote-frac` | 0.05 | minimum vote fraction for the top element |
| `validate --max-assign` | 3 | drop k-mers with more than this many element/locus assignments |
| `validate --t3` | `1.5,5,25` | T3 near-miss gate row (`rv,mc,ms`) |
| `validate --tiers` | none | extra OR-gate rows appended to T1+T2 |
| `validate --rescue-panels` | off | per-element panels for SW rescue |
| `validate --rescue-tau` | — | normalized SW threshold; required with panels |
| `validate --window-size` / `--step-size` | 150 / 25 | genome sliding window (~6X) |
| `assemble -k` | 19 | single-k mode; reads mode uses `--multi-k` by default |
| `assemble --multi-k` | 27 | reads-mode k list |
| `assemble --min-count` | 1 | minimum (k+1)-mer occurrence |
| `assemble --trim` | on | trim to confident k-mer region; use `--no-trim` to disable |
| `assemble --trim-margin` | 10 | bp kept beyond the first/last element k-mer |
| `assemble --max-assign` | 3 | must match validate |
| `assemble --max-reads` | 200 | max reads per element/locus |
| `assemble --max-loci-gap` | 5000 | max gap for genome-mode locus clustering |
| `assemble --snp` | off | scan candidate SNPs from validated reads |

### Input formats

`--cne` expects a `###`-delimited block MSA (each block a FASTA alignment,
headers `>Species.contig[+-]:start-end`). Whole-genome alignments in **MAF**
format (UCSC Multiz/TBA, or Cactus exported with `hal2maf`) are also
supported: convert them first with

```bash
python3 utils/maf_convert.py alignment.maf          # -> most-cons-cne.fa (.maf.gz ok)
python3 src/00.filter_msa.py most-cons-cne.fa       # -> most-cons-cne.filtered.fa
cnex pipeline genomes/ --cne most-cons-cne.filtered.fa ...
```

`maf_convert.py` rewrites each MAF block as a `###`-delimited FASTA block
(coordinates preserved); use `-o FILE` to change the output path.

## Phylogeny Pipeline

### Workflow

```
per-element FASTA
    → MAFFT              (one MSA per element,  results/aln/all/)
    → trimal             (drop gap-rich columns → *.trimmed.aln)
    → concat_msa         (supermatrix per tag: all / intergenic / intron)
    → IQ-TREE 3 + RAXML-NG
```

`--species-whitelist FILE` (one species per line) restricts tree building to a
subset of species; empty = use all species.

### Per-Tag Supermatrix

One supermatrix per tag (all/intergenic/intron), built without quantile filtering.
Per-element species-coverage filtering (`--concat-cov`, default 75) keeps only
elements present in ≥75% of species; `--concat-cov 0` disables the filter. Each
tag emits a `full/` supermatrix and (when the filter is on) a `cov75/` one:

- `--concat-tool` (default `both`): write IQ-TREE 3 and/or RAXML-NG commands into
  `run.sh`. IQ-TREE keeps `-p` partitions with model selection (`MFP`, or
  `GTR+F+R4` for large sets); RAXML-NG runs one unpartitioned model
  (`--raxml-model`, default `GTR+R4`) with `--raxml-bs` bootstrap replicates
  (default 200).
- Tree outputs are prefixed per tool: `iqtree.*` / `raxml.*`.

```
results/iqtree/{type}/
├── all/full/run.sh      # full supermatrix
├── all/cov75/run.sh     # coverage-filtered supermatrix
└── run_all.sh
```

### Element Tags (Optional)

Classify elements as intron vs intergenic (elements overlapping exons are dropped):

```bash
# --gff-dir is optional: GFF3 files named by species (*.gff3, .gz ok)
python3 utils/classify_elements.py \
    --msa blocks_10k.fa \
    --gff-dir gff/ \
    -o element_tags.tsv
```

Without `--gff-dir` (or if no GFF matches), every element is labeled
`intergenic`. Then set `ELEMENT_TAGS_FILE = "element_tags.tsv"` in `config.py`,
or pass `--element-tags element_tags.tsv` on the CLI.

With tags, each type (`all`, `intergenic`, `intron`) gets its own subdirectory
and an independent supermatrix.

### Run Modes

| Mode | Behavior |
|------|----------|
| **dry-run** (default) | Generate `run.sh` + `run_all.sh` only; execute via `bash run_all.sh` |
| **`--submit`** | Generate and execute immediately |

## Output

| File | Description |
|------|-------------|
| `results/aln/*.trimmed.aln` | Per-element trimmed alignments (trimal) |
| `results/iqtree/{type}/all/full/iqtree.treefile` | IQ-TREE tree (prefix `iqtree.`) |
| `results/iqtree/{type}/all/full/raxml.raxml.bestTree` | RAXML-NG tree (prefix `raxml.`) |
| `variants.tsv` | SNP/INDEL candidates (from `--snp`) |
| `snp_elements.gfa` | GFA 1.0 graph with bubble paths |

## Dependencies

| Tool | Method | Default path |
|------|--------|-------------|
| MAFFT | phylogeny | `mafft` (aligner; `--mafft-mode`: ginsi high-precision / auto default / fftn2 fastest) |
| IQ-TREE 3 | phylogeny | `iqtree3` |
| RAXML-NG | phylogeny | `raxml-ng` (single unpartitioned model; `--raxml-model` GTR+R4, `--raxml-bs` 200; https://github.com/amkozlov/raxml-ng) |
| trimal | phylogeny | `trimal` (alignment trimming) |
| bedtools | classification | `bedtools` (for `classify_elements.py`) |

## Benchmark

Genome-mode validation takes ~33 s for the human genome at 6X with 16 threads
in our benchmark; the full validate+assemble pipeline is ~3 min. The full
benchmark (accuracy and timing against MMseqs2 / BLASTN) and the reads-depth
benchmark against Megahit+MMseqs2 live in `benchmark/` — see
`benchmark/README.md`.

## License

MIT
