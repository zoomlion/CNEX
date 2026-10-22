# CNEX read validation, assembly, and false-positive control

## Read validation

CNEX validates each read, or each sliding genome window, against a k-mer table of conserved noncoding elements (CNEs) built from the study species. In reads mode, each 150 bp read is scanned with a rolling hash, and every 13-mer found in the table casts a vote for the element that owns it, giving at most 138 votes per read. The vote already acts as a first filter: a read that shares k-mers with no element, or with several at once, is not useful for assembly.

A read is assigned to the element with the largest vote count. It is accepted only if that element satisfies four conditions, which together reject matches that are mere noise, shared by cross-element k-mers (paralogous or repeated conserved sequences), or scattered across the genome.

1) The element receives n matching k-mers, so a read with only a handful of incidental matches is excluded.

2) Its vote count c₁ must exceed 10% of all matches and be at least three times the second-largest count c₂ (c₁ > 0.1n and c₁ ≥ 3c₂). The first inequality ensures the element accounts for a substantial part of the read, and the second that no competing element owns a large share of the same k-mers; a split such as 8 votes for one element and 6 for another therefore fails.

3) The matched k-mers must be collinear with the element sequence. K-mers are scanned in read order and their loci are kept in that order, never sorted, and at least 7 of the c₁−1 consecutive matched pairs must preserve their spacing, |xᵢ₊₁ − xᵢ − (lᵢ₊₁ − lᵢ)| ≤ 2 bp, where xᵢ are the positions of the matched k-mers in the read and lᵢ their loci on the element. This screens out reads whose matches scatter across the genome or fall on rearranged sequence.

4) The matched positions must span at least 25 bp, so acceptance cannot rest on a few matches clustered at one end of the read.

Both strands of each read are processed independently through the same chain, and the strand with the most matches is kept. Reads that pass these checks are assigned to the element and provide the material for sequence assembly.

## T3 rescue

In the production reads pipeline, near-miss reads can be recovered by
Smith-Waterman verification against per-element panels derived from the
leave-one-out block MSA (`--rescue-panels blocks.fa --rescue-tau 0.90`).
The near-miss gate is controlled by `--t3` (default `1.5,5,25`); additional
OR-gate rows can be supplied with `--tiers`. Rescued reads are written with
their normalized SW score and handed to assembly, where they are kept only if
they pass the assembly-level rescue checks.

## De Bruijn assembly

For each element, the validated reads are assembled into a De Bruijn graph. In single-k genome mode the graph uses 19-mers; the reads pipeline uses multi-k assembly with k=27 by default. An edge connects a node to its overlapping successor when the longer context that spans them appears in the reads; edges seen fewer than a minimum number of times (default 1) are removed. Each node is scored by the number of element-specific 13-mers from the CNE k-mer table that it contains, so the graph carries a prior over which paths are likely to belong to the element. A greedy traversal starts from the pair of nodes with the largest in-degree and out-degree imbalances, walks through the BFS distances between candidate starts and ends, and extends forward and backward while choosing the neighbor with the longest remaining distance, breaking ties by the element-specific score. The contig is reconstructed by reading off overlapping k-mers along the path, so no alignment step is required. Because edges supported by fewer than the configured minimum (default 1) are removed, sporadic sequencing errors carried by single reads do not persist into the assembled sequence; the validated reads fold into one consensus sequence per locus rather than propagating their individual errors. Over the five-species benchmark, the combined validation and assembly pipeline reached a precision of 98.0-99.3% and a recall of 84.3-96.5%, consistent with the per-read false-positive bound derived below.

## False-positive probability of a random read

We estimate the probability that a read drawn from genomic DNA unrelated to any CNE is falsely accepted. For the human reference library the table contains Nₜ = 1.27 × 10⁷ distinct k-mers contributed by Nₑₗₑ = 9,991 elements. Each of the 138 k-mers of a 150 bp read hits the table with probability ρ = Nₜ/4¹³ = 0.19; treating the hits as independent, the number of hits on an element with mₑ distinct k-mers follows a binomial distribution.

Under this model the probability that a random read is falsely accepted is bounded by

P_FP ≤ 2 · C(138,8) · P_coll · Σₑ (mₑ/4¹³)⁸ ≈ 5.2 × 10⁻²³ per read    (1)

where P_coll = 7.5 × 10⁻⁶ is the probability that a specific set of eight matched k-mers is collinear and spans at least 25 bp, estimated by Monte Carlo over the real k-mer-to-locus mapping (95% upper confidence bound, 4 × 10⁵ draws). The factor of 2 accounts for the two strands. The term C(138,8) counts the ways eight hits can be placed among the read k-mers, and Σₑ (mₑ/4¹³)⁸ weights the elements by their k-mer content. Within this expression the factors that depend on ρ and on the total library size Nₜ cancel exactly, so the bound contains neither ρ nor Nₜ: enlarging or restructuring the library cannot, by itself, inflate the random false-positive probability, which is set only by the k-mer content of the elements.

## Expected versus observed false positives

With the production parameters, the expected number of false acceptances is fewer than 1.1 × 10⁻¹⁵ per haploid human genome (3.1 × 10⁹ bp). The human benchmark produced 162 false acceptances at a precision of 97-99%, roughly 1.5 × 10¹⁷-fold more than this expectation. The observed false positives therefore cannot originate from random genomic DNA; they must arise from regions that genuinely share k-mer content with an element, such as paralogous or repeated conserved loci. The full derivation, the relaxation audit of the bounds, and the Monte-Carlo protocol are provided in Supplementary Note S1.