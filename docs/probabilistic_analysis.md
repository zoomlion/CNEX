# CNEX: False-Positive Probability of a Random 150 bp Read — and Why It Is Invariant to $\rho$ and Library Scale

This document estimates the probability that the CNEX validator accepts a
**random** read (i.e. a read drawn from genomic DNA unrelated to any CNE
element) as a valid hit for some element. It combines the full acceptance
chain with its production preset parameters and gives the resulting
symbolic upper bound together with its scaling behaviour.

A key finding is that, under the model, this probability is **independent
of the table occupancy $\rho = N_t/4^k$ and of the library size $N_t$**: in
the dominant term the factors $\rho^8$ and $\sum_e p_e^8$ cancel exactly
into $\sum_e (m_e/4^k)^8$ (§3). The estimate is therefore robust to how
many, or how closely related, the species in the panel are: enlarging or
restructuring the library cannot, by itself, inflate the random-read
false-positive probability — only the k-mer content of the elements
themselves matters. This robustness is not incidental: the validation
chain is a single left-to-right pass of rolling-hash lookups
(validator_core.hpp) that **probes the k-mer table but never scans it** —
expected O(1) per read position, independent of $N_t$ and $N_{ele}$. The
probabilistic result and the running time are two sides of the same
property: the library's presence enters neither the bound (§3) nor the
cost.

---

**Contents:** [1. Acceptance chain](#1-acceptance-chain-and-preset-parameters) ·
[2. Model](#2-model) · [3. Symbolic derivation](#3-symbolic-derivation) ·
[4. Symbolic bound](#4-symbolic-bound) ·
[5. Which filter dominates](#5-which-filter-dominates-and-sensitivity) ·
[6. Interpretation](#6-interpretation)

## 1. Acceptance chain and preset parameters

In reads mode, each read of length $R$ is processed independently by
`validate_read` (validator_core.hpp), which runs the forward read and its
reverse complement through `gate_pass()` and keeps the strand with more
matches. A read is accepted for element $e$ iff the following four filters
all pass for the winning element:

1. **Hit count.** $n$ = number of the read's k-mers present in the k-mer
   table. (Rolling-hash lookup; mismatching k-mers contribute nothing.)
2. **Vote majority.** The per-element vote histogram $c_1 \ge c_2 \ge \cdots$
   must satisfy $c_1 > f\, n$ and $c_1 \ge r\, c_2$, with the fraction
   threshold $f = 0.1$ and ratio $r = 3$.
3. **Collinearity.** The $m$ matched k-mers of the winning element are
   scanned in read order (their loci are *not* sorted); among the $m-1$
   consecutive pairs, at least $\mu = 7$ must satisfy
   $\left|\,(\text{read pos}_{i+1} - \text{read pos}_{i}) - (\text{loci}_{i+1} - \text{loci}_i)\right| \le 2$,
   so read order and locus order must roughly agree
   (`gate_pass()`, `hip/validator_core.hpp`).
4. **Span.** The matched positions of the winning element must span at
   least $\ell = 25$ bp.

**Algorithmic cost.** Each read is processed once in a single pass: a
rolling hash (expected O(1) per position) probes the k-mer table, a hash
map (`VoteTally`, `hip/validator_core.hpp`), so validating
a read costs O(R) with O(1) expected per probe — independent of $N_t$ and
$N_{ele}$. The vote step (`VoteTally`, `hip/validator_core.hpp`) sorts at
most $R - k + 1 = 138$ hit ids, O(R log R), again constant in the library
size. The chain never builds an index over the library, never aligns, and
never scans other elements: the table is only probed.

**Table statistics used.** $k = 13$, $R = 150$, so $n_w = R - k + 1 = 138$
k-mers per read. The table contains $N_t$ distinct k-mers over $N_{ele}$
elements; element $e$ owns $m_e$ of them, $p_e = m_e/N_t$. The occupancy
of the table in k-mer space is

$$\rho = \frac{N_t}{4^{k}} .$$

Note that $\rho$ and the shares $\{p_e\}$ are **not independent**:
$p_e = m_e/N_t$ and $\rho = N_t/4^k$ imply the identity
$\rho^{m} S_m = \sum_e (m_e/4^{k})^{m}$ for every $m$. The only free
quantities are the elements' own k-mer contents $\{m_e\}$, the element
count $N_{ele}$, and $k$ — the library size $N_t$ drops out everywhere.
Shared k-mers between elements make $\sum_e p_e \ge 1$; the union bound
below tolerates this as an over-count, so the resulting bound stays
conservative.

Table occupancy is $\rho = N_t/4^k$ and element $e$ owns a share
$p_e = m_e/N_t$. Ownership shares are typically skewed: a few long elements
carry a large fraction of the table, the rest carry a long tail.

---

## 2. Model

For a random read (unrelated to any element), the accepted modelling is
assumption **A4**: each of the $n_w$ read k-mers independently hits the
table with probability $\rho$, and a hit lands on element $e$ with
probability $p_e$ (independently across hits). No element is a "true"
origin, so **any** accepted read is a false positive. The independence of
overlapping k-mers (neighbouring 13-mers share 12 bp) is an idealization;
the resulting caveats are stated in Lemma 3 and §4.

---

## 3. Symbolic derivation

### Lemma 1 — hit count

$$n \sim \mathrm{Binomial}(n_w,\, \rho).$$

### Lemma 2 — vote concentration (voting filters)

Because the collinearity step needs at least $\mu = 7$ passing consecutive
pairs, the winning element must have at least $m = \mu + 1 = 8$ matched
k-mers; we bound the corresponding vote event (winning element holding
$m_0 = 8$ of the $n$ hits) by the union bound, ignoring the ratio
condition on the remaining votes (a relaxation, hence an upper bound):

$$\Pr(\exists\, e : c_e \ge m_0 \mid n) \;\le\; \binom{n}{m_0}\; S_{m_0}, \qquad S_{a} = \sum_{e=1}^{N_{ele}} p_e^{\,a}.$$

This is a **relaxation** of the true vote rule
(`VoteTally`, `hip/validator_core.hpp`), which also
requires $c_1 > 0.1\, n$ and $c_1 \ge 3\, c_2$: with $c_1 = 8$ votes the
fraction condition holds only for $n \le 79$, and the ratio condition
further restricts the runner-up. Dropping both conditions can only enlarge
the probability, so the result is an upper bound, not an exact
probability — it certifies only the *necessary condition* $c_1 \ge 8$.
The full acceptance event for element $e$ is

$$c_1 \ge 8 \;\wedge\; c_1 > 0.1\, n \;\wedge\; c_1 \ge 3\, c_2 \;\wedge\; \text{collinear} \;\wedge\; \text{span},$$

with $c_1, c_2$ the two largest per-element vote counts; every term except
$c_1 \ge 8$ is deliberately relaxed in the Theorem. The ratio rule
$c_1 \ge 3\, c_2$ deserves emphasis: it is what stops a giant element
from winning on hit count alone when a second element shares many of its
k-mers ($p_A \approx p_B$); a split like $c_A = 8$, $c_B = 6$ fails
($8 < 3 \cdot 6$). $S_a$ is the $a$-th power sum of the ownership shares; for
$a \ge 3$, $S_a \le p_{\max}^{\,a-2}\, S_2$, so the tail collapses once
$m_0$ exceeds a few votes.

### Lemma 3 — collinearity and span

No closed-form gap model is used: the validator evaluates the actual
match set. For the winning element with $m$ matched k-mers, the read
positions $x_1 < \cdots < x_m$ are ordered by construction (the read is
scanned left to right), while the loci $l_1, \dots, l_m$ are **kept in
read order and never sorted** (`gate_pass()`) —
consecutive differences $l_{i+1}-l_i$ may be negative or arbitrarily
large. The collinearity indicators are

$$I_i \;=\; \mathbf{1}\left\{\left|(x_{i+1}-x_i) - (l_{i+1}-l_i)\right| \le 2 \right\}, \qquad i = 1,\dots,m-1,$$

and the chain passes iff $\sum_i I_i \ge \mu = 7$ and the span
$x_m - x_1 + k \ge 25$ (see `gate_pass()`). The effective requirement is
therefore that read order and locus order roughly agree; an
exponential-gap approximation (used in an earlier draft) does not
correspond to this implementation and is abandoned. The span
half-condition alone is exact and simple: with $k = 13$ it is
$x_8 - x_1 \ge 12$ for 8 read order statistics of $[0,\, R-k]$.

**Validator-matched Monte-Carlo protocol (subset level).** For element
$e$, the sampling frame is its $m_e$ table rows — pairs (k-mer, locus),
with repeated k-mers (multiple rows) allowed. One trial: draw 8 rows
uniformly **without replacement**; draw the read positions $x_1 < \cdots < x_8$ as order statistics of $[0, 137]$; evaluate the exact checks above
with the loci exactly as drawn, i.e. in read order, unsorted. Repeat $N$
times per element; with $X$ passes report the one-sided 95% upper
confidence bound ($3/N$ for $X = 0$; Clopper–Pearson otherwise).
$P_{\mathrm{coll}}$ is a **subset-level** probability: the pass
probability of one *specific* 8-match subset, which is exactly the object
multiplied by $\binom{n}{8}$ in the Theorem. Sorting the loci before the
checks would silently reintroduce the favourable permutation the
implementation never performs.

### Theorem — false-positive probability of a random read

Applying the chain decomposition,

$$P_{\mathrm{FP}}(R) \;\le\; \sum_{n=0}^{n_w} \binom{n_w}{n}\, \rho^{n}\, (1-\rho)^{n_w-n} \cdot \binom{n}{8}\, S_{8} \cdot P_{\mathrm{coll}} .$$

The sum factors through the binomial identity
$\sum_n \binom{n_w}{n} \rho^n (1-\rho)^{n_w-n} \binom{n}{8} = \binom{n_w}{8} \rho^8$, so with
$\rho^8 S_8 = \sum_e (m_e/4^k)^8$,

$$P_{\mathrm{FP}}(R) \;\le\; \binom{n_w}{8}\; \sum_{e} \left(\frac{m_e}{4^{k}}\right)^{\!8} \; P_{\mathrm{coll}} .$$

In this last form **$N_t$ and $\rho$ have cancelled completely** — the
bound depends only on the elements' own k-mer contents $\{m_e\}$, the
k-mer space size $4^k$, and the collinearity/span factor.
$P_{\mathrm{coll}}$ is the subset-level probability of Lemma 3 (one
specific 8-match subset passing collinearity and span; $m = 8$ is the
dominant least-count term); each additional matched k-mer contributes a
further factor of order $\rho\, p_e$, so the $m_0 = 8$ term dominates. Both strands run the same chain independently
(`validate_read()`, `hip/validator_core.hpp`) and are identically distributed for a
random read; the read-level probability is exactly $2p - p^2 \le 2p$, and
this is the quantity used in §4. Because the fraction and ratio rules are
relaxed (Lemma 2), the whole bound is a **necessary-condition bound**: it
certifies that some element receives 8 collinear-spanning hits, not that
it wins the vote majority.

---

## 4. Symbolic bound

Combining Lemmas 1–3, the probability that a specific element collects at
least $m_0 = 8$ hits is bounded by the vote-concentration term, and the
collinearity/span filter multiplies it by the subset probability
$P_{\mathrm{coll}}$. Summing over elements and over both strands gives

$$P_{\mathrm{FP}}(R) \;\le\; 2\, \binom{n_w}{8}\, \rho^{8} S_{8}\, P_{\mathrm{coll}} \;=\; 2\, \binom{n_w}{8}\, \sum_{e} \left(\frac{m_e}{4^{k}}\right)^{8} P_{\mathrm{coll}} ,$$

with $n_w = R-k+1$ read k-mers, $\rho = N_t/4^k$, $S_8 = \sum_e p_e^8$ and
$p_e = m_e/N_t$. The first form shows where library occupancy enters; the
second uses $\rho^8 S_8 = \sum_e (m_e/4^k)^8$ and shows that it cancels.
The expected count for a genome of length $G$ follows as

$$\mathbb{E}[\mathrm{FP}] \;\le\; \frac{G}{R}\, P_{\mathrm{FP}}(R) .$$

Both expressions are upper bounds: the vote term relaxes the fraction and
ratio rules of Lemma 2, and $P_{\mathrm{coll}}$ bounds a single $m_0$-hit
subset (Lemma 3). They need no table measurements — only the per-element
k-mer counts $\{m_e\}$ and, for $P_{\mathrm{coll}}$, the read length and the
spacing tolerance.

**Scaling.** At fixed element content the bound grows linearly in
$\sum_e (m_e/4^k)^8$, falls as $4^{-8k}$ with the k-mer length, and is
independent of $N_t$ and $\rho$: enlarging the table dilutes every $p_e$ at
the same rate as it raises $\rho^8$. The table therefore enters only through
how much k-mer content individual elements contribute, not through its size.

## 5. Which filter dominates, and sensitivity

1. **Collinearity is the decisive suppressor.** The vote-concentration term
   already requires several matched k-mers on one element; the collinearity
   filter additionally demands that their spacing against the element's loci
   agree within 2 bp. That second requirement is what makes the bound small:
   a random read must both land repeatedly on the same element and preserve
   that element's coordinate ordering.
2. **The vote-ratio rule ($c_1 \ge 3\, c_2$) breaks the "longest element
   dominates" inference.** The $\sum_e (m_e/4^k)^8$ term ranks which
   elements are most likely to collect the required hits, but acceptance also
   demands a 3:1 separation from the runner-up (Lemma 2), a condition the
   bound relaxes. A long element therefore wins disproportionately often only
   when its k-mers are mostly private ($p_A \gg p_B$); two largely shared
   elements ($p_A \approx p_B$, e.g. paralogs) disqualify each other even when
   both collect many hits. The element-length tail predicts where
   false-positive *candidates* concentrate, not which elements are accepted.
3. **Invariance to $\rho$ and library scale.** The bound is
   $O(\sum_e (m_e/4^k)^8)$ and contains no $\rho$ and no $N_t$: the
   $\rho^8 S_8$ combination cancels exactly (§3). Doubling the table occupancy
   cannot multiply $P_{\mathrm{FP}}$ by $2^8$, because enlarging $N_t$ at
   fixed element content dilutes every $p_e$ at the same rate. Even a fully
   occupied table ($\rho = 1$) does not raise the bound: hits are then spread
   over all elements rather than concentrated on one. The effective knobs are
   the element k-mer counts, the k-mer length (entering as $4^{-8k}$), and the
   vote/collinearity thresholds.
4. **Sensitivity to $\mu$ (min-c).** min-c counts *passing consecutive
   pairs* (`gate_pass()`); it is not a vote threshold. Requiring one more
   passing pair multiplies the collinearity factor by roughly the per-pair
   pass probability and, since $m$ matches yield only $m-1$ pairs, forces one
   extra matched k-mer as well. min-c is therefore the coarse control knob for
   the random-read false-positive probability.

## 6. Interpretation

The bound separates two sources of false positives. Random genomic DNA can
contribute only the amount governed by the expression above, and that amount
is invariant to how large or how diverse the k-mer library is. Any observed
false positives beyond it must therefore come from sequence that genuinely
shares k-mers with an element — paralogous loci, duplicated or repeated
conserved elements, and cross-element shared k-mers. The measured specificity
of the pipeline thus reflects how much real sequence sharing exists among the
panel species, not the size of the library.
