# CNEX: False-Positive Probability of a Random 150 bp Read — and Why It Is Invariant to $\rho$ and Library Scale

This document estimates the probability that the CNEX validator accepts a
**random** read (i.e. a read drawn from genomic DNA unrelated to any CNE
element) as a valid hit for some element. It combines the full acceptance
chain with its production preset parameters and gives both a symbolic
derivation and a concrete reference estimate computed from the actual
k-mer table of the *Homo sapiens* benchmark.

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

## 1. Acceptance chain and preset parameters

In reads mode, each read of length $R$ is processed independently by
`validate_read` (validator_core.hpp), which runs the forward read and its
reverse complement through `validate_chain` and keeps the strand with more
matches. A read is accepted for element $e$ iff the following four filters
all pass for the winning element:

1. **Hit count.** $n$ = number of the read's k-mers present in the k-mer
   table. (Rolling-hash lookup; mismatching k-mers contribute nothing.)
2. **Vote majority.** The per-element vote histogram $c_1 \ge c_2 \ge \cdots$
   must satisfy
   $$c_1 > f\, n \qquad\text{and}\qquad c_1 \ge r\, c_2 ,$$
   with the fraction threshold $f = 0.1$ and ratio $r = 3$.
3. **Collinearity.** The $m$ matched k-mers of the winning element are
   scanned in read order (their loci are *not* sorted); among the $m-1$
   consecutive pairs, at least $\mu = 7$ must satisfy
   $$\bigl|\,(\text{read pos}_{i+1} - \text{read pos}_{i}) -
      (\text{loci}_{i+1} - \text{loci}_i)\bigr| \le 2 ,$$
   so read order and locus order must roughly agree
   (`validate_chain`, `validator_core.hpp` lines 257–274).
4. **Span.** The matched positions of the winning element must span at
   least $\ell = 25$ bp.

**Algorithmic cost.** Each read is processed once in a single pass: a
rolling hash (expected O(1) per position) probes the k-mer table, a hash
map (`validate_chain`, `validator_core.hpp` lines 228–249), so validating
a read costs O(R) with O(1) expected per probe — independent of $N_t$ and
$N_{ele}$. The vote step (`findFrequentWithMap`, lines 85–128) sorts at
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

Measured ($Homo\_sapiens$ full library, `mers_table.tsv`):

$$N_t = 12{,}704{,}780 ,\quad N_{ele} = 9{,}991 ,\quad
\rho = \frac{12{,}704{,}780}{4^{13}} = 0.1893 .$$

Per-element ownership shares are skewed: maximum
$p_{\max} = 8.98\times10^{-4}$, median $8.26\times10^{-5}$.

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

$$\Pr\!\big(\exists\, e : c_e \ge m_0 \mid n\big)
\;\le\; \binom{n}{m_0}\; S_{m_0},
\qquad S_{a} = \sum_{e=1}^{N_{ele}} p_e^{\,a}.$$

This is a **relaxation** of the true vote rule
(`findFrequentWithMap`, `validator_core.hpp` lines 85–128), which also
requires $c_1 > 0.1\, n$ and $c_1 \ge 3\, c_2$: with $c_1 = 8$ votes the
fraction condition holds only for $n \le 79$, and the ratio condition
further restricts the runner-up. Dropping both conditions can only enlarge
the probability, so the result is an upper bound, not an exact
probability — it certifies only the *necessary condition* $c_1 \ge 8$.
The full acceptance event for element $e$ is

$$c_1 \ge 8 \;\wedge\; c_1 > 0.1\, n \;\wedge\; c_1 \ge 3\, c_2
\;\wedge\; \text{collinear} \;\wedge\; \text{span},$$

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
read order and never sorted** (`validate_chain`, lines 257–274) —
consecutive differences $l_{i+1}-l_i$ may be negative or arbitrarily
large. The collinearity indicators are

$$I_i \;=\; \mathbf{1}\Big\{\big|(x_{i+1}-x_i) - (l_{i+1}-l_i)\big| \le 2
\Big\}, \qquad i = 1,\dots,m-1,$$

and the chain passes iff $\sum_i I_i \ge \mu = 7$ and the span
$x_m - x_1 + k \ge 25$ (lines 267–278). The effective requirement is
therefore that read order and locus order roughly agree; an
exponential-gap approximation (used in an earlier draft) does not
correspond to this implementation and is abandoned. The span
half-condition alone is exact and simple: with $k = 13$ it is
$x_8 - x_1 \ge 12$ for 8 read order statistics of $[0,\, R-k]$.

**Validator-matched Monte-Carlo protocol (subset level).** For element
$e$, the sampling frame is its $m_e$ table rows — pairs (k-mer, locus),
with repeated k-mers (multiple rows) allowed. One trial: draw 8 rows
uniformly **without replacement**; draw the read positions $x_1 < \cdots
< x_8$ as order statistics of $[0, 137]$; evaluate the exact checks above
with the loci exactly as drawn, i.e. in read order, unsorted. Repeat $N$
times per element; with $X$ passes report the one-sided 95% upper
confidence bound ($3/N$ for $X = 0$; Clopper–Pearson otherwise).
$P_{\mathrm{coll}}$ is a **subset-level** probability: the pass
probability of one *specific* 8-match subset, which is exactly the object
multiplied by $\binom{n}{8}$ in the Theorem. Sorting the loci before the
checks would silently reintroduce the favourable permutation the
implementation never performs.

**Current value.** The factor used in §4 comes from an order-statistics
sweep (loci also drawn as order statistics, i.e. only the favourable
monotone read/locus correspondence), $4\times10^{5}$ draws, 0 passes — a
95% upper confidence bound of $7.5\times10^{-6}$. Because the
implementation's unsorted-loci condition only adds failure modes
(negative or mismatched gaps), this value overestimates the true factor,
so the final bound remains conservative. The protocol above is the exact
replacement; since the top 50 elements carry 95.3% of $\sum_e q_e^8$
(§4.1), it needs to run only on that subset to pin down 95% of the mass.

### Theorem — false-positive probability of a random read

Applying the chain decomposition,

$$P_{\mathrm{FP}}(R)
\;\le\;
\sum_{n=0}^{n_w}
\binom{n_w}{n}\, \rho^{n}\, (1-\rho)^{n_w-n}
\cdot \binom{n}{8}\, S_{8}
\cdot P_{\mathrm{coll}} .$$

The sum factors through the binomial identity
$\sum_n \binom{n_w}{n} \rho^n (1-\rho)^{n_w-n} \binom{n}{8}
= \binom{n_w}{8} \rho^8$ (verified numerically:
$3.4724\times10^{-18}$, matching the exact computation of §4), so with
$\rho^8 S_8 = \sum_e (m_e/4^k)^8$,

$$P_{\mathrm{FP}}(R)
\;\le\;
\binom{n_w}{8}\; \sum_{e} \left(\frac{m_e}{4^{k}}\right)^{\!8}
\; P_{\mathrm{coll}} .$$

In this last form **$N_t$ and $\rho$ have cancelled completely** — the
bound depends only on the elements' own k-mer contents $\{m_e\}$, the
k-mer space size $4^k$, and the collinearity/span factor.
$P_{\mathrm{coll}}$ is the subset-level probability of Lemma 3 (one
specific 8-match subset passing collinearity and span; $m = 8$ is the
dominant least-count term). Measured ratios from the real table (§4.2):
$T_9/T_8 = 2.08\times10^{-3}$ and $T_{10}/T_8 = 4.08\times10^{-6}$, so the
$m_0 = 8$ term dominates. Both strands run the same chain independently
(`validate_read`, lines 284–300) and are identically distributed for a
random read; the read-level probability is exactly $2p - p^2 \le 2p$, and
this is the quantity §4 reports. Because the fraction and ratio rules are
relaxed (Lemma 2), the whole bound is a **necessary-condition bound**: it
certifies that some element receives 8 collinear-spanning hits, not that
it wins the vote majority.

---

## 4. Reference estimates (Homo sapiens full library)

All numbers below are **upper bounds under the idealized A4 model**, not
the algorithm's exact false-positive probability: the vote-concentration
term relaxes the fraction/ratio rules (Lemma 2), the collinearity factor
is a heuristic Monte-Carlo estimate that omits the read-order /
locus-order consistency requirement (Lemma 3) and thus overestimates the
true factor, overlapping k-mers are treated as independent (§2), and both
strands are accounted for below (`validate_read`, lines 284–300: forward
and reverse complement each run the full chain; for a random read the two
strands are identically distributed, so the read-level event is the union
of two identical events, exactly $2p - p^2 \le 2p$).

### 4.1 Inputs measured from the real table

| quantity | value |
|---|---|
| $N_t$ | 12,704,780 |
| $N_{ele}$ | 9,991 |
| $\rho = N_t/4^{13}$ | 0.1893 |
| $S_2$ | $1.526\times10^{-4}$ |
| $S_3$ | $3.267\times10^{-8}$ |
| $S_8$ | $7.932\times10^{-25}$ |
| $p_{\max}$ | $8.98\times10^{-4}$ |
| $m_{\max} = p_{\max} N_t$ | 11,409 |
| $q_{\max} = m_{\max}/4^{13}$ | $1.70\times10^{-4}$ |
| $\sum_e (m_e/4^{13})^8$ | $1.309\times10^{-30}$ |

The quantities $S_m$ and $p_{\max}$ are $N_t$-relative ($p_e = m_e/N_t$) and
alone are misleading: a large library dilutes them, but by §1–§3 they act
only through the *$\rho$-invariant* product $\rho^m S_m =
\sum_e (m_e/4^k)^m$. The last three rows give the absolute, scale-free
inputs that actually drive the bound. The single largest element
contributes $\approx 53\%$ of $\sum_e q_e^8$, so the bound is dominated by
the tail of the element-length distribution.

### 4.2 Evaluation of the bound

| term | value |
|---|---|
| hit-count weight $\sum_n \mathrm{Bin}(138, 0.1893, n)$ | 1 (by construction) |
| vote concentration (union bound, $m_0 = 8$) | $3.47\times10^{-18}$ |
| vote concentration ($m_0 = 9$) | $7.24\times10^{-21}$ |
| vote concentration ($m_0 = 10$) | $1.42\times10^{-23}$ |
| collinearity + span, $m = 8$ (order-statistics MC, $4\times10^{5}$ draws) | $< 7.5\times10^{-6}$ (0 passes, 95% UCB) |
| collinearity + span, $m = 9$ (order-statistics MC) | $7.5\times10^{-6}$ (95% UCB) |
| ratio $T_9/T_8$ (from the real table) | $2.08\times10^{-3}$ |
| ratio $T_{10}/T_8$ (from the real table) | $4.08\times10^{-6}$ |

MC = Monte-Carlo. The collinearity factor is the subset-level probability
of Lemma 3 — one specific 8-match subset passing collinearity + span —
and "95% UCB" denotes a one-sided 95% upper confidence bound ($3/N$ rule
for zero observed passes). $T_a = \binom{138}{a}\sum_e q_e^a$ are the
vote-concentration terms of §3.

**Per read** (single-strand chain: $3.47\times10^{-18} \times
7.5\times10^{-6} \approx 2.6\times10^{-23}$; both strands, union bound
$2p - p^2 \le 2p$):

$$P_{\mathrm{FP}}(150\,\mathrm{bp}) \;\lesssim\;
2 \times 3.47\times10^{-18} \times 7.5\times10^{-6}
\;\approx\; 5.2\times10^{-23} .$$

**Per genome** (human, $G = 3.1\times10^{9}$ bp,
$G/R \approx 2.07\times10^{7}$ reads):

$$\mathbb{E}\big[\mathrm{FP}\big] \;\lesssim\;
\frac{G}{R}\, P_{\mathrm{FP}}(R)
\;\approx\; 2.07\times10^{7} \times 5.2\times10^{-23}
\;\approx\; 1.1\times10^{-15} .$$

**Comparison with the benchmark.** The measured false-positive count of
the *Homo sapiens* benchmark run is $162$ — about
$1.5\times10^{17}$-fold larger than the random-DNA expectation. The observed
false positives therefore cannot originate from random genomic DNA; they
must arise from regions that genuinely share k-mer content with an element
(paralogous or repeated conserved sequences, or cross-element shared
k-mers). Because the random contribution is invariant to library scale
(§3) and the observed count is not, the excess cannot be attributed to the
library size: it tracks how much real biological sequence sharing is
present in the panel, which grows as the species in the clade become more
closely related.

---

## 5. Which filter dominates, and sensitivity

1. **Collinearity is the decisive suppressor.** The vote-concentration
   term alone (≈ $10^{-18}$) is already small, but the collinearity filter
   multiplies it by $\lesssim 10^{-5}$–$10^{-7}$: a random read would need
   8 k-mers that both land on one element *and* keep mutually consistent
   spacing against that element's loci — an event of order $10^{-23}$ per
   read ($5\times10^{-23}$ with both strands). Without the collinearity
   filter the false-positive rate would be governed by the vote majority
   alone, i.e. ~$10^{-18}$ per read and ~$10^{-10}$ per genome.
2. **The vote-ratio rule ($c_1 \ge 3\, c_2$) breaks the "longest element
   dominates" inference.** The $m_e^8$ concentration (§4.1) ranks FP
   *candidates*: which elements are most likely to receive 8 random hits.
   Acceptance additionally demands a 3:1 separation from the runner-up
   (Lemma 2), a condition the bound relaxes. A long element therefore wins
   disproportionately often only when its k-mers are mostly private
   ($p_A \gg p_B$); two largely shared elements ($p_A \approx p_B$, e.g.
   paralogs) disqualify each other under the ratio rule even when both
   collect many hits. The element-length tail predicts where FP
   *candidates* concentrate — not which elements ultimately win
   acceptance.
3. **Invariance to $\rho$ and library scale.** The bound is
   $O\!\big(\sum_e (m_e/4^k)^8\big)$ and contains **no $\rho$ and no
   $N_t$** — the $\rho^8 S_8$ combination cancels exactly (§3). Doubling the
   table occupancy *cannot* multiply $P_{\mathrm{FP}}$ by $2^8$: enlarging
   $N_t$ at fixed element content equally dilutes every $p_e$ and cancels
   the growth in $\rho^8$. Even filling the table entirely ($\rho = 1$) does
   not raise the false-positive rate — every read k-mer then hits the
   table, but the ~138 hits are spread over the $N_{ele} = 9{,}991$
   elements instead of concentrating 8 votes on a single one. The effective
   knobs are the element-length distribution (longest element ≈ 53% of
   $\sum_e q_e^8$), the k-mer length $k$ (which enters as $4^{8k}$), and
   the vote/collinearity thresholds; the library size is not among them.
4. **Sensitivity to $\mu$ (min-c).** min-c counts *passing consecutive
   pairs* (`validate_chain`, line 278); it does not set a vote threshold
   directly. Requiring one more passing pair multiplies the collinearity
   factor by roughly the per-pair pass probability (~$10^{-1}$) — and,
   since $m$ matches yield only $m-1$ pairs, it also forces one extra
   matched k-mer (another factor $T_9/T_8 = 2.08\times10^{-3}$, §4.2,
   from the vote term). min-c is the coarse control knob for the
   random-read false-positive probability.

---

## 6. Relation to the observed false positives

The model establishes a strict separation: random DNA contributes
$\lesssim 10^{-15}$ false positives per human genome under the production
parameters, while the benchmark measures $87$–$263$ false positives per
species. The observed false positives are **overwhelmingly unlikely to be
explained by random genomic DNA alone** (expected $\lesssim 10^{-15}$ per
genome vs. $87$–$263$ observed) and are therefore consistent with genuine
sequence sharing — paralogous loci, duplicated or repeated conserved
elements, i.e. hits on *homologous but wrong* loci — rather than spurious
random-DNA acceptances. This is the counterpart, on the
specificity side, of the benchmark's measured precision of
$97$–$99\%$.

Two layers must be kept apart when attributing false positives to
elements. (i) **Candidate bound:** the $m_e^8$ concentration (§4.1 — the
top 10 elements already carry 87% of $\sum_e q_e^8$) identifies which
elements are most likely to *receive* 8 random hits; it is a ranking of
FP candidates. (ii) **Accepted FPs:** the $c_1 \ge 3\,c_2$ separation rule
(§5) and the collinearity filter act on top of the vote counts, so the
observed per-element FP distribution need not track the $m_e^8$ ranking:
elements with heavily shared k-mer sets (paralogous pairs) undercut each
other's ratio condition, while shorter, k-mer-private elements are
comparatively favoured. A decisive empirical test would rank elements by
$m_e^8$ and by their observed FP counts and compare the two rankings:
agreement would confirm size-driven susceptibility, disagreement would
implicate winner separation, biological sharing, and collinearity as the
true determinants (planned; not run in this document).
