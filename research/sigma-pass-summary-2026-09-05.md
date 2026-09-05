# Sigma pass: six parallel directions on the two-set Paley conjecture

**Status: the two-set Paley conjecture remains open. This pass proves no
new cancellation exponent. It does establish, with exact verification, where
every known route stops, corrects two claims of earlier passes, refutes four
structural hypotheses about biased rectangles, machine-checks the elementary
layer in Lean 4, and prepares (but does not upload) a prove2me mission.**

Seven notes, seven verifiers, all re-run by the orchestrator on 2026-09-05
(every run exited 0 with zero failures; the referee count differs by 30
from the worker's own run because of its sampled cases):

| direction | note | verifier checks | outcome |
|---|---|---|---|
| barriers (root) | [sigma-barriers](sigma-barriers-2026-09-05.md) | 12,143 exact | seven barriers, subfield diagnostic, corrected constant |
| sumproduct | [sigma-sumproduct](sigma-sumproduct-2026-09-05.md) | 737 exact | no published θ<1/2 for two arbitrary sets; no-go for shift chains |
| structured | [sigma-structured](sigma-structured-2026-09-05.md) | 220,917 exact | even the subgroup case is open; exact equivalence proved |
| crux | [sigma-crux](sigma-crux-2026-09-05.md) | 59,738 exact | dilation-invariance obstruction; H1–H4 refuted; Conjecture SI |
| referee | [sigma-referee](sigma-referee-2026-09-05.md) | 1,417,118 exact (re-run) | passes 9–10 plausible, citation-level gaps, no witness |
| subgroup | [sigma-subgroup](sigma-subgroup-2026-09-05.md) | 307,645 exact | pass-7 bound is Konyagin 2002 in weaker form; class = E_3 hypothesis |
| lean | [sigma-lean](sigma-lean-2026-09-05.md) | 36,820 exact + Lean | four theorems machine-checked; proposal package ready |

The pass ran under one interruption: all six workers were killed by an
account usage limit and relaunched from their partial artifacts.

## 1. What is now proved (PROVED, complete arguments in the notes)

1. **Elementary layer, machine-checked.** In Lean 4.30.0 / Mathlib
   `c5ea003`, with axioms exactly `propext, Classical.choice, Quot.sound`:
   shift orthogonality `Σ_x χ(x)χ(x+c) = −1`; the exact second moment
   `Σ_x (Σ_{b∈B}χ(x−b))² = |B|(p−|B|)`; the Chung bound
   `(Σ_{A,B} χ(a−b))² ≤ |A||B|(p−|B|)`; and the implication from
   `S([1,N],[1,N]) < N²` to a quadratic non-residue in `[2,2N]`.
   Files: `~/prove2me_workspace/Theorems/Thm_paley_*.lean`,
   `Solutions/Sol_paley_*.lean`.
2. **Reductions.** The conjecture implies Vinogradov's least-non-residue
   conjecture `n_p ≤ 2⌊p^ε⌋`, and more strongly that no run of
   `2⌊p^ε⌋−1` consecutive residues or non-residues exists in `[1,p−1]`
   for large `p` (crux note, via the exact triangle-weight identity
   `S([1,N],[1,N]) = Σ_{s=2}^{2N} min(s−1, 2N+1−s) χ(s)`).
3. **The dilation crux.** With `g(t) = S(tA,B)`, `0 ∉ A`, `k ≥ 1`,
   `η ∈ (0,1]`:
   `#{t : |g(t)| ≥ η|A||B|} ≤ η^{−2k}[(2k−1)!!·p/max(|A|,|B|)^k + (2k−1)√p]`.
   The right side never drops below `(2k−1)√p·η^{−2k}` for `k ≥ 2`, and
   never below `√p` for any `k` when `|A||B| ≤ √p`. Every termwise-Weil
   dilation-moment bound depends on `(A,B)` only through the multiset of
   signed ratio-class sums, which is invariant under `A ↦ uA`; hence every
   such bound assigns the same value to all `p−1` dilates and cannot single
   out `t = 1`. The conjecture `C(ε,δ)` is equivalent to the level-set
   count being zero. The even-multiplicity tuple count is
   `N_k ≤ (2k−1)!!·R^k` with `N_2 = 3n²−2n` exactly; the constant `k!`
   used informally elsewhere is wrong.
4. **Subgroup case is exactly the shifted-subgroup problem.** For a
   multiplicative subgroup `H` and `M_H = max_{c≠0}|Σ_{h∈H}χ(h+c)|`: for
   every `A` and every union `B` of `k` cosets of `H`,
   `|S(A,B)| ≤ k(M_H|A| + |H|)`, and conversely the set
   `A_c = c·{h ∈ H : χ(h)=1}` attains `|S(A_c,H)| = M_H|A_c|`. So the
   two-set conjecture for `B = H` holds if and only if the uniform bound
   `M_H ≤ p^{−δ}|H|` holds (structured note, Theorem 2.2).
5. **Explicit Burgess–Chang chain and its no-go.** Theorem E of the
   sum-product note gives, for any `A, B`, shift set `Z`, interval `I`,
   `r ≥ 1`, with `λ = |B − Z·I|/|B|` and `Φ` the ratio-energy deviation,
   `|S(A,B)| ≤ λ|A||B|·Φ^{1/4r}·[((2r−1)!!)^{1/2r}H^{−1/2} + (2r−1)^{1/r}p^{−1/2r}]^{1/2}`,
   and Proposition F shows the right side is always
   `≥ |A||B|·|B|/E⁺(B)^{1/2}`. For additively unstructured `B`
   (`E⁺(B) ≤ C|B|²`) the chain saves at most `√C`, whatever energy or
   incidence bound is inserted. The exponent bookkeeping
   `θ = 1/(8−a−b)` locates Chang's `4/9` exactly at her energy exponent
   `11/4`, and would give `2/5` with the point–plane exponent `5/2` and
   `1/3` with the optimal exponent `2`, all for one-sided structured sets.
6. **Konyagin's inequality and the subgroup line.** For every subgroup
   `H ≤ F_p^*` of order `n`, all `k,l ≥ 1`, and every `a ≠ 0`:
   `|η(a)|^{kl} ≤ n^{kl−k−l}·√(p·E_k(H)·E_l(H))`
   (self-contained proof; cited as Shkredov arXiv:1311.5726, Lemma 7,
   attributed to Konyagin 2002). With `k=l=3`, `p ≤ n⁴`, `E_3 ≤ Bn³` this
   gives `M ≤ B^{1/9}n^{8/9}`, strictly better than pass 7's
   `(B+2)^{1/9}n^{8/9}`. Every exponent in passes 5–7 follows from it with
   the same energy inputs; its optimum over all `(k,l)` with ideal energies
   is exactly `7/8`.
7. **Barrier map.** Seven barriers (Fourier, Hölder–Weil, polynomial
   method, subfield, Burgess, degree-4 SOS, RIP bottleneck), each pinned to
   its lossy step, with the `F_{p²}` subfield test as a decisive
   diagnostic: any step of a proposed proof that holds verbatim over
   `F_{p²}` cannot pass `ε = 1/2`.

## 2. What was refuted (REFUTED, with witnesses in the results files)

- **Chang's `4/9` is not an arbitrary-two-set theorem.** Its hypothesis is
  `|B+B| < K|B|`. No published result gives `θ < 1/2` for two arbitrary
  sets; Fouvry–Shparlinski–Xi (2024) state that Vinogradov's and
  Karatsuba's bounds "still stand". The brief of this pass and earlier
  workspace notes that presumed otherwise are corrected.
- **The shifted-subgroup bound at `|H| = p^ε` is not in the literature.**
  Chang's survey lists it as Bourgain's Problem 5 even at `|H| ≈ √p`; the
  best both-structured result found is Chang–Shparlinski Theorem 4
  (`A = a+G`, `B = G`, `|G| ≥ p^{13/33+ε}`). Consequently the subgroup
  case of the conjecture is open, and the multiplicatively structured case
  (`|B·B| ≤ K|B|`) contains it at `K = 1`.
- **Container-to-subset inheritance fails.** `p = 433`, `|H| = 24`,
  `M_H = 8 = |H|/3`, but the index-3 subgroup `H'` of order 8 has
  `M_{H'} = 8 = |H'|`, so a shifted-subgroup bound for `H` says nothing
  about dense subsets of `H`.
- **Coset bootstrapping of a generic `B` reduces nothing.** `p = 1000003`,
  random `|B| = 251`, `|H| = 6`: all 251 cosets are singletons.
- **Shift chains cannot be rescued by energy.** `p = 2003`, an explicit
  `|B| = 31 = p^{0.45}` with `E⁺(B) = 2.39|B|²`: exhaustively over all
  shifts, `min |B − z{1,2}|/|B| = 1.903`, so Theorem E's bound is
  `≥ 1.60|A||B|` and Proposition F gives `≥ 0.65|A||B|` for every
  `(Z,I,r)`.
- **Structural hypotheses about biased rectangles.** Over 1,500 rectangles
  with `|A||B| ≥ 64` and bias `≥ 1/2`: H1 (small additive doubling) and H2
  (small multiplicative doubling) fail, with paired doubling ratios to
  random sets of `1.00–1.01`; H3 (biased dilate set is a union of subgroup
  cosets) fails in all but 6 cases; H4 (fixed-`k` maximum bias decays with
  `p`) fails because for `(m,n)` up to `(6,6)` the exact maximum bias is
  `1` at every prime `101 ≤ p ≤ 1499`, and the counting identity
  `Σ_U |N(U)| = p·C((p−1)/2, k)` proves this for all large `p`.
- **Konyagin plus known Stepanov energies is trivial at `p = n⁴`.** An
  exact-rational ledger shows every pair `k,l ≤ 40` is trivial with
  `E_2 ≪ n^{49/20}`, `E_3 ≪ n⁴ log n`, `E_d ≪ n^{2d−2+2^{1−d}}`; the
  recurrence-only descent gives `E_2 ≤ (3/8)N⁴`, weaker than trivial.
- **Pass 7's prime class is not vacuous.** 14 of the 24,379 quartic primes
  at `N = 64` violate `E_3 ≤ (15+log N)N³`, e.g. `p = 7204033` with
  `E_3 = 20.93·64³`.
- **A tacit uniformity in pass 9.** The fixed-variable sum `I(y)` is not
  `O(p^{(k−1)/2})` uniformly in `y`: for the word
  `({0,1},{1,2},{0,2},{0,1},{1,2},{0,2})`, `k = 6`, at the surviving anchor
  the exact value is `p³(1+O(1/p))` at `p = 10009` and `100049`. Pass 9's
  separate treatment of anchor fibres is therefore necessary, not optional.

## 3. Referee verdict on passes 9 and 10

Scope: the concurrent session wrote passes 11–20 while this pass ran;
those are not covered by this review.

Neither is refuted: 1,417,148 exact checks (dense necklace traces to
`p = 3389`, one-anchor traces to `p = 2·10⁶`, fixed-variable inner sums and
Gram matrices to `p = 1,000,033`) found no witness against the inequality
`(T)`, the inner-sum bound, or the word-aggregate Gram bound. Both are
classified **plausible with unreviewed gaps**: the rank-growth lemma, the
pseudoreflection, the word recovery, the raw-rank recurrence and the
dimension obstruction are re-derived from verbatim-quoted Katz theorems
(Rigid Local Systems 2.9.7, 3.3.6–3.3.7); the remaining gaps are
citation-level (arithmetic Frobenius weights of the constant and punctual
kernels asserted via Weil II 1.8.4 without argument; the BBD scan has no
text layer). The referee signs the conditional statement: on Katz
2.9.4(3), 2.9.7, 3.3.3, 3.3.6–7, Deligne Weil II 1.8.4/3.3.1, BBD
5.1.8/5.1.14/5.3.1 and base change, inequality `(T)` and Kunisky's
Conjecture 1.14 hold at every fixed anchor count. Two omissions in pass 9
are recorded: the one-anchor case is Kunisky's Theorem 1.18 with the better
constant `k`, and the constant two-anchor word is his Theorem B.1. The
bound `(T)` is non-vacuous only for `p > (3a(2a+2)^{k−2})²`; the exponent
`(k+1)/2` is sharp for one anchor and its Möbius conjugates and a method
artefact for the other words.

## 4. What the data say (HEURISTIC, exact numbers)

- Biased dilate sets are tiny. Among 1,243 stored rectangles with bias
  `≥ 1/2`, the set `D_{1/2} = {t : |g(t)| ≥ ½|g(1)|}` is `{1}` in 95.4 %,
  and every case with `|D| > log₂ p` is a geometric progression, the one
  mechanism proved to produce `≥ ¼log₂ p − O(log log p)` biased dilates.
  **Conjecture SI** (crux note): for `A,B ⊆ F_p^*`, `|A||B| ≥ 64`,
  `|S| ≥ ½|A||B|`, `|D_{1/2}| ≤ 8 log₂ p`, and `≤ 2 log₂ p` unless `A` or
  `B` lies in a geometric progression of at most twice its length. Zero
  violations; it implies a weak subgroup case, so it needs new input.
- The Burgess–Chang energy deviation `Φ` is `1–26` on every tested family
  at `p ≤ 2003`; the expansion factor `λ` is `7–21` for random or coset
  `B` and `1.5–2.6` for intervals. Energy is never the bottleneck;
  shift-invariance is.
- Dyadic Gauss periods: over 51,780 scanned quartic primes,
  `E_3/N³ ≤ 20.93`; the ratio `M/√(n log(p/n))` at `N = 64` has mean
  `1.29` and maximum `1.70`, both below an iid-Gaussian extreme-value null
  model (`1.34`, `1.74`). No growing constant through `N = 64`.

## 5. prove2me

Four theorems compile in the platform's `c5ea003` environment; the OPEN
goal and three literature statements (Karatsuba, Hanson–Petridis
Corollary 1.5, two sentences) are faithful `sorry` statements. A complete
proposal package (`~/prove2me_workspace/proposals/paley-mission-proposal.json`,
description, README with the exact `curl` sequence) is prepared. A
read-only authenticated search of the platform found no existing
Paley, Legendre-symbol or quadratic-character theorems and no related
mission. **Nothing was uploaded**: submission publishes content under the
user's account and was left for the user. A `credentials.json` written by
the concurrent Proximity-Prize session exists and was used only for the
read-only environment and search calls.

## 6. Open obligations, in order of leverage

1. Any estimate for `S(A,B)` that depends on the representation function
   `r_{A+B}` rather than the dilation-invariant ratio-class sums. Every
   moment bound in this workspace is dilation-invariant and therefore
   cannot distinguish `t = 1`; this is the exact form of the obstruction.
2. The uniform shifted-subgroup bound `M_H ≤ p^{−δ}|H|` for `|H| = p^ε`
   (Bourgain's Problem 5). It is the subgroup case of the conjecture and
   the `K = 1` case of the multiplicatively structured case.
3. The single-prime hypothesis `E_3(μ_N) ≤ N^{3+o(1)}` at every quartic
   prime, which would make the subgroup exponent `8/9` unconditional; even
   `Ẽ_4 = O(N⁴)` caps the moment method at `7/8`, so the target
   `C√(n log(p/n))` needs the logarithmic-depth estimate (SG).
4. Closing the citation-level gaps in passes 9–10 from a text source of
   Weil II and BBD, and acknowledging Kunisky 1.18 and B.1.
5. Proving or refuting Conjecture SI.
6. Rewriting Chang's Freiman step with `K` tracked, to put the one-sided
   `2/5` exponent in print.

Every bound in this workspace that has been written down is field-agnostic
except the sum-product inputs; by the subfield test, the missing step must
enter exactly where a prime-field-only fact is used.

## 7. Second wave: three reformulations built on the dilation obstruction

After the first wave, three further workers attacked the obstruction
itself. Each produced a theorem, a refutation with witnesses, and a
verifier re-run by the orchestrator.

| direction | note | verifier checks | outcome |
|---|---|---|---|
| dual | [sigma-dual](sigma-dual-2026-09-05.md) | 310,324 exact | coset expansion exact; "cancellation across cosets" is a relabelling |
| tuple | [sigma-tuple](sigma-tuple-2026-09-05.md) | 375,169 exact | Paley ⟺ tuple-moment domination (t=0 excluded); Conjecture T(k) |
| si | [sigma-si](sigma-si-2026-09-05.md) | 47,163 exact | Conjecture SI REFUTED; corrected SI* with signal-to-noise hypothesis |
| frobenius (root) | [sigma-frobenius](sigma-frobenius-2026-09-05.md) | 1,258 exact | tuple sum = signed Frobenius traces; subgroup loci give Jacobi sums, the same open input as the dual note |

**Fourier dual (PROVED / REFUTED).** For a subgroup `H`, with Gauss sum
`G` and periods `Ĥ(ξ) = Σ_{h∈H} e(ξh/p)`, `Ĥ_χ(η) = Σ_{h∈H} χ(h)e(ηh/p)`,
the expansion `G·Σ_{h∈H}χ(h+c) = Σ_{cosets ξH} χ(ξ)Ĥ(ξ)Ĥ_χ(cξ)` holds
exactly with no boundary term (checked in `Z[ζ_p]` for all 353 subgroups
of all `p ≤ 200`). Norm-only information about the two period families
gives at best `√p(1−o(1))`, and any pointwise Bourgain–Glibichuk–Konyagin
saving `p^{−β}` on `max|Ĥ|` has `β ≤ ε/2`, which yields nothing for
`|H| ≤ √p`. The coset-cancellation ratio `ρ` satisfies an exact sandwich
around `|T_H(c)|/√|H|`, so square-root cancellation across cosets is the
original problem relabelled; `ρ = O(1)` is false (Weil plus Dirichlet
give `ρ ≥ √m` for order-`m` subgroups at infinitely many `p`; witnesses
`p = 2437, |H| = 29, ρ = 3.80` and `p = 2081, |H| = 13, ρ = 3.76`). The
genuine open input is a Jacobi-sum statement over a character subgroup of
order `p^{1−ε}`.

**Tuple level (PROVED / REFUTED / OPEN).** With `T_ns(A,B;k)` the sum of
the non-square Weil terms `Σ_{t∈F_p^*} χ(Π_i(ta_i+b_i))` over
`(A×B)^{2k}`: for `0<ε<1`, the conjecture holds for `ε` with some `δ₀`
iff there are `k` and `δ'` with `k(ε−2δ') > 1` such that
`T_ns ≤ p^{−2kδ'}(mn)^{2k}` for all large `p` and all `|A|,|B| > p^ε`
(Theorem 1.5; the dilate `t = 0` must be excluded, since `g(0) = mn`
when `B ⊆ Q`). The square tuples are automatically negligible. The
conjecture therefore asks for a saving of `(2k−1)p^{1/2+2kδ'}` over the
triangle bound, and this is not the spectral line's "cancellation
between words". The natural strong form `|T_ns| ≤ C_k√p·R_ν^k` is refuted
(`A = {1}`, `B = Q`, `k = 2`: ratio `≈ 3√p` at `p = 61, 101, 151`). On
663 instances the non-square sum for random rectangles sits at the
pattern-level random-sign scale (`|T_ns|/R_D` median `0.6–1.0`, no
growth in `p`), while for complete rectangles it equals
`(#biased dilates)·(mn)^{2k}` with no cancellation, as the identity
forces. Conjecture T(k): `|T_ns| ≤ C_k√p(mn)^{2k}/min(m,n)^k`, which
holds with `C₂ ≤ 7.8`, `C₃ ≤ 36` on every instance and would imply the
conjecture for `δ < ε/2 − 1/(2k)`. It is HEURISTIC and untested where it
should be tight.

**Conjecture SI (REFUTED).** Both clauses fail: `p = 97`, `H` the
subgroup of order 8, `A = 2H`, `B = H` has bias exactly `½` and
`|D_{1/2}| = 88 > 8 log₂ 97`; 966 of 1,009 subgroup rectangles
`(cH⁺,H)` with bias `≥ ½` at `p ≤ 4000` violate clause 1; random-like
`8×8` rectangles of bias `½` violate both clauses from `p = 941` on
(`p = 1009`: `|D| = 67 > 2 log₂ p`), and complete `8×8` rectangles at
`p = 4000037` have `|D| = 327`. Exact formulas show `|D_{1/2}|` is linear
in `p` at fixed sizes. PROVED: for every subgroup `H`,
`g_{A,H}(th) = χ(h)g_{A,H}(t)`, so `D_{1/2}(A,H)` is a union of
`H`-cosets containing `H`, a biased `A` exists iff `M_H ≥ |H|/2`, and
clause 1 restricted to `B = H`, `|H| > 8 log₂ p` is equivalent to
`M_H < |H|/2` (true for all 3,890 subgroups at `11 ≤ p ≤ 4000`, open in
general). The corrected Conjecture SI*: if bias `≥ ½` and
`S² ≥ 8·R·ln p` (signal-to-noise), then `|D_{1/2}| ≤ 2 log₂ p`; zero
violations on about 371,000 hypothesis-satisfying rectangles, observed
maximum `1.07 log₂ p`, all near-extremal cases of geometric-progression
type. SI* is OPEN and does not bear on `t = 1`.

**Net effect on the obligations in §6.** Item 1 (a non-invariant
estimate) is now sharpened: the only non-invariant data in the tuple
decomposition are the signed ratio-class sums `σ_r`, entering through
`c(D)`. Item 5 is replaced by SI*. Item 2 gains an exact dual form
(Jacobi-sum cancellation over a character subgroup).

## 8. Third wave: stress tests and server verification

| direction | note | verifier | outcome |
|---|---|---|---|
| stress | [sigma-stress](sigma-stress-2026-09-05.md) | 4,626 exact | Conjecture T(k) REFUTED for every constant; Conjecture SI* REFUTED; corrected T′(k) survives |
| p2m | [sigma-p2m](sigma-p2m-2026-09-05.md) | 6,241 exact | four proofs ACCEPTED by prove2me's server; everything private |

**Conjecture T(k) is false for every `k ≥ 2` and every `C_k`.** With
`A = B = Q` (the nonzero squares), `h = (p−1)/2`, the non-square tuple sum
has the closed form `T_ns = h^{2k}·h[1 − N_k(h) − N_k(h−1)]`, so the ratio
`|T_ns|·min(m,n)^k/(√p(mn)^{2k})` equals `2(h−1)(3h−2)/(h√p)`, which is
`2.9968√p` at `p = 4001`, `k = 2`, and tends to `(2k−1)!!·√p`. Subgroups
`A = B = H ⊆ Q` do worse in absolute terms: ratio `152.3` at `p = 2833`,
`|H| = 59`, `k = 2` (and `2753` at `k = 3`), `785.9` at `p = 13183`,
`|H| = 169`, against the observed constants `C₂ ≤ 7.8`, `C₃ ≤ 36` of the
tuple note; the median ratio over subgroups grows like `|H|^{0.65}`, flat
in `p` at fixed `|H|`. So the "pattern random-sign" scale is the wrong
null for structured pairs, exactly as the Frobenius note predicts (their
Weil terms are correlated Jacobi sums): the profile is constant on
`H`-cosets, so the moment is `|H|` times a short sum and coherence costs a
factor `√min(m,n)`. Theorem 1.5's equivalence stands; only this candidate
sufficient condition dies. The corrected statement
T′(k): `|T_ns| ≤ C_{k,η} p^{1/2+η}(mn)^{2k}/min(m,n)^{k−1/2}` has the sharp
exponent `k−½` (attained by both refuting families), survives every family
tested (normalized ratio rms `5.7` at `k=2`, one spike of `60.5` at
`p = 13183` where `M_H = 65 = 5.0√|H|`), and would still imply the
conjecture for `δ < ε/2 − 1/(2k)`. It is OPEN and, for subgroup pairs, is
a Gaussian-moment statement for `T_H` across cosets, i.e. again the
shifted-subgroup problem.

**Conjecture SI* is false.** Designed `1×n` rectangles `A = {1}`, `B` = the
`n` best shifts for a weighted random target satisfy the signal-to-noise
hypothesis with margin (`z² = 10 ln p`) and have `|D_{1/2}| = 32 > 23.9`
at `p = 4001`, `47 = 3.54 log₂ p` at `p = 10007`, `70 = 4.21 log₂ p` at
`p = 100003`, `124 = 6.22 log₂ p` at `p = 1000003`: the ratio to `log₂ p`
grows. A capacity heuristic suggests only constants near `22 log₂ p`
(relative threshold) or `5.5 log₂ p` (absolute threshold) could survive;
no logarithmic statement with a small constant survives this pass.

**prove2me (server-verified, private).** In the platform's environment
Mathlib `c5ea003` / Lean v4.30.0, the four proofs (shift orthogonality,
second moment, Chung bound, interval-to-non-residue) were submitted and
returned `ACCEPTED` on first submission (submission ids in the p2m note);
the eight statements exist as private theorems, the goal and the three
literature statements as `Open`; a private draft mission proposal with
seven milestones is ready for the user to audit and launch. Nothing was
made public, launched, commented, voted or rated.

## 9. Where this leaves the programme

Every candidate sufficient condition formulated in this pass has now
been either shown equivalent to the conjecture (tuple-moment domination,
level-set form), shown equivalent to an older open problem (the subgroup
case, coset cancellation, the Jacobi-sum input), or refuted with exact
witnesses (SI, SI*, T(k), the shift-chain rescue, four structural
hypotheses). The obstruction is stated exactly in three equivalent
languages: dilation invariance of every termwise-Weil bound, the need for
a saving of `(2k−1)p^{1/2+2kδ'}` over the triangle bound on the
non-square tuple sum, and the need for cancellation in
`Σ_{ψ∈Ψ} ψ(u)J(ψ,χ)` over character subgroups of order `p^{1−ε}`. A proof
needs an input that is none of Weil, orthogonality, moments, spectral
identities or the polynomial method, and that fails over `F_{p²}`.

## 10. Fourth wave: directions the workspace had not touched

| direction | note | verifier checks | outcome |
|---|---|---|---|
| lit2026 | [sigma-lit2026](sigma-lit2026-2026-09-05.md) | 14.3 million exact | no claimed proof anywhere; Sárközy's conjecture now proved (Kalmynin 2025); nothing below √p |
| extractor | [sigma-extractor](sigma-extractor-2026-09-05.md) | 168,410 exact | exact extractor equivalence; condenser composition provably cannot reach Paley |
| biclique | [sigma-biclique](sigma-biclique-2026-09-05.md) | see results file | max complete-biclique product is (p+3)/2 for 43≤p≤1000, attained at |A|=2; balanced bicliques stay near ω; sum-product route capped at (2/9)p |

**Literature 2024–2026 (CITED).** Thirty-five arXiv queries, a citation
sweep of the core papers and web searches found no source claiming the
two-set conjecture, a clique bound `p^{1/2−c}`, or a sub-`√p` biclique
bound. The best exponent for two arbitrary sets is still Karatsuba's
`|A| > p^{1/2+η}` (Fouvry–Shparlinski–Xi, arXiv:2404.09295, §1.2:
"the classical inequalities of Vinogradov and Karatsuba still stand");
the best clique bound is Hanson–Petridis `ω(G_p) ≤ (√(2p−1)+1)/2`.
Genuinely new: Kalmynin (arXiv:2504.10202) proves Sárközy's conjecture,
that the quadratic residues are never `A+B` with `|A|,|B|>1`, and
`μ_d = A+B` forces `|A|=|B|=√d`; Rudnev–Tyrrell (arXiv:2607.24270) and
Yip–Yoo (arXiv:2608.02568, 2607.25711) classify additive decompositions
of subgroups. All are Stepanov-method complete-biclique statements whose
prime-field input (binomials nonvanishing mod `p`, degree below `p`) is
exactly what caps them at `√(p/2)`. Wang–Shen–Kobzar (2024) show the
Lovász-type localization hierarchy has value `≥ √p/2^{t−1}` at every fixed
level, another excluded proof shape. Satake (arXiv:2405.08608, Theorem 18)
derives the extractor statement from the Paley-ETF RIP, which is itself
conditional. A small elementary theorem was added: if `Φ_p` is
`(K,δ)`-RIP and `ω(G_p) ≤ K` then `ω(G_p) ≤ δ√p+1`.

**Extractor formulation (PROVED / CITED).** The two-set conjecture is
equivalent, up to `ε' < ε`, to: `½(1+χ(x+y))` is an `(ε log₂ p, p^{−δ})`
two-source extractor, for flat sources, for all sources of that
min-entropy (constructive flat decomposition, checked in exact rationals),
and in the strong sense. No published work analyses Paley itself below
rate `1/2`: Bourgain 2005 and Lewko 2019 use Hadamard on `(x,x²)` and
paraboloids, Barak–Impagliazzo–Wigderson condense `ab+c` on three or more
sources, and Chattopadhyay–Zuckerman, Cohen, Li and Ben-Aroya–Doron–Ta-Shma
are non-algebraic; every working sub-1/2 argument fails the `F_{p²}` test
at its incidence step. Exact obstruction to composition: no single-source
condenser exists; the Fourier/energy route satisfies
`S(A,B)² ≤ p·E⁺(A,B)` with `E⁺ ≥ |A||B|`, so it is capped at Chung; and the
Cauchy–Schwarz squaring that drives Bourgain's analysis merely relabels
Paley's sources bijectively (`χ(x₁+y)χ(x₂+y) = χ(d)χ(d^{−1}+(x₁+y)^{−1})`,
6,300 exact checks) whereas for Hadamard it produces `X−X`. Refuted: the
orbit-average hypothesis `|S(A,B)| ≤ 2·mean_c|S(A,cB)|` (witness
`p = 8009`, `|A|=|B|=57`, bias `0.376` against orbit mean `0.0139`).

**Complete bicliques (PROVED / REFUTED, see the note for the final table).**
For `A+B ⊆ Q∪{0}` the Hanson–Petridis bound refines to
`|A||B| ≤ (p−1)/2 + min(|A|,|B|,ω(G_p))` because `A ∩ (−B)` is a clique, so
the balanced biclique number and the sum-clique number obey the same
`(1+√(2p−1))/2` bound as `ω(G_p)`. With `|A| = 2` one gets exactly
`|B| = (p+3)/4`, so the maximum product `P(p) ≥ (p+3)/2`, and an exact
finite computation gives `P(p) = (p+3)/2` for every prime `43 ≤ p ≤ 1000`,
`p ≡ 1 (mod 4)` (the small primes `13, 37, 41` are clique-type). So the
unbalanced complete case is completely understood up to `p = 1000` and sits
at `p/2`, while the balanced case is governed by the clique number. Two
no-go lemmas show the sum-product route on the complete case (ratio set of
the sumset inside `Q`, weighted multiplicative energy of the sumset) cannot
push below `|A||B| ≤ (2/9+o(1))p`. Refuted with witnesses: `b(p) = ω(G_p)`
(`p = 89`) and `b(p) ≤ ω(G_p)+1` (`p = 233`). The worker's final table and
verifier count were still being assembled under heavy machine load when
this summary was written; the note carries the final numbers.

