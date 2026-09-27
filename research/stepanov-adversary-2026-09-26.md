# Stepanov wave, worker `adversary`: is robust Hanson–Petridis true, and what does f look like?

**Status: PARTIAL — the worker was stopped by an account usage limit. Sections 1–4 are final (PROVED for the stated ranges). Sections 5–8 were never written; their data are in `results/stepanov_adversary_2026_09_26_search.json` and summarised in the root addendum at the end. The verifier re-checks every stored witness (61,052 exact checks, 0 failures after a root fix on 2026-09-27: part H had applied the Legendre second-moment identity to the random ±1 null model, where it does not hold).**

Worker slug `adversary`. Verifier `experiments/stepanov_adversary_2026_09_26.py`
(writes `results/stepanov_adversary_2026_09_26.json`), C++ search helper
`experiments/stepanov_adversary_2026_09_26.cpp`, driver
`experiments/stepanov_adversary_2026_09_26_driver.py`, witness data in
`results/stepanov_adversary_2026_09_26_search.json`.

## 0. Definitions

`p` odd prime, `χ` Legendre, `χ(0)=0`, `d=(p−1)/2`. For `A ⊆ F_p`, `m=|A|`,
`e_b = #{a∈A : χ(a+b) = −1}`, `B_e(A) = {b : e_b ≤ e}`,
`r_e(A) = |B_e(A) ∩ (−A)|`, `F_A(b) = Σ_{a∈A} χ(a+b)`,

    R_e(A) = (m·|B_e(A)| − r_e(A)) / d,        η = e/m.

HP (Hanson–Petridis, arXiv:1905.09134 Thm 1.2 with `d=(p−1)/2`, local text
`sources/sigma-hanson-petridis-1905.09134.txt` l. 120–123, CITED): "Let p be a
prime and suppose A, B ⊆ F_p satisfy A + B ⊆ Z_d ∪ {0} for some d properly
dividing p − 1. Then |A||B| ≤ d + |B ∩ (−A)|." With `B = B_0(A)` this is
exactly `R_0(A) ≤ 1` for every `A`. RHP(f,g) is the statement of
`research/stepanov-brief-2026-09-26.md`; taking `B = B_e(A)` it says
`R_e(A) ≤ f(e/m) + g(m)/d` for `|A| = m ≤ √p`.

Symmetries: `A ↦ qA + t` with `q ∈ Q` (and `B ↦ qB − t`) preserves every
`χ(a+b)`, hence `R_e`; a non-residue multiplier swaps `+1` and `−1` and is
not a symmetry (for `p ≡ 3 mod 4` this includes `x ↦ −x`). Every `A` with
`|A| ≥ 2` is equivalent to one containing `{0,1}` or `{0,ν}` (`ν` a fixed
non-residue); the exhaustive search uses this normalisation.

Notation: `C(m,≤e) = Σ_{j≤e} C(m,j)`, `g(m,e) := 2m·C(m,≤e)/2^m`, and the
**second-moment share** `s_e(A) := Σ_{b∈B_e(A)} F_A(b)² / (m(p−m))`; since
`Σ_b F_A(b)² = m(p−m)` exactly, `0 ≤ s_e ≤ 1`.

## 1. The Weil range: RHP holds there, and it forces the shape of f (PROVED, using Weil)

**Proposition 1.1 (fixed-m expansion).** For every `A ⊆ F_p` with `|A| = m`
and every `0 ≤ e ≤ m`, with `C = C(m,≤e)` and

    W_m(p) = (m·2^{m−1} − m) + C(m,2) + √p·((m−2)·2^{m−1} + 1 − C(m,2)),

one has `| |B_e(A)| − p·C/2^m | ≤ m + C·W_m(p)/2^m`, and consequently

    | R_e(A) − g(m,e) | ≤ g(m,e)/(p−1) + 2m(m + 1 + C·W_m(p)/2^m)/(p−1).

*Proof.* For `b ∉ −A` every `χ(a+b) = ±1`, so
`[χ(a+b) = −1] = (1−χ(a+b))/2`, `[χ(a+b) = 1] = (1+χ(a+b))/2`, and
`#{b ∉ −A : e_b ≤ e} = Σ_{b∉−A} Σ_{J⊆A,|J|≤e} Π_{a∈J}(1−χ(a+b))/2 · Π_{a∉J}(1+χ(a+b))/2
= 2^{−m} Σ_{|J|≤e} Σ_{I⊆A} (−1)^{|I∩J|} T(I)`, where
`T(I) = Σ_{b∉−A} χ(Π_{a∈I}(a+b))`. `T(∅) = p − m`. For `I ≠ ∅` write
`T(I) = Σ_{b∈F_p} χ(P_I(b)) − Σ_{b∈−A} χ(P_I(b))`; the second sum has at most
`m − |I|` nonzero terms. The complete sum is `0` for `|I| = 1`, `−1` for
`|I| = 2` (`Σ_b χ(b)χ(b+c) = −1`, `c ≠ 0`), and at most `(|I|−1)√p` in absolute
value for `|I| ≥ 3` by Weil's bound for the squarefree polynomial `P_I`
(`research/sigma-brief-2026-09-05.md`, background facts; CITED). Summing,
`Σ_{I≠∅} |T(I)| ≤ W_m(p)`, and `|Σ_{|J|≤e}(−1)^{|I∩J|}| ≤ C`. So
`#{b∉−A : e_b ≤ e} = C(p−m)/2^m + err`, `|err| ≤ C·W_m(p)/2^m`. Finally
`|B_e| = #{…} + r_e` and `|r_e − mC/2^m| ≤ m` (both lie in `[0,m]`). For the
second statement write `δ = |B_e| − pC/2^m`; then
`R_e = (m|B_e| − r_e)/d = g·p/(p−1) + 2(mδ − r_e)/(p−1)`. ∎

The verifier checks the exact character-sum identity (420 cases, `p ≤ 43`,
`m ≤ 5`, all `e`), the first bound on 228 random sets (`p ∈ {1009, 10007, 100003}`,
`m ≤ 8`; largest deviation/bound ratio `0.79`), and the second on 36 more.

**Corollary 1.2 (lower bounds on any admissible f).** If RHP(f,g) holds, then
`f(e/m) ≥ g(m,e)` for every `0 ≤ e < m/2`. In particular
`f(1/6) ≥ 21/16`, `f(1/5) ≥ 15/8`, `f(1/4) ≥ 5/2`, `f(2/7) ≥ 203/64`,
`f(1/3) ≥ 1191/256 ≈ 4.652` (at `(m,e) = (12,4)`), `f(2/5) ≥ 10.87`
(at `(35,14)`), and `g(m,(m−1)/2) = m` for odd `m`, so no admissible `f` is
bounded on `[0,1/2)`.

*Proof.* Fix `A_0 ⊂ Z` with `m` elements; for primes `p` larger than
`max(m², diam A_0)` its reduction `A` has `m ≤ √p` elements. RHP with
`B = B_e(A)` gives `R_e(A) ≤ f(e/m) + g(m)/d`, Proposition 1.1 gives
`R_e(A) ≥ g(m,e) − o(1)`, and `g(m) = o(p)`. ∎

**Lemma 1.3.** If `1 ≤ e < m/6` then `g(m,e) ≤ 7/8`, with equality only at
`(m,e) = (7,1)`. *Proof.* Exact rational check for `7 ≤ m ≤ 400` (verifier
A1). For `m > 400`: `C(m,≤e) ≤ 2^{mH(e/m)} ≤ 2^{mH(1/6)}` (standard entropy
bound, `e/m < 1/6 < 1/2`), `H(1/6) = 0.65002… < 0.6501`, and
`2m·2^{−0.3499m} < 10^{−38}` for `m ≥ 401`. ∎

**Theorem 1.4 (RHP in the Weil range).** Let `c < c_0 := 1/(2H(1/6)) = 0.7692…`.
There is `p_0(c)` such that for all primes `p ≥ p_0`, all `A ⊆ F_p` with
`|A| = m ≤ c·log₂ p` and all integers `0 ≤ e < m/6`: `R_e(A) ≤ 1`, and in
fact `R_e(A) ≤ 7/8 + o(1)` when `e ≥ 1`. Thus RHP(f, 0) holds for these `m`
with `f ≡ 1` on `[0, 1/6)`.

*Proof.* `e = 0` is HP. For `e ≥ 1` we have `m ≥ 7` and `g(m,e) ≤ 7/8`
(Lemma 1.3). In Proposition 1.1, `2^{−m}W_m(p) ≤ m(1+√p)` and
`C(m,≤e) ≤ 2^{mH(1/6)} ≤ p^{cH(1/6)} = p^{1/2−δ}` with
`δ = (c_0 − c)H(1/6) > 0`, so the error is
`O(m²p^{−δ} + m²/p) = O((log p)² p^{−δ}) → 0`. ∎

*Remarks.* (i) The range is where the Paley problem itself is easy, so
Theorem 1.4 is a consistency check, not progress; the explicit error bound
only bites for large `p` (verifier A6: the largest `M` such that
`g(m,e) + err ≤ 1` for all `7 ≤ m ≤ M`, `1 ≤ e < m/6`, is `M = 6` (vacuous)
at `p = 10^6`, `15` at `10^9`, `22` at `10^{12}`, `41` at `10^{20}`, `240`
at `10^{100}`; `M/log₂p` rises from `0.50` to `0.72`, towards `c_0`). (ii) Together with Corollary 1.2 it
pins down the smallest `f` that the Weil range allows: the step function

    F(η) := sup_{e/m ≤ η} g(m,e) = 1 on [0,1/6),  21/16 at 1/6,  15/8 on [1/5,1/4),  5/2 at 1/4, …

(`F(2/7) = 203/64`, `F(3/10) = 55/16`, `F(1/3) = 1191/256`, `F(2/5) ≈ 10.87`,
`F(0.45) ≈ 38.06`; verifier A3). `F` jumps at `η = 1/6` because
`g(6,1) = 21/16 > 1 > 7/8 = g(7,1)`. **So "f(η) → 1 as η → 0" is consistent
with the Weil range, and in that range the limit is attained exactly:
f = 1 for all η < 1/6.** The content of RHP is entirely in
`log p ≪ m ≤ √p`.

## 2. The second-moment share: what an RHP counterexample must look like (PROVED)

**Proposition 2.1.** For every `A` with `|A| = m < p` and every `e` with `m > 2e+1`:

(i) `R_e(A) ≤ 2·(p−m)/(p−1) · s_e(A) / (1 − 2η − 1/m)²`;

(ii) (HP ⇒) `s_0(A) ≤ (d − r_0 + r_0/m)/(p−m) ≤ d/(p−m)`;

(iii) for `A = {0,1}` (`p ≡ 1 mod 4`), `s_0 = (p−3)/(2(p−2))`, so (ii) is sharp up to `O(1/p)`.

*Proof.* (i) For `b ∈ B_e ∖ (−A)`, `F_A(b) = m − 2e_b ≥ m − 2e`; for `b ∈ B_e ∩ (−A)`
exactly one sum is `0` and `F_A(b) = m − 1 − 2e_b ≥ m − 2e − 1 > 0`. Hence
`|B_e|(m−2e−1)² ≤ s_e·m(p−m)`, and `R_e ≤ m|B_e|/d`. (ii) On `B_0`,
`F = m` off `−A` and `m−1` on `−A`, so
`s_0·m(p−m) = (|B_0|−r)m² + r(m−1)² = m(m|B_0|) − 2rm + r ≤ m(d+r) − 2rm + r`
by HP. (iii) `|B_0({0,1})| = (p+3)/4` with `r = 2`
(`research/sigma-biclique-2026-09-05.md` Prop. 2.3), `F = 2` on `(p−5)/4`
points and `1` on two. ∎

So HP says: *complete positive rows carry at most half of the second
moment.* Proposition 2.1(i) says an RHP counterexample (`R_e ≥ 1 + c`,
`η → 0`, `m → ∞`) must have `s_e(A) ≥ (1+c)/2 − o(1)`: **its near-complete
positive rows must carry more than half of `Σ_b F_A(b)²`**. The second-moment
ceiling `R_e ≤ 2/(1−2η)² + O(1/m)` of the brief is the case `s_e = 1`. Over
`F_{p²}` the subfield has `s_0 = (p−1)/p → 1` (section 3). This is the
precise sense of "approximate subfield" used below; the verifier checks (i)
on every stored witness and (ii) on every exhaustive `e = 0` witness.

## 3. The F_{p²} contrast (PROVED, exact)

`F_{p²} = F_p[x]/(x² − ν)` with `ν` the least non-residue mod `p`, `q = p²`,
`d_q = (q−1)/2`, `χ_q(z) = z^{(q−1)/2}` (computed by exponentiation in the
verifier and checked equal to `χ_p(N(z))`, `N(u+vx) = u² − νv²`, for all `z`,
`p ≤ 31`).

**Proposition 3.1.** (i) `A = B = F_p` satisfies `A + B ⊆ Q_q ∪ {0}` with
`|A||B| − r = p² − p > d_q`, so HP with `r` retained fails for every `p ≥ 3`.
(ii) `F_{F_p}(b) = p − 1` for `b ∈ F_p` and `−1` for `b ∉ F_p`; hence
`R_0(F_p) = 2p/(p+1)`, `s_0(F_p) = (p−1)/p`: the subfield puts all but
`1/p` of the second moment on complete rows and attains the ceiling at `η = 0`.
(iii) For any `X` with `|X| = k ≤ e`, `B_e(F_p ∪ X) ⊇ F_p`, so
`R_e(F_p ∪ X) ≥ ((p+k)p − p)/d_q ≥ 2p/(p+1)`, with `η = e/(p+k) → 0`.
Hence over `F_{p²}` every admissible `f` has `f(η) ≥ 2 − o(1)` for all
`η ∈ [0, 1/2)`: **RHP fails over `F_{p²}` at every `η`, including `η = 0`.**

*Proof.* Every element of `F_p^*` is a square in `F_{p²}` (`F_p^* ⊆ F_{q}^{*2}`
since `(q−1)/(p−1) = p+1` is even). For `b = u + vx`, `v ≠ 0`:
`Σ_{a∈F_p} χ_q(a+b) = Σ_{a} χ_p((a+u)² − νv²) = −1` (the leading coefficient
`1` is a square and the discriminant `4νv² ≠ 0`). The rest is counting. ∎

**Exhaustive and search data over `F_q`.** Exact maxima (C++ `exhq`, all `A`
normalised to contain `{0,1}` or `{0,ν_q}`): `max R_0 = 3/2, 5/3, 7/4` at
`m = √q` for `q = 9, 25, 49`, attained by `F_p` — above 1, while over `F_p`
the exhaustive `max R_0 ≤ 1` for every `p ≤ 61` and every tabulated `m`.
Positive control for the search: started from random sets, the same
annealer that is used over `F_p` below finds `R_0 = 2p/(p+1)` exactly at
`m = p` for every `p ∈ {5,7,11,13,17,19}` (the witnesses are affine
`F_p`-lines; verifier part C), and `R_1, R_2 ≈ 1.87–1.92` at `m = p`.
Near-subfield constructions (`p ≤ 31`): `F_p` plus `k ≤ 3` random points has
`R_k = 2.00–2.13` at `η = k/(p+k)`; `F_p` with `k` points swapped out has
`R_k ≈ 1.87–1.94`, share `0.79–0.91` at `p = 31`; a random `p`-subset of
`F_{p²}` has `R_e ≈ 0`.

**Proposition 3.2 (the diagonal: near-cliques).** If every `a ∈ A` has at
most `e` elements `a' ∈ A` with `χ(a' − a) = −1`, then `−A ⊆ B_e(A)` and
`R_e(A) ≥ (m² − m)/d`. Hence RHP(f,g) implies the *robust clique bound*
`m² − m ≤ f(η)d + g(m)`, i.e. `m ≤ (√(f(η)/2) + o(1))√p`, whereas Chung's
bound gives only `m(1 − 2η) ≤ √p + 1` (from
`m(m − 1 − 2e) ≤ Σ_{a,a'} χ(a'−a) ≤ m√p`). *Proof.* `b = −a` has
`e_b = #{a' : χ(a' − a) = −1} ≤ e`, and `r_e = m`. ∎
In `F_{p²}` the subfield is an exact clique with `m = √q`, `R_0 → 2`; in
`F_p` the extremal `(m,e) = (7,1)` set at `p = 37` (below) is exactly such a
near-clique (defect `≤ 1` per vertex, `m = 7 > √37`).

## 4. Exhaustive data over F_p, p ≤ 61 (PROVED for the stated ranges; C++ + independent Python)

C++ mode `exh` enumerated every `A` of size `m ≤ m_max(p)` containing `{0,1}`
or `{0,ν}` (773,809,516 sets in total) and recorded, for every `(m, e)`, the
exact maximum of `m|B_e| − r_e`, a lexicographically first witness and the
number of normalised sets attaining it. `m_max = p` for `p ≤ 23`; `15, 13,
10, 9, 9, 9, 8, 8, 8` for `p = 29, 31, 37, 41, 43, 47, 53, 59, 61` (so
`m_max ≥ ⌊√p⌋ + 1` everywhere). The verifier recomputes all 1,521 witnesses,
and re-derives the maxima independently in numpy for every `p ≤ 47` and
`m ≤ ⌊√p⌋` (`m ≤ ⌊√p⌋ + 3` for `p ≤ 23`): 356 maxima compared, all equal.

Findings (exact):

1. **HP is tight only at `m ≤ 2` and at small primes.** `max R_0 ≤ 1` in
   every entry (as it must be). For `3 ≤ m ≤ √p` equality occurs only at
   `p ∈ {11,13,17,19,23,29,31,37,41}` (`(p,m) = (11,3), (13,3), (17,3),
   (17,4), (19,3), (23,3), (23,4), (29,3), (29,5), (31,4), (37,3), (41,5)`),
   never for `p ≥ 43`; `m = 1, 2` give equality at every `p`
   (`B_0({0,1})` has `(p+3)/4` elements). For `m > √p` equality also occurs
   at `m = d, d+1` (`A = Q − b`, `|B| = 1`) and a few other small cases.
2. **Robust ratios at small p exceed 1 at η < 1/6.** Row `p = 61`:
   `max R_1 = 2.07` at `m = 7` (`η = 1/7`) and `1.77` at `m = 8`;
   over all `p ≤ 61`, `max R_1` is `2.33` (`m = 7`, `p = 37`), `2.00`
   (`m = 8`, `p = 29`), `1.70` (`m = 9`), `1.50` (`m = 10`), `1.11`
   (`m = 11`), `1.09` (`m = 12, 13`), `1.00` (`m = 14`, `p = 29`). The
   witnesses have tiny `B` (`|B_1| = 4–9` against `d = 14–30`), e.g.
   `p = 61`, `A = {0,1,2,6,20,43,57}`, `|B_1| = 9`, `r = 1`, `R_1 = 31/15`.
   These are small-number effects: the Weil mean is `g(7,1) = 7/8` and the
   Weil error at `p = 61` exceeds the main term.
3. `max R_e` at fixed `(m,e)` is almost flat in `p` over `29 ≤ p ≤ 61`
   (e.g. `(4,1)`: `3.10–3.39`, limit `2.5`; `(6,1)`: `2.24–2.44`, limit
   `1.31`; `(8,2)`: `4.04–4.29`, limit `2.31`): at `p ≤ 61` nothing is in
   the asymptotic regime.

## Root addendum (2026-09-27; sections 5–8 were not written by the worker)

HEURISTIC, from the verifier's recomputation of the stored search witnesses:

- **Empirical f.** Best `R_e` found by annealing with structured seeds, `e = ⌊ηm⌋`:
  at `η = 1/16`, `0.70` (`p = 503`, `m = 16`) and `0.54` (`p = 1009`); at `η = 1/8`,
  `1.18` (`p = 503`, `m = 8`, i.e. `e = 1`) and `1.00` (`p = 1009`); at `η = 1/4`,
  `3.04` and `2.88`, all attained at `m = 8`, where `e ≥ 1` forces `η ≥ 1/8`. Values
  above 1 occur only at small `m`, where the Weil error dominates. These are far
  below the proved `f(η) = (1−√(2η))^{−2}` of `research/stepanov-robust-2026-09-26.md`
  (`f(1/16) ≈ 2.39`), consistent with that theorem and with room to spare.
- **Rectangles.** The most biased rectangles found with `|A||B| ≥ p/2` at
  `p = 1009` have bias `0.72–0.78` (e.g. `16 × 32`, `|A||B| = 0.507p`, bias `0.775`).
  All 165 stored `F_p` rectangles satisfy Theorem 2.5 of the robust note (root check,
  `results/stepanov_star_exhaustive_2026_09_27.json`).
- Obligations left by the worker: the full per-`(p, m, η)` tables, the null-model
  comparison write-up, and the `F_{p²}` near-subfield constructions beyond §3.
