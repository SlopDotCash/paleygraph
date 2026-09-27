# The Hanson–Petridis pencil and the ramification budget

**Status: PROVED (elementary, with exact checks) — the Hanson–Petridis
auxiliary polynomials of all subsets `A∖E` form one explicit linear system
(a pencil when `|E| = 1`); points with at most one bad partner are
ramification points of the pencil's rational map; and a single inequality
bounds complete and one-bad points together:
`(2M−2)·N₀ + (M−2)·N₁⁼ ≤ 2d−2`. This reproves, and slightly sharpens, the
constant 2 of the second moment and of the `algebra` worker's Hankel bound.
It does NOT prove robust Hanson–Petridis (RHP). It reduces RHP at `e = 1`
exactly to a statement about where the ramification of one rational map of
degree `d` lies. OPEN — that statement.**

Root note of the Stepanov wave (brief: `stepanov-brief-2026-09-26.md`).
Verifier: `experiments/stepanov_pencil_2026_09_26.py` (standard library,
18,191 exact checks at 13 primes `29 ≤ p ≤ 149`, 0 failures, 8 s) →
`results/stepanov_pencil_2026_09_26.json`.

## Setting

`p` odd prime, `d = (p−1)/2`, `A = {a_1,…,a_M}` distinct, `3 ≤ M ≤ (p+1)/2`,
`D = d + M − 1 ≤ p − 1`, `c_k = 1/Π_{l≠k}(a_k − a_l)`. The HP polynomial is

    F_A(x) = −1 + Σ_k c_k (x + a_k)^D,        deg F_A = d,

normalised so that `G(x) = Σ_k c_k (x+a_k)^{M−1} ≡ 1` (Lagrange:
`Σ_k c_k a_k^s = [s = M−1]` for `0 ≤ s ≤ M−1`). For `b ∈ F_p ∖ (−A)` let
`E(b) = {k : χ(b+a_k) = −1}`, `y_k = b + a_k`, `T_t(b) = Σ_{k∈E(b)} c_k y_k^t`,
and `N₀`, `N₁⁼` the numbers of such `b` with `|E(b)| = 0`, `= 1`.

## 1. Rational values of `F_A` and its derivatives (PROVED)

For `b ∈ F_p ∖ (−A)` and `1 ≤ j ≤ M−1`,

    F_A(b) = −2 T_{M−1}(b),     F_A^{(j)}(b) = −2 (D)_j T_{M−1−j}(b),

because `(b+a_k)^{D−j} = χ(b+a_k)(b+a_k)^{M−1−j}` and
`Σ_k c_k (x+a_k)^t = [t = M−1]` identically for `t ≤ M−1`. In particular
complete points (`E(b) = ∅`) are roots of order `≥ M`: this is HP's proof.

## 2. All subset polynomials lie in one linear system (PROVED)

**Lemma 2.1 (pencil).** For every `k`, `F_{A∖{a_k}} = F_A − (x + a_k)F_A′/D`
exactly (with HP's normalisation for `A∖{a_k}`).

*Proof.* The coefficients for `A∖{a_k}` are `c_l′ = c_l (a_l − a_k)`, and
`D′ = D − 1`. So `F_{A∖{a_k}} = −1 + Σ_l c_l[(x+a_l) − (x+a_k)](x+a_l)^{D−1}`,
the `l = k` term being zero, which is `F_A − (x+a_k)F_A′/D`. The
normalisation `G′ = G − (x+a_k)G′/(M−1) = 1` is inherited. ∎

**Lemma 2.2 (linear system).** For every `E ⊆ A` with `|E| = e ≤ M−2`,

    F_{A∖E} = Σ_{i=0}^{e} (−1)^i σ_i({x + a_k : k ∈ E}) · F_A^{(i)} / (D)_i,

where `σ_i` is the `i`-th elementary symmetric polynomial. Expanding
`σ_i(x + a_E)` in powers of `x`, every `F_{A∖E}` is
`Σ_j σ_j(a_E)·Y_j` for fixed polynomials `Y_0,…,Y_e` of degree `≤ d`
depending only on `A`: the HP polynomials of all `(M−e)`-subsets lie in one
`(e+1)`-dimensional linear system, with coordinates the coefficients of the
split polynomial `Π_{k∈E}(T + a_k)`. *Proof.* Induction on `e` using
Lemma 2.1; checked exactly for `e = 1, 2, 3`. ∎

Write `W = F_A′/D` and `V = F_A − xW`, so `F_{A∖{a_k}} = V − a_k W`.

## 3. The ramification inequality (PROVED)

Let `N = V′W − VW′`. Then

    N = ((D−1)·F_A′² − D·F_A·F_A″)/D²,        deg N = 2d − 2,

the leading coefficient being `d(M−1)λ²/D²` for the leading coefficient
`λ` of `F_A`, nonzero because `d, M−1 < p`. At rational points,
by §1, `N(b) = 4(D−1)·(T_{M−2}² − T_{M−1}T_{M−3})(b)`.

* If `|E(b)| = 0`, `F_A` vanishes to order `≥ M`, so `N` vanishes to order
  `≥ 2M − 2`.
* If `E(b) = {k}`, then `V − a_kW = F_{A∖{a_k}}` vanishes to order `≥ M−1`
  (HP for `A∖{a_k}`), and `N = (V − a_kW)′W − (V − a_kW)W′` vanishes to
  order `≥ M − 2`. Also `W(b) ≠ 0`, since otherwise `F_A(b) = 0` against
  `F_A(b) = −2c_k y_k^{M−1} ≠ 0`.
* If `|E(b)| = 2`, `E(b) = {k,l}`, then
  `N(b) = −4(D−1)c_kc_l(y_ky_l)^{M−3}(y_k − y_l)² ≠ 0`.

Summing multiplicities of the nonzero polynomial `N`:

**Theorem 3.1.** For every `A ⊆ F_p` with `3 ≤ |A| = M ≤ (p+1)/2`,

    (2M − 2)·N₀(A) + (M − 2)·N₁⁼(A)  ≤  2d − 2.

Consequences, with `N₁ = N₀ + N₁⁼` (points with at most one bad partner):

1. `M·N₁ ≤ M(2d − 2 − M·N₀)/(M − 2)`; in particular `M·N₁ ≤ (2+O(1/M))·d`.
   This is the constant 2 of the second moment and of the `algebra`
   worker's Hankel bound `N₁ ≤ (2d−2)/(M−2)`, sharpened by the `−M·N₀` term.
2. **Trade-off.** Complete points and one-bad points compete for one budget.
   If `A` is Hanson–Petridis extremal (`M·N₀ ≈ d`), then
   `(M−2)·N₁⁼ ≤ 2d − 2 − (2M−2)N₀ ≈ 0`: there are essentially no one-bad
   points at all.
3. **Riemann–Hurwitz form.** `R = V/W` is a rational map `P¹ → P¹` of degree
   at most `d < p`, hence tamely ramified, and after cancelling the common
   factor `Π_{b: E(b)=∅}(x−b)^{M−1}` its ramification divisor is the divisor
   of `N` there. One-bad points are ramification points of index `≥ M−1`
   lying over the `M` branch values `a_1,…,a_M`. Robust HP at `e = 1`,
   namely `M·N₁ ≤ (1+o(1))d`, holds **if and only if** at most about half of
   the total ramification `2 deg R − 2` of `R` sits at `F_p`-rational points
   over `A`. The inequality above is the statement that all of it could.
4. **Prime sensitivity.** The argument uses `D ≤ p−1` (so `(D)_i ≠ 0`),
   `d, M−1 < p` (leading coefficient of `N`) and tameness (`deg R < p`).
   Over `F_{p²}` all three fail, consistent with the `F_{p²}` test.

## 4. General `e` and the Plücker barrier (PROVED, reproduces known constant)

For the `(e+1)`-dimensional system `⟨Y_0,…,Y_e⟩` of Lemma 2.2, the Wronskian
has degree `≤ (e+1)(d − e)`. At a point with bad set `E(b)`, `|E(b)| ≤ e`,
some member (`F_{A∖E}` for any `E ⊇ E(b)`) vanishes to order `≥ M − e`, so
the Wronskian vanishes to order `≥ M − 2e`. When the Wronskian is nonzero,
`(M − 2e)·N_e ≤ (e+1)(d − e)`: constant `e+1`, the Plücker bound for a
rational curve of degree `d` in `P^e` (the same constant as the `algebra`
worker's Proposition 2.7). RHP asks for constant `1 + o(1)`, i.e. that the
osculating hyperplanes of this specific curve at `F_p`-points do not come
from the special family `{Σ_j σ_j(a_E) y_j = 0}` for many different `E`.

## 5. What would finish RHP at `e = 1`, and what does not

- Subset-HP (HP applied to every `A∖K`) gives `N₁ ≤ 4d/M` by splitting `A`
  into two halves, and no better from those constraints alone: the linear
  programme in the counts `x_∅, x_k` with every constraint
  `x_∅ + Σ_{k∈K} x_k ≤ d/(M − |K|)` has optimum `4d/M` (for even `M`, at
  `x_∅ = 0`, `x_k = 4d/M²`). This is an elementary LP computation, not in the
  verifier; the `robust` worker treats the set-system version. The pencil
  gives `(2+O(1/M))d/M`. Neither reaches `(1+o(1))d/M`.
- Needed: an argument that the polynomial `N ∝ (D−1)F′² − D·F·F″`, which is
  `F^{2−1/D}·(F^{1/D})″` up to a constant, cannot have more than about half
  of its `2d−2` roots (with multiplicity) at rational points of type `≤ 1`.
  Its rational values are the `2×2` Hankel determinants
  `T_{M−2}² − T_{M−1}T_{M−3}` of the bad-set sequences.
- Data (exact, small `p`, random `A`): the budget fraction
  `((2M−2)N₀ + (M−2)N₁⁼)/(2d−2)` is recorded per instance in the results
  file. The largest `M·N₁/d` observed is `3.5` at `p = 37`, `M = 3`, where
  `e = 1` means `η = 1/3` (a regime RHP does not constrain).

## References

Hanson–Petridis, arXiv:1905.09134, Theorem 1.2 and §2 (local text
`sources/sigma-hanson-petridis-1905.09134.txt`). Riemann–Hurwitz for tame
maps of the projective line and the Plücker formula for linear systems on
`P¹` are standard; only the elementary polynomial identities above are used
in the proofs, and each is checked in the verifier.
