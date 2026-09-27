# Stepanov wave, worker `sharpen`: the optimal use of the Hankel-minor inequality (2026-09-27)

**Status:** PROVED, taking as given the Hankel-minor inequality (★) (Thm 2.1 of
`research/stepanov-robust-2026-09-26.md`, cited below as [R]), its order-`k` version
([R, Prop. 2.9]) and, where stated, Weil's bound for the Legendre symbol: (1) a sharpened
constant saving `|S(A,B)| ≤ (u(w) + O(p^{−1/2}))|A||B|` for `|A||B| ≥ (1/2+κ)p`, where
`w = 1/(1+2κ)` and `u(w) = (√(12w−3w²) − w)/2`. The saving is
`η★(κ) = 1 − u(w) = (4/3)κ² − (40/9)κ³ + O(κ⁴)`, against `κ² − 3κ³ + …` in [R, Thm 2.5]
(Theorem 1.5). (2) `η★` is exactly the best saving that the family (★) can certify. The
linear programme in the counts `n_j` has value `u(w) − Θ(1/m)` in the cases computed, and
explicit feasible profiles show that it tends to `u(w)` (Prop. 2.1). So the (★)-saving is
quadratic in `κ`, not linear. (3) Using Weil, the bounded-`|A|` examples show exactly the
following. The upper bound `2κ/(1+2κ)` on the true saving ([R, Remark 2.7]) is the best such
example for `κ ≤ 11/48`. For `κ > 11/48` it improves strictly: `|A| = 5` examples give a
saving `≤ 0.3625` at `κ = 1/2` (Cor. 3.3). (4) Balanced
sets: with `min(|A|,|B|) ≥ k`, the threshold is exactly `k/2^k` for all pairs outside a
"balanced window" `M(k,κ) < |A| ≤ |B| < 12√p` (Thm 4.3 and Prop. 4.1, sharp). Inside the
window only the threshold `1/2` is proved. (5) Characters of order `k`: the threshold is
`|A||B| ≈ (p−1)/k` (sharp), with saving at least `θ(1−θ)(1−cos(2π/k))`, where
`θ = (1−u(1/(1+λ)))/2 ≈ λ²/6` (Thm 5.1, Cor. 5.2). REFUTED: "optimising the weights and
combining (★) over several `e` gives a saving linear in `κ`". The witnesses are the explicit
LP-feasible profiles of Prop. 2.1 (verifier §D3). OPEN: whether the true saving is linear in
`κ` (it lies between `η★(κ)` and `1 − sup_m β_m(w)`). Also open: any threshold below `1/2`
inside the balanced window, including balanced sets of size `c√p` with `c < 1/√2`, and hence
whether `τ_k → 0` without restriction.

Verifier: `experiments/stepanov_sharpen_2026_09_27.py` → `results/stepanov_sharpen_2026_09_27.json`
(numpy, scipy, sympy; 5,911,513 checks, 0 failures, 45 s CPU; about 25 min wall time
because the machine's load average was 120–260). Every inequality used is re-checked there
in exact integer or rational arithmetic unless marked "float". Nothing in this note has been
independently refereed.

---

## 0. Notation and inputs

`p` is an odd prime, `χ` the Legendre symbol with `χ(0)=0`, and `d = (p−1)/2`.
`S = S(A,B) = Σ_{a∈A,b∈B} χ(a+b)`, `m = |A|`, `n = |B|`, `r = |B∩(−A)|`, `δ_b = [b∈−A]`,
`e_b = #{a∈A : χ(a+b) = −1}`, `N_− = Σ_{b∈B} e_b`, `t̄ = N_−/(mn)`, `x̄ = N_−/n = m t̄`.
Always `S = mn − r − 2N_−`, so `S/(mn) ≤ 1 − 2t̄`.

For `0 < w ≤ 1` put

    u(w) = (√(12w − 3w²) − w)/2,        η★(κ) = 1 − u(1/(1+2κ)).

`u` is increasing and concave on `(0,1]`, and `u(1) = 1`: `u'(w) ≥ 0 ⇔ (w−1)(w−3) ≥ 0`, and
`u''(w) = −2√3/(w(4−w))^{3/2}` (sympy). Also `u'(1/4) = 1.0652…`.

**Input (★)** ([R, Thm 2.1], taken as given). For `1 ≤ m ≤ (p+1)/2` and every integer
`0 ≤ e ≤ (m−1)/2`:
`Σ_{b∈F_p, e_b≤e} (e+1−e_b)(m − (3e+e_b)/2 − δ_b) ≤ (e+1)(d−e)`.
**Input (★)_k** ([R, Prop. 2.9]). The same holds with `d` replaced by `d_k = (p−1)/k` and
`e_b = #{a : a+b ≠ 0, (a+b)^{d_k} ≠ 1}`, for `k | p−1`, `m + d_k − 1 ≤ p−1` and
`2e ≤ min(m−1, d_k)`.
*Re-check (verifier §A).* (★) holds for every `A ∋ 0` with `p ∈ {7,11,13,17,19}` (all
admissible `m, e`; 872,312 pairs `(A,e)`). (★)_k holds for every `A ∋ 0` and
`(p,k) ∈ {(7,3),(11,5),(13,3),(13,4),(13,6),(17,4),(17,8),(19,3),(19,6),(19,9)}`
(2,386,844 pairs `(A,e)`), and for 400 random `A` for each of 13 pairs `(p,k)` with
`p ∈ {31,37,61}` (two of them with `k = 2`). There are no violations. Equality holds (e.g. `A = {0}`, `e = 0`).
**Weil** (brief, background facts; CITED there): `|Σ_x χ(f(x))| ≤ (deg f − 1)√p` for `f`
not a constant times a square.

## 1. The mean form of (★) and the optimal bias bound (PROVED)

**Lemma 1.1 (mean form).** Let `1 ≤ m ≤ (p+1)/2`, `B ⊆ F_p`, `n ≥ 1`. For every integer `E`
with `1 ≤ E ≤ (m+1)/2`, let `W_E(x) = (E−x)(m + 3/2 − (3E+x)/2)` for `x ≤ E` and
`W_E(x) = 0` for `x ≥ E`. Then

    Σ_{b∈B} W_E(e_b) ≤ E(d − E + 1 + r),   and hence   n·W_E(x̄) ≤ E(d − E + 1 + r).

*Proof.* Take `e = E−1` in (★). Since `m − (3e+e_b)/2 = m + 3/2 − (3E+e_b)/2`, the term of `b`
equals `W_E(e_b) − δ_b(E−e_b)`. Every term is `≥ 0`, because
`m − (3e+e_b)/2 − δ_b ≥ m − 2e − 1 ≥ 0`. So restricting to `b ∈ B` gives
`Σ_{b∈B, e_b≤e} W_E(e_b) ≤ E(d−e) + Σ_{b∈B∩−A, e_b≤e}(E−e_b) ≤ E(d−e) + Er`. For `e_b ≥ E`,
`W_E(e_b) = 0`.

`W_E` is convex on `[0,∞)`. On `[0,E]` it is a quadratic with leading coefficient `+1/2`.
Its left derivative at `E` is `−(m + 3/2 − 2E) ≤ 0`, because `2E ≤ m+1`. It is `0` on
`[E,∞)`. Jensen's inequality gives the second claim. ∎

**Normalisation.** Put `s = E/m`, `β = 3/(2m)`, `a = 1 − 2s + β` (so `a ≥ 1/(2m) > 0`), and
`w' = (d+r)/(mn)`. For `t̄ < s`, `W_E(x̄) = m²(s−t̄)(a + (s−t̄)/2)`, and
`E(d−E+1+r)/n ≤ E(d+r)/n = m²·s·w'`. With `y = s − t̄`, Lemma 1.1 therefore reads
`y(a + y/2) ≤ s w'`.

**Corollary 1.2 (one constraint at a time).** With the notation above,

    S/(mn) ≤ 1 − 2t̄ ≤ h_m(E,w') := −a − β + 2√(a² − a w' + (1+β) w')   for every admissible E,
    S/(mn) ≤ 1 − 2t̄ ≤ 1 − 2(1−w')/m.

*Proof.* If `t̄ ≥ s`, then `1 − 2t̄ ≤ 1 − 2s = a − β ≤ h_m`, because
`a² − aw' + (1+β)w' = a² + 2s w' ≥ a²`. If `t̄ < s`, then `y ↦ y(a+y/2)` is increasing on
`y ≥ 0`, so `y ≤ −a + √(a² + 2sw')`. Hence `1 − 2t̄ = 1 − 2s + 2y ≤ h_m` (using
`2s = 1 + β − a`).

For the second line, (★) at `e = 0` gives `m·n_0 − r ≤ d`, where
`n_0 = #{b∈B : e_b = 0}`. Every other `b` has `e_b ≥ 1`, so
`N_− ≥ n − n_0 ≥ n − (d+r)/m`. ∎

Both right-hand sides are increasing in `w'`. So each bound stays valid with `w'` replaced
by any `w ≥ w'`.

**Theorem 1.3 (the one-variable inequality).** For every integer `m ≥ 1` and every
`w ∈ [1/4, 1)`:

    min{ 1 − 2(1−w)/m ,  min_{1≤E≤(m+1)/2} h_m(E,w) }  ≤  u(w).

*Proof.* Write `ε = 1−w`, `c = w(4−w)/4`, `P(a) = a² − aw + w = (a − w/2)² + c`,
`h₀(a) = −a + 2√P(a)`, and `v = (3w + √(12w−3w²))/6`.

(F1) `3v² − 3vw + w² − w = 0`, so `√P(v) = 2v − w`, `h₀'(v) = 0` and
`h₀(v) = 3v − 2w = u(w)`.
(F2) `h₀''(a) = 2c/P(a)^{3/2} ≤ 2/√c`.
(F3) `√P` is 1-Lipschitz, since `(2a−w)² = 4P − 4c ≤ 4P`.
(F4) `q := 2v − w = √P(v) ≤ 1`, because `P(v) = v² + w(1−v) ≤ v² + 1 − v ≤ 1`. Here
`0 ≤ v ≤ 1`, where `v ≤ 1 ⇔ (w−1)(w−3) ≥ 0`.
(F5) `v − w = 2wε/(√(12w−3w²) + 3w) ≥ wε/(3−2ε)`, using
`√(12w−3w²) = √(9−6ε−3ε²) ≤ 3−ε`.
(F6) `1 − u(w) ≤ ε²/2` for `ε ≤ √3−1`. Indeed `1 − u = (3 − ε − 3R)/2` with
`R = √(1 − 2ε/3 − ε²/3)`, and `R ≥ 1 − ε/3 − ε²/3`, since squaring reduces this to
`ε² + 2ε − 2 ≤ 0`.

*Jensen step.* Let `M_J(w) = (1 + 2/(3√c))(3−2ε)/(2wε)` and suppose `m ≥ max(2, M_J(w))`.
The admissible `a_E = 1 − (2E − 3/2)/m` step by `2/m` from `1 − 1/(2m)` down to a value
`≤ 3/(2m)`. Since `v ≥ 0.40` for `w ≥ 1/4`, some `E` has `|a_E − v| ≤ 1/m`. Using
`√(X+δ) ≤ √X + δ/(2√X)`, then (F2), (F1) and (F3):

    h_m(E,w) ≤ h₀(a_E) − β(1 − w/√P(a_E)) ≤ u + 1/(√c m²) − (3/(2m))(1 − w/(q − 1/m)).

By (F5), `q − w = 2(v−w) ≥ 2wε/(3−2ε) ≥ (1 + 2/(3√c))/m > 1/m`. With (F4) this gives
`0 < q − 1/m ≤ 1` and `1 − w/(q−1/m) ≥ q − w − 1/m`. Hence
`h_m − u ≤ (3/(2m))[(1 + 2/(3√c))/m − (q−w)] ≤ 0`.

*HP step.* If `m ≤ 4/ε` and `ε ≤ √3−1`, then `1 − 2ε/m ≤ 1 − ε²/2 ≤ u` by (F6).

*Case `w ∈ [11/20, 1)`.* `ε·M_J = (1 + 2/(3√c))(3−2ε)/(2(1−ε))` is increasing in `ε`, and
at `ε = 9/20` it is `≤ 3.757 < 4` (verifier C1, exact with a rational lower bound for `√c`).
So every `m` either has `m ≤ 4/ε` (HP step) or `m > 4/ε ≥ max(8, M_J)` (Jensen step).

*Case `w ∈ [1/4, 11/20]`.* The first factor of `M_J` is at most its value at `ε = 3/4`, and
`(3−2ε)/(ε(1−ε))` has its only critical point, a minimum, at `(6−√12)/4`. Hence
`M_J ≤ 10.085` (C1), and every `m ≥ 11` is covered by the Jensen step. For `1 ≤ m ≤ 10`
the claim is a finite check. Both sides are increasing in `w`: for `h_m`,
`∂_w(a² − aw + (1+β)w) = 2s > 0`. So it suffices, on each of 64 subintervals
`[w₁,w₂]` of `[1/4, 11/20]`, that some option evaluated at `w₂` is `≤ u(w₁)`. The verifier
(C2) does this in rational arithmetic, with square roots enclosed by integer square roots.
All 640 boxes pass on the first subdivision, with margin at least `0.046`. ∎

(Float sanity scan C3: `min(HP, Jensen) − u(w) ≤ −3.3·10⁻⁷` for all `m ≤ 3000` on a
151-point grid in `[1/4, 0.999]`. The worst case is at the largest `m`, where the gap
closes like `1/m`.)

**Theorem 1.4 (pair form).** If `1 ≤ m ≤ (p+1)/2` and `w' = (d+r)/(mn) < 1`, then
`|S(A,B)| ≤ u(max(w', 1/4))·|A||B|`.

*Proof.* Corollary 1.2 and Theorem 1.3 at `w = max(w',1/4)` bound `S`. For `−S`, apply
them to `(νA, νB)` with `ν` a non-residue. This gives `S(νA,νB) = −S(A,B)` with the same
`m`, `n` and `r`. ∎

The verifier (§B) checks Lemma 1.1 (683,145 instances) and Theorem 1.4 (1,190,240 instances)
exactly. It covers every `A ∋ 0` for `p ≤ 13` and 20,000 random `A` for `p = 17`, with `B`
equal to the top-`n` and bottom-`n` sets of `F_A`, all level sets `{e_b ≤ e'}`, and random
sets. The largest observed `S/(u(w')mn)` is `0.70, 0.81, 0.81, 0.85` at `p = 11, 13, 17, 19`.

**Theorem 1.5 (sharpened constant cancellation).** Let `p ≥ 11`, `0 < κ ≤ 3/2`,
`w = 1/(1+2κ)`, `m₁ = ⌈√((1/2+κ)p)⌉`, `w₁ = (p − 1 + 2m₁)/((1+2κ)p)`. For all
`A, B ⊆ F_p` with `|A||B| ≥ (1/2+κ)p`:

    |S(A,B)| ≤ u(min(1, w₁))·|A||B| ≤ [ 1 − η★(κ) + 2.14·(√(2p)+1)/p ]·|A||B|.

*Proof.* `S` is symmetric in `(A,B)`, so assume `m ≤ n`. Put `m₀ = ⌈(1/2+κ)p/n⌉`. Then
`m₀ ≤ m`, and `m₀ ≤ m₁ ≤ √(2p)+1 ≤ (p+1)/2`, since `n ≥ √(mn)` and `p ≥ 11`. Every
`A' ⊆ A` with `|A'| = m₀` has `m₀n ≥ (1/2+κ)p` and `r(A',B) ≤ m₀ ≤ m₁`, so
`w'(A',B) ≤ (d+m₁)/((1/2+κ)p) = w₁`. If `w₁ < 1`, Theorem 1.4 and monotonicity of `u` give
`|S(A',B)| ≤ u(w₁)m₀n`. Note `w₁ ≥ w ≥ 1/4`. Averaging,
`S(A,B) = (m/m₀)·avg_{A'} S(A',B)`, finishes the case `w₁ < 1`. If `w₁ ≥ 1` the bound is
trivial.

For the second inequality: `u` is concave with `u' ≤ u'(1/4) < 1.066` on `[1/4,1]`, and
`w₁ − w = (2m₁−1)/((1+2κ)p) ≤ 2(√(2p)+1)/p`. When `w₁ > 1`, the right side is `≥ 1` because
`u(w) + 1.066(1−w) ≥ 1`. ∎

| κ | η★(κ) (this note) | (1−√w)² ([R] Thm 2.5) | η★/κ² | upper bound on the true saving from bounded `|A|` (Cor. 3.3) | 2κ/(1+2κ) |
|---|---|---|---|---|---|
| 0.01 | 1.290e-4 | 0.971e-4 | 1.290 | 0.01961 (m=2) | 0.01961 |
| 0.05 | 2.844e-3 | 2.166e-3 | 1.137 | 0.0909 (m=2) | 0.0909 |
| 0.1 | 9.838e-3 | 7.591e-3 | 0.984 | 0.1667 (m=2) | 0.1667 |
| 0.25 | 0.04234 | 0.03367 | 0.677 | 0.3167 (m=5) | 0.3333 |
| 0.5 | 0.10436 | 0.08579 | 0.417 | 0.3625 (m=5) | 0.5 |
| 1 | 0.20924 | 0.17863 | 0.209 | 0.4777 (m=7) | 0.6667 |
| 1.5 | 0.28647 | 0.25 | 0.127 | 0.5417 (m=6) | 0.75 |

`η★(κ) > (1−√w)²` on a 5000-point grid of `w ∈ (0,1)` (float; the ratio tends to `4/3`).
The two methods differ as follows. [R] bounds every weight crudely by `m − 2e` and uses
one `e` chosen for the worst case. Lemma 1.1 keeps the exact weight. Its convexity turns
the whole distribution of the `e_b` into a statement about the mean. The optimal `e` is
`e ≈ s*m` with `s* = (1−v)/2 ≈ 2κ/3`.

## 2. Optimality within (★): the linear programme (PROVED)

For given `(p, m, n)` with `m ≤ (p+1)/2`, the task's LP has variables `n_j ≥ 0` (the number
of `b ∉ −A` with `e_b = j`, `0 ≤ j ≤ m`) and `n'_j ≥ 0` (those `b ∈ −A` with `e_b = j`).
Its constraints are `Σ n_j + Σ n'_j = n`, `Σ n'_j ≤ min(m,n)`, and (★)_e for every
`0 ≤ e ≤ (m−1)/2`, with weights `(e+1−j)(m − (3e+j)/2) − [δ](e+1−j)`. Its objective is
`S = Σ(m−2j)n_j + Σ(m−1−2j)n'_j`. Let `V(p,m,n)` be its value. Any nonnegative combination
of the inequalities (★) for `A` bounds `S` by at least `V`.

**Proposition 2.1.** (i) If `w'' := (d + min(m,n))/(mn) < 1`, then
`V/(mn) ≤ u(max(w'',1/4))`.

(ii) Let `w ∈ [1/4,1)` and `t* = (1−u(w))/2`, and let `m ≥ 2`, `d ≥ 1` and `J` be integers
with `t = J/m ≥ t*` and

    (t − t*)(1/2 − 1/(2m) − t) ≥ 3/(2m) + m/(2d).

With `n = ⌊d/(mw)⌋`, the point `n_J = n` (all other variables `0`) is feasible. So
`V/(mn) ≥ 1 − 2J/m` while `d/(mn) ≥ w`.

(iii) Consequently, for fixed `w` and `m → ∞` with `m/d → 0`, taking
`J = ⌈(t* + δ)m⌉` with `δ = 5(3/(2m) + m/(2d))` gives `V/(mn) = u(w) − O(1/m + m/d)`.
Combining (★) over any set of `e`, with the exact weights, certifies no saving larger
than `η★ + o(1)`. Balanced pairs with `m, n ≍ √p` and `mn ≈ (1/2+κ)p` have `m → ∞` and
`m/d → 0`. So `inf (1 − V/(mn))` over `m ≤ min(n, (p+1)/2)` and `mn ≥ (1/2+κ)p` tends to
exactly `η★(κ)` as `p → ∞`. The lower bound comes from (i): `n ≥ √(mn) ≥ √(p/2)` gives
`w'' ≤ d/(mn) + 1/n ≤ w + O(p^{−1/2})`. This is the best bound on the saving obtainable
from (★).

*Proof.* (i) Lemma 1.1 and Corollary 1.2 use only the LP constraints. The `δ`-terms enter
as `Σ n'_j ≤ min(m,n)`. Then apply Theorem 1.3.

(ii) Put `K(s,t) = (s−t)(1 − (3s+t)/2)` and `τ(s) = 1 − s − √((1−2s)² + 2sw)`. With
`a = 1−2s`, `1 − 2τ(s) = h₀(a)`, so `max_{s∈ℝ} τ(s) = (1 − min h₀)/2 = t*` by (F1).

For `0 < s ≤ 1/2` we have `K(s,t*) ≤ sw`. If `τ(s) ≥ 0`, then `K(s,τ(s)) = sw`, and
`K(s,·)` is decreasing on `[0,s]` (`∂_tK = −(1−s−t)`) with `t* ≥ τ(s)`. If `τ(s) < 0`, then
`K(s,t*) ≤ K(s,0) < sw`. For `1/2 < s ≤ 1/2 + 1/(2m)`, directly
`K(s,t*) ≤ (s−t*)(1/4 − t*/2) ≤ s/4 ≤ sw`.

The constraints with `e < J` have zero left side. For `J ≤ e ≤ (m−1)/2` put `E = e+1` and
`s = E/m ∈ (t, 1/2 + 1/(2m)]`. The left side is
`n m²[K(s,t) + β(s−t)] ≤ (dm/w)[sw − X + βs]`, where
`X = (t−t*)(1/2 − 1/(2m) − t) ≤ K(s,t*) − K(s,t)` by integrating `∂_tK`. The right side is
`(e+1)(d−e) ≥ smd − sm²/2`. So the constraint follows from `X ≥ βs + swm/(2d)`, which the
hypothesis implies (`s ≤ 1`, `w ≤ 1`). The objective is `(m−2J)n`.

(iii) `t* ≤ 0.144` for `w ≥ 1/4`. With `J = ⌈(t*+δ)m⌉` we have
`t ∈ [t*+δ, t*+δ+1/m]`. For `m ≥ 40` and `δ ≤ 0.05`,
`X ≥ δ(1/2 − 3/(2m) − t* − δ) ≥ 0.2δ = 3/(2m) + m/(2d)`, so (ii) applies. The profile
bias is then `1 − 2t ≥ u(w) − 2δ − 2/m`, and `n ≥ d/(mw) − 1`. ∎

*Checks.* Exact rational LPs (sympy simplex), primal and dual solved separately, with
feasibility and equal objectives verified exactly, for 180 triples: `p ∈ {101, 1009, 10007}`,
`m ∈ {2,…,30}`, `n = ⌈(1/2+κ)p/m⌉`, `κ ∈ {1/20, 1/10, 1/4, 1/2, 1}` (D1). All satisfy (i).
The largest `V/(mn)/u(w'')` is `0.9964`. The optimum uses `r > 0` in 9 of the 180 cases.

Float LPs with `r = 0` and `d = ∞` (D2) give `m·(u(w) − V/(mn)) ≈ 0.52, 0.20, 0.10, 0.05`
for `w = 0.5, 0.8, 0.9, 0.95` and `m` from 50 to 800. So `V/(mn) = u(w) − Θ(1/m)` in these
cases, and the full LP improves on the best single `E` of Corollary 1.2 only at the third
decimal for `m ≤ 25`, and not at all asymptotically.

The profiles of (ii) are checked constraint by constraint in exact integers (D3; 5,714
constraints). Examples: `w = 9/10`, `m = 4000`, `d = 10^{10}` gives profile bias `0.9950`
against `u = 0.99655`. `w = 2/3`, `m = 2000` gives `0.954` against `0.9577`.

**Answer to task 1.** No, the (★)-saving is not linear in `κ`. The extremal (★)-feasible
configuration is a *uniformly slightly defective near-biclique*. Every `b ∈ B` has
`e_b ≈ t*m` with `t* ≈ (2/3)κ²`, and `mn ≈ (1/2+κ)p`. The family (★), for every `e` and
with the exact weights, cannot exclude it. The same holds with the roles of `A` and `B`
swapped, since the bound depends only on `t̄` and `mn/d`. Applied to both sides, the doubly
regular profile `e_b ≡ tm`, `e_a ≡ tn` satisfies both families when `m, n → ∞` and
`m, n = o(p)`: this is the argument of (ii) applied on each side. This last point is a
statement about the counts only; realisability by a 0/1 matrix is not claimed. A linear
saving needs an input that excludes this profile. (★) applied to subsets `A' ⊆ A` is not
covered by (iii) (HEURISTIC: the same profile, restricted to large subsets, satisfies those
constraints too).

## 3. Comparison with the truth

**Lemma 3.1 (Weil counts; PROVED given Weil).** For every `A` with `|A| = m` and every
`0 ≤ j ≤ m`, `N_j(A) := #{b∈F_p : e_b = j}` satisfies

    |N_j(A) − C(m,j)p/2^m| ≤ C(m,j)(m/2)(√p + 1) + 2m.

*Proof.* For `b ∉ −A`,
`[e_b = j] = 2^{−m} Σ_{|J|=j} Σ_{I⊆A} (−1)^{|I∩J|} χ(Π_{a∈I}(a+b))`. Put
`T(I) = Σ_{b∉−A} χ(P_I(b))`. Then `T(∅) = p−m`. For `I ≠ ∅`, the complete sum is `0`
(`|I| = 1`), `−1` (`|I| = 2`) or at most `(|I|−1)√p` in absolute value (Weil; `P_I` is
squarefree). At most `m − |I|` terms are removed. So `|T(I)| ≤ (|I|−1)√p + m − |I|`, and
`2^{−m} Σ_{I≠∅}|T(I)| ≤ (m/2)(√p+1)`. Finally account for `b ∈ −A` (at most `m`) and
`C(m,j)m/2^m ≤ m`. ∎
(F1: 302 exact checks at `p ∈ {10007, 100003}`, `m ≤ 8`; the largest ratio of the deviation
to the bound is `0.15`.)

**Proposition 3.2 (fixed `|A|`; PROVED given Weil).** For `0 < q ≤ 1` let `β_m(q)` be the
mean of `1 − 2j/m` over the lowest mass `q` of Binomial`(m,1/2)`. For fixed `m`, `w` with
`q = 1/(2mw) ≤ 1`, and `n = ⌈p/(2mw)⌉`:
`max_{|A|=m,|B|=n} S(A,B)/(mn) = β_m(q) + O_m(p^{−1/2})`.
*Proof.* `S(A,B) ≤ Σ_{b∈B}(m − 2e_b)`, and the best `B` takes the smallest `e_b`. The counts
are given by Lemma 3.1 for every `A` of size `m`, which gives both the upper bound and a
matching construction. ∎
(F3: at `p = 100003`, random `A` with `m ∈ {2,3,4,6}` realise `β_m` to within `2·10⁻³`.)

**Corollary 3.3 (upper bounds on the true saving).** As `p → ∞`, the optimal saving `η(κ)`
over all `|A||B| ≥ (1/2+κ)p` satisfies `η(κ) ≤ 1 − sup_{m≥1} β_m(1/(2mw))`. For
`w ≥ 1/3`, `β_2 = w`, which recovers [R, Remark 2.7]'s `2κ/(1+2κ)`. For `w ∈ [8/15, 1]`,
`β_5 = 3/5 + w/8`, which exceeds `w` exactly when `w < 24/35`, i.e. `κ > 11/48`.
Conversely, `β_m(1/(2mw)) ≤ w` for all `m ≥ 3` and `w ∈ [24/35, 1]`. For `m ≥ 200`,
Hoeffding gives `β_m ≤ 1/2 + 2m e^{−m/8} < 24/35`. For `3 ≤ m < 200`, `β_m(1/(2mw)) − w`
is piecewise linear in `w`, and the verifier (F4) checks it exactly at every breakpoint.
`β_1 = 2w − 1 < w`. So `2κ/(1+2κ)` is the best bounded-`|A|` upper bound exactly for
`κ ≤ 11/48`. Beyond that it improves; see the table in §1 (e.g. `0.3625` at `κ = 1/2`,
attained at `m = 5`). ∎

**Exhaustive data (§E; PROVED for the stated `p`).** For every `p ≤ 23`, every `A ∋ 0`
(translation invariance; `|S|` is covered by the non-residue dilation) and every `n`, the
exact `max_B S(A,B)` is computed. All `1,052` pairs `(m,n)` with `w'' < 1` satisfy
Theorem 1.4 and Corollary 1.2 exactly. All 558 with `mn > d` satisfy `S ≤ V` (float LP).
The maximal bias over `|A||B| ≥ (1/2+κ)p` is, at `p = 23`: `0.917` (`κ = 0`, `m = 1`),
`0.8125` (`κ ∈ [0.05, 0.15]`, `4×4`), `0.75` (`κ ∈ [0.2, 0.35]`), `0.64` (`κ = 1/2`) and
`0.528` (`κ = 1`). At `p = 19`, `κ = 0.05`, it is `0.833`. These are far below
`u(w) ≈ 0.997`. At such small `p` the extremisers are `3×3`, `3×4` and `4×4` rectangles.
This is a finite-size effect: the `m = 2` example is not yet dominant.

## 4. Balanced sets (task 2)

Let `τ_k` be the infimum of `τ` such that for every `κ > 0` there are `c > 0` and `p₀` with
`|S| ≤ (1−c)|A||B|` whenever `p ≥ p₀`, `min(|A|,|B|) ≥ k` and `|A||B| ≥ (τ+κ)p`.
Theorem 1.5 gives `τ_k ≤ 1/2` for every `k`. `τ_1 = τ_2 = 1/2` by the `m = 1, 2` examples.

**Proposition 4.1 (lower bound; PROVED given Weil).** `τ_k ≥ k/2^k` for every `k ≥ 1`.
*Proof.* By Lemma 3.1 every `A` of size `k` has `N_0(A) ≥ p/2^k − (k/2)(√p+1) − 2k`.
`B = B_0(A)` is a complete biclique with `|A||B| = kp/2^k − O(k²√p)`, `|B| ≥ k` for large
`p`, and `S = |A||B| − r ≥ |A||B| − k`. ∎

**Lemma 4.2 (fourth moment; PROVED given Weil).**
`Σ_{b∈F_p} F_A(b)⁴ ≤ 3m²p + 3m⁴√p`, where `F_A(b) = Σ_{a∈A} χ(a+b)`.
*Proof.* Expand over `(a₁,…,a₄) ∈ A⁴` and write `Π(b+a_i) = Q(b)²R(b)` with `R` squarefree.
There are at most `3m²` tuples with `deg R = 0`; each contributes at most `p`. If
`deg R = 2`, the sum is `−1` minus at most one term, so at most `2`. If `deg R = 4`, Weil
gives at most `3√p`. ∎ (G2: 60 exact checks; the largest ratio to the bound is `0.70`.)

**Theorem 4.3 (threshold `k/2^k` outside the balanced window; PROVED given Weil).**
Let `k ≥ 2`, `0 < κ ≤ 1`, `τ = k/2^k + κ`, `M = ⌈12/τ⌉`, and `p ≥ max(25, (2M²/κ)²)`.
Suppose `k ≤ |A| ≤ |B|` and `|A||B| ≥ τp`. If `|A| ≤ M` then
`|S| ≤ (1 − κ/(τM))|A||B|`. If `|A| > M` and `|B| ≥ 12√p` then `|S| ≤ 0.841|A||B|`.

*Proof.* (a) `S ≤ (m−2)n + 2N_0(A)`, since every `b` with `e_b ≥ 1` has `F_A(b) ≤ m−2`.
The same holds for `−S` via `νA`. For `m ≥ k ≥ 2`, `m/2^m ≤ k/2^k`. Lemma 3.1 gives
`N_0/n ≤ (m/2^m)/τ + M²(√p+5)/(2τp) ≤ 1 − κ/τ + M²/(τ√p) ≤ 1 − κ/(2τ)`. So the saving is
`(2/m)(1 − N_0/n) ≥ κ/(τM)`.

(b) Hölder and Lemma 4.2 give
`|S|/(mn) ≤ (3p/(m²n) + 3√p/n)^{1/4} ≤ (3/(τm) + 1/4)^{1/4} ≤ 2^{−1/4} < 0.841`. ∎
(G2: 180 Hölder checks.)

So `τ_k = k/2^k` exactly (sharp by Prop. 4.1) for pairs outside the **balanced window**
`W = {M(k,κ) < |A| ≤ |B| < 12√p}`, and `k/2^k → 0`. Inside `W` the only tool is
Theorem 1.5, with threshold `1/2`. The second moment (Chung) gives nothing below
`|A||B| ≈ p`, and Weil counting says nothing about sets of size `O(√p)`: the error term
`(m/2)√p` already exceeds `|B|`. (★) is LP-consistent with complete bicliques of every
shape up to `mn = d`, since `n_0 = d/m` satisfies every (★)_e whenever `m ≤ 3d/2`.

**Balanced sets of size `c√p`, `c < 1/√2`.** These lie in `W` for large `p`. Nothing is
proved for them (OPEN). A constant saving there would imply, in particular, that the
bipartite graph `χ(a+b) = 1` has no `K_{c√p, c√p}`. That is a constant-factor improvement
of Hanson–Petridis for balanced bicliques. For cliques the best prime-order bound stated in
Yip, arXiv:2304.13213 (v5, Aug 2025; introduction, fetched via a summarising tool; not
checked verbatim) is Hanson–Petridis's `√(p/d)+1`. So whether `τ_k → 0` without the window
restriction is OPEN.

*Exhaustive biclique data (G1, `p ≤ 23`).* `M_m(p) = max_{|A|=m} N_0(A)` is
`[12,6,4,3,2,2,1,…]` at `p = 23`. The largest complete biclique with `min ≥ 3` has
`mn/p = 0.52, 0.47, 0.53, 0.69` at `p = 23, 19, 17, 13`, against `3/8`. These are pure
small-`p` effects: the Weil error exceeds the main term. HP (`m·M_m ≤ d + m`) holds in all
44 cases.

## 5. Characters of order `k` (task 3; PROVED given [R, Prop. 2.9])

Let `k ≥ 2`, `k | p−1`, `ψ` a character of order `k` (`ψ(0) = 0`), `d_k = (p−1)/k`,
`H_k = ker ψ`, `ζ = e^{2πi/k}`, `S_ψ = Σ_{a,b} ψ(a+b)`. For `ω ∈ μ_k`, call `(a,b)`
**ω-bad** if `a+b ≠ 0` and `ψ(a+b) ≠ ω`. Let `N_ω = #{(a,b) : ψ(a+b) = ω}` and
`M = mn − r = Σ_ω N_ω`.

**Theorem 5.1.** Let `1 ≤ m ≤ d_k + 1`, `w_k = (d_k + r)/(mn) < 1`, and
`θ = (1 − u(max(w_k, 1/4)))/2 ∈ [0, 1/2]`. Then for every `ω ∈ μ_k` the number of ω-bad
pairs is at least `θmn`, so `N_ω ≤ (1−θ)M`. Moreover

    |S_ψ(A,B)| ≤ M·|1 − θ + θζ| = M·√(1 − 2θ(1−θ)(1 − cos(2π/k))).

*Proof.* Choose `g` with `ψ(g) = ω` and put `A' = g^{−1}A`, `B' = g^{−1}B`. Then
`ψ(a'+b') = ω̄ψ(a+b)`, so the ω-bad pairs of `(A,B)` are the pairs of `(A',B')` counted by
the `e_{b'}` of [R, Prop. 2.9]. The hypotheses of (★)_k hold for every `e ≤ (m−1)/2`,
because `m − 1 ≤ d_k` and `m + d_k − 1 ≤ p − 1`. Lemma 1.1, Corollary 1.2 (whose HP part is
(★)_k at `e = 0`) and Theorem 1.3 use nothing else, so the mean number of ω-bad partners
satisfies `t̄_ω ≥ θ`. Then `N_ω = M − t̄_ω mn ≤ (1−θ)M`.

For the modulus: `|Σ_ω ωN_ω| = max_{|z|=1} Σ_ω N_ω Re(z̄ω)`. For fixed `z`, the linear form
is maximised on `{N ≥ 0, ΣN = M, N_ω ≤ (1−θ)M}` by putting `(1−θ)M` on the largest
coefficient and `θM` on the second (as `θ ≤ 1/2`). If the nearest root is at angle
`α ∈ [0, π/k]` from `z`, the second nearest is at `2π/k − α`. The value is
`M·Re(e^{−iα}(1 − θ + θζ)) ≤ M|1 − θ + θζ|`. ∎

**Corollary 5.2 (bias statement).** Let `0 < λ ≤ 3` and `p ≥ (1+λ)k + 1`. Put
`m₁ = ⌈√((1+λ)d_k)⌉`, `w₁ = (d_k + m₁)/((1+λ)d_k)` and
`θ₁ = (1 − u(min(1, w₁)))/2`. If `|A||B| ≥ (1+λ)(p−1)/k`, then

    |S_ψ(A,B)| ≤ √(1 − 2θ₁(1−θ₁)(1 − cos(2π/k))) · |A||B|.

As `p → ∞`, `θ₁ → θ(λ) = (1 − u(1/(1+λ)))/2 = λ²/6 + O(λ³)`. The saving is therefore at
least `θ(1−θ)(1−cos(2π/k)) ≈ (λ²/6)(1 − cos(2π/k))`. For `k = 2` this is Theorem 1.5 with
`λ = 2κ`, because `|1 − 2θ| = u`.

*Proof.* Subsample as in Theorem 1.5. Use `m₀ = ⌈(1+λ)d_k/n⌉ ≤ m₁ ≤ d_k + 1` and the
triangle inequality for the average. Then `w_k(A') ≤ w₁`, and the bound is decreasing in
`θ` and increasing in `w`. ∎

**Threshold and sharpness.** `A = {0}`, `B = H_k` gives `|A||B| = (p−1)/k` and
`S_ψ = |A||B|`, so the threshold `1/k` is sharp. Adding `λd_k` points of `gH_k`
(`ψ(g) = ζ`) gives `|S_ψ|/(|A||B|) = |1 + λζ|/(1+λ)`, so the saving is at most
`1 − √(1 + 2λcos(2π/k) + λ²)/(1+λ) ≈ λ(1 − cos(2π/k))` (H3: exact values at
`p ∈ {1009, 10009}`, `k ∈ {3,4,6}`). So the order-`k` problem has the same `λ²` versus `λ`
gap as `k = 2`.

*Checks (§H).* Lemma 1.1 per value (91,199 checks), `N_ω ≤ (1−θ)M` (456,706 checks) and
`|S_ψ|² ≤ M²(1 − 2θ(1−θ)(1−cos(2π/k)))` (112,417 checks). For `k ∈ {3,4,6}` the last is
exact, since `|S_ψ|² = Σ N_iN_j cos(2π(i−j)/k)` is rational, with `θ` replaced by a
rational lower bound. Sets used: exhaustive or 6,000 sampled `A ∋ 0` at `p = 13, 19`, and
random `A` at `p = 31, 37, 61`, with greedy and random `B`. The vertex lemma is checked by
brute force for `k ∈ {3,4,6}`, `M ≤ 12` (36 cases).

## 6. Open obligations

1. Referee (★) and (★)_k ([R, Thm 2.1, Prop. 2.9]). Everything here is conditional on them.
2. **Linear saving (OPEN).** `η★(κ) ≤ η(κ) ≤ 1 − sup_m β_m`. The gap is exactly the
   uniformly defective near-biclique of §2 with `m → ∞`, `m = O(√p)`. Excluding it needs
   information beyond (★): Weil does this for bounded `m` only.
3. **Balanced window (OPEN).** Any threshold `< 1/2` for `M < |A| ≤ |B| < 12√p` would improve
   Hanson–Petridis for balanced bicliques. `τ_k → 0` is proved only outside the window.
4. The order-`k` statement inherits obligations 1–2. The Weil-range analogue of §3–4 for
   `ψ` (thresholds `k'/k^{k'}`) was not written, since only the quadratic Weil bound is
   among the brief's background facts.
5. Theorem 1.3's finite range `m ≤ 10`, `w ∈ [1/4, 11/20]` is a computer-assisted
   (rational interval) verification. An analytic proof there would be cleaner.
