# Referee report `referee2` on the robust Hanson–Petridis theorem (2026-09-27)

**Status:** PROVED (independent re-derivation, §§1–7): Lemma 1.1, Theorem 2.1 (all five steps),
Lemma 2.2, Corollary 2.3, Theorem 2.4, Theorem 2.5 and Corollary E of the Stepanov wave are
**correct as stated in `research/stepanov-robust-2026-09-26.md`**; I found no mathematical error.
Verdicts: all CORRECT. Presentation fixes (§8): the pass summary's Theorem B omits the hypothesis
`|A| ≤ (p+1)/2` (its Theorem C omits it too, but I prove C holds without it); the robust note's
remark "`≈ 0.92` at `e ≥ 1`" is wrong, since `(★)` is attained with equality at `e = 1`; the
summary's `(p+3)/2` example needs `p ≡ 1 (mod 4)`. REFUTED: nothing in the claim; no
counterexample among 257,753,636 exact checks (§9). OPEN: whether the threshold `1/2` is sharp
when both `|A|,|B| → ∞` (the sharpness witness has `|A| = 2`); the gap `κ² ≲ η(κ) ≲ 2κ`; novelty
(negative literature search only, §10).

Verifier (independent; the worker's verifier was neither read nor imported):
`experiments/stepanov_referee2_2026_09_27.py` → `results/stepanov_referee2_2026_09_27.json`
(`/opt/miniconda3/bin/python3`, numpy only; 60–90 s CPU, 70 s wall when the machine is idle;
257,753,636 checks including 4,180 for the count in §8.2; 0 failures).

Object under review: `research/stepanov-robust-2026-09-26.md` (Lemma 1.1, Thm 2.1, Lemma 2.2,
Cor 2.3, Thm 2.4, Thm 2.5) and Corollary E of `research/stepanov-pass-summary-2026-09-27.md`.

---

## 0. Notation

`p` odd prime, `d = (p−1)/2`, `χ` Legendre with `χ(0)=0`. `A = {a_1,…,a_m} ⊆ F_p` with
`1 ≤ m ≤ (p+1)/2`; `D = d+m−1`. Note `D ≤ p−1 ⇔ m ≤ (p+1)/2 ⇔ m−1 ≤ d`; this is the only role of the size
hypothesis on `A`. `c_k = Π_{l≠k}(a_k − a_l)^{−1}`, `F(x) = −1 + Σ_k c_k (x+a_k)^D`.
For `b ∈ F_p`: `y_k = b + a_k`, `E(b) = {k : χ(y_k) = −1}`, `e_b = |E(b)|`, `δ_b = [b ∈ −A]`.
For `B`: `n = |B|`, `r = |B ∩ (−A)|`, `N_− = Σ_{b∈B} e_b`. Since each `b ∈ B ∩ (−A)` has exactly one
`a` with `a + b = 0`, the pairs split as `r` zeros, `N_−` non-residues, `mn − r − N_−` residues, so
`S(A,B) = mn − r − 2N_−` (PROVED).

**Fact 0 (PROVED).** `Σ_k c_k a_k^j = [j = m−1]` for `0 ≤ j ≤ m−1` (it is the coefficient of
`x^{m−1}` in the Lagrange interpolant of `x^j` on `A`). By the binomial theorem the same holds with
`a_k` replaced by `b + a_k`, for every `b`.

**Order criterion (PROVED).** For `P ∈ F_p[x]` and `M ≤ p`: `ord_b P ≥ M` iff `P^{(j)}(b) = 0` for
`j < M` (the ordinary derivative is `j!` times the Hasse derivative and `j! ≢ 0` for `j < p`). This
is used with `M ≤ m ≤ (p+1)/2`.

## 1. Lemma 1.1 (local form of `F`) — CORRECT

(a) Coefficient of `x^{D−l}` in `F + 1` is `C(D,l)·Σ_k c_k a_k^l`; by Fact 0 it is `0` for
`l ≤ m−2` and `C(D, m−1)` for `l = m−1`. As `D ≤ p−1`, `C(D,m−1) ≢ 0 (mod p)` (Lucas). Since
`d ≥ 1` the constant `−1` does not interfere, so `deg F = d`, leading coefficient `C(D,m−1)`.

(b) `b ∉ −A`, `0 ≤ j ≤ m−1`. For `j ≥ 1`, `F^{(j)}(b) = (D)_j Σ_k c_k y_k^{D−j}`, and
`y^{D−j} = y^d·y^{m−1−j} = χ(y) y^{m−1−j}` because `D − j = d + (m−1−j)` with `m−1−j ≥ 0` and
`y^d = χ(y)` for `y ≠ 0`. Write `χ(y_k) = 1 − 2[k∈E(b)]`. The unweighted sum
`Σ_k c_k y_k^{m−1−j}` is `[j=0]` by Fact 0. Hence `F^{(j)}(b) = −2(D)_j Σ_{k∈E} c_k y_k^{m−1−j}` for
`j ≥ 1`; for `j = 0` the `−1` cancels the `1`, giving the same formula with `(D)_0 = 1`.

(c) `b = −a_{k₀}`, so `y_{k₀} = 0` and `k₀ ∉ E(b)`. For `j ≤ m−1 < D` the `k₀` term of
`Σ_k c_k y_k^{D−j}` is `0`. The remaining unweighted sum is
`Σ_{k≠k₀} c_k y_k^{m−1−j} = [j=0] − c_{k₀}·0^{m−1−j}`, i.e. `[j=0]` for `j ≤ m−2` (using `m ≥ 2` at
`j = 0`) and `−c_{k₀}` for `j = m−1`. This gives (c), including
`F^{(m−1)}(b) = (D)_{m−1}(−c_{k₀} − 2Σ_E c_k)`. (For `m = 1`, `F(−a) = −1 = −c_1`: consistent.)

(d) `R_E(x) = Σ_{k∈E} c_k (x+a_k)^D` has `R_E^{(j)}(b) = (D)_j Σ_E c_k χ(y_k) y_k^{m−1−j}
= −(D)_j Σ_E c_k y_k^{m−1−j}` for `j ≤ m−1`. Comparing with (b),(c):
`(F − 2R_{E(b)})^{(j)}(b) = 0` for `0 ≤ j ≤ m−1−δ_b`, so by the order criterion
`ord_b(F − 2R_{E(b)}) ≥ m − δ_b`. ∎

Numerics (Part A-local): (a)–(c) at every `b` and every `j ≤ m−1`, and (d) in the form of Step 1
below, for **every** `A ∋ 0` with `m ≤ (p+1)/2`, `p ∈ {5,7,11,13}`, plus 96 random `A` at
`17 ≤ p ≤ 43`: 0 failures (counts `L11a`, `L11bc`, `L11d_step1` in the JSON).

## 2. Theorem 2.1 (Hankel-minor Stepanov inequality) — CORRECT

Put `u_s = F^{(s)}/(D)_s` (`(D)_s` is a product of integers in `[1, p−1]`, a unit). For `s ≥ 1`,
`u_s = Σ_k c_k (x+a_k)^{D−s}`, and `u_0 = F`. `H := H_{e+1} = det[u_{i+j}]_{0≤i,j≤e}`,
`0 ≤ e ≤ (m−1)/2`. Crucially `H` does **not** depend on `b`; only the decomposition below does.

**Step 1.** Fix `b` with `e_b ≤ e`; `E = E(b)`, `δ = δ_b`. By Lemma 1.1(d),
`F = 2R_E + (x−b)^{m−δ}Q` for some `Q ∈ F_p[x]`. Differentiate `s ≤ 2e ≤ m−1` times and divide by
`(D)_s`: `u_s = 2ρ_s + ε_s` with `ρ_s = Σ_{k∈E} c_k (x+a_k)^{D−s}` (because
`R_E^{(s)} = (D)_s ρ_s`) and `ε_s = ((x−b)^{m−δ}Q)^{(s)}/(D)_s`. Leibniz: each term of the
`s`-th derivative contains `((x−b)^{m−δ})^{(t)} = (m−δ)_t (x−b)^{m−δ−t}` with `t ≤ s ≤ m−δ`
(`s ≤ m−1 ≤ m−δ`), so `(x−b)^{m−δ−s} | ε_s`. ✓

**Step 2.** Since `D − 2e ≥ d ≥ 0`, for `0 ≤ i,j ≤ e`:
`ρ_{i+j} = Σ_{k∈E} [c_k (x+a_k)^{D−2e}]·(x+a_k)^{e−i}·(x+a_k)^{e−j}`, i.e.
`[ρ_{i+j}] = Uᵀ·diag(c_k (x+a_k)^{D−2e})_{k∈E}·U` with `U_{k,j} = (x+a_k)^{e−j}` an `e_b × (e+1)`
matrix of **polynomials**. So `rank_{F_p(x)}[ρ_{i+j}] ≤ e_b` (a sum of `e_b` rank-one matrices;
`= 0` if `e_b = 0`). No division by `x + a_k` is needed. ✓

**Step 3.** Row `i` of `[u_{i+j}]` is `2·(ρ-row i) + (ε-row i)`. The determinant is multilinear in
rows, so `H = Σ_{S⊆{0,…,e}} det M_S`, where row `i` of `M_S` is the ε-row if `i ∈ S`, the 2ρ-row
otherwise. If `e+1−|S| > e_b`, the 2ρ-rows of `M_S` are more than `e_b` rows of a matrix of rank
`≤ e_b` over `F_p(x)`, hence linearly dependent over `F_p(x)`, so `det M_S = 0` in `F_p(x)`, hence as
a polynomial. Otherwise every Leibniz term of `det M_S` contains, for each `i ∈ S`, one entry
`ε_{i+σ(i)}` with `ord_b ≥ m − δ − i − σ(i) ≥ m − δ − i − e` (as `σ(i) ≤ e`), and polynomial
entries elsewhere. Each such number is `≥ m − 1 − 2e ≥ 0`, so
`ord_b det M_S ≥ Σ_{i∈S}(m − δ − i − e)`. Because every summand is nonnegative, the minimum over
admissible `S` (`|S| ≥ e+1−e_b`) is attained by the `e+1−e_b` **largest** indices,
`S = {e_b,…,e}`, and `ord_b H ≥ min_S ord_b det M_S`:

    ord_b H ≥ Σ_{i=e_b}^{e} (m − e − i − δ_b) = (e+1−e_b)·(m − (3e+e_b)/2 − δ_b) ≥ 0.

`δ_b` is handled uniformly (it only lowers the order in Step 1 by one). ✓

**Step 4.** For `s ≤ 2e ≤ m−1 ≤ d`, the coefficient of `x^{D−s−l}` in `u_s` is
`C(D−s,l)Σ_k c_k a_k^l`, zero for `l ≤ m−2`, `C(D−s, m−1)` for `l = m−1` (the `−1` in `u_0` sits in
degree `0 ≤ d`). So `deg u_s ≤ d − s` and the coefficient of `x^{d−s}` is `C(D−s, m−1)` read mod
`p` (possibly zero a priori). Every Leibniz term `Π_i u_{i+σ(i)}` has degree
`≤ Σ_i (d − i − σ(i)) = (e+1)(d−e)`, and the coefficient of `x^{(e+1)(d−e)}` in `H` is
`Λ := det[C(D−i−j, m−1)]_{0≤i,j≤e}` mod `p`. With Lemma 2.2 (`p ∤ Λ`): `deg H = (e+1)(d−e)` exactly,
in particular `H ≠ 0`. (Check: `C(D,m−1)(d)_s/(D)_s = C(D−s,m−1)` as rationals, so the worker's
normalisation agrees.) ✓

**Step 5.** For `H ≠ 0`, `Σ_{b∈F_p} ord_b H ≤ deg H`. Summing Step 3 over `b` with `e_b ≤ e`:

    (★)  Σ_{b: e_b ≤ e} (e+1−e_b)(m − (3e+e_b)/2 − δ_b) ≤ (e+1)(d−e).   ∎

`e = 0` gives `Σ_{e_b=0}(m − δ_b) ≤ d`, i.e. Hanson–Petridis Theorem 1.2 at `d=(p−1)/2`.

**Sharpness observation (PROVED by the hand computation below; JSON `star_equality_e1`).** `(★)` holds with **equality** at `e = 1`
for, e.g., `p = 13`, `A = {0,2,3,5}`: the `b` with `e_b ≤ 1` are `b = 1,7,9,12` (`e_b = 1`,
`δ = 0`, weight 2 each) and `b = 10, 11` (`e_b = 1`, `δ = 1`, weight 1 each), total
`10 = 2·(6−1)` (hand check; also exhaustive search: equality at `e = 1` occurs for `p = 5,7,11,13`;
maximal `e ≥ 1` ratios `13/14` at `p = 17`, `15/16` at `p = 19`, `19/20` at `p = 23`).
The robust note's sentence "at `e ≥ 1` the maximum observed is `≈ 0.92`" is therefore inaccurate
(harmless); the pass summary's "tight at `e = 1` for `p ≤ 13`" is right.

Numerics (Part A): the full polynomial `H_{e+1}` was computed over `F_p` (memoised Laplace
expansion of the polynomial matrix, with `u_s` built from power sums, **not** from the claimed
formulas) for 2509 triples `(p, A, e)`, `11 ≤ p ≤ 293`, `e ≤ 4` (`p ≤ 100`), `≤ 3` (`p ≤ 200`),
`≤ 2` (`p ≤ 293`), `A` = intervals, multiplicative subgroups and shifted cosets, `Q`, `Q∪{0}`
(`m = (p+1)/2`, the extreme `D = p−1`), greedy Paley cliques, greedy "sum-cliques"
(`A+A ⊆ Q∪{0}`), common neighbourhoods `N(0)∩N(1)`, `N(0)`, and random sets. For each: the
top coefficients of `u_s` vanish (Fact 0), `lc(u_s) = C(D−s,m−1)`, `deg H = (e+1)(d−e)` and
`lc(H) = Λ` (Λ computed separately by elimination), and the **exact** multiplicity of every
`b ∈ F_p` (via Hasse derivatives, valid for `deg H ≥ p`) is `≥` the Step-3 bound — 25,122
root checks, 0 failures, minimum slack `0` (the bound is attained at 14,190 points with `e ≥ 1`).
Step 1 (`ord_b ε_s ≥ m−δ−s`) was checked directly on `ε_s = u_s − 2ρ_s` at every `b`, every
`s ≤ m−1` (Part A-local). `(★)` was also checked for **every** `A ∋ 0`, `p ≤ 23`
(`2^22` sets at `p = 23`; translation invariance of `(★)` makes `A ∋ 0` exhaustive), every
admissible `e` (Part C).

## 3. Lemma 2.2 (leading coefficient) — CORRECT

Write `K = m−1`. Vandermonde convolution with `x = e−i ≥ 0`, `y = D−e−j ≥ 0`:
`C(D−i−j, K) = Σ_{l=0}^{e} C(e−i,l)·C(D−e−j, K−l)`, i.e. `[C(D−i−j,K)] = U·W`, `U_{il} = C(e−i,l)`
(zero when `i + l > e`, equal to `1` on `i + l = e`), so `det U = sgn(i ↦ e−i) = (−1)^{e(e+1)/2}`.
Again `C(D−e−j, K−l) = Σ_t C(e−j,t) C(D−2e, K−l−t)`, so `W = Z·Uᵀ` with
`Z_{lt} = C(N', K−l−t)`, `N' = D−2e`. Hence `Λ = det U · det Z · det U = det Z`. Reversing the
columns of `Z` (`t ↦ e−t`, sign `(−1)^{e(e+1)/2}`) gives `[C(N', c − l + t)]_{0≤l,t≤e}`,
`c = K − e = m−1−e ≥ e` (since `2e ≤ m−1`). This is Krattenthaler's (3.12) at `q = 1` with
`A = N'`, `L_i = c − i`, `n = e+1`. I checked the statement of (3.12) in the local text
`sources/stepanov-robust-krattenthaler-math9902004.txt` (Theorem 26, lines 2177–2221): at `q = 1`,
`det_{1≤i,j≤n} C(A, L_i+j) = Π_{i<j}(L_i−L_j)·Π_i (A+i−1)! / (Π_i (L_i+n)! · Π_i (A−L_i−1)!)` —
the worker's transcription is exact. `Π_{i<j}(j−i) = Π_{k=1}^{n−1} k!`.

**Nonvanishing mod p.** All factorials in the product formula have arguments in `[0, D−e]`:
numerator `k! (k ≤ e)`, `(N'+i−1)! ≤ (D−e)!`; denominator `(c−i+n)! ≤ (m−1)!`,
`(N'−c+i−1)! ≤ d!` (and `N'−c = d−e ≥ 0`). Since `D−e ≤ p−1`, both sides of
`Λ·Π(den) = ±Π(num)` are products of units mod `p`, and `Λ ∈ ℤ`; so `p ∤ Λ`. ✓

Numerics (Part B): (i) `Λ ≢ 0 (mod p)` by batched modular Gaussian elimination, directly from the
binomial matrix (no product formula), for **every** prime `p ≤ 1500`, every `m ≤ (p+1)/2` and every
`e ≤ min((m−1)/2, 12)`: 1,037,922 determinants, 0 zero. (ii) The closed formula equals the integer
determinant (fraction-free Bareiss) for every `D ≤ 60`, `1 ≤ m ≤ D+1`, `0 ≤ e ≤ (m−1)/2` with
`d = D−m+1 ≥ e`: 13,881 cases, 0 mismatches.

## 4. Corollary 2.3 — CORRECT

Restrict `(★)` to `b ∈ B` (dropping nonnegative terms). For `e_b ≤ e`: `(3e+e_b)/2 ≤ 2e`, so each
term is `≥ (e+1−e_b)(m−2e) − (e+1−e_b)δ_b ≥ (e+1−e_b)(m−2e) − (e+1)δ_b`. For `b ∈ B` with `e_b > e`
the quantity `(e+1−e_b)(m−2e)` is `< 0` (`m − 2e ≥ 1`), so it may be added to the left. Summing:
`(m−2e)Σ_{b∈B}(e+1−e_b) ≤ (e+1)(d−e) + (e+1)r`, i.e.
`(m−2e)((e+1)n − N_−) ≤ (e+1)(d−e+r)`; with `S = mn − r − 2N_−` this is the stated bias bound. ✓
(The sign of the `e_b > e` terms is exactly what lets the inequality be extended to all of `B`.)

Numerics: for every `A` and `e`, the worst `B` is `{b : w_b > 0}` with
`w_b = (m−2e)(e+1−e_b) − (e+1)δ_b`, so the check `Σ_b max(w_b,0) ≤ (e+1)(d−e)` covers **all** `B`
at once; done for all `A ∋ 0`, `p ≤ 23` (14.2M `(A,e)` pairs): 0 failures; the ratio reaches `1`.

## 5. Theorem 2.4 (robust HP) — CORRECT

First display: `(★)` restricted to `B`; for `b ∈ B`, `e+1−e_b ≥ e+1−e' > 0` and
`m − (3e+e_b)/2 − δ_b ≥ m − (3e+e')/2 − δ_b ≥ m − 2e − 1 ≥ 0`, so each term is
`≥ (e+1−e')(m − (3e+e')/2 − δ_b)`; sum. Then `(m−2e)n − r ≤ (e+1)d/(e+1−e')` and the identity
`mn − r = (m/(m−2e))((m−2e)n − r) + r·2e/(m−2e)` give the second display. η-form: if `e' = 0` use
HP. Else `m = e'/η ≥ 1/η ≥ 8`; `x = √(η/2) ≤ 1/4`; `e = ⌈xm⌉ − 1` has `xm ≤ e+1`, `e < xm ≤ m/4`,
so `2e ≤ m−1`; `η ≤ 1/8 ⇔ η ≤ x/2`, whence `m(x−η) ≥ mx/2 ≥ 1/(2√(2η)) ≥ 1` and `e ≥ xm−1 ≥ ηm = e'`.
Then `(e+1)/(e+1−e') ≤ 1/(1−η/x)`, `η/x = √(2η)`; `m/(m−2e) < 1/(1−2x) = 1/(1−√(2η))`;
`2er/(m−2e) < m·2x/(1−2x) ≤ m` (`r ≤ m`, `x ≤ 1/4`). ✓

Numerics: first display for all `e' ≤ e` and the η-form (with the tightest `η = max_{b∈B} e_b/m`)
on level sets `B = {e_b ≤ e'}` (the worst `B`), all `A ∋ 0`, `p ≤ 23`: 0 failures.

## 6. Theorem 2.5 (constant cancellation above `p/2`) — CORRECT

**Core lemma (PROVED; slightly more general than used).** If `1 ≤ m ≤ (p+1)/2`,
`mn ≥ (1/2+κ)p`, `0 < κ ≤ 3/2`, `u = (1+2κ)^{−1/2}`, then
`S(A,B) ≤ [1 − (1−u)² + u(1−u)m/d]·mn`.
*Proof.* `(1/2+κ)p > (1+2κ)d`, so `mn > d/u²`. `u ∈ [1/2,1)`, `θ = (1−u)/2 ∈ (0,1/4]`,
`e = ⌊θm⌋`: `e ≤ m/4 ≤ (m−1)/2` for `m ≥ 2`, `e = 0` for `m = 1`; `e+1 > θm`; `m−2e ≥ m(1−2θ) = um`.
Cor 2.3 with `d−e+r ≤ d+m` and `m−2e ≥ um`: `N_− ≥ (e+1)Y`, `Y := n − (d+m)/(um)`.
If `Y < 0` then `umn < d+m`, and with `umn > d/u` this gives `(1−u)/u < m/d`, so the bracket exceeds
`1 − (1−u)² + (1−u)² = 1` and the claim is trivial (`S ≤ mn`). If `Y ≥ 0`, `N_− ≥ θmY` and
`S ≤ mn − 2N_− ≤ mn[1 − (1−u)(1 − (d+m)/(umn))]`; `(d+m)/(umn) < u(1+m/d)` (from `mn > d/u²`) and
`1−u ≥ 0` give `S ≤ mn[1 − (1−u)² + u(1−u)m/d]`. ∎

**Reduction.** `S` is symmetric; WLOG `m = |A| ≤ n = |B|`, so `n ≥ √(mn) ≥ √((1/2+κ)p)`.
`m₀ := ⌈(1/2+κ)p/n⌉ ≤ m` (as `m ≥ (1/2+κ)p/n` and `m ∈ ℤ`), `m₀ ≤ √((1/2+κ)p)+1 ≤ √(2p)+1 ≤ (p+1)/2`
for `p ≥ 11` (at `p = 11`: `5.69 ≤ 6`; false at `p = 7`), and `m₀n ≥ (1/2+κ)p`.
*Sub-sampling identity (PROVED):* with `s(a) = Σ_{b∈B}χ(a+b)`, averaging over all `m₀`-subsets
`A'` of `A`, `E[S(A',B)] = Σ_a s(a)·Pr[a∈A'] = (m₀/m)S(A,B)`. Hence the core lemma applied to every
`(A',B)` gives `S(A,B) ≤ [1 − (1−u)² + u(1−u)m₀/d]mn`. *Dilation:* for a non-residue `ν`,
`S(νA,νB) = −S(A,B)` with the same sizes, so the same bound holds for `−S`. *Constants:*
`u(1−u) ≤ 1/4` and `m₀/d ≤ 2(√((1/2+κ)p)+1)/(p−1)` give the stated error term
`(√((1/2+κ)p)+1)/(2(p−1))`. `κ > 3/2`: the hypothesis with `κ = 3/2` holds. Expansion
`(1−u)² = κ² − 3κ³ + O(κ⁴)` checked. ∎

**Remarks.** (1) The bound is vacuous (factor `≥ 1`) for `p ≤ 13` and, for small `κ`, until `p`
is large: e.g. at `κ = 1/2` the factor is `< 1` only for `p ≥ 47`, at `κ = 0.1` only for
`p ≥ 2741`; at `κ = 3/2` from `p = 17`. It is an asymptotic statement with explicit constants; that is fine.
(2) The sharpness in Remark 2.7 uses `|A| = 2` (for `p ≡ 1 mod 4`, `A = {0,1}`,
`B = {b : χ(b),χ(b+1) ≠ −1}`, `|B| = (p+3)/4`; for `p ≡ 3 mod 4` the same set has `(p+1)/4`
elements). Whether `1/2` is still the right threshold when **both** `|A|, |B| → ∞` is OPEN (the
Paley conjecture would imply it is not; HEURISTIC). The upper cap `η(κ) ≤ 2κ/(1+2κ)+o(1)` is correct for
`κ ≤ 1` by the same `|A| = 2` family padded with `b` having `e_b = 1`.

Numerics (Part C): for every `A` (all subsets at `p ≤ 13`; all `A ∋ 0` at `p = 17, 19, 23`) and
every `n`, the extreme values of `S(A,B)` over `|B| = n` are the sums of the `n` largest / smallest
row sums `s_A(b)`; the final bound was tested at every admissible `κ` on a 400-point grid plus the
endpoint `K = min(3/2, mn/p − 1/2)`, and the core lemma at its exact minimiser `κ = K` (the core
bracket is concave in `u`, so its minimum over `κ ∈ (0,K]` is at `κ = K` or is `1`). 0 failures in
1.0·10⁸ `(A, n)` cases for the final bound and 5.8·10⁷ for the core lemma; non-vacuous cases (factor `< 1`) exist from `p = 17` on (≈9·10⁷ at
`p = 23`), minimum absolute slack `2.20` at `p = 23`. Local search over `A` (sizes `2…√(2p)+1`)
at `p = 101, 211, 409, 1009, 2003` maximising `max_n [S_max(A,n)/(mn) − bound]`: every gap
negative (JSON `C_adversarial`); the closest cases are the `|A| = 2` family. No instance had bias
above even the error-free `1 − (1−u)²` at `p ≤ 23`.

## 7. Corollary E (Paley graph) — CORRECT (with an explicit error term)

`p ≡ 1 (mod 4)`, `|A| = m`, `m² ≥ (1/2+κ)p`, `B = −A` (`r = m`). `χ(−1) = 1`, so
`S(A,−A) = Σ_{a≠a'} χ(a−a') = 4E(A) − m(m−1)` and, with `ρ = E(A)/C(m,2)`,
`S = m(m−1)(2ρ−1)`. Theorem 2.5: `|2ρ − 1| ≤ Φ·m/(m−1)`, `Φ = 1 − c(κ) + ε_p`,
`ε_p = (√((1/2+κ)p)+1)/(2(p−1))`. Hence

    c(κ)/2 − ε_p/2 − Φ/(2(m−1)) ≤ ρ ≤ 1 − c(κ)/2 + ε_p/2 + Φ/(2(m−1)),

and `ε_p, 1/(m−1) = O(p^{−1/2})` since `m ≥ √(p/2)`. This is the summary's statement with the
`o(1)` made explicit. ✓ Numerics (Part D): all 38 primes `p ≡ 1 (mod 4)`, `13 ≤ p ≤ 409`, plus
`p = 1009, 2017`, `κ ∈ {0.05, 0.25, 0.5, 1, 1.5}`, `m = ⌈√((1/2+κ)p)⌉`, tabu local search for the
densest `m`-set (min density is `1 −` max by self-complementarity); exact maxima for `p = 13, 17`
(`m ≤ 7`) and `p = 29` (`m ≤ 6`). All satisfy the explicit bound; the densest sets found are far
inside it (e.g. `p = 2017`, `κ = 1.5`, `m = 64`: density `0.729` vs proven `≤ 0.889`;
`κ = 0.05`, `m = 34`: `0.840` vs proven `1.017`, vacuous). The comparison claims in the summary
("HP gave only 'not a clique'", "Chung needs `|A| ≥ (1+κ)√p`") are accurate: from HP alone one
gets only `#non-edges ≥ m − √(p/2) − O(1)`, i.e. density `1 − O(1/m)`.

## 8. Required fixes (presentation only; none affects the mathematics of the robust note)

1. **Summary, Theorem B** ("For all `A, B`"): add `|A| ≤ (p+1)/2`. Cor 2.3 is proved only under
   `D ≤ p−1`. (Numerically no violation was found for `|A| > (p+1)/2` at `p ≤ 23` — counts
   `C_hypfree_cor23` — but that range is unproved.)
2. **Summary, Theorem C** omits `|A| ≤ (p+1)/2`, but the statement is nevertheless true
   (PROVED here). Suppose `m > (p+1)/2` and `n ≥ 2`; take `b₁ ≠ b₂` in `B`, `c = b₂ − b₁`. Each
   `b_i` has `≥ 7m/8` partners `a` with `χ(a+b_i) ≠ −1`, so `≥ 3m/4` elements `a` have both.
   With `y = a + b₁`: `#{y : χ(y) ≠ −1, χ(y+c) ≠ −1} = (p−3−χ(c)−χ(−c))/4 + [χ(c)=1] + [χ(−c)=1]
   ≤ (p+3)/4` (Jacobsthal: `Σ_{y≠0,−c}(1+χ(y))(1+χ(y+c)) = p−3−χ(c)−χ(−c)`). So
   `m ≤ (p+3)/3 ≤ (p+1)/2` for `p ≥ 3`, a contradiction. Hence `n ≤ 1` and
   `mn − r ≤ m ≤ f(η)d + m`. (Exhaustive check at `p ≤ 23` agrees: count `C_hypfree_thmC`.)
   Still, the summary should state the hypothesis or cite this remark.
3. **Robust note, end of §2 Lemma 2.2 verifier paragraph**: "at `e ≥ 1` the maximum observed is
   `≈ 0.92`" is wrong; `(★)` is attained with equality at `e = 1` (witness in §2).
4. **Summary, Theorem D paragraph**: "`A = {0,1}` … `|A||B| = (p+3)/2`" needs `p ≡ 1 (mod 4)`
   (for `p ≡ 3 (mod 4)` one gets `(p+1)/2`); and "the threshold `1/2` is sharp" should say "sharp
   for the statement over all set sizes (witness `|A| = 2`)".
5. Optional: state in Theorem 2.5 that it is vacuous for small `p` (Remark (1) of §6).

## 9. Verification summary

| part | what | cases | failures |
|---|---|---|---|
| A-local | Lemma 1.1(a)–(c) at every `b, j`; Step 1 at every `b, s` | ≈5·10⁵ | 0 |
| A | exact `H_{e+1}`, degree, `lc = Λ`, exact root orders vs bound | 2509 polys / 25,122 roots | 0 |
| B | `Λ ≢ 0 mod p`, `p ≤ 1500`, `e ≤ 12` | 1,037,922 | 0 |
| B | product formula over ℤ, `D ≤ 60` | 13,881 | 0 |
| C | `(★)`, Cor 2.3 (all `B`), Thm 2.4, core lemma, Thm 2.5, exhaustive `p ≤ 23` | 2.38·10⁸ | 0 |
| C | hypothesis-free variants (`|A| > (p+1)/2`) of summary Thms B, C, `p ≤ 23` | 1.78·10⁷ | 0 violations |
| C | adversarial local search `p ≤ 2003` | ≈2.7·10⁴ evals | 0 |
| D | Corollary E, 40 primes, exact for `p ≤ 29` | 200 + 1.4·10⁵ sets | 0 |

Exact counts are in `results/stepanov_referee2_2026_09_27.json` (`checks`, `notes`).

## 10. Literature (searched 2026-09-27)

- Hanson–Petridis, arXiv:1905.09134 (local `sources/sigma-hanson-petridis-1905.09134.txt`,
  lines 81–95): after Vinogradov's bound `|S_χ(A,B)| ≤ √(p|A||B|)` they write that it
  "is still the best known estimate for Sχ(A,B)", nontrivial only for `|A||B| > p`. So as of
  2019–21 no constant saving below `p` was known to them. Their Theorem 1.2 is exact containment.
- Shkredov, *Sumsets in quadratic residues*, arXiv:1305.4093 (Acta Arith. 164 (2014)), fetched:
  Theorem 3.2 treats **approximate** equality `A+A ≈ R` (lower bound `(1/6−o(1))|A|` on
  `max{|R∖(A+A)|, Σ_{x∈(A+A)∖R}(A∗A)(x)}`) by Cauchy–Schwarz/second moment; it is a sumset
  statement at `|A|² ≈ p`, not a bias bound for arbitrary pairs in `p/2 < |A||B| < p`. Its
  inequality (24) (`A+B ⊆ R ⇒ p + ab/p ≥ ab + a + b`) is second-moment.
- Heath-Brown–Konyagin, Q. J. Math. 51 (2000) (ORA preprint fetched): Lemma 5 bounds
  `#E ≪ (hT)^{2/3}` for `E = ∪_{u∈U}{y : uy ∈ μ_h, uy − u ∈ μ_h}` (`h⁴T < p³`): Stepanov for several
  shifts simultaneously, exact membership in a subgroup; no exceptions, no arbitrary sets.
- Kalmynin arXiv:2504.10202; Rudnev–Tyrrell arXiv:2607.24270; Yip–Yoo arXiv:2608.02568 (v2,
  9 Sep 2026, "global differential identities"); Cochrane arXiv:2607.28559; Kim–Yip–Yoo
  arXiv:2602.20919, 2607.24370: abstracts (and local texts for the first two) concern exact
  decompositions/containments only. Local grep of Kalmynin and Rudnev–Tyrrell for
  "robust/exception/stability": nothing relevant.
- Volostnov arXiv:1712.09355 (local): for `|A| ∼ |B| ∼ p^{1/2}` "inequality (2) is unknown" (power
  saving); results need small doubling. Schoen–Shkredov arXiv:2004.01885 (fetched): character-sum
  bounds near `√p` need bounded doubling.
- Yip, PhD thesis *Topics in arithmetic combinatorics* (UBC, 2024): UBC server blocked the fetch
  (**UNVERIFIED** content); per its catalogue record it consists of four articles, and Yip's
  publication list (checked by the worker) has no robust HP.
- Searches for robust/defect/Hankel/syndrome/Wronskian versions of Stepanov or HP and for
  Paley-graph induced-subgraph density at `|A| ≍ √p` found nothing. The Randomstrasse101 Paley
  problems (25–29) contain no density statement.

Conclusion: I did **not** find Theorem 2.5, `(★)`, or any constant-bias bound for arbitrary sets
with `p/2 < |A||B| < p` in the literature. This is a negative search, not a proof of novelty.

## 11. Verdicts

| result | verdict |
|---|---|
| Lemma 1.1 | CORRECT (§1) |
| Thm 2.1, Steps 1–5 | CORRECT (§2) |
| Lemma 2.2 | CORRECT (§3) |
| Cor 2.3 | CORRECT (§4); summary's Thm B needs `|A| ≤ (p+1)/2` |
| Thm 2.4 | CORRECT (§5); summary's Thm C true even without `|A| ≤ (p+1)/2` (§8.2) |
| Thm 2.5 | CORRECT (§6) |
| Corollary E | CORRECT (§7), explicit error term given |

I believe Theorem 2.5 is true (I have a complete proof above, and it survives every exact test).

## 12. Open obligations

1. Threshold for `|A|, |B| → ∞` (is any constant saving available below `p/2` for large sets?).
2. `κ² ≲ η(κ) ≲ 2κ`.
3. A human specialist's reading, and a Lean formalization of Steps 1–5 (the argument is short
   and entirely algebraic; the only external input is Krattenthaler (3.12), which could be replaced
   by the direct modular check only for bounded `p`).
