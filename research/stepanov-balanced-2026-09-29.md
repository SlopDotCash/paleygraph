# Stepanov wave, worker `balanced`: near-tight Hanson–Petridis bicliques (2026-09-29)

**Status:** PROVED: (i) an exact reformulation of complete bicliques through the
Hanson–Petridis polynomial (Theorem 1.1): `A + B ⊆ Q ∪ {0}` iff
`P₁^m P₀^{m−1} | F_A`, and then `F_A = λ P₁^m P₀^{m−1} G_A` with **exact** multiplicities
`m`, `m−1` at `B`, `gcd(G_A, P_B) = 1`, `deg G_A = g := d − mn + r`, the same defect `g`
for `F_B`, and every root of `F_A` of multiplicity `≤ m`; (ii) residue identities fixing
`G_A` on `B` (Lemma 1.2), hence `χ(G_A(b)) = χ(2m)` on `B ∖ (−A)` for even `m`; (iii) the
forced factorisations of all derivatives, Hankel minors and the pencil map (Proposition 1.3:
the pencil map has degree exactly the number of distinct roots of `F_A`, `≤ n + g`); (iv) a
"derivative Hanson–Petridis" inequality `(m−1)(N₀ + N_m) ≤ d − 1 + r₀ + r_m` for complete and
anti-complete points together, sharp for `A = {0,1}` (Proposition 1.4); (v) precise
obstructions for the four proof routes of the brief (§3): every Wronskian/abc
(Brownawell–Masser–Voloch type) inequality for `λP^mG + 1 = Σ_k c_k(x+a_k)^D` reduces to
"multiplicities `≤ m`" and is blind to `B` (Proposition 3.3); two-sided multiplicity
Schwartz–Zippel with the available polynomials returns at most Hanson–Petridis (§3.1); and an
exact rational linear-programming certificate (Theorem 3.5) shows that **no nonnegative
combination** of the known count constraints — the Hankel-minor inequalities (★) for `A` and
for `νA`, every averaged subset Hanson–Petridis inequality, the pencil inequality, the
derivative inequality of (iv), the exact first and second moment identities, the Weil moment
bounds of orders 4, 6, 8 and the Weil counts — excludes a Hanson–Petridis-tight balanced
biclique profile (`m = n ≈ √(p/2)`, `mn ≥ d − m`) at `p = 1009, 10009, 40009, 100049, 1000033`,
even one in which every point off `B` has between `0.4m` and `0.6m` bad partners; such a
profile forces `f_A = Σ_{a∈A}χ(·+a)` to put half its `L²` mass on `B`, and excluding that for
`|A| ≍ √p` is a square-root-barrier statement.
REFUTED: "near-tightness forces extra zeros of `F_A` outside `B` or extra vanishing of Hankel
minors on `B`" (witnesses: the tight configurations `A = {0,1}` for every `p`, and
`p = 13, 37, 41`, where `G_A = 1`; exact Hankel orders in all tested cases); "the count
constraints above exclude balanced near-tight bicliques" (witness: the rational profiles of
Theorem 3.5). HEURISTIC (exact data, §2): `P_K(p) := max{|A||B| : min(|A|,|B|) ≥ K}`
satisfies `P_K(p)/p ≈ 0.41, 0.35, 0.29, 0.25` for `K = 3,4,5,6` at `p ≈ 500` (exact for all
`p ≡ 1 mod 4`, `p ≤ 509`) and keeps decreasing towards `K/2^K` (best values found at
`p ≈ 3000`: `0.389, 0.290, 0.22` for `K = 3, 4, 5`); the extremal pairs have smaller side
exactly `K` in most rows and at most 11 in all, and the largest `k × k` bicliques have
`k²/p ≤ 0.354` for `200 ≤ p ≤ 509` (`0.43` at `p = 113`) — the data never reach the balanced
window `m, n ≍ √p`. OPEN: the target
`|A||B| ≤ (1/2 − c_K)p + o(p)` for `min(|A|,|B|) ≥ K`; in particular any improvement of the
Hanson–Petridis constant for prime Paley cliques (none is known: §4).

Verifier: `experiments/stepanov_balanced_2026_09_29.py` → `results/stepanov_balanced_2026_09_29.json`
(run with `/opt/miniconda3/bin/python3`; about 32 s; 1,624,403 exact checks, 0 failures).
Data driver: `experiments/stepanov_balanced_2026_09_29_data.py` with the C helper
`experiments/stepanov_balanced_2026_09_29.c` → `results/stepanov_balanced_2026_09_29_data.json`
(27 min with 3 processes; 106 primes); the verifier re-checks the data file.

---

## 0. Notation

`p` odd prime (for cliques `p ≡ 1 mod 4`), `χ` the Legendre symbol, `χ(0) = 0`, `d = (p−1)/2`,
`Q` the nonzero squares, `ν` a fixed non-residue. `A = {a_1,…,a_m}` with `2 ≤ m ≤ (p+1)/2`,
`D = d + m − 1 ≤ p − 1`, `c_k = Π_{l≠k}(a_k − a_l)^{−1} = 1/P_A'(a_k)`,

    F_A(x) = −1 + Σ_k c_k (x + a_k)^D,     deg F_A = d,  leading coefficient λ := C(D, m−1)

(robust note, Lemma 1.1(a)). For `b ∈ F_p`: `y_k = b + a_k`, `E(b) = {k : χ(y_k) = −1}`,
`e_b = |E(b)|`. For a set `B`: `n = |B|`, `B₀ = B ∩ (−A)`, `B₁ = B ∖ (−A)`, `r = |B₀|`,
`P₁ = Π_{b∈B₁}(x−b)`, `P₀ = Π_{b∈B₀}(x−b)`, `P = P_B = P₀P₁`, and the **defect**

    g := d − mn + r.

Hanson–Petridis (CITED: B. Hanson, G. Petridis, *Refined estimates concerning sumsets contained
in the roots of unity*, Proc. LMS 122 (2021) 353–358, arXiv:1905.09134, Theorem 1.2; local text
`sources/sigma-hanson-petridis-1905.09134.txt`) is `g ≥ 0` for complete bicliques. The target
of this note is `g ≥ (2c_K − o(1))d` when `min(m,n) ≥ K`.

`(D)_j = D(D−1)⋯(D−j+1)`. Every integer in `[1, p−1]` is a unit mod `p`; since `D ≤ p−1`, all
`(D)_j`, `j ≤ D`, and all factorials below are units. This is the only place the prime field
enters, exactly as in Hanson–Petridis.

## 1. The exact reformulation (task 1)

### 1.1 Statement

**Theorem 1.1 (PROVED).** Let `2 ≤ m ≤ (p+1)/2`.

1. For `b ∉ −A`: `ord_b F_A ≥ m` iff `b + A ⊆ Q`. In that case `ord_b F_A = m` exactly and
   `F_A^{(m)}(b) = (−1)^{m−1}(D)_m / Π_k (b + a_k)`.
2. For `b = −a_{k₀}`: `ord_b F_A ≥ m−1` iff `b + A ⊆ Q ∪ {0}`. In that case `ord_b F_A = m−1`
   exactly and `F_A^{(m−1)}(b) = −(D)_{m−1} c_{k₀}`. If `b + A ⊄ Q ∪ {0}` then `ord_b F_A ≤ m−2`.
3. Hence for every `B ⊆ F_p`: `A + B ⊆ Q ∪ {0}` iff `P₁^m P₀^{m−1}` divides `F_A`, and then

       F_A = λ · P₁^m · P₀^{m−1} · G_A,    G_A monic,  deg G_A = g,  gcd(G_A, P_B) = 1.

4. If moreover `n ≥ 2`, then `n ≤ (p+1)/2`, `|A ∩ (−B)| = r`, and
   `F_B = λ_B · Π_{a∈A∖(−B)}(x−a)^n · Π_{a∈A∩(−B)}(x−a)^{n−1} · G_B` with `deg G_B = g`:
   **both Hanson–Petridis polynomials of a complete biclique have the same defect.**
5. Every root of `F_A` in `F̄_p` has multiplicity `≤ m` (`≤ m − 1` on `−A`), and at a
   rational point `x ∉ −A ∪ B(A)` the multiplicity is at most `min(e_x − 1, m − e_x)`; in
   particular points with exactly one bad partner are never roots.

*Proof.* Robust note Lemma 1.1(b),(c) (re-proved in the referee note): for `b ∉ −A` and
`0 ≤ j ≤ m−1`, `F^{(j)}(b) = −2(D)_j Σ_{k∈E(b)} c_k y_k^{m−1−j}`; for `b = −a_{k₀}` the same holds
for `j ≤ m−2` (with `k₀ ∉ E(b)`), and `F^{(m−1)}(b) = (D)_{m−1}(−c_{k₀} − 2Σ_{E(b)} c_k)`.

(1) If `E(b) = ∅` all `F^{(j)}(b)`, `j < m`, vanish. Conversely, if they vanish then
`Σ_{k∈E} (c_k)·y_k^t = 0` for `t = 0,…,m−1`; the `y_k` (`k ∈ E`) are distinct and `|E| ≤ m`, so the
Vandermonde matrix `(y_k^t)` has full column rank and all `c_k = 0` for `k ∈ E`, impossible
unless `E = ∅`. For the exact order: `F^{(m)}(b) = (D)_m Σ_k c_k y_k^{D−m}` and
`y^{D−m} = y^{d−1} = χ(y)/y`; with `E = ∅` this is `(D)_m Σ_k c_k/(b+a_k)`, and the partial
fraction identity `Σ_k c_k/(t + a_k) = (−1)^{m−1}/Π_k(t + a_k)` (both sides have simple poles at
`−a_k` with residues `c_k` and `(−1)^{m−1}/Π_{l≠k}(a_l − a_k) = c_k`, and vanish at `∞`) gives the
formula; `(D)_m` is a unit.

(2) Vanishing of `F^{(j)}(b)` for `j ≤ m−2` means `Σ_{k∈E}(c_k y_k)·y_k^{t−1} = 0` for
`t = 1,…,m−1`; `|E| ≤ m−1`, `c_k y_k ≠ 0`, the `y_k` distinct: `E = ∅`. Then
`F^{(m−1)}(b) = −(D)_{m−1}c_{k₀} ≠ 0`. If `b + A ⊄ Q ∪ {0}` some `F^{(j)}(b)`, `j ≤ m−2`, is
nonzero.

(3) The divisibility is (1)+(2). `deg F_A = d` with leading coefficient `λ`, a unit, so
`G_A := F_A/(λP₁^mP₀^{m−1})` is monic of degree `d − m(n−r) − (m−1)r = g`. The exact orders give
`G_A(b) ≠ 0` for `b ∈ B`.

(4) `mn ≤ d + r ≤ d + m` (Hanson–Petridis, i.e. `g ≥ 0`) gives
`n ≤ d/m + 1 ≤ d/2 + 1 ≤ (p+1)/2`. `x ↦ −x` maps `B ∩ (−A)` onto `A ∩ (−B)`. Apply (3) with the
roles exchanged: the defect is `d − nm + r = g`.

(5) If `β ∉ −A` and `F^{(j)}(β) = 0` for `j = 1,…,m`, then
`Σ_k [c_k(β+a_k)^{D−m}]·(β+a_k)^{m−j} = 0` for `m−j = 0,…,m−1`, a Vandermonde system in the
distinct `β + a_k ≠ 0`; so `c_k(β+a_k)^{D−m} = 0`, impossible. Hence `ord_β F ≤ m`. At `−A`,
(2) gives `ord ≤ m−1`. At rational `x ∉ −A` with `e_x = m`, `F(x) = −2 ≠ 0`. At rational
`x ∉ −A` with `1 ≤ e_x ≤ m−1` (so `μ := ord_x F ≤ m`): `ord_x F ≥ μ` forces
`Σ_{E}c_k y_k^t = 0` for `t = m−μ,…,m−1` (so `μ < e_x` by Vandermonde), and, writing
`Σ_E = [t=m−1] − Σ_{Ē}` for `t ≤ m−1` (Fact 0 of the robust note), `Σ_{Ē}c_k y_k^t = 0` for
`t = m−μ,…,m−2` (so `μ − 1 < m − e_x`). ∎

So near-tightness (`g` small) means exactly: **`F_A` is an `m`-th power up to the factor
`P₀^{−1}` and a cofactor of degree `g`, and it is simultaneously a sum of `m` pure `D`-th powers
minus 1.** The equivalence in (3) is exact: every constraint on complete bicliques is a
constraint on this identity, and conversely.

*Tight examples (checked).* `A = {0,1}`: `F_A = −1 − x^{d+1} + (x+1)^{d+1} = (d+1)·x(x+1)·
Π_{b∈B₁}(x−b)²` with `|B| = (p+3)/4`, `r = 2`, `g = 0`, for every `p ≡ 1 mod 4` (biclique note
Prop. 2.3). With `min(m,n) ≥ 3`: `p = 13, A = {0,1,4}, B = −A`; `p = 37, A = {0,1,11}`,
`|B| = 7`; `p = 41`, the 5-clique `{0,1,2,10,33}`, `B = −A`; all have `g = 0`, `G_A = 1`
(verifier §A; Kalmynin, arXiv:2504.10202v2, Lemma 4 and Theorem 4 classify the clique case).

### 1.2 The values of the defect polynomial on `B`

**Lemma 1.2 (PROVED).** Under Theorem 1.1(3):

- for `b ∈ B₁`: `P'(b)^m · G_A(b) · Π_{a∈A}(a+b) / P₀(b) = (−1)^m/(2m)`;
- for `b = −a_{k₀} ∈ B₀`: `P₁(b)^m · P₀'(b)^{m−1} · G_A(b) = −c_{k₀}`.

Consequently, if `m` is even, `χ(G_A(b)) = χ(2m)` for every `b ∈ B₁`; if `m` is odd,
`χ(G_A(b)) = χ(−2m)·χ(P_B'(b))`. If `r = 0` the first identity reads
`P_B'(b)^{−m} = 2m·G_A(b)·P_A(−b)`, and symmetrically `P_A'(a)^{−n} = 2n·G_B(a)·P_B(−a)`
for `a ∈ A`.

*Proof.* Near `b ∈ B₁`, `F_A = λ(x−b)^m P₁'(b)^m P₀(b)^{m−1}G_A(b)(1 + O(x−b))` and
`P₁'(b)P₀(b) = P'(b)`; compare `F_A^{(m)}(b)/m!` with Theorem 1.1(1) and use
`(D)_m/(m!·λ) = (D−m+1)/m = d/m ≡ −1/(2m)`. Near `b ∈ B₀` compare `F_A^{(m−1)}(b)/(m−1)!` with
Theorem 1.1(2) and use `λ = (D)_{m−1}/(m−1)!`. For the characters: `P'(b)^m` is a square for even
`m`, `a + b ∈ Q`, and `P₀(b) = Π_{a'∈A∩(−B)}(b + a') ∈ Q`. For `r = 0`,
`Π_a(a+b) = (−1)^m P_A(−b)`. ∎

**Remark 1.2′ (the reach of coefficient comparison; PROVED).** If `g ≤ n − 2`, Lemma 1.2
prescribes `G_A` on `n > g + 1` points, so the prescribed values satisfy the `n − 1 − g`
linear relations `Σ_{b∈B} G_A(b) b^j / P_B'(b) = 0` (`0 ≤ j ≤ n−2−g`), polynomial identities in
`A` and `B`; symmetrically `G_B` on `A` gives `m − 1 − g` relations if `g ≤ m − 2`. These are
the relations Kalmynin exploits when `g = 0` (his Lemmas 4–6 compare the top coefficients of
`F_A = C·Π(x−b)^{α−ε(b)}`). They exist only while `g ≤ max(m,n) − 2`, that is within `O(√p)`
of the Hanson–Petridis bound in the balanced window; for `g ≥ max(m,n) − 1` any values on `A` and
on `B` are interpolated by polynomials of degree `g`, and Lemma 1.2 imposes nothing beyond
`G_A(b) ≠ 0`, `G_B(a) ≠ 0`. **Coefficient comparison of this kind therefore cannot reach a
constant-factor improvement.** (Verifier: the relations hold in all 378 instances with
`g ≤ n − 2` among the exhaustive cases.)

### 1.3 Forced factors of derivatives, Hankel minors and the pencil

**Proposition 1.3 (PROVED).** Under Theorem 1.1(3), with `u_s = F_A^{(s)}/(D)_s` and
`H_{e+1} = det[u_{i+j}]_{0≤i,j≤e}` (robust note, §2):

1. `F_A^{(s)} = P₁^{m−s}P₀^{m−1−s}Q_s` with `deg Q_s = g + s(n−1)` for `0 ≤ s ≤ m−2`, and
   `F_A^{(m−1)} = P₁·Q_{m−1}`, `deg Q_{m−1} = d − m − n + r + 1`.
2. For `0 ≤ e ≤ (m−1)/2`: `H_{e+1} = P₁^{(e+1)(m−e)}P₀^{(e+1)(m−1−e)}K_e` with
   `deg K_e = (e+1)(g + e(n−1))`. At `b ∈ B₁` the order is exactly `(e+1)(m−e)` iff
   `Δ_e := det[(m)_{i+j}/(D)_{i+j}]_{0≤i,j≤e} ≢ 0 (mod p)`; `Δ₁ = m(m−D)/(D²(D−1)) ≢ 0`.
3. The pencil map `R = V/W`, `V = F_A − xF_A'/D`, `W = F_A'/D` (pencil note), satisfies
   `F_{A∖{a}} = W·(R − a)` and has degree exactly `n₀(F_A)`, the number of distinct roots of `F_A`
   in `F̄_p`; for a complete biclique `n₀(F_A) = n + n₀(G_A) ≤ n + g`.

*Proof.* (1) Leibniz: `(P₁^mP₀^{m−1}G)^{(s)}` is divisible by `P₁^{m−s}P₀^{m−1−s}`, and
`deg F^{(s)} = d − s` (leading coefficient `λ(d)_s`, a unit). (2) Each Leibniz term
`Π_i u_{i,σ(i)}` has order `≥ Σ_i(m − i − σ(i)) = (e+1)(m−e)` at `b ∈ B₁` and
`≥ (e+1)(m−1−e)` at `b ∈ B₀`; `deg H_{e+1} = (e+1)(d−e)` is robust Lemma 2.2. At `b ∈ B₁`, write
`F = φ(x−b)^m(1 + O(x−b))`, `φ ≠ 0`; then `u_s = φ·(m)_s/(D)_s·(x−b)^{m−s}(1 + O(x−b))`; removing
`(x−b)^{m−i}` from row `i` and `(x−b)^{−j}` from column `j` leaves `φ^{e+1}Δ_e` as leading
coefficient. `Δ₁` is a two-line computation, nonzero as `1 ≤ D − m = d − 1 < p` for `p ≥ 5`.
(3) `gcd(V, W) = gcd(F, F') =: Γ`. Put `F = ΓF₁`, `F' = ΓF₂`; since all multiplicities are
`< p`, `F₁` is the radical of `F` up to a constant (degree `n₀`) and `F₂` has degree `n₀ − 1` and
no common root with `F₁`. `R = (DF₁ − xF₂)/F₂` in lowest terms, and the numerator has leading
coefficient `(D − d)·lc(F₁) = (m−1)lc(F₁) ≠ 0`, so `deg R = n₀`. The identity for `F_{A∖{a}}`
is pencil note Lemma 2.1. ∎

**What this gives.** Every object built from `F_A` inherits the factor `P_B` to a power, but
never with better "efficiency" than `F_A` itself: `F_A` has degree `d` and order `m` on `B₁`
(ratio `d/m`), `H_{e+1}` has degree `(e+1)(d−e)` and order `(e+1)(m−e)` (ratio
`(d−e)/(m−e) > d/m`), `F^{(s)}` has ratio `(d−s)/(m−s)`. Root counting with any of them is
weaker than Hanson–Petridis. The pencil map has degree `n + n₀(G_A)` instead of the generic
`≈ d`: near-tightness makes `R` a map of degree `≈ n` (`≍ √p` in the balanced window), and
Riemann–Hurwitz for `R` gives `(m−2)·#{x ∉ −A : e_x = 1} ≤ 2(n + g) − 2`, which is the pencil
note's Theorem 3.1 again. All these are upper bounds on further special points; none bounds `n`.

### 1.4 Complete and anti-complete points together

**Proposition 1.4 (derivative Hanson–Petridis; PROVED).** Let `2 ≤ m ≤ (p+1)/2`,
`N₀ = #{x : x + A ⊆ Q ∪ {0}}`, `N_m = #{x : x + A ⊆ νQ ∪ {0}}`, and `r₀, r_m` the numbers of
those `x` lying in `−A`. Then

    (m − 1)(N₀ + N_m) ≤ d − 1 + r₀ + r_m,

with equality for `A = {0,1}` and every `p ≡ 1 mod 4` (`N₀ = (p+3)/4`, `N₂ = (p−1)/4`, `r₀ = 2`,
`r₂ = 0`).

*Proof.* `F_A'` has degree `d−1` and leading coefficient `dλ ≠ 0`. For `x ∉ −A` with all
`χ(x + a_k) = ε ∈ {±1}` and `1 ≤ j ≤ m−1`,
`F^{(j)}(x) = (D)_j Σ_k c_k χ(y_k)y_k^{m−1−j} = ε(D)_j Σ_k c_k y_k^{m−1−j} = 0` (Fact 0, as
`m−1−j ≤ m−2`), so `ord_x F' ≥ m−1`. For `x = −a_{k₀}` with all other `χ(x+a_k) = ε` and
`1 ≤ j ≤ m−2`: `F^{(j)}(x) = ε(D)_j(Σ_k c_k y_k^{m−1−j} − c_{k₀}·0^{m−1−j}) = 0`, so
`ord_x F' ≥ m−2`.
The two kinds of points are disjoint (`m ≥ 2`). Count roots of `F'`. For `A = {0,1}` the counts
are those of the biclique note, Prop. 2.3, and `(p+3)/4 + (p−1)/4 = d + 1 = d − 1 + 2`. ∎

For a complete biclique with `B = B(A)` this says `N_m ≤ (g + n − 1 + r_m)/(m−1)`:
near-tight pairs have almost no anti-complete points. It is again an upper bound.

## 2. Exact data (task 2)

**Definitions.** `P_K(p) := max{|A||B| : A + B ⊆ Q ∪ {0}, min(|A|,|B|) ≥ K}` and
`M_m(p) := max_{|A|=m}|B(A)|`. If the smaller side of an extremal pair has size `m`, then
`P_K = m·M_m`; hence `P_K(p) = max_{m ≥ K, M_m ≥ m} m·M_m(p)`.

**Method (exact where flagged).** `M_m(p)` for `m ≤ 8` (`p ≤ 769`) and `m ≤ 6` (`p ≤ 997`) are
the exact values of `results/sigma_biclique_2026_09_05_search.json`; our C helper
(`profile` mode) re-derived `M_3,…,M_6` for every `p ≤ 997` (all 80 rows agree with the table)
and computed `M_3, M_4, M_5` for 26 primes in `[1009, 2969]`. With `T_K = max_{K≤m≤m_k} m·M_m`
(attained, witness recorded), the helper's `cert` mode searches exhaustively (branch and bound
over normalised smaller sides, pruning with `M_m ≤ M_{m_k}`) for bicliques whose smaller side
has size `m_k + 1 ≤ m ≤ b_max` and product `> T_K`, where `b_max` is the exact balanced number
`b(p)` of the biclique note (`p ≤ 509`) or the Hanson–Petridis bound `⌊(1+√(2p−1))/2⌋`
otherwise. When the search finishes, `P_K(p)` is exact (PROVED for that `p`, by computation);
otherwise the entry is a lower bound `*` with witness. The verifier recomputes `|B(A)|` for all
383 witnesses and `M_3` in pure Python for `p ≤ 400`. For `p ≤ 509` every entry is exact
(`K = 3,…,6`); in nine cases the search found a larger product with smaller side `≥ 9`
(e.g. `p = 257`, `K = 6`: `9 × 9`). Between 510 and 999 the searches mostly timed out (30 s;
`K = 5, 6` not attempted), so those entries are the exact maxima over `|A| ≤ 8` (`p ≤ 769`) or
`|A| ≤ 6`, i.e. lower bounds for `P_K`.

**Selected rows** (`P_K/p`, extremal shape `|A| × |B|`; `*` = lower bound):

| p | P_3/p | P_4/p | P_5/p | P_6/p |
|---|---|---|---|---|
| 13 | 0.692 (3×3) | — | — | — |
| 37 | 0.568 (7×3) | 0.432 (4×4) | — | — |
| 41 | 0.610 (5×5) | 0.610 (5×5) | 0.610 (5×5) | — |
| 61 | 0.443 (3×9) | 0.410 (5×5) | 0.410 (5×5) | — |
| 101 | 0.446 (3×15) | 0.396 (8×5) | 0.396 (8×5) | 0.356 (6×6) |
| 149 | 0.423 (3×21) | 0.403 (4×15) | 0.369 (5×11) | 0.362 (6×9) |
| 197 | 0.442 (3×29) | 0.365 (4×18) | 0.335 (6×11) | 0.335 (6×11) |
| 257 | 0.432 (3×37) | 0.358 (4×23) | 0.331 (5×17) | 0.315 (9×9) |
| 293 | 0.420 (3×41) | 0.358 (5×21) | 0.358 (5×21) | 0.311 (7×13) |
| 349 | 0.421 (3×49) | 0.344 (4×30) | 0.330 (5×23) | 0.281 (7×14) |
| 401 | 0.411 (3×55) | 0.349 (4×35) | 0.312 (5×25) | 0.269 (6×18) |
| 449 | 0.408 (3×61) | 0.339 (4×38) | 0.290 (5×26) | 0.281 (9×14) |
| 509 | 0.407 (3×69) | 0.346 (4×44) | 0.285 (5×29) | 0.251 (8×16) |
| 601 | 0.404* (3×81) | 0.326* (4×49) | 0.275* (5×33) | 0.250* (6×25) |
| 701 | 0.407* (3×95) | 0.325* (4×57) | 0.278* (5×39) | 0.231* (6×27) |
| 797 | 0.403* (3×107) | 0.316* (4×63) | 0.263* (5×42) | 0.226* (6×30) |
| 997 | 0.400* (3×133) | 0.313* (4×78) | 0.256* (5×51) | 0.223* (6×37) |
| 1009 | 0.401* (3×135) | 0.317* (4×80) | 0.253* (5×51) | — |
| 1061 | 0.399* (3×141) | 0.324* (4×86) | 0.254* (5×54) | — |
| 1117 | 0.400* (3×149) | 0.312* (4×87) | 0.255* (5×57) | — |
| 1201 | 0.397* (3×159) | 0.313* (4×94) | 0.254* (5×61) | — |
| 1249 | 0.396* (3×165) | 0.307* (4×96) | 0.248* (5×62) | — |
| 1321 | 0.397* (3×175) | 0.306* (4×101) | 0.246* (5×65) | — |
| 1429 | 0.397* (3×189) | 0.308* (4×110) | 0.241* (5×69) | — |
| 1493 | 0.396* (3×197) | 0.305* (4×114) | 0.244* (5×73) | — |
| 1609 | 0.393* (3×211) | 0.306* (4×123) | 0.245* (5×79) | — |
| 1669 | 0.394* (3×219) | 0.300* (4×125) | 0.237* (5×79) | — |
| 1733 | 0.393* (3×227) | 0.300* (4×130) | 0.234* (5×81) | — |
| 1801 | 0.391* (3×235) | 0.300* (4×135) | 0.236* (5×85) | — |
| 1901 | 0.393* (3×249) | 0.301* (4×143) | 0.242* (5×92) | — |
| 1993 | 0.393* (3×261) | 0.295* (4×147) | 0.233* (5×93) | — |
| 2069 | 0.393* (3×271) | 0.296* (4×153) | 0.242* (5×100) | — |
| 2137 | 0.392* (3×279) | 0.294* (4×157) | 0.227* (5×97) | — |
| 2221 | 0.390* (3×289) | 0.294* (4×163) | 0.232* (5×103) | — |
| 2293 | 0.391* (3×299) | 0.297* (4×170) | 0.233* (5×107) | — |
| 2357 | 0.391* (3×307) | 0.294* (4×173) | 0.231* (5×109) | — |
| 2417 | 0.391* (3×315) | 0.293* (4×177) | 0.228* (5×110) | — |
| 2521 | 0.389* (3×327) | 0.290* (4×183) | 0.236* (5×119) | — |
| 2617 | 0.391* (3×341) | 0.295* (4×193) | 0.239* (5×125) | — |
| 2689 | 0.389* (3×349) | 0.292* (4×196) | 0.225* (5×121) | — |
| 2749 | 0.390* (3×357) | 0.291* (4×200) | 0.224* (5×123) | — |
| 2801 | 0.389* (3×363) | 0.290* (4×203) | 0.223* (5×125) | — |
| 2897 | 0.388* (3×375) | 0.291* (4×211) | 0.219* (5×127) | — |
| 2969 | 0.389* (3×385) | 0.290* (4×215) | 0.221* (5×131) | — |


**Ranges of `P_K/p`** (all `p ≡ 1 mod 4` up to 999; 26 sampled primes in `[1009, 2969]`):

| p range | K=3 | K=4 | K=5 | K=6 |
|---|---|---|---|---|
| 100–299 | 0.413–0.468 (18/18 exact) | 0.357–0.440 (18/18 exact) | 0.331–0.434 (18/18 exact) | 0.299–0.434 (18/18 exact) |
| 300–509 | 0.406–0.426 (16/16 exact) | 0.333–0.363 (16/16 exact) | 0.282–0.363 (16/16 exact) | 0.251–0.331 (16/16 exact) |
| 510–999 | 0.397–0.411 (6/35 exact) | 0.311–0.340 (1/35 exact) | 0.253–0.305 (0/35 exact) | 0.215–0.253 (0/35 exact) |
| 1000–1999 | 0.391–0.401 (0/14 exact) | 0.295–0.324 (0/14 exact) | 0.233–0.255 (0/14 exact) | — |
| 2000–3000 | 0.388–0.393 (0/13 exact) | 0.290–0.297 (0/13 exact) | 0.219–0.242 (0/13 exact) | — |

For comparison, the Weil limits (sharpen note Prop. 4.1 and Thm 4.3, outside the balanced
window) are `K/2^K = 0.375, 0.25, 0.156, 0.094`.

**HEURISTIC reading.** (1) `P_K(p)/p` is below `0.47` for `p ≥ 100` and below `0.43` for
`p ≥ 300`, and decreases with `p` towards `K/2^K`: at `p ≈ 3000` the best values found are
`0.389` (`K = 3`), `0.290` (`K = 4`), `0.22` (`K = 5`).
The only values above `1/2` are the Hanson–Petridis-tight configurations at `p = 13, 37, 41`
(REFUTED hypothesis "`P_K(p) ≤ p/2` with a uniform gap for all `p`": `P_3(13) = 9 = 0.69p`,
`P_3(37) = 21 = 0.57p`, `P_5(41) = 25 = 0.61p`). (2) The extremal shapes are unbalanced:
the smaller side has size exactly `K` in 103/106 (`K = 3`), 93/106 (`K = 4`), 92/106 (`K = 5`)
and 45/71 (`K = 6`) rows (exact maxima for `p ≤ 509`, best found above), and otherwise at most
11. (3) Balanced pairs are tiny: `b(p) ≈ log₂ p`, and
`b(p)²/p ≤ 0.354` for `200 ≤ p ≤ 509`. **So the data say nothing about the balanced window
`m, n ≍ √p`; what they show is that the product maximum is governed by the Weil range, where
`K/2^K` is already proved (sharpen note §4).**

**Factorisation type of `F_A` for extremal configurations.** By Theorem 1.1,
`F_A = λP₁^mP₀^{m−1}G_A` exactly; the table gives the defect and the factorisation of `G_A` over
`F_p` (squarefree part; distinct-degree factorisation), for the profile extremals, the largest
balanced pair and a maximum clique:

| p | kind | m×n | r | g | g/d | distinct roots of G | roots in F_p | factor degrees of G (degree: count) |
|---|---|---|---|---|---|---|---|---|
| 233 | profile | 3×33 | 3 | 20 | 0.172 | 20 | 0 | 5:4 |
| 233 | profile | 4×21 | 1 | 33 | 0.284 | 33 | 1 | 1:1, 4:1, 7:1, 21:1 |
| 233 | profile | 5×16 | 3 | 39 | 0.336 | 39 | 1 | 1:1, 14:1, 24:1 |
| 233 | profile | 6×12 | 1 | 45 | 0.388 | 45 | 5 | 1:5, 2:2, 3:2, 4:2, 22:1 |
| 233 | profile | 7×11 | 1 | 40 | 0.345 | 40 | 0 | 40:1 |
| 233 | profile | 8×9 | 1 | 45 | 0.388 | 45 | 1 | 1:1, 44:1 |
| 233 | balanced | 9×9 | 1 | 36 | 0.310 | 36 | 0 | 18:2 |
| 233 | clique | 7×7 | 7 | 74 | 0.638 | 74 | 2 | 1:2, 8:2, 14:4 |
| 509 | profile | 3×69 | 3 | 50 | 0.197 | 50 | 2 | 1:2, 3:8, 6:4 |
| 509 | profile | 4×44 | 4 | 82 | 0.323 | 82 | 0 | 2:16, 10:5 |
| 509 | profile | 5×29 | 1 | 110 | 0.433 | 110 | 2 | 1:2, 9:4, 14:2, 22:2 |
| 509 | profile | 6×21 | 0 | 128 | 0.504 | 128 | 0 | 2:1, 4:2, 58:1, 60:1 |
| 509 | profile | 7×18 | 1 | 129 | 0.508 | 129 | 0 | 2:2, 11:1, 30:1, 84:1 |
| 509 | profile | 8×16 | 1 | 127 | 0.500 | 127 | 1 | 1:1, 11:1, 20:1, 95:1 |
| 509 | balanced | 11×11 | 1 | 134 | 0.528 | 134 | 2 | — |
| 509 | clique | 9×9 | 9 | 182 | 0.717 | 182 | 1 | — |
| 997 | profile | 3×133 | 3 | 102 | 0.205 | 102 | 2 | 1:2, 2:2, 3:8, 6:12 |
| 997 | profile | 4×78 | 1 | 187 | 0.376 | 187 | 0 | 2:1, 5:1, 15:1, 165:1 |
| 997 | profile | 5×51 | 5 | 248 | 0.498 | 248 | 0 | 2:2, 53:2, 69:2 |
| 997 | profile | 6×37 | 1 | 277 | 0.556 | 277 | 0 | 2:1, 6:1, 10:1, 30:1, 229:1 |
| 997 | balanced | 12×12 | 0 | 354 | 0.711 | 354 | 1 | — |
| 997 | clique | 13×13 | 13 | 342 | 0.687 | 342 | 0 | — |

The defect is large: among the 642 table configurations, `g/d ≥ 0.11` whenever
`min(m,n) ≥ 3` and `p ≥ 100` (minimum at `p = 109`, `3 × 17`); `g/d ≈ 0.2` for the `K = 3` extremals, `0.5–0.56`
for `K = 5, 6` at `p ≈ 1000`, and `≈ 0.7` for balanced pairs and maximum cliques. `F_A` is an
`m`-th power only along `B`; the cofactor `G_A` has as many distinct roots as its degree in all
but 23 of 642 configurations (the exceptions, mostly maximum cliques, are listed in the results
file), and its factor degrees are often strongly non-generic (§3.4).

## 3. Proof attempts and where they stop (task 3)

Throughout, "a gain" means a proof of `mn ≤ (1/2 − c)p + o(p)` for complete bicliques with
`m, n ≍ √p`. For each route we state what it proves, test it, and name the exact point where it
stops. None of the four gives a gain; the obstructions are PROVED where so labelled.

### 3.1 Route (a): two-sided arguments

**Lemma 3.1 (multiplicity Schwartz–Zippel on a grid; PROVED).** If `Φ ∈ F_p[x,y]` is nonzero of
total degree `t`, then `Σ_{(a,b)∈A×B} mult_{(a,b)}Φ ≤ t·max(|A|,|B|)`.

*Proof.* Write `Φ = Π_{b∈B}(y−b)^{s_b}·Ψ` with `Ψ(x,b) ≢ 0` for every `b ∈ B`. Then
`mult_{(a,b)}Φ = s_b + mult_{(a,b)}Ψ ≤ s_b + ord_{x=a}Ψ(x,b)`, and
`Σ_a ord_{x=a}Ψ(x,b) ≤ deg Ψ(x,b) ≤ deg Ψ = t − Σ s_b`. Summing,
`Σ mult ≤ |A|Σ_b s_b + |B|(t − Σ_b s_b) ≤ t·max(|A|,|B|)`. ∎ (Verifier: 300 random products of
linear forms over `F_7,…,F_17`.)

So a two-sided certificate for `m ≤ n` needs a nonzero `Φ` of degree `t` vanishing to order
`≥ M` on `A × B` with **`t/M < m`**; Hanson–Petridis corresponds to `t/M = d/m`. Every
polynomial we can write down from the data of the problem has `t/M ≥ d/max(m,n)`:

| polynomial | degree `t` | order `M` on `A × B₁` | `t/M` |
|---|---|---|---|
| `F_A(y)` | `d` | `m` | `d/m` |
| `F_B(x)` | `d` | `n` | `d/n` |
| `(x+y)^d − 1` | `d` | `1` | `d` |
| products `Π Φ_i^{α_i}` of these | `Σα_i t_i` | `Σα_i M_i` | `≥ min_i t_i/M_i = d/max(m,n)` |
| `Δ_AΔ_B[(x+y+s+t)^{d+m+n−2}] − C(m+n−2, m−1)` (a polynomial in `x+y`) | `≤ d` | `m+n−1` at `x+y = 0` only | — |

The last line is the natural "two-sided Hanson–Petridis polynomial" (double divided difference
over `A` in `s` and over `B` in `t`). It is a polynomial `ψ(x+y)` of degree `≤ d` that vanishes to
order `m+n−1` at every `z` with `z + A + B ⊆ Q ∪ {0}`, i.e. on `B(A+B)`; it bounds
`(m+n−1)|B(A+B)| ≤ d`, which is weaker than Hanson–Petridis for the pair `(A+B, B(A+B))`
(`|A+B| ≥ m+n−1`), and says nothing about `mn`. With `t/M ≥ d/max(m,n)` the grid lemma returns
`m·n ≤ d + O(m)` at best, i.e. Hanson–Petridis.

The two Hanson–Petridis polynomials also carry **no joint constraint**: by Theorem 1.1(3), each
of "`F_A` has the factor `P₁^mP₀^{m−1}`" and "`F_B` has the factor `Q₁^nQ₀^{n−1}`" is equivalent
to `A + B ⊆ Q ∪ {0}`, so any contradiction derivable from both is derivable from either. The only
two-sided identities with content are those of Lemma 1.2 (the values of `G_A` on `B` and of `G_B`
on `A`), and they constrain `(A,B)` only when `g ≤ max(m,n) − 2` (Remark 1.2′).

**Obstruction (a).** A gain needs an auxiliary polynomial vanishing on `A × B` with
`t/M ≤ (1−c)d/m` for `m ≈ n`, i.e. vanishing *more efficiently than any product of one-sided
objects*. Its vanishing would have to come from a cancellation that couples `a` and `b` beyond
`a + b ∈ Q`; the Lagrange trick of Hanson–Petridis produces cancellation in one variable only.
OPEN.

### 3.2 Route (b): abc, Mason–Stothers, Brownawell–Masser, Voloch

Put `f_k = c_k(x+a_k)^D`, so `F_A = −1 + Σ_k f_k` is an `(m+2)`-term identity
`F_A − f_1 − ⋯ − f_m + 1 = 0`, all degrees `≤ D < p`.

**Proposition 3.2 (Wronskian identity; PROVED).**
`W(F_A, f_1, …, f_m) = κ·Π_k (x+a_k)^{D−m}` with
`κ = −Π_k c_k · Π_{i=1}^{m}(D)_i · det[(x+a_k)^{m−i}]_{i,k} ∈ F_p^*` (the determinant is a
Vandermonde in the `x + a_k`, a nonzero constant).

*Proof.* Column operations replace `F_A` by `−1`; expanding along that column leaves
`−W(f_1', …, f_m') = −det[(D)_i c_k (x+a_k)^{D−i}]_{1≤i≤m, k}`; take `c_k(x+a_k)^{D−m}` out of
column `k` and `(D)_i` out of row `i`. All `(D)_i` are units because `D < p`. ∎
(Verifier: exact polynomial determinant for every `A ∋ 0` normalised, `m ≤ 3`, `p ∈ {13, 17}`;
53 cases.)

**Proposition 3.3 (the Wronskian route is blind to `B`; PROVED).**
1. The linear system `L = ⟨1, (x+a_1)^D, …, (x+a_m)^D⟩` (dimension `m+1`) has no inflection
   points off `−A`: at every `β ∈ F̄_p ∖ (−A)` the vanishing orders of the members of `L` are
   exactly `{0, 1, …, m}`. In particular for **every** `β ∉ −A` some member of `L` vanishes to
   order exactly `m` at `β`, and none to higher order.
2. `b ∉ −A` is complete (`b + A ⊆ Q`) iff the *fixed* member `F_A ∈ L` is that member at `b`
   (Theorem 1.1(1)). So complete points are not distinguished by any local invariant of `L`:
   they are the points where the osculating hyperplane of the curve `φ_L : P¹ → P^m` is one fixed
   hyperplane `H_{F_A}`. Counting them through `H_{F_A} ∩ φ_L(P¹)` is Bézout,
   `Σ_b ord_b F_A ≤ deg F_A = d`, which is Hanson–Petridis.
3. For an `N`-term identity, the Wronskian `W` of any `N − 1` of the terms satisfies
   `ord_β W ≥ Σ_i max(0, ord_β f_i − (N−2))` (divide column `i` by
   `(x−β)^{ord_β f_i − (N−2)}`); this is the engine of the Brownawell–Masser and Voloch
   inequalities (UNVERIFIED CITATION for their exact published forms: W. D. Brownawell,
   D. W. Masser, Math. Proc. Camb. Phil. Soc. 100 (1986) 427–434; J. F. Voloch, Bol. Soc. Brasil.
   Mat. 16 (1985) 29–39). Here `N − 2 = m` and `deg W = m(D − m)`, and the pure powers alone
   contribute `Σ_k (D − m) = m(D−m)` at the points `−a_k`. So the inequality says exactly
   `Σ_{β ∉ −A} max(0, ord_β F_A − m) ≤ 0`, i.e. Theorem 1.1(5), and nothing about how many roots
   of multiplicity `m` there are. Eliminating terms with the pencil (`F_{A∖E}` has `m − |E| + 2`
   terms and order `m − |E|` on `B`) keeps "truncation level = order on `B`", so no
   sub-identity helps either.
4. Mason–Stothers for the three-term identity `(F_A + 2) − F_A − 2 = 0` (equivalently
   Riemann–Hurwitz for the map `F_A : P¹ → P¹`, tame since `d < p`) gives
   `d + 1 ≤ n₀(F_A) + n₀(F_A + 2)`, which is Proposition 1.4: complete points enter with weight
   `m − 1`, **less** than Hanson–Petridis's `m`. Riemann–Hurwitz for the pencil map `R`
   (Proposition 1.3(3), degree `n + n₀(G_A)`) sees only points with one bad partner.

*Proof.* (1) `W(β) ≠ 0` (Proposition 3.2) means the `(m+1)×(m+1)` matrix of the derivatives of
orders `0..m` of the basis at `β` is invertible; so for each `s ≤ m` there is a member with order
exactly `s`, and a member with order `≥ m+1` would give a nonzero kernel vector. (2) is
Theorem 1.1(1) plus Bézout. (3): the column bound is immediate; the degree count is
Proposition 3.2. (4): Mason–Stothers, `max deg ≤ n₀(uvw) − 1` for coprime `u + v = w` with
`u', v'` not both zero, holds in characteristic `p` (standard; here `u = F_A + 2`, `v = −F_A`,
`w = 2`, `u' = F_A' ≠ 0` as `0 < d < p`). With `n₀(F_A) ≤ d − (m−1)n + r` (Theorem 1.1(3)) and,
since `F_A + 2` vanishes to order `≥ m` (`≥ m−1` on `−A`) at anti-complete points by the same
computation, `n₀(F_A + 2) ≤ d − (m−1)n⁻ + r⁻`, it becomes Proposition 1.4. ∎

**Obstruction (b).** The abc/Wronskian machinery measures ramification; complete points are
unramified for the relevant linear system (Proposition 3.3(1)). Any inequality built from local
contact orders of `L` alone is satisfied by every configuration with `mn ≤ d + r`, including the
tight ones (`A = {0,1}`, `p = 13, 37, 41`). A gain would need a *global* input about the
hyperplane `H_{F_A}` beyond Bézout, for instance the rationality of its contact points
(`B ⊆ F_p`), which the Wronskian does not see.

### 3.3 Route (c): Hankel minors with `e ≥ 1`, exchanged roles, extra vanishing

**Lemma 3.4 (PROVED).** Let `A + B ⊆ Q ∪ {0}` with `B ⊆ {b : e_b = 0}`, `m ≤ (p+1)/2`,
`0 ≤ e ≤ (m−1)/2`. Then the robust note's inequality (★)_e reads

    Σ_{x ∉ B, e_x ≤ e} (e+1−e_x)(m − (3e+e_x)/2 − δ_x)  ≤  (e+1)(g − e + 3en/2),

so (★)_e restricted to `B` is implied by Hanson–Petridis (`e ≤ g + 3en/2` always), for `A` and,
with the roles exchanged, for `B`. At `B₁` the minor `H_{e+1}` vanishes to order exactly
`(e+1)(m−e)` whenever `Δ_e ≢ 0` (Proposition 1.3(2)).

*Proof.* The points of `B` contribute `(e+1)(n(m − 3e/2) − r)`; subtract from `(e+1)(d−e)` and use
`d = mn − r + g`. ∎

*Checks.* `Δ_e ≢ 0 (mod p)` for all `3 ≤ m ≤ 40`, `1 ≤ e ≤ min(5,(m−1)/2)` and
`p ∈ {13, 29, 101, 509, 997, 1009, 10009}` (904 cases), and `Δ₁ = m(m−D)/(D²(D−1))`
exactly. The cofactor `K_e` is nonzero at every `b ∈ B₁` in all 8,108 computed cases
(`e = 1, 2`, `p ≤ 29` exhaustive, and the table configurations with `p ≤ 200`): **no extra
vanishing at `B`**.

So `e ≥ 1` can only help through points `x ∉ B` with few bad partners. Quantitatively, every
such `x` with `1 ≤ e_x ≤ e = ηm` has weight `≥ m − 2e` in (★)_e, so if
`X_η := #{x ∉ −A ∪ B : 1 ≤ e_x ≤ ηm}` then
`n(m − 3e/2) ≤ d + r − e − X_η(m − 2e)/(e+1)`, and a constant gain requires

    X_η  ≳  (3/4)·η²p/(1 − 2η).                                                  (3.1)

(3.1) is an anti-concentration statement for `f_A(x) = Σ_{a∈A}χ(x+a)`: at least `≈ η²p` points
with `f_A(x) ≥ (1−2η)m`. HEURISTIC: for `|A| = m ≍ √p` and pseudo-random `A` the expected count
is `p·P(Bin(m,1/2) ≤ ηm) = p·2^{−(1−H(η))m + o(m)}`, astronomically smaller. In the extremal
configurations (verifier, `few_bad_partner_counts`): for the largest balanced pairs,
`m = 11` at `p = 509` and `m = 12` at `p = 997`, `X_η = 4, 40, 100` and `14, 45, 154` for
`η = 0.2, 0.3, 0.4`, against the required `25, 86, 305` and `50, 168, 598`. For the `K = 5, 6`
extremals (`m = 5, 6`) `X_{0.2}` does exceed the requirement (`56, 29` at `p = 509`; `125, 58` at
`p = 997`): these are in the Weil range, where `m` is bounded and `P_K/p → K/2^K` is already
known; there (3.1) holds and (★) is consistent with that gain.

**Theorem 3.5 (the count constraints do not exclude a tight balanced profile; PROVED by exact
rational certificates).** For each `(p, m)` in the table below, with `n₀ = ⌊d/m⌋` (so
`m n₀ ≥ d − m + 1`), there are nonnegative rationals `n_j` (`0 ≤ j ≤ m`; the number of points
`x ∉ −A` with `e_x = j`) and `n'_j` (`0 ≤ j ≤ m−1`; points of `−A`), with `n_0 = n₀`, `n'_0 = 0`
(`r = 0`), **all other mass in the band `0.4m ≤ e_x ≤ 0.6m`** (so `|f_A(x)| ≤ 0.2m + 1` off `B`,
where `f_A(x) = Σ_{a∈A}χ(x+a) = m − δ_x − 2e_x`), satisfying every one of:

- the exact identities `Σ(n_j + n'_j) = p`, `Σ n'_j = m`, `Σ_x f_A(x) = 0`,
  `Σ_x f_A(x)² = m(p−m)`;
- (★)_e for every `0 ≤ e ≤ (m−1)/2`, for `A` and for `νA` (robust note Thm 2.1);
- the averaged subset Hanson–Petridis inequalities for every `1 ≤ t ≤ m`, for `A` and `νA`
  (robust note Prop. 1.2);
- the pencil inequality `(2m−2)N₀ + (m−2)N₁ ≤ 2d − 2` for `A` and `νA` (pencil note Thm 3.1);
- Proposition 1.4;
- the Weil moment bounds `Σ_x f_A(x)^{2k} ≤ (2k−1)!!m^kp + (2k−1)m^{2k}√p`, `k = 2,3,4`, and the
  Weil counts `|N_j − C(m,j)p/2^m| ≤ C(m,j)(m/2)(√p+1) + 2m` (sharpen note Lemma 3.1), with `√p`
  replaced by `⌊√p⌋` (which only makes them stricter).

| `p` | `m` | `n₀` | `m n₀/p` | constraints (theorem rows) | support of the certificate | mean square of `f_A` off `B` |
|---|---|---|---|---|---|---|
| 1009 | 22 | 22 | 0.4797 | 158 (124) | `n₀=22`, `n₉≈281.9`, `n₁₁≈307.8`, `n₁₃≈375.4`, `n'₁₃=22` | 11.2 |
| 10009 | 71 | 70 | 0.4966 | 484 (370) | `n₀=70`, `n₃₃≈4978`, `n₃₇≈2524`, `n₄₀≈2366`, `n'₄₂=71` | 35.5 |
| 40009 | 141 | 141 | 0.4969 | 946 (720) | `n₀=141`, `n₆₆≈13207`, `n₆₈≈7138`, `n₇₅≈19382`, `n'₇₀=141` | 70.7 |
| 100049 | 224 | 223 | 0.4993 | 1492 (1134) | `n₀=223`, `n₁₀₉≈53369`, `n₁₁₄≈41365`, `n₁₃₃≈4869`, `n'₁₁₂=224` | 111.9 |
| 1000033 | 707 | 707 | 0.4998 | 4680 (3550) | `n₀=707`, `n₃₄₇≈416627`, `n₃₅₇≈567137`, `n₄₁₉≈14855`, `n'₃₅₃=707` | 353.4 |

The certificate is found by a floating LP (HiGHS interior point, maximising the common relative
slack of the theorem rows), rounded to denominator `1000`, repaired exactly on the four
equalities (one `n'_j` and three `n_j` re-solved in rationals), and then every constraint is
checked in exact rational arithmetic (7,740 inequalities and 20 equalities in total, all
satisfied; verifier `lp_certificates`). By weak LP duality, **no nonnegative combination of these
constraints implies `m·n₀ ≤ (1/2 − c)p` for any `c > m/p`** at these parameters. Since the
constraints are one-sided, the same profile serves for `B` (`m = n`); the only coupling of the
two sides that is linear in these counts is `N₀(A) ≥ n`, `N₀(B) ≥ m`. The band restriction
shows more: the certificate has **no** point off `B` with fewer than `0.4m` bad partners (or
fewer than `0.4m` good ones), so every inequality that only counts points with `e_x < 0.4m` or
`e_x > 0.6m` — any conceivable robust or joint (complete + anti-complete) Hankel-type
inequality in that range — sees only `B`, where Hanson–Petridis-tightness is consistent.

*A remark on the unrestricted LP.* Without the band, the LP optimum at large `p` puts the `m`
points of `−A` (or about `m` other points) at `e_x ≈ m − O(1)` — `A` nearly independent in the
Paley graph — so that they carry half of the second moment. Such a profile is not excluded by any
known bound either (a nearly independent `m`-set with `m ≈ √(p/2)` is itself at the
Hanson–Petridis limit for `νA`), but the band certificate avoids the question.

**What a gain would have to show (PROVED reformulation).** For a complete biclique with
`B ⊆ B(A)`, the points of `B` carry `n m² = m(mn)` of the total `L²` mass
`Σ_x f_A(x)² = m(p−m)`; Hanson–Petridis says this share is at most `(d + r)/(p − m) ≈ 1/2`. So the
target is: *for `|A| = m ≍ √p`, the function `f_A(x) = Σ_{a∈A}χ(x+a)` cannot carry more than
`(1/2 − c)` of its `L²` mass on `≍ √p` points where it attains its maximum `m`.* In the
certificates `f_A` has mean square `11.2, 35.5, 70.7, 111.9, 353.4` off `B`, against
`m/2 = 11, 35.5, 70.5, 112, 353.5` (a random `A` would give `≈ m`): sub-binomial, exactly as this
forces. The
Weil moment bounds cannot see this: `B` contributes `n m^{2k} ≈ (p/2)m^{2k−1}` to the `2k`-th
moment against a Weil error `(2k−1)m^{2k}√p`, a ratio `√p/(2(2k−1)m) = 1/(√2(2k−1)) + o(1) < 1`
for `m ≈ √(p/2)` (verifier `weil_moment_ratio`). This is the square-root barrier in the form it
takes here.

**Obstruction (c).** Hankel minors are less efficient than `F_A` on `B` (ratio `(m−e)/m`), have
no extra vanishing at `B`, and can gain only through (3.1), an anti-concentration input at the
square-root scale that no available constraint supplies (Theorem 3.5).

### 3.4 Route (d): zeros of `F_A` outside `B`

The premise fails as stated: near-tightness does not force zeros outside `B`.

- **REFUTED** ("near-tightness forces `F_A` to vanish at further points"): in the tight
  configurations `A = {0,1}` (every `p ≡ 1 mod 4`) and `p = 13, 37, 41` one has `G_A = 1`, so
  `F_A` has no zeros at all outside `B`.
- Rational zeros of `G_A` are constrained: at `x ∉ −A ∪ B(A)` the multiplicity is at
  most `min(e_x − 1, m − e_x)` (Theorem 1.1(5); checked at every such `x` for all normalised
  `A`, `p ≤ 29`), so points with one bad partner are never zeros. In the 642 table
  configurations (`29 ≤ p ≤ 997`; profile extremals, balanced witnesses, maximum cliques,
  sum-cliques) the number of roots of `G_A` in `F_p` has histogram (roots: configurations)
  `0:267, 1:125, 2:149, 3:32, 4:24, 5:3, 6:25, 7:1, 8:8, 10:4, 14:2, 18:1, 20:1`, and `G_A` is
  squarefree in 619 of them (profile extremals 428/434, balanced 51/54, cliques 68/77,
  sum-cliques 72/77). The exceptions are mostly maximum cliques and
  clique-like pairs (`r = m`): at `p = 113` the 7-clique has `g = 14` and `G_A` has 10 rational
  roots, one of multiplicity 3; at `p = 349` the 9-clique has `g = 102` and 20 rational roots of
  `G_A`. So `F_A` *can* vanish at many rational points outside `B` (HEURISTIC: far more often for
  cliques than for random polynomials, which have about one), but these configurations are far
  from tight (`g/d ≈ 0.6–0.7`), and the extra zeros cost defect rather than create a
  contradiction: every zero outside `B` is paid for inside `deg G_A = g`.
- What near-tightness does force is the opposite: few distinct roots (`n₀(F_A) ≤ n + g`),
  few anti-complete points (Proposition 1.4), few one-bad points (pencil), all upper bounds.
- The formal `m`-th root. Since
  `F_A/λ = Σ_{j=0}^{d} [(−1)^j (1/2)^{(j)}/m^{(j)}]·h_j(A)·x^{d−j} − 1/λ`
  (rising factorials; from `C(D, m−1+j)/C(D, m−1) = (d)_j/m^{(j)}` and `d ≡ −1/2`), the
  statement `F_A = λP₁^mP₀^{m−1}G_A` is a system of `d` polynomial equations in the `m + n + g`
  unknowns `(A, B, G_A)` (coefficient formula checked in the verifier). It is overdetermined by
  `d − m − n − g = (m−1)(n−1) − 1 − r` equations (before quotienting by the two-dimensional
  symmetry group), but a count of equations is not an inequality: complete bicliques with
  `mn ≤ d + r` exist for every `p`, so the system is solvable in that range; Kalmynin's method extracts
  contradictions from its top coefficients only when `g = 0` (and Remark 1.2′ shows the
  interpolation version works only for `g ≤ max(m,n) − 2`).
- HEURISTIC (structure of extremal configurations). The factorisation patterns of `G_A` over
  `F_p` are often far from those of random polynomials: e.g. `p = 233, m = 3`: four quintics;
  `p = 509, m = 4`, `A = {0,1,121,122}`: sixteen quadratics and five decics;
  `p = 997, m = 3`: degrees `1,1,2,2`, eight cubics and twelve sextics. For
  `A = {0,1,121,122}` this is explained by the affine involution `x ↦ 122 − x` of `A`, which forces
  `F_A(−x − 122) = F_A(x)` (checked), so `F_A` and `G_A` are polynomials in `(x+61)²`. The other
  cases have no affine stabiliser; the patterns suggest that extremal sets are special
  (symmetric) configurations. Not pursued.

**Obstruction (d).** Zeros of `F_A` outside `B` are not forced by near-tightness (they are
absent in the tight examples) and, when present, are paid for inside the defect `g`; the
structure that near-tightness does force is "few distinct roots", which is exactly what
Hanson–Petridis counts.

## 4. Literature check

Searches and fetches made on 2026-09-29 (WebSearch/WebFetch); earlier pins in
`research/sigma-biclique-2026-09-05.md` §1 were not repeated.

- **Randomstrasse101** (A. S. Bandeira, D. Dmitriev, K. Lucca, P. Nizić-Nikolac, A. Rödder,
  *Randomstrasse101: Open Problems of 2025*, arXiv:2603.29571, 31 Mar 2026; the Paley post
  https://randomstrasse101.math.ethz.ch/posts/PaleyGraph/, fetched through a summarising tool, so
  the quotations are close paraphrases, not verified verbatim): Problem 25 asks for
  `ω(G_p) = O(polylog p)`; the post records Hanson–Petridis' `ω(G_p) ≤ (1+o(1))√(p/2)` as the best
  bound, and its Conjecture 27 (`ϑ(Ḡ_{p,2}) ≤ (2/3)√p` for the 2-localised theta function) is
  presented as something that *would* give "a concrete improvement over the Hanson–Petridis
  bound", with Kunisky's empirical observation `≈ (√(1/2) − ε)√p` as motivation. So, as of
  March 2026, the constant `1/√2` for prime `p` was open.
- **Yip.** *On the clique number of Paley graphs of prime power order*, FFA 77 (2022) 101930
  (arXiv:2004.01175): `q = p^{2s+1}` only. *Exact values and improved bounds on the clique number
  of cyclotomic graphs*, DCC 93 (2025) (arXiv:2304.13213): `ω ≤ √|S/S| + √(q/p)`, which for
  `q = p`, `S = Q` is `√((p−1)/2) + 1`, Hanson–Petridis strength (abstract-level check, as in
  the biclique note).
- **The critical case.** Kalmynin (arXiv:2504.10202v2, local copy; Lemma 4, Theorems 2–4),
  Rudnev–Tyrrell (arXiv:2607.24270, local copy; §3 "reciprocal transfer identities") and
  Yip–Yoo (arXiv:2608.02568) treat `g = 0`: `A + B = μ_d` (then `r = 0`, `g = 0`) and
  `A − A = μ_d ∪ {0}` (cliques with `g = 0`). Kalmynin and Rudnev–Tyrrell compare coefficients
  of the exact factorisation `F_A = C·Π_{b∈B}(x−b)^{α−ε(b)}` (Kalmynin's Lemma 4 is our
  Theorem 1.1(3) with `g = 0`); Yip–Yoo's later proof (described in Rudnev–Tyrrell's
  introduction as based on "differential identities for the root polynomials"; not re-read here)
  also starts from the critical case. Their identities use `g = 0`. Remark 1.2′ shows that the
  value-interpolation relations of Lemma 1.2 exist only for `g ≤ max(m,n) − 2`; whether other
  coefficient identities survive a defect of order `p` is not settled here (it would contradict
  nothing we prove, but we found no mechanism).
- A web search for improvements of the Hanson–Petridis clique bound for prime `p`
  ("improvement Hanson-Petridis bound clique number Paley graph prime order", and a 2025–2026
  variant) returned only the items above. **No improvement of the constant for prime `p` was
  found; completeness of the search is not claimed.**

Nothing in this note improves the constant either.

## 5. Open obligations

1. **The target (OPEN).** `mn ≤ (1/2 − c_K)p + o(p)` for complete bicliques with
   `min(m,n) ≥ K`; only the window `M(K,κ) < m ≤ n < 12√p` (sharpen note §4) is open. For
   `B = B(A) ∖ (−A)` it is equivalent (§3.3) to: *for `|A| = m ≍ √p`,
   `f_A = Σ_{a∈A}χ(· + a)` puts at most `(1/2 − c)` of its `L²` mass `m(p−m)` on the points where
   `f_A = m`.*
2. **A two-sided auxiliary polynomial** (route (a)): a nonzero `Φ(x,y)` of degree `t` vanishing
   to order `M` on `A × B` with `t/M ≤ (1−c)d/m` (Lemma 3.1). Nothing of this kind is known;
   Rudnev–Tyrrell's reciprocal transfer is two-sided but needs `g = 0`.
3. **An anti-concentration input** (route (c)): `#{x : 1 ≤ e_x ≤ ηm} ≳ η²p` for the `A`-side of a
   near-tight balanced biclique would suffice through (★) (inequality (3.1)); Theorem 3.5 shows
   no current count constraint implies it, and for pseudo-random `A` it is false, so it would
   have to come from near-tightness itself.
4. **Extend Theorem 3.5** to a two-sided LP with genuinely joint constraints (e.g. counts of
   `x` by the pair `(e_x(A), e_x(B))`, or the Möbius/reciprocal maps of the critical-case papers)
   and test whether it becomes infeasible below `1/2`.
5. **Near-critical Kalmynin.** Remark 1.2′ gives `n − 1 − g` explicit relations when
   `g ≤ n − 2` (and `m − 1 − g` when `g ≤ m − 2`). Working them out could at best give additive
   improvements `mn ≤ d + r − c·max(m,n)`, not multiplicative ones. Not attempted.
6. **Exact data beyond `p = 509`.** For `510 ≤ p ≤ 999` only 6/35 (`K = 3`) and 1/35 (`K = 4`)
   entries were certified; `K = 5, 6` were not attempted there, nor any certification above
   `p = 1000`. The helper's search needs better caps for `M_m` beyond the computed profile
   (`m ≤ 8` or `m ≤ 6`). The lower bounds are recorded with witnesses.
7. **Why are the `G_A` of extremal and clique configurations non-generic** (repeated roots,
   many rational roots, factor degrees in orbits)? Only the affine-involution case is explained.
8. Independent review of Theorem 1.1, Lemma 1.2, Propositions 1.3–1.4, 3.2–3.3 and of the LP
   certificate script.
