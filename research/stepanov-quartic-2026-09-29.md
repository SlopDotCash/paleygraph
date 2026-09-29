# Stepanov wave, worker `quartic`: the quartic sign matrix of a Paley biclique (2026-09-29)

**Status: the target `|A||B| ≤ (1/2 − c)p` for balanced complete bicliques stays OPEN; the quartic idea does not reach it. PROVED: (i) the exact second-moment ("Jacobi") identities for every pair of characters of 2-power order, one- and two-sided (Lemma 1.2), including the forced phase coherence `|Σ_x f_ψ(x)²| = √p·m(m−1)` for a clique of size `m`; (ii) the sign-class rows `(★)_k`, a pattern Hanson–Petridis inequality, and the relation between the quadratic and the order-`k` Hanson–Petridis polynomials: the product of the `k/2` class polynomials has degree exactly `d` and vanishes (to order `≥ m`) only at constant columns, so the quartic budget `d/2` per sign is the quadratic budget split in two (Proposition 1.7); (iii) Theorem 2.1, the exact solution of the optimisation "quadratic HP + all order-`k` class rows + moments, restricted to the biclique": the class rows cap the biclique weight `mn − r` at `d·min(1, Θ)`, with `Θ` an explicit function of the row and column statistics of the sign matrix. `Θ ≥ 1` (no gain) for every half-weight statistic, in particular every flip-invariant one with `r = 0`. A constant gain occurs only for lopsided sign matrices, with an explicit bound via the robust inequality. PROVED by exact computation (Theorem 3.1): rational certificates (two-sided joint profiles plus an explicit sign matrix) satisfy every row of the family: quadratic, quartic and (for `p ≡ 1 mod 8`) octic class rows, subset HP, all one- and two-sided Jacobi identities, Weil and Krawtchouk–Weil rows, union rows for `A ∪ B`, and pattern and monochromatic-rectangle HP. They reach `mn/p = 0.4952–0.4998` at eight primes up to `1,000,033`. One-sided clique certificates reach `m(m−1)/d` up to `0.990` (and `1` in the critical case `g = 0`) and realise the forced coherence; the middle Krawtchouk–Weil rows are checked in floating point with an error allowance (§3.3). REFUTED: "the quartic constraints force `mn ≤ (1/2 − c)p`" and "a balanced sign matrix conflicts with near-tightness through the quartic second moment". Witnesses: the certificates, in which every `(★)_4` and `(★)_8` row vanishes identically, and the actual exactly-tight cliques at `p = 13, 41`, whose sign matrices are exactly balanced (`Θ = ∞`). Also REFUTED: "the quartic rows never see actual bicliques"; witness: the `11 × 11` sum-clique at `p = 593`, which has `Θ = 0.95`. CONDITIONAL (Proposition 5.1): if the order-`k` analogue of the two-set conjecture holds, then `Θ ≥ 1` (rows and columns) for every complete biclique in the balanced window, so this route cannot give a gain. HEURISTIC (813 actual bicliques, `p ≤ 1500`): the median of `Θ` is 4–14 by type, and `Θ < 1` occurs only at `p = 17` and `p = 593`.**

Verifier: `experiments/stepanov_quartic_2026_09_29.py` → `results/stepanov_quartic_2026_09_29.json`
(run with `/opt/miniconda3/bin/python3`; about 100 s on this machine; **2,157,238 checks, 0 failures**;
standard library + numpy + scipy (HiGHS LP, floating point; every certificate is re-checked exactly)
+ sympy (primitive roots)). The JSON contains every certificate (profiles, joint profiles, the sign
matrices for `mn ≤ 3000`), the exact Jacobi constants `K(φ^i, φ^j)` at each prime, and all
statistics quoted below. Scratch prototypes were kept outside the workspace.

---

## 0. Notation

`p ≡ 1 (mod 4)`, `d = (p−1)/2`, `Q` the nonzero squares, `g` a fixed primitive root, `ind` the
discrete logarithm. For `k | p−1` put `d_k = (p−1)/k`; the *class* of `y ≠ 0` is `ind(y) mod k`, and
`φ_k(y) = e^{2πi·ind(y)/k}` (with `φ_k(0) = 0`) is a character of order `k`. `ψ = φ_4`, `χ = φ_2 = ψ²`.
For `k = 2^s`, `φ_k^{k/2} = χ`. For `X ⊆ F_p`, `|X| = m`, and `x ∈ F_p`:

- `δ_x = [x ∈ −X]`; the *type* of `x` (order `K`) is `(δ_x; n_0(x), …, n_{K−1}(x))`,
  `n_c(x) = #{a ∈ X : a + x ≠ 0, ind(a+x) ≡ c mod K}`;
- `f^X_α(x) = Σ_{a∈X} α(x + a)`; so `f^X_{φ^j}(x) = Σ_c n_c(x) ζ^{jc}`, `ζ = e^{2πi/K}`;
- the class-`(k,c)` bad count `e^{(k,c)}_x = #{a : a + x ≠ 0, ind(a+x) ≢ c mod k}` (for `k = 2`,
  `c = 0` this is the usual `e_x`).

A complete biclique is `A + B ⊆ Q ∪ {0}`, `|A| = m`, `|B| = n`, `r = |B ∩ (−A)|`; Hanson–Petridis (HP;
CITED: arXiv:1905.09134, Thm 1.2, local `sources/sigma-hanson-petridis-1905.09134.txt`) is
`mn − r ≤ d`. Its order-`K` **sign matrix** is `M_{ab} = ind(a+b) mod K` (`a + b ≠ 0`); for
`K = 4` we also write `M_{ab} = ψ(a+b) ∈ {±1}`. For `b ∈ B` and an even class `c`,
`t^c_b := e^{(k,c)}_b` is the number of entries of column `b` of `M^{(k)}` that differ from `c`.
"The robust note" is `research/stepanov-robust-2026-09-26.md`, "the balanced note" is
`research/stepanov-balanced-2026-09-29.md`, "the sharpen note" is `research/stepanov-sharpen-2026-09-27.md`,
"the hybrid note" is `research/stepanov-hybrid-2026-09-29.md`.

## 1. Identities and inequalities

### 1.1 The sign classes

**Lemma 1.1 (PROVED).** Let `k = 2^s ≥ 4` divide `p − 1`. If `a + b ∈ Q` then `ind(a+b) mod k` is even,
i.e. `φ_k(a+b) ∈ μ_{k/2}`. In particular, on a complete biclique every nonzero entry of the quartic
sign matrix is `±1`, and every nonzero entry of the octic one (`p ≡ 1 mod 8`) lies in `{±1, ±i}`.
*Proof.* `Q = {g^{2j}}`, and `k` is even. ∎

Multiplying `A`, `B` (and `x`) by a square `h` preserves completeness and rotates classes by
`ind h`, an even number: this is the **rotation symmetry** `c ↦ c + 2` used throughout. With
`ψ(h) = −1` it swaps the two quartic sign classes.

### 1.2 Jacobi second moments

For characters `α, β` of `F_p^*` put `K(α,β) := Σ_{v∈F_p} α(v)β(1+v) = α(−1)J(α,β)`, where
`J(α,β) = Σ_w α(w)β(1−w)`.

**Lemma 1.2 (PROVED).** For nontrivial `α, β` and `X, Y ⊆ F_p`,

    Σ_{x∈F_p} f^X_α(x) f^Y_β(x) = K(α,β) · Σ_{a∈X, b∈Y, a≠b} (αβ)(b − a)  +  [αβ = 1]·(p − 1)|X ∩ Y|,

where, if `αβ = 1`, `(αβ)(u) := 1` for `u ≠ 0`; then `K(α, ᾱ) = −1`. Note
`Σ_{a≠b}(αβ)(b−a) = Σ_{x∈−X} f^Y_{αβ}(x)`, so the right side is itself a sum of profile data over `−X`.

*Proof.* For `a = b` the inner sum is `Σ_{u≠0}(αβ)(u) = (p−1)[αβ = 1]`. For `a ≠ b` put `c = b − a`
and `u = x + a`. Then `Σ_u α(u)β(u+c)`, and the substitution `u = cv` (`c ≠ 0`) gives
`α(c)β(c)Σ_v α(v)β(v+1) = (αβ)(c)K(α,β)`. If `αβ = 1`, the map `v ↦ v/(1+v)` is a bijection
`F_p∖{0,−1} → F_p∖{0,1}`, so `K(α,ᾱ) = Σ_{t≠0,1} α(t) = −1`. ∎

(Verifier A1: every pair `(φ^i, φ^j)`, `K ∈ {4, 8}`, random and overlapping `X, Y` at 14 primes,
exact in `Z[ζ_K]`, 2,814 checks. A2: `|J(φ^i,φ^j)|² = p` whenever `i, j, i+j ≢ 0`, and
`K(φ^i, φ^{−i}) = −1`, at those primes and at every certificate prime. That `|J| = √p` is standard
(UNVERIFIED CITATION: Ireland–Rosen, *A Classical Introduction to Modern Number Theory*, Ch. 8,
Thm 1; it is used only in the discussion, never in a certificate, where the exact values are
computed.)

**Corollary 1.3 (PROVED).** Let `|X| = m`.
1. `Σ_x |f^X_α(x)|² = m(p − m)` for every nontrivial `α`.
2. `Σ_x f^X_ψ(x)² = K(ψ,ψ) · σ_X` with `σ_X := Σ_{a≠a'} χ(a' − a) = Σ_{x∈−X} f^X_χ(x)`, and
   `Σ_x f^X_ψ f^X_χ = K(ψ,χ) · conj(τ_X)` with `τ_X := Σ_{x∈−X} f^X_ψ(x)`.
3. **Cliques.** If `X = C` is a clique of the Paley graph, then `σ_C = m(m−1)` and
   `|Σ_x f_ψ(x)²| = √p·m(m−1)`. Since `|Σ f_ψ²| ≤ Σ|f_ψ|² = m(p−m)`, this gives
   `√p(m−1) ≤ p − m` (the Delsarte–Hoffman / ϑ bound `m ≲ √p` again). At `m ≈ √(p/2)` the numbers
   `f_ψ(x)²` must have **phase coherence** `|Σ f_ψ²|/Σ|f_ψ|² = √p(m−1)/(p−m) ≈ 1/√2`.
4. **Two-sided.** `Σ_x f^A_α(x)·conj f^B_α(x) = p|A∩B| − mn`, and
   `Σ_x f^A_ψ(x) f^B_ψ(x) = K(ψ,ψ)·D_χ(A,B)` with `D_χ(A,B) = Σ_{a∈A,b∈B} χ(b−a)`. The *difference*
   structure of `A` against `B` enters; it is not constrained by `A + B ⊆ Q ∪ {0}`.

These are the new prime-specific relations the quartic (and octic) character brings: exact,
bilinear, and involving the Gauss/Jacobi constants `K(φ^i,φ^j)` of the prime. For a biclique the
right sides of items 2 and 4 contain the free quantities `σ_A`, `τ_A`, `D_χ(A,B)`. Only for cliques
is the right side forced (item 3).

### 1.3 Weil-type bounds

**Lemma 1.4 (PROVED given Weil's bound).** For `α` of order `k ≥ 3` and `|X| = m`:
`Σ_x |f^X_α(x)|⁴ ≤ (2m² − m)p + 3(m⁴ − 2m² + m)√p`.

*Proof.* Expand over `(a_1,…,a_4) ∈ X⁴`. The summand is `α((x+a_1)(x+a_2))·ᾱ((x+a_3)(x+a_4)) = α(P(x))`
with `P = (x+a_1)(x+a_2)(x+a_3)^{k−1}(x+a_4)^{k−1}`; both sides vanish at the roots. The multiplicity of a
root `−c` is `≡ #{i ≤ 2 : a_i = c} − #{i ≥ 3 : a_i = c} (mod k)`, a number in `[−2, 2]`. Since
`k ≥ 3`, `P` is a constant times a `k`-th power iff `{a_1,a_2} = {a_3,a_4}` as multisets; there are
`2m² − m` such tuples, each contributing at most `p`. For the others Weil's bound gives
`|Σ_x α(P(x))| ≤ (r−1)√p ≤ 3√p`. ∎

Citation status: Weil's bound in the form "`|Σ_x α(P(x))| ≤ (r−1)√p` for `α` of order `k`, `P`
not of the form `c·h^k`, `r` distinct roots" is UNVERIFIED CITATION (Iwaniec–Kowalski, Thm 11.23, as
referenced in the local copy of Fouvry–Shparlinski–Xi, Rev. Mat. Iberoam. 41 (2025), p. 1935,
`sources/bilinear-2025.txt`). The weaker bound actually stated there (Lemma 3.1, any nontrivial
character: `Σ_x |Σ_a α(x+a)|^{2r} ≤ (2r)^r(A^r p + A^{2r}p^{1/2})`) is CITED. Every certificate
below satisfies the stronger row of Lemma 1.4, hence the cited one a fortiori. (Verifier A6: 450
exact checks.) The **Krawtchouk–Weil rows** of the hybrid note (Prop. 2.3) extend verbatim to every
`φ^j`: `Σ_x e_t((φ^j(x+a))_{a∈X}) = Σ_{|U|=t} Σ_x φ^j(Π_{a∈U}(x+a))`, so
`|Σ_x e_t| ≤ C(m,t)(t−1)√p` for `t ≥ 3` (same citation). Here `e_t` is the `t`-th elementary
symmetric function and zero entries are omitted. At `t = m` this is the "parity" row
`|Σ_x φ^j(P_X(x))| ≤ (m−1)√p`.

### 1.4 Class rows

**Proposition 1.5 (PROVED given robust note Prop. 2.9).** Let `k | p − 1`, `c ∈ Z/k`, `|X| = m`,
`m + d_k − 1 ≤ p − 1`. For every integer `e` with `0 ≤ 2e ≤ min(m−1, d_k)`,

    (★)_{k,c,e}    Σ_{x: e^{(k,c)}_x ≤ e} (e + 1 − e^{(k,c)}_x)(m − (3e + e^{(k,c)}_x)/2 − δ_x) ≤ (e+1)(d_k − e),

and for `1 ≤ t ≤ m`: `Σ_x C(m−1−e_x, t−1)(m − e_x − δ_x) ≤ C(m,t)d_k` with `e_x = e^{(k,c)}_x`
(terms with a negative argument omitted).

*Proof.* Choose `h` with `ind h ≡ c (mod k)` and put `X' = h^{−1}X`, `x' = h^{−1}x`. Then
`(a' + x')^{d_k} = 1` iff `ind(a+x) ≡ c`, so `e^{(k,c)}_x` is the bad count of `x'` for `X'` in the
robust note's Prop. 2.9, which gives `(★)`. The subset rows follow from the robust note's Prop. 1.2
with `d` replaced by `d_k`, using the case `e = 0` of Prop. 2.9 (HP of order `k`) for every
`T ⊆ X`. ∎

For `k = 2` these are `(★)` for `X` and `νX`. (Verifier A3: exhaustively every `A ∋ 0` with
`|A| ≤ 6` at `p = 13, 17`, `|A| ≤ 4` at `p = 29, 37, 41`, for `k ∈ {4, 8}`, plus random sets up to
`p = 401`; 37,691 sets, 1,477,972 exact row checks.)

### 1.5 Pattern Hanson–Petridis

**Proposition 1.6 (PROVED).** Let `k | p−1`, `X = {a_1,…,a_m}`, `D_k = d_k + m − 1 ≤ p − 1`,
`c_i = Π_{l≠i}(a_i − a_l)^{−1}`, and `σ ∈ μ_k(F_p)^m`. Put
`N_σ = #{x ∉ −X : (x + a_i)^{d_k} = σ_i ∀i}` and
`l(σ) = min{l ≥ 0 : Σ_i c_i σ_i^{−1} a_i^l ≠ 0}`. Then `l(σ) ≤ m − 1`, with equality iff `σ` is
constant, and

    m·N_σ ≤ d_k + m − 1 − l(σ).

*Proof.* Let `F_σ = −1 + Σ_i c_iσ_i^{−1}(x + a_i)^{D_k}`. At a point `x_0` with pattern `σ`,
`(x_0+a_i)^{D_k−j} = σ_i(x_0+a_i)^{m−1−j}`. Hence
`F_σ^{(j)}(x_0)/(D_k)_j = Σ_i c_i(x_0+a_i)^{m−1−j} = 0` for `1 ≤ j ≤ m−1`, and `F_σ(x_0) = 0`, by
the partial-fraction identity `Σ_i c_i(x_0 + a_i)^s = [s = m−1]`, `s ≤ m−1` (robust note, Fact 0).
So `ord_{x_0} F_σ ≥ m`, using that `j!` is a unit because `D_k < p`. The coefficient of
`x^{D_k − l}` in `F_σ + 1` is `C(D_k,l)Σ_i c_iσ_i^{−1}a_i^l`, so
`deg F_σ = D_k − l(σ) ≥ d_k ≥ 1`. The linear map `y ↦ (Σ_i c_iy_ia_i^l)_{l ≤ m−2}` has rank `m−1`
and contains the constant vectors in its kernel (Fact 0), so its kernel is exactly the constants.
Count zeros. ∎

For constant `σ` this is HP of order `k`; a non-constant pattern costs at most `m − 1` in the
budget. (Verifier A4: 6,122 patterns, exact degrees and orders.)

### 1.6 The quadratic and the order-`k` Hanson–Petridis polynomials

**Proposition 1.7 (PROVED).** Let `k = 2^s | p−1`, `|X| = m ≥ 2`, `D_k = d_k + m − 1 ≤ p − 1`,
`W_k(x) = Σ_i c_i(x + a_i)^{D_k}`. Then:
1. `deg W_k = d_k` exactly (leading coefficient `C(D_k, m−1) ≢ 0`), and `Φ_k := W_k^{k/2} − 1 =
   Π_{ω∈μ_{k/2}}(W_k − ω)` has degree exactly `d = (k/2)d_k`. For `k = 2`, `Φ_2 = F_X` is the HP
   polynomial.
2. For `b ∉ −X` with `b + X ⊆ Q`, `ord_b Φ_k ≥ m` **iff** the column `(ind(b+a) mod k)_a` is
   constant; at `b ∈ −X` with a constant column, `ord_b Φ_k ≥ m − 1`.

*Proof.* (1) The coefficient of `x^{D_k−l}` is `C(D_k,l)Σc_ia_i^l`, which is `0` for `l < m−1` and
`C(D_k,m−1)` for `l = m−1` (Fact 0; `D_k < p`). (2) With `(b+a_i)^{d_k} = ω_i`,
`W_k^{(j)}(b)/(D_k)_j = Σ_i c_iω_i(b+a_i)^{m−1−j}` for `0 ≤ j ≤ m−1`. If all `ω_i = ω`, this is
`ω[j = 0]`, so `ord_b(W_k − ω) ≥ m`. Conversely `ord_b(W_k − ω) ≥ m` gives
`Σ_i c_i(ω_i − ω)(b+a_i)^t = 0` for `0 ≤ t ≤ m−1`, a Vandermonde system in the distinct
`b + a_i`, so `ω_i = ω` for all `i`. At most one factor `W_k − ω` vanishes at `b`. The case
`b ∈ −X` is the same computation with the zero term removed. ∎

**Consequence (relation (iv) of the brief).** For a complete biclique, HP of order `k` for the class
`ω` at `e = 0` reads `m·N_ω − r_ω ≤ d_k`, where `N_ω` is the number of columns constant `ω`.
Summed over the `k/2` classes this is `m·N_const − r_const ≤ d`, which is exactly the quadratic HP
inequality restricted to the constant columns. The halved degree `d_4 = d/2` is the quadratic
budget split between the two signs, and the product of the two quartic polynomials is a polynomial
of the same degree `d` as `F_A` that vanishes on a subset of `B`. The quartic family can beat HP at
`e = 0` only if more than `2/k` of the constant-column weight sits in one class. Theorem 2.1 makes
this precise for all `e`. (Verifier A5: exact degrees; orders at every complete point; 1,312 checks.)

## 2. The row and column statistics of the sign matrix (the optimisation of task 2)

Fix a complete biclique `(A, B)` and `k = 2^s ≥ 4` dividing `p − 1`. For an even class `c` and
`0 ≤ e ≤ min((m−1)/2, d_k/2)` put

    L^B_{k,c,e} := Σ_{b∈B, t^c_b ≤ e} (e + 1 − t^c_b)(m − (3e + t^c_b)/2 − δ_b),
    W^c_e := Σ_{b∈B, t^c_b ≤ e} (m − δ_b)            (so Σ_{b∈B}(m − δ_b) = mn − r),

and define the **scaling factor**
`Θ(A,B) := min_{k,c,e : L > 0} (e+1)(d_k − e)(mn − r)/(d·L^B_{k,c,e})` (`+∞` if every `L` vanishes).
`Θ` depends only on the column statistics of the sign matrices `M^{(k)}` and on `m, p`. `Θ^T`, from
the rows (roles of `A` and `B` exchanged), is defined the same way.

**Theorem 2.1 (PROVED; given robust note Thm 2.4 for (iv)).**
1. *(The biclique-restricted system.)* Every term of every class row `(★)_{k,c,e}` is `≥ 0` (as
   `m ≥ 2e+1`), so each row implies `L^B_{k,c,e} ≤ (e+1)(d_k − e)`. The moment identities imply
   nothing new on `B`: `|f_{φ^j}(b)| ≤ f_χ(b) = m − δ_b` on `B`, so every `B`-restricted second-moment
   bound follows from the quadratic identity `Σ_x f_χ² = m(p−m)`. The `B`-parts of the Weil
   fourth-moment rows (Lemma 1.4) are at most `n m⁴ ≤ 3m⁴√p` whenever `n ≤ 3√p`, which holds in the
   balanced window.
2. *(Exact value.)* Replace `B` by `λ` copies of each of its columns (`λ > 0` rational: "the same
   statistic, scaled"). HP holds iff `λ(mn − r) ≤ d`, and all class rows restricted to the scaled
   biclique hold iff `λ(mn−r) ≤ d·Θ`. So the largest biclique weight compatible with this sign
   statistic under HP and all class rows is **`d·min(1, Θ)`**, and the same holds with rows and `Θ^T`.
   The combined constraint is `mn − r ≤ d·min(1, Θ, Θ^T)`.
3. *(No gain for half-weight statistics.)* If `W^c_e ≤ (2/k)(mn − r)` for every even `c` and every
   `e`, then `Θ ≥ 1`. This holds in particular if `r = 0` and the multiset of columns of `M^{(k)}` is
   invariant under the rotation `c ↦ c + 2` (then the sets `{t^c_b ≤ e}`, `c` even, are disjoint and
   carry equal weight). For exactly balanced columns (`t^c_b ≥ m/2` for all `c`) every `L` vanishes
   and `Θ = ∞`.
4. *(Lopsided gain.)* If `η ≤ 1/8` and `W^c_{⌊ηm⌋} ≥ φ(mn − r)` for some class `c`, then

       mn − r ≤ (2/(kφ))·(1 − √(2η))^{−2}·d + m/φ,

   which is below `d` by a constant factor when `φ > (2/k)(1 − √(2η))^{−2}`.

*Proof.* (1) is stated in the text. (2) Scaling multiplies both `L` and `mn − r` by `λ`, while the
right sides `(e+1)(d_k−e)` and `d` are fixed. (3) Take `λ = d/(mn−r) ≥ 1` (HP for `(A,B)`). If
`L > 0` then `N_e := #{b : t^c_b ≤ e} ≥ 1`, and `(e+1−t)(m − (3e+t)/2 − δ) ≤ (e+1)(m − δ − 3e/2)`
for `t ≤ e`. So
`λL ≤ (e+1)[λW^c_e − (3e/2)λN_e] ≤ (e+1)[(2/k)d − (3e/2)] ≤ (e+1)(d_k − e)`, using
`(2/k)d = d_k` and `λN_e ≥ 1 ≥ 2/3`. For the rotation-invariant case: `t^c_b + t^{c'}_b ≥ m − δ_b`
for `c ≠ c'`, so with `r = 0` and `2e ≤ m−1` no column lies in two of the sets, and each carries at
most a `2/k` share. (4) Put `B' = {b : t^c_b ≤ ηm}`. Every `b ∈ B'` has at most `ηm` class-`c` bad
partners, and the robust note's Theorem 2.4 holds verbatim for order `k`, because its proof uses
only `(★)` restricted to `B'` (here `(★)_{k,c,e}`, valid for `2e ≤ min(m−1, d_k)`, which covers the
`e ≤ m/4` it uses when `m ≤ 2d_k`). It gives `m|B'| − r' ≤ (1−√(2η))^{−2}d_k + m`. The left side is
`W^c_{⌊ηm⌋} ≥ φ(mn−r)`. ∎

(Verifier A7: part 3 on 3,000 random abstract column statistics, 134,933 exact row checks. Section B
computes `Θ` exactly for 813 actual bicliques and checks the `B`-restricted rows (51,116) and the
bound of part 4 (308).)

**Answer to task 2.** In terms of the sign statistics, the strongest consequence of quadratic HP,
quartic `(★)_4` (and octic `(★)_8`), and the quartic moments is exactly
`mn − r ≤ d·min(1, Θ, Θ^T)`. It forces `mn ≤ (1/2 − c)p` **only if the sign matrix is lopsided**:
more than a `2/k` share of the weight in columns (or rows) close to one constant class. It forces
nothing for half-weight, flip-invariant or balanced sign matrices. Whether near-tight balanced
bicliques must have lopsided sign matrices is not decided by the family. Theorem 3.1 shows the family
permits exactly balanced ones at Hanson–Petridis tightness, and the actual tight cliques at
`p = 13, 41` have them (§4).

The robust machinery also pushes the sign matrix *away* from being constant, not towards it. The
order-`k` bias theorem (sharpen note, Thm 5.1, PROVED there given robust Prop. 2.9) gives, at
`mn − r ≈ d`, that every quartic sign class misses at least a fraction `θ ≈ 0.052` of the entries.

## 3. Exact certificates for the full constraint family

### 3.1 The family `𝓕_K(p; m, n)`

`K = 8` if `p ≡ 1 (mod 8)`, else `K = 4`. The unknowns are a **joint profile**: nonnegative rational
weights on pairs (order-`K` type w.r.t. `A`, order-`K` type w.r.t. `B`), one pair per point `x`,
together with an explicit `m × n` matrix `M` with entries in the even classes. The rows are:

- **One-sided, for `X = A` (size `m`) and `X = B` (size `n`)** on the marginal profile:
  total weight `p`; weight `|X|` on `−X`; the complete points of `X` (`δ = 0`, `e^{(2,0)} = 0`) are
  exactly the partner cell and `r = 0`. All first moments `Σ_x f_{φ^j} = 0` and all Jacobi identities
  of Lemma 1.2 for every pair `1 ≤ i ≤ j ≤ K−1` (this includes `Σ f_χ² = m(p−m)`). The class rows
  `(★)_{k,c,e}` and subset HP for every `k | K` (`k = 2, 4, 8`) and every class (Prop. 1.5).
  Pencil (pencil note Thm 3.1) for `X` and `νX`; derivative HP (balanced note Prop. 1.4). Weil:
  quadratic moments of orders 4, 6, 8 and quadratic Weil counts (sharpen note Lemmas 3.1, 4.2,
  balanced note Thm 3.5), the fourth moment of Lemma 1.4 for every `φ^j` of order `≥ 3`, and the
  Krawtchouk–Weil rows for every `φ^j`, `1 ≤ j ≤ K/2`, and every `3 ≤ t ≤ |X|`. Everywhere `√p` is
  replaced by `⌊√p⌋` except in the KW rows, which are compared with `√p` exactly (squared).
- **Two-sided:** all two-sided Jacobi identities (Lemma 1.2 with `X = A`, `Y = B`, all `1 ≤ i, j ≤
  K−1`), with `D_l := Σ_{x∈−A} f^B_{φ^l}(x)` and the consistency
  `Σ_{x∈−B} f^A_{φ^l}(x) = φ^l(−1)D_l`. **Union rows** for `U = A ∪ B` (size `m + n`, disjoint):
  `(★)_{k,c,e}` for every `k | K`, class and `e`; averaged subset HP over `T_A ∪ T_B` for
  `(t_A, t_B)` on a grid (all pairs when `m, n ≤ 30`, otherwise
  `{0,…,6, 8, 10, m/8, m/4, m/2, 3m/4, m}` × the same for `n`), namely
  `Σ_x[(t_A+t_B)C(m−e^A_x,t_A)C(n−e^B_x,t_B) − δ^A_x C(m−1−e^A_x,t_A−1)C(n−e^B_x,t_B) −
  δ^B_x C(m−e^A_x,t_A)C(n−1−e^B_x,t_B−1)] ≤ C(m,t_A)C(n,t_B)d_k`; pencil and derivative HP for `U`;
  Krawtchouk–Weil rows for `U` (every `φ^j`, `3 ≤ t ≤ m+n`).
- **Sign matrix:** the columns of `M` (as types) are exactly the `A`-types of the points of `B`, and
  its rows are exactly the `B`-types of the points of `A`. Every column and row pattern satisfies
  pattern HP (Prop. 1.6) for every `k | K`, `k ≥ 4`. Every submatrix `A' × B'` constant in one class
  of order `k` satisfies `|A'||B'| ≤ d_k` (HP of order `k` for the sub-biclique), certified by
  (max class count in a column) × (max class count in a row) `≤ d_k`.

Two restrictions that are *not* theorems make the certificates stronger: every non-complete point
has quadratic bad count in `[0.4m, 0.6m]`, and `|A ∩ B| = 0`, `A ∩ (−A) = B ∩ (−B) = ∅`.

**Not included** (explicitly outside `𝓕`): identities coupling the types at `x` and `−x` (for
example `Σ_x f^A_ψ(x)f^A_ψ(−x)`, which involves `A + A`); Kalmynin/Yip–Yoo relations, which need
`g = 0` or `g ≤ max(m,n) − 2` (balanced note Remark 1.2′); robust (`e ≥ 1`) pattern rows for
non-constant patterns (whether the Hankel minor of `F_σ` is nonzero is not established) and robust
rows for individual sub-rectangles; Weil counts of refined types; characters of order `16` and
higher.

### 3.2 Construction

Three linear programmes (HiGHS). (a) Quadratic levels for each side: the balanced-note rows plus
the quadratic KW rows, which were added as cutting planes and are essential: without them the
parity row `t = m` fails. (b) Refinement of each level into order-`K` types, using seeded
multinomial class splits as candidates and rotation orbits as variables (so every row that is odd
under rotation vanishes identically), with the Jacobi identities as equalities and the complex KW
rows as cutting planes. (c) The joint coupling of the two sides over the cells
`B, A, −A, −B, rest`, with the two-sided identities and the union rows. Each LP maximises a common
relative slack. Its solution is rounded to denominator `10⁶`, and every equality is then repaired
exactly by rational Gauss–Jordan on a pivot set. The sign matrix starts from
`M_{ab} = 2((a+b) mod K/2)` and is randomised by seeded `2 × 2` switches, which preserve every row
and column class count. Only the final objects matter: **every row of §3.1 is re-evaluated on the
final joint profile and matrix**, in exact rational / `Z[ζ_K]` / Python-integer arithmetic.

### 3.3 Result

**Theorem 3.1 (PROVED by exact computation, except as stated for the middle KW rows).** For each row
of Table 3.1 there is a joint profile and an explicit sign matrix satisfying every row of
`𝓕_K(p; m, n)`. Consequently no nonnegative combination of these rows implies
`mn ≤ (1/2 − c)p` for any `c > 1/2 − mn/p` at these parameters. Precision statement: the KW rows with
`12 < t < M − 12` (`M = |X|` or `m+n`) are evaluated in double precision by a numerically stable
normalised recurrence (every step is a convex combination), with an additive allowance of
`p·10⁻⁹` on the normalised value. All other rows, including every KW row with `t ≤ 12` or
`t ≥ M − 12`, are exact. The largest normalised ratio among the floating rows is `≤ 0.07`.

**Table 3.1** (`g = d − mn`; mean squares are over the points off `B`, side `A`; slack = the smallest
relative slack over sides `A`, `B` and the union, where `1` means the row is identically `0` on the
support).

| p | K | m | n | g | mn/p | types A/B/joint | mean `f_χ²` off B | mean `|f_ψ|²` off B | slack `(★)_4, (★)_8` | slack union `(★)_2` | max KW ratio (exact) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1009 | 8 | 20 | 25 | 4 | 0.4955 | 45/56/156 | 9.94 | 20.10 | 0.997, 1 | 0.52 | 0.75 |
| 1013 | 4 | 24 | 21 | 2 | 0.4975 | 17/20/48 | 11.73 | 23.93 | 0.991, — | 0.60 | 0.99 |
| 10009 | 8 | 72 | 69 | 36 | 0.4964 | 45/48/148 | 35.99 | 71.98 | 0.999, 1 | 0.56 | 0.86 |
| 10037 | 4 | 70 | 71 | 48 | 0.4952 | 19/16/46 | 35.10 | 70.01 | 0.9996, — | 0.52 | 0.98 |
| 40009 | 8 | 136 | 147 | 12 | 0.4997 | 49/52/156 | 67.83 | 136.04 | 0.9999, 1 | 0.49 | 0.99 |
| 100049 | 8 | 224 | 223 | 72 | 0.4993 | 49/48/152 | 111.91 | 224.00 | 1, 1 | 0.49 | 0.99 |
| 100069 | 4 | 224 | 223 | 82 | 0.4992 | 17/16/44 | 111.93 | 224.00 | 1, — | 0.53 | 0.99 |
| 1000033 | 8 | 708 | 706 | 168 | 0.4998 | 49/52/156 | 353.87 | 708.00 | 1, 1 | 0.58 | 0.99 |

The binding rows are HP itself (slack `g/d`) and the one-sided parity-type KW rows (ratio up to
`0.99`). The quartic and octic `(★)` rows are *identically zero* on side `A` and for the union:
every column of the sign matrix is exactly balanced (`m/2` of each sign; `m/4` of each octic class),
and every other point has all class counts below `(m+1)/2`. On side `B`, when `n` is odd, each row of
`M` has a one-entry majority, which uses at most 1% of a `(★)_4` budget. The subset rows of order 4 and 8 sit at about `1/2` and `1/4` of their
budgets, as for random data. The quartic second moment is met with `|f_ψ|²` of mean exactly `m`
off `B` (a binomial-like spread), while `f_χ²` has mean `≈ m/2` there (the sub-binomial spread that
HP-tightness forces, balanced note §3.3). All column patterns and row patterns of the sign matrices
are distinct.

**Cliques (one-sided; PROVED by exact computation, same precision statement).** For a clique
`C = −A` (`r = m`), the complete cell is `−A` itself, and the Jacobi identity has the forced right
side `K(ψ,ψ)m(m−1)`. With a circulant sign matrix (a regular tournament when `ψ(−1) = −1`, i.e.
`p ≡ 5 mod 8`; a symmetric `±1` matrix with balanced rows when `ψ(−1) = 1`), there are rational
one-sided profiles satisfying every one-sided row of §3.1 for `K = 4`, including the forced
coherence, which the certificates realise exactly:

| p | ψ(−1) | m | m(m−1)/d | m²/p | coherence `|Σf_ψ²|/Σ|f_ψ|²` = `√p(m−1)/(p−m)` | max KW ratio |
|---|---|---|---|---|---|---|
| 1013 | −1 | 21 | 0.830 | 0.435 | 0.642 | 0.52 |
| 1009 | +1 | 21 | 0.833 | 0.437 | 0.643 | 0.52 |
| 10037 | −1 | 71 | 0.990 | 0.502 | 0.704 | 0.98 |
| 10009 | +1 | 69 | 0.938 | 0.476 | 0.684 | 0.93 |
| 100069 | −1 | 223 | 0.989 | 0.497 | 0.703 | 0.98 |
| 100049 | +1 | 221 | 0.972 | 0.488 | 0.697 | 0.96 |
| 1013 | −1 | 23 | 1.000 | 0.522 | 0.707 | 0.77 |

The last row has `m(m−1) = d`, i.e. `g = 0`: an exactly Hanson–Petridis-tight 23-clique profile at
`p = 1013` (actual `ω(G_1013)` is far smaller). `g = 0` is the critical case, where further identities
(Kalmynin; balanced note Remark 1.2′) are known and lie outside `𝓕`, so this row is recorded only as
evidence about `𝓕`. The
two-sided coupling of a clique is the reflection `x ↔ −x`, which is not imposed. So the clique
certificates are one-sided only.

**REFUTED:** "if the sign matrix is balanced, the quartic second moment conflicts with near-tightness
of the quadratic HP". Witness: every certificate of Table 3.1, where every column of the sign matrix
is exactly balanced, `Σ_x |f_ψ|² = m(p−m)` holds exactly, and `mn/p` reaches `0.4998`.

## 4. Actual bicliques (task 2, test)

Section B of the verifier takes every complete biclique with `p ≡ 1 (mod 4)`, `p ≤ 1500` from
`results/sigma_biclique_2026_09_05_search.json` (maximum cliques `(C, −C)`, balanced witnesses
`(A, B(A))`, sum-cliques `(A, A)`, profile witnesses `k = 3..8`, greedy pairs) and from
`results/stepanov_balanced_2026_09_29_data.json` (the `P_K` witnesses): 813 configurations after
deduplication. On each it re-checks exactly, as theorems: every class row and subset row (subset
rows for `t ≤ 8`) for `A` and for `B` at every order `k | K`; all one- and two-sided Jacobi identities;
Weil fourth moments; KW rows for `t ≤ 5` and `t ≥ M − 5` (vectorised Newton identities in `Z[ζ_K]`,
cross-checked against brute force in A8); pattern HP; and Theorem 2.1(1), (4). In total 380,018
checks, 0 failures. It also computes `Θ` exactly.

**Table 4.1** (`Θ = min(Θ, Θ^T)`; "col²·m" is `m·Σ_b(Σ_a M_{ab})²/(nm²)`, which is `≈ 1` for
binomial-like columns and `= m` for constant ones).

| kind | configs | min Θ | median Θ | Θ < 1 | max `|S_ψ|/(mn−r)` | mean col²·m |
|---|---|---|---|---|---|---|
| maximum cliques | 116 | 0.75 | 13.1 | 1 | 0.53 | 0.81 |
| balanced witnesses | 57 | 0.75 | 10.9 | 1 | 0.60 | 0.77 |
| sum-cliques | 116 | 0.75 | 13.9 | 2 | 0.85 | 0.91 |
| profile witnesses | 421 | 0.75 | 5.2 | 1 | 0.52 | 0.92 |
| greedy pairs | 79 | 1.73 | 10.4 | 0 | 0.55 | 0.89 |
| `P_K` witnesses | 24 | 3.27 | 4.4 | 0 | 0.14 | 0.96 |

- The five configurations with `Θ < 1` are three `3 × 3` configurations and a `4 × 2` profile at
  `p = 17` (`Θ = 0.75`), and the sum-clique `A = B`, `|A| = 11`, at `p = 593` (`Θ = 0.950`, at
  `k = 4`, `e = 1`). The latter is quartic-structured: `A = −A`, `A ∖ {0} ⊆ Q_4`, and 102 of its
  110 nonzero sums are fourth powers, so `S_ψ/(mn−r) = 0.855`. It sits near the *quartic* HP bound
  (`m² = 121` against `d_4 + m = 159`) and far from the quadratic one (`mn/p = 0.20`). This is the
  **witness** that the quartic rows do see some actual bicliques. The one non-tiny example lives
  inside `Q_4` and is capped by the quartic budget `d/2`; whether this is typical is HEURISTIC.
- **The exactly tight configurations** (`g = 0`; balanced note §1.1): the `3 × 3` clique at `p = 13`
  and the `5 × 5` clique `{0,1,2,10,33}` at `p = 41` have **exactly balanced** quartic sign matrices
  (each column: one `+1`, one `−1`, one `0` at `p = 13`; two `+1`, two `−1`, one `0` at `p = 41`).
  Every quartic and octic class row restricted to `B` vanishes, so `Θ = ∞`. Actual
  Hanson–Petridis-tight bicliques therefore coexist with balanced sign matrices, and no argument of
  the form "tight ⇒ lopsided" is possible.
- For `p ≡ 5 (mod 8)` every clique has `S_ψ(C, −C) = 0` exactly: `ψ(−1) = −1` makes `M` a tournament.
- The crude half-weight condition of Theorem 2.1(3) fails in 642 of the 813 configurations, but only
  at `e` close to `(m−1)/2`, where the weights are tiny. `Θ`, not the crude condition, is the right
  measure.

**The tight family `A = {0,1}`** (`n = (p+3)/4`, `g = 0` for every `p ≡ 1 mod 4`; Section D, all
`p ≤ 3000`). The four quartic column patterns `(ψ(b), ψ(b+1))` occur about `n/4` times each, for
example `(68, 61, 61, 61)` at `p = 1009` and `(171, 190, 190, 190)` at `p = 2969`. The constant-column
weight is about half of the threshold `(2/k)(mn−r)` for `k = 4` and about a quarter for `k = 8`. It
reaches the threshold only at `p = 5, 17` (826 exact checks).

## 5. The obstruction

1. **On the biclique, the order-`k` family is a lopsidedness test and nothing more** (Theorem 2.1).
   HP of order `k` summed over the classes is HP restricted to constant columns (Prop. 1.7). The
   robust class rows see only columns within Hamming distance `e < m/2` of a constant column. Every
   second-moment and Weil fourth-moment row restricted to `B` is implied by the quadratic ones in the
   balanced window.
2. **Off the biclique, the order-`k` refinement is free** (Theorem 3.1). The only exact couplings the
   higher characters add are the Jacobi identities (Lemma 1.2). Their right sides for bicliques are
   the free internal and cross sums `σ_A, τ_A, D_l(A,B)`. In the certificates, multinomial class
   splits of the quadratic profiles meet them with every `(★)_4, (★)_8` row (essentially) zero.
3. **The expected behaviour of `ψ` is exactly the invisible case.**

**Proposition 5.1 (CONDITIONAL; complete proof).** Let `k = 2^s ≥ 4`, `0 < η < 1/2`, `n_0 ≥ 1`,
`η' = η + 1/m`, and let `(A, B)` be a complete biclique with `m ≥ 2` such that
`|Σ_{a∈A,b∈B'} φ_k(a+b)| ≤ ηm|B'|` for every `B' ⊆ B` with `|B'| ≥ n_0`. If

    m n_0/(mn − r) + η'(η'm + 1)/((m − 1)(1 − η')) + m/(2d) ≤ 2/k,

then `Θ ≥ 1` for the columns of `M^{(k)}`: the class rows of order `k` never cut the column statistic
of `(A, B)` below Hanson–Petridis tightness. In particular, assume the order-`k` analogue of the
two-set conjecture: for every `ε` there are `δ, p_0` with `|Σ_{A'×B'}φ_k(a+b)| ≤ p^{−δ}|A'||B'|`
whenever `|A'|, |B'| ≥ p^ε`. Take `η = p^{−δ}` and `n_0 = p^ε`. Then for all large `p`, every complete
biclique with `p^ε ≤ m ≤ p^{1−ε}` and `n ≥ p^{2ε}` has `Θ ≥ 1`, and likewise `Θ^T ≥ 1` with the roles
exchanged.

*Proof.* Let `c` be an even class with value `ω`. A column with `t^c_b ≤ e` has
`Re(ω̄Σ_aφ_k(a+b)) ≥ (m − δ_b − t_b) − t_b ≥ m − 1 − 2e`, because entries of class `c` contribute `1`,
the others at least `−1`, and a zero entry `0`. So if `B'_e := {b : t^c_b ≤ e}` has `|B'_e| ≥ n_0`, then
`m − 1 − 2e ≤ ηm`. Hence `|B'_e| < n_0` for `e < e_0 := ((1−η)m − 1)/2`. Put `N = mn − r` and scale by
`λ = d/N ≥ 1` as in Theorem 2.1(2); we must show `λL^B_{k,c,e} ≤ (e+1)(d_k − e)`.

If `e < e_0`, every weight is at most `(e+1)m`, so `λL ≤ (e+1)mλn_0`.

If `e_0 ≤ e ≤ (m−1)/2`, the fewer than `n_0` columns with `t < e_0` contribute at most `(e+1)m` each.
A column with `e_0 ≤ t ≤ e` contributes at most `((η'm+1)/2)·η'm`, since `e + 1 − t ≤ (m+1)/2 − e_0
= (η'm + 1)/2` and `m − (3e+t)/2 − δ ≤ m − 2t ≤ m − 2e_0 = η'm`. There are at most `n` of them, and
`λn = dn/N ≤ d/(m−1)` because `r ≤ n`. With `e + 1 ≥ (1 − η')m/2`,

    λL/(e+1) ≤ mλn_0 + (d/(m−1))·η'(η'm+1)/(1−η').

In both cases `λL ≤ (e+1)(d_k − e)` follows from `mλn_0 + dη'(η'm+1)/((m−1)(1−η')) + m/2 ≤ d_k`.
Divided by `d` (with `d_k = 2d/k`, `λ = d/N`, `e ≤ m/2`), this is the displayed hypothesis. For the
asymptotic statement: `mn_0/N ≤ n_0/(n(1 − 1/m)) → 0` when `n ≥ p^{2ε}`; `η' → 0` since `m → ∞`;
the middle term is `O(η'² + η'/m) → 0`; and `m/(2d) → 0` since `m ≤ p^{1−ε}`. ∎

So a proof of a constant-factor improvement of Hanson–Petridis *through the quartic sign matrix*
needs an input forcing `ψ(a+b)` to be **biased** on `A × B'` for large `B' ⊆ B`. That is the negation
of the (expected) order-4 analogue of the conjecture this programme aims at. By Proposition 5.1,
combined with the free off-`B` refinement of Theorem 3.1, the order-4 and order-8 characters cannot
supply the missing prime-specific input for balanced bicliques unless the Paley-type conjecture for
`ψ` fails. What they do supply is Lemma 1.2's identities. These are genuinely prime-specific, but
for bicliques they are linear in free profile data, and for cliques the forced coherence is
satisfiable.

## 6. Open obligations

1. **The target** `mn ≤ (1/2 − c)p` for balanced complete bicliques, and the Paley clique constant
   `1/√2`: OPEN. Nothing here improves Hanson–Petridis.
2. **Reflection identities.** Identities coupling the types at `x` and `−x` (for bicliques, via
   `A + A`, `B + B`; for cliques, the entire two-sided structure) are not in `𝓕`. An LP over
   `(type(x), type(−x))` pairs would test them. Not attempted.
3. **Robust pattern rows.** For a non-constant pattern `σ`, is the `(e+1) × (e+1)` Hankel minor of the
   derivatives of `F_σ` (Prop. 1.6) nonzero, and what is its exact degree? If so, it gives
   `(★)`-type rows for column patterns, which random-like sign matrices would still satisfy.
4. **Characters of order 16+** (`p ≡ 1 mod 16`): the same construction (multinomial class splits)
   should extend. Not run.
5. **Floating rows.** The middle KW rows (`12 < t < M−12`) are verified in floating point with an
   allowance (§3.3). An exact check needs `O(M²)` big-integer work per type; feasible for `p ≤ 10⁵`,
   not done.
6. **Citations.** The exact form `(r−1)√p` of Weil's bound for characters of order `k`, and `|J| = √p`,
   are UNVERIFIED CITATIONS (the weaker Fouvry–Shparlinski–Xi Lemma 3.1 is CITED and suffices for the
   certificates, which satisfy the stronger rows).
7. **Independent review** of Lemma 1.2, Propositions 1.5–1.7, Theorem 2.1 and Proposition 5.1, and of
   the certificate verifier (`verify_side`, `verify_joint`, `verify_matrix` and the clique section).

## 7. Verification summary

| section | content | checks | failures |
|---|---|---|---|
| A1–A2 | Jacobi identities; `|J|² = p`, `K(α,ᾱ) = −1` | 3,283 | 0 |
| A3 | class rows `(★)_k`, subset HP_k (exhaustive small `p`, random) | 1,477,972 | 0 |
| A4 | pattern HP: exact degrees and orders | 26,265 | 0 |
| A5 | `W_k^{k/2} − 1`: degrees, orders iff constant column | 1,312 | 0 |
| A6 | Weil fourth moments (order 4, 8) | 450 | 0 |
| A7 | Theorem 2.1(3) on random statistics | 134,933 | 0 |
| A8 | exact KW engine vs brute force over subsets; float vs exact | 206 | 0 |
| B | 813 actual bicliques: all rows as theorems, `Θ` | 380,018 | 0 |
| C | 8 biclique and 7 clique certificates, every row | 130,936 | 0 |
| D | `A = {0,1}`, all `p ≤ 3000` | 1,863 | 0 |
| **total** | | **2,157,238** | **0** |
