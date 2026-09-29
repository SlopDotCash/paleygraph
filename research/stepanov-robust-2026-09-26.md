# Stepanov wave, worker `robust`: Robust Hanson–Petridis via Hankel minors (2026-09-26)

**Status:** PROVED: (i) a weighted Stepanov inequality (Theorem 2.1) obtained from the
`(e+1)×(e+1)` Hankel minor of the normalised derivatives of the Hanson–Petridis polynomial;
(ii) Robust Hanson–Petridis RHP(f,g) with `f(η) = (1−√(2η))^{−2}` for `0 ≤ η ≤ 1/8` and
`g(m) = m`, for every `m ≤ (p+1)/2` (Theorem 2.4); (iii) the target consequence
**unconditionally**: if `|A||B| ≥ (1/2+κ)p`, `0<κ≤3/2`, `p ≥ 11`, then
`|S(A,B)| ≤ [1 − (1−(1+2κ)^{−1/2})² + (√((1/2+κ)p)+1)/(2(p−1))]·|A||B|` (Theorem 2.5);
(iv) the baselines (HP with the exact derivative formula; the cheap subset-averaging bound and a
lower bound `≥ e·e^{1−e/m}` on its constant; the general implication RHP ⇒ bias saving);
(v) a no-go: the constraints "HP for every subset `T ⊆ A`" plus set sizes and universe size do
**not** imply RHP — exact set systems with `m·#{≥ m−1 coverage} = 2(m/(m−1))²·d` (Theorem 4.1).
REFUTED: "subset-HP constraints imply RHP" (witness: Theorem 4.1, `m = 43`, `2d+1` prime; see
the results file). Literature (§3, CITED sources): nothing we could fetch states a robust/defect
form of Hanson–Petridis or a constant-bias bound in the window `p/2 < |A||B| < p`; novelty is
**not** asserted beyond that search (web search hit a session limit before Shkredov 2014 and
Heath-Brown–Konyagin 2000 could be fetched). OPEN: the optimal `f` (we get `1+O(√η)`, not `1+O(η)`), the
optimal `η(κ)` (we get `≍ κ²`; two-element examples cap it at `≤ 2κ/(1+2κ)`), and any power
saving. The threshold `1/2` in Theorem 2.5 is sharp (Remark 2.7). This note has **not** been
independently refereed; the verifier re-checks every step exactly.

Verifier: `experiments/stepanov_robust_2026_09_26.py` → `results/stepanov_robust_2026_09_26.json`
(standard library + numpy + sympy; about 30 s CPU on an idle machine). Local copies of the
fetched sources: `sources/stepanov-robust-krattenthaler-math9902004.txt`,
`sources/stepanov-robust-kalmynin-2504.10202.txt`,
`sources/stepanov-robust-rudnev-tyrrell-2607.24270.txt`, and the pre-existing
`sources/sigma-hanson-petridis-1905.09134.txt`.

---

## 0. Notation and one standard fact

`p` odd prime, `χ` Legendre with `χ(0)=0`, `d=(p−1)/2`, `Q` the nonzero squares.
`A={a_1,…,a_m} ⊆ F_p` with `1 ≤ m ≤ (p+1)/2`; `D = d+m−1 ≤ p−1`.
`c_k = Π_{l≠k}(a_k−a_l)^{−1}`, `G(x)=Σ_k c_k(x+a_k)^{m−1}`, `F(x) = −1+Σ_k c_k (x+a_k)^D`.
For `b ∈ F_p`: `y_k = b+a_k`, `E(b) = {k : χ(y_k) = −1}`, `e_b = |E(b)|`, `δ_b = [b ∈ −A]`.
For `B ⊆ F_p`: `n=|B|`, `r = |B∩(−A)|`, `N_− = Σ_{b∈B} e_b` (number of pairs with `χ(a+b)=−1`),
`S = S(A,B) = Σ_{a,b} χ(a+b) = mn − r − 2N_−`.

**Fact 0 (partial fractions).** `Σ_k c_k a_k^j = 0` for `0 ≤ j ≤ m−2` and `= 1` for `j = m−1`
(the leading coefficient of the Lagrange interpolant of `x^j` on `A`). Hence `G ≡ 1`, and
translating, `Σ_k c_k (b+a_k)^j = [j=m−1]` for every `b` and `0 ≤ j ≤ m−1`.

All derivatives below are of polynomials of degree `≤ D ≤ p−1`, so `j!` and `(D)_j`
(`j ≤ D`) are units mod `p`, Taylor's formula `P(x+h) = Σ_j P^{(j)}(x)h^j/j!` holds, and
`ord_b P ≥ M` iff `P^{(j)}(b)=0` for `j<M`. **This is where `D ≤ p−1` (i.e. the prime field)
enters; over `F_{p²}` with `A = F_p` one has `D ≥ p` and `(D)_j ≡ 0` already for some
`j ≤ m−1` (verifier §F).**

## 1. Baselines

### 1.1 HP and the local form of `F` (PROVED)

**Lemma 1.1.** (a) `deg F = d` exactly, with leading coefficient `C(D,m−1) ≢ 0`.
(b) For `b ∉ −A` and `0 ≤ j ≤ m−1`:
`F^{(j)}(b) = −2(D)_j Σ_{k∈E(b)} c_k (b+a_k)^{m−1−j}`.
(c) For `b = −a_{k_0}` the same holds for `0 ≤ j ≤ m−2`, and
`F^{(m−1)}(b) = (D)_{m−1}(−c_{k_0} − 2Σ_{k∈E(b)} c_k)` (for `m ≥ 2`).
(d) Let `R_E(x) = Σ_{k∈E} c_k (x+a_k)^D`. Then `F − 2R_{E(b)}` vanishes at `b` to order
`≥ m − δ_b`.

*Proof.* (a) The coefficient of `x^{D−l}` in `F+1` is `C(D,l)Σ_k c_k a_k^l`, zero for `l<m−1` and
`C(D,m−1)` for `l = m−1` (Fact 0), and `C(D,m−1) ≢ 0` as `D<p`. (b) For `j ≥ 1`,
`F^{(j)}(b) = (D)_j Σ_k c_k y_k^{D−j} = (D)_j Σ_k c_k χ(y_k) y_k^{m−1−j}` (as `y^{D−j} = y^d y^{m−1−j}`
and `y^d = χ(y)`); write `χ(y_k) = 1 − 2[k∈E]` and use Fact 0 (`Σ_k c_k y_k^{m−1−j} = 0` for
`1 ≤ j ≤ m−1`). For `j = 0`, `F(b) = −1 + Σ c_k y_k^{m−1} − 2Σ_E(…) = −2Σ_E c_k y_k^{m−1}`.
(c) Same computation with the `k_0` term deleted (`y_{k_0}^{D−j} = 0` for `j<D`); the deleted term
`c_{k_0} y_{k_0}^{m−1−j}` is `0` except at `j = m−1`, where it is `c_{k_0}`. (d) For `k∈E(b)`,
`R_E^{(j)}(b) = (D)_j c_k χ(y_k) y_k^{m−1−j}` summed, i.e. `2R_E^{(j)}(b) = −2(D)_jΣ_E c_k y_k^{m−1−j}`,
which equals `F^{(j)}(b)` for `j ≤ m−1−δ_b` by (b),(c). ∎

HP follows: for `B ⊆ {b: e_b = 0}`, `F` has order `≥ m − δ_b` at each `b ∈ B`, so
`mn − r ≤ d`. By (b), a single bad partner already gives `F(b) = −2c_k y_k^{m−1} ≠ 0`.
Verifier §A: `deg F = d` and (b),(c) at every `b` and every `j ≤ m−1` for all `A ∋ 0` with
`p ≤ 13`, and for structured/random `A` with `p ≤ 80`.

### 1.2 The cheap robust bound from HP on subsets (PROVED)

**Proposition 1.2.** For `1 ≤ t ≤ m`:
`Σ_{b∈F_p} C(m−1−e_b, t−1)·(m − e_b − δ_b) ≤ C(m,t)·d`. Consequently, with
`N_e = #{b : e_b ≤ e}` and `r_e = #{b ∈ −A : e_b ≤ e}`,

    (m−e)·N_e − r_e ≤ K(m,e)·d,     K(m,e) := min_{1≤t≤m−e} C(m,t)/C(m−1−e,t−1).

*Proof.* For `T ⊆ A`, `|T| = t`, let `B_T = {b : T ∩ E(b) = ∅}`; `T + B_T ⊆ Q∪{0}`, so HP gives
`t|B_T| − |B_T∩(−T)| ≤ d`. Sum over the `C(m,t)` sets `T`: `b` lies in `C(m−e_b,t)` of the `B_T`,
and `b = −a_{k_0}` lies in `B_T ∩ (−T)` for `C(m−1−e_b,t−1)` of them (`T ∋ a_{k_0}`, `T` avoids
`E(b) ∌ k_0`). Use `t·C(m−e_b,t) = (m−e_b)C(m−1−e_b,t−1)`. For the consequence bound each
term with `e_b ≤ e` below by `C(m−1−e,t−1)(m−e−δ_b)`. ∎

**Proposition 1.3 (the cheap constant grows linearly in `e`).** `K(m,e) ≥ e·exp(1−e/m)` for
`1 ≤ e < m`. *Proof.* For `1 ≤ t ≤ m−e`,
`C(m,t)/C(m−1−e,t−1) = (m/t)·C(m−1,t−1)/C(m−1−e,t−1)` and
`C(m−1,t−1)/C(m−1−e,t−1) = Π_{i<e}(m−1−i)/(m−t−i) ≥ ((m−1)/(m−t))^e` (each factor is
`≥ (m−1)/(m−t)` because `t ≥ 1`). With `y = (t−1)/(m−1)`: `(m−1)/(m−t) = 1/(1−y)` and
`m/t ≥ 1/(y+1/m)`, so the ratio is `≥ exp(ey)/(y+1/m)`; putting `z = ey` this is
`e·exp(z)/(z+e/m) ≥ e·exp(1−e/m)` (minimum over `z ≥ 0` at `z = 1−e/m`). ∎
For fixed `e` and `m→∞`, `K(m,e) → (e+1)(1+1/e)^e ∈ [2(e+1), e_{Euler}(e+1))`
(exact values: `K(40,1) = 3.9`, `K(100,1) = 3.96`, `K(100,5) = 14.48`, `K(100,10) = 26.93`,
`K(1000,10) = 28.37`; results file, `cheap_constants`). So subset averaging gives
`m N_e ≲ C·e·d`: it depends on `e`, not on `e/m`, and gives nothing for the bias application,
where `e = ηm → ∞`. Verifier §B checks the identity-level inequality for every `t` and the
lower bound for all `2 ≤ m < 120`.

### 1.3 The conditional consequence (PROVED, conditional on RHP)

**RHP(f,g).** For all `p`, all `A` with `m ≤ √p`, and all `B` with `e_b ≤ e` (`b∈B`):
`mn − r ≤ f(e/m)·d + g(m)`.

**Proposition 1.4.** Assume RHP(f,g) where `f*(η) := sup_{η'≤η} f(η') → 1` as `η → 0` and
`G(p) := max_{m≤√p} g(m) = o(p)`. Let `0<κ≤1/4`, pick `η' ∈ (0,1/2)` with
`f*(η') ≤ 1+κ/2`, and put `η(κ) = κη'/(1+2κ)`. Then for all `p` with
`κd > G(p) + √p + 1` and `p ≥ 60`, and all `A,B` with `|A||B| ≥ (1/2+κ)p`:
`|S(A,B)| ≤ (1−η(κ))|A||B|`. (For `κ > 1/4` apply the case `κ = 1/4`.)

*Proof.* `S` is symmetric in `(A,B)`; let `m = |A| ≤ n = |B|`. *Sub-sampling.* Put
`m₀ = ⌈(1/2+κ)p/n⌉ ≤ m`; since `n ≥ √(mn) ≥ √((1/2+κ)p)`, `m₀ ≤ √(3p/4)+1 ≤ √p` (`p ≥ 60`).
`S(A,B) = (m/m₀)·avg_{A'⊆A,|A'|=m₀} S(A',B)`, so it suffices to treat `(A',B)` with
`|A'| = m₀ ≤ √p` and `m₀n ≥ (1/2+κ)p`. *Sign.* For a non-residue `ν`, `S(νA',νB) = −S(A',B)`
with the same sizes and `r`, so it suffices to prove `S(A',B) ≤ (1−η₀)m₀n` with `η₀ = η(κ)`.
Suppose instead `S > (1−η₀)m₀n`. Then `N_− = (m₀n − r − S)/2 ≤ η₀m₀n/2`.
With `λ = 2(1+2κ)/κ`, the set `B' = {b∈B : e_b ≤ λη₀m₀/2}` has `|B'| > (1−1/λ)n` (Markov).
RHP with `e = ⌊λη₀m₀/2⌋`, `e/m₀ ≤ λη₀/2 = η'`: `m₀|B'| ≤ f*(η')d + g(m₀) + m₀`. Since
`m₀n > (1+2κ)d` and `(1−1/λ)(1+2κ) = 1+3κ/2`, we get `(1+3κ/2)d < (1+κ/2)d + G(p) + √p + 1`,
contradicting `κd > G(p)+√p+1`. ∎

*Remark (on the brief's example `g = O(m²)`).* The only regime not already covered
unconditionally is `m ≍ n ≍ √p` (for `m ≤ 24` use HP as in §2.5; for `n ≥ 24√p`, `m ≥ 24` the
fourth-moment/Weil bound `|S| ≤ (6/m + 3√p/n)^{1/4}mn` (with `mn ≥ p/2`) already saves a
constant). There `g(m) = Cm²` is of order `p`, so RHP with `g = O(m²)` would give the
consequence only for `κ ≳ C`; the hypothesis `G(p) = o(p)` is genuinely needed. The RHP proved
below has `g(m) = m`, so this issue does not arise.

## 2. The Hankel-minor Stepanov inequality (PROVED)

For `0 ≤ s ≤ D` let `u_s(x) = F^{(s)}(x)/(D)_s`, a polynomial of degree `d−s` (for `s ≤ d`) with
leading coefficient `C(D,m−1)(d)_s/(D)_s = C(D−s, m−1)`. For `0 ≤ e ≤ (m−1)/2` put

    H_{e+1}(x) = det[ u_{i+j}(x) ]_{0 ≤ i,j ≤ e}.

**Theorem 2.1.** For every odd prime `p`, every `A ⊆ F_p` with `1 ≤ m ≤ (p+1)/2` and every
integer `0 ≤ e ≤ (m−1)/2`, `H_{e+1}` has degree exactly `(e+1)(d−e)`, and for every `b` with
`e_b ≤ e`, `ord_b H_{e+1} ≥ Σ_{i=e_b}^{e} (m−e−i−δ_b)`. Consequently

    (★)   Σ_{b∈F_p, e_b ≤ e} (e+1−e_b)·( m − (3e+e_b)/2 − δ_b )  ≤  (e+1)(d−e).

For `e = 0` this is exactly HP (`m·#{e_b=0} − r_0 ≤ d`).

*Proof.* **Step 1 (local structure).** Fix `b` with `e_b ≤ e` and write `E = E(b)`,
`δ = δ_b`. By Lemma 1.1(d), `F(z) − 2R_E(z) = (z−b)^{m−δ}Q_b(z)` for a polynomial `Q_b`.
Differentiating `s` times and dividing by `(D)_s`:

    u_s(x) = 2ρ_s(x) + ε_s(x),  ρ_s(x) := Σ_{k∈E} c_k (x+a_k)^{D−s},  ord_b ε_s ≥ m − δ − s,

because `R_E^{(s)} = (D)_s Σ_E c_k (x+a_k)^{D−s}` and `((z−b)^{m−δ}Q_b)^{(s)}` is divisible by
`(z−b)^{m−δ−s}` (Leibniz). Here `s ≤ 2e ≤ m−1`.

**Step 2 (rank).** In the field `F_p(x)`, `ρ_{i+j} = Σ_{k∈E} [c_k(x+a_k)^D]·w_k^i·w_k^j` with
`w_k = (x+a_k)^{−1}`. So the matrix `[ρ_{i+j}]_{0≤i,j≤e}` is `V^T·diag·V` with `V` of size
`e_b × (e+1)`: its rank over `F_p(x)` is `≤ e_b`.

**Step 3 (order).** Expand `det[2ρ_{i+j} + ε_{i+j}]` multilinearly in the rows:
`H_{e+1} = Σ_{S⊆{0..e}} det M_S`, where `M_S` has row `i` equal to the `ε`-row for `i∈S` and to the
`2ρ`-row for `i∉S`. If `|S| < e+1−e_b`, the `2ρ`-rows of `M_S` are more than `e_b` rows of a
rank-`≤e_b` matrix, hence dependent, and `det M_S = 0`. If `|S| ≥ e+1−e_b`, every term of the
Leibniz expansion of `det M_S` contains one entry `ε_{i+j}` from each row `i∈S`, so
`ord_b det M_S ≥ Σ_{i∈S}(m−δ−i−e)`. The minimum over admissible `S` is at `S = {e_b,…,e}`:
`ord_b H_{e+1} ≥ Σ_{i=e_b}^{e}(m−e−i−δ) = (e+1−e_b)(m − (3e+e_b)/2 − δ)`; all terms are `≥ 0`
because `m−2e−1 ≥ 0`.

**Step 4 (degree and non-vanishing).** Each Leibniz term `Π_i u_{i,σ(i)}` has degree exactly
`Σ_i(d−i−σ(i)) = (e+1)(d−e)`, so the coefficient of `x^{(e+1)(d−e)}` in `H_{e+1}` is
`Λ := det[C(D−i−j, m−1)]_{0≤i,j≤e}`, and `deg H_{e+1} = (e+1)(d−e)` iff `Λ ≢ 0 (mod p)`,
which is Lemma 2.2.

**Step 5 (count).** A nonzero polynomial has `Σ_b ord_b ≤ deg`. ∎

**Lemma 2.2 (leading coefficient).** With `N' = D−2e`, `c = m−1−e`, `n = e+1`,

    Λ = (−1)^{e(e+1)/2} · det[C(N', c−i+j)]_{1≤i,j≤n}
      = (−1)^{e(e+1)/2} · Π_{k<n} k! · Π_{i=1}^{n} (N'+i−1)! / [ (c−i+n)!·(N'−c+i−1)! ],

and `Λ ≢ 0 (mod p)`.

*Proof.* Vandermonde's convolution (`x = e−i ≥ 0`, `y = D−e−j ≥ 0`) gives
`C(D−i−j, K) = Σ_l C(e−i,l) C(D−e−j, K−l)` (`K = m−1`), i.e. `[C(D−i−j,K)] = U·W` with
`U_{il} = C(e−i,l)`, `W_{lj} = C(D−e−j,K−l)`. `U` vanishes below the anti-diagonal `l = e−i` and
has `1`'s on it, so `det U = (−1)^{e(e+1)/2}`. Again `C(D−e−j,K−l) = Σ_t C(e−j,t)C(D−2e,K−l−t)`,
so `W = Z·U^T` with `Z_{lt} = C(N', K−l−t)`, and `Λ = det Z` (the two signs cancel). Reversing the
columns of `Z` (`t ↦ e−t`, sign `(−1)^{e(e+1)/2}`) gives the Toeplitz matrix
`[C(N', c−l+t)]`, `c = K−e`. Its determinant is Krattenthaler's (3.12) at `q=1` with `A = N'`,
`L_i = c−i` (CITED: C. Krattenthaler, *Advanced determinant calculus*, Sém. Lothar. Combin. 42
(1999) B42q, Theorem 26, eq. (3.12):
`det_{1≤i,j≤n} C(A, L_i+j) = Π_{i<j}(L_i−L_j)·Π_i (A+i−1)! / (Π_i (L_i+n)!·Π_i (A−L_i−1)!)`;
local text `sources/stepanov-robust-krattenthaler-math9902004.txt`), with
`Π_{i<j}(L_i−L_j) = Π_{k<n}k!`. All factorial arguments are nonnegative (`c ≥ e` because
`2e ≤ m−1`; `N'−c = d−e ≥ 0`). Non-vanishing: `Λ` is an integer and
`Λ·Π_i(c−i+n)!(N'−c+i−1)! = ±Π_{k<n}k!·Π_i(N'+i−1)!`; the right side is a product of integers
in `[1, N'+n−1] = [1, D−e] ⊆ [1, p−1]`, so `p ∤ Λ`. ∎

(Equivalently, by the dual Jacobi–Trudi identity `det Z = ±s_{((e+1)^{m−1−e})}(1^{D−2e})` and
the hook-content formula, `Λ = ±Π_{cells}(D−2e+content)/hook` with all factors in `[1,p−1]`; we
did not fetch a source for this form, so the proof above uses only (3.12).)
Verifier §C1: the integer identity for all `D ≤ 40`, `m ≤ 25`, admissible `e`, and an
independent modular elimination showing `Λ ≢ 0` for all `p ≤ 300`, `m ≤ min(40,(p+1)/2)`.
Verifier §C2 computes the full polynomial `H_{e+1}` for a sample and confirms the exact degree
and orders; §C3–C4 confirm the order bounds and (★) for **every** `A ∋ 0` with `p ≤ 13`
(`(★)` is translation invariant), every `A ∋ 0` with `p = 17`, `m ≤ 6`, and random/structured
`A` (intervals, subgroups, Paley cliques, squares) for `p ≤ 211`. No failure. The tightest
instances reach `LHS/RHS = 1` only at `e = 0` (HP equality); at `e ≥ 1` the maximum observed
is `≈ 0.92` (results file, `theorem1_tightest`). **Erratum (root, 2026-09-27):** this is
wrong in general; (★) holds with equality at `e = 1`, e.g. `p = 13`, `A = {0,2,3,5}` (both
sides equal 10; `research/stepanov-referee2-2026-09-27.md`), and the root's exhaustive check
(`results/stepanov_star_exhaustive_2026_09_27.json`) finds ratio 1 at `e = 1` for every
`p ≤ 13`. The inequality itself is unaffected.

**Corollary 2.3 (bias inequality).** For `A` as above, any `B ⊆ F_p` and `0 ≤ e ≤ (m−1)/2`:

    (m−2e)·((e+1)n − N_−) ≤ (e+1)(d−e+r),   hence
    S(A,B) ≤ mn − r − 2(e+1)·[ n − (d−e+r)/(m−2e) ].

*Proof.* In (★) restrict to `b ∈ B` (terms are `≥ 0`), use `(3e+e_b)/2 ≤ 2e`, and add the
nonpositive terms `(e+1−e_b)(m−2e)` for `b∈B` with `e_b > e`:
`Σ_{b∈B}(e+1−e_b)(m−2e) ≤ (e+1)(d−e) + Σ_{b∈B∩(−A)}(e+1−e_b) ≤ (e+1)(d−e+r)`. Then
`S = mn − r − 2N_−`. ∎ (Verifier §D: every `e`, on level sets `B = {e_b ≤ e'}` and random `B`.)

**Theorem 2.4 (RHP, PROVED).** Let `m ≤ (p+1)/2`, `B` with `e_b ≤ e'` for all `b∈B`. For every
integer `e` with `e' ≤ e ≤ (m−1)/2`:
`(e+1−e')[(m − (3e+e')/2)n − r] ≤ (e+1)(d−e)`, hence

    mn − r ≤ (e+1)m·d / ((e+1−e')(m−2e)) + 2e·r/(m−2e).

In particular, if `η = e'/m ≤ 1/8`: `mn − r ≤ (1−√(2η))^{−2}·d + m`.
So RHP(f,g) holds with `f(η) = (1−√(2η))^{−2}` on `[0,1/8]` and `g(m) = m`, for all
`m ≤ (p+1)/2` (not only `m ≤ √p`). (On `(1/8, 1/2−δ]` the second moment
`n(m−2e'−1)² ≤ m(p−m)` supplies a bounded `f`; it is not needed.)

*Proof.* The first display is (★) restricted to `B` with each weight bounded below using
`e_b ≤ e'`. Divide by `e+1−e'`, bound `m−(3e+e')/2 ≥ m−2e`, and write
`mn − r = (m/(m−2e))((m−2e)n − r) + r(m/(m−2e) − 1)`. For the `η`-form: if `e' = 0` take `e = 0`
(HP). Otherwise `m ≥ 8`; let `x = √(η/2) ≤ 1/4` and `e = ⌈xm⌉−1`, so `e+1 ≥ xm`, `e < xm ≤ m/4`.
Since `m ≥ 1/η` and `x − η ≥ x/2` (as `η ≤ 1/8`), `m(x−η) ≥ 1/(2√(2η)) ≥ 1`, so `e ≥ xm−1 ≥ e'`.
Then `(e+1)/(e+1−e') ≤ 1/(1−η/x) = 1/(1−√(2η))`, `m/(m−2e) ≤ 1/(1−2x) = 1/(1−√(2η))`, and
`2er/(m−2e) ≤ m·2x/(1−2x) ≤ m`. ∎

Comparison: the second moment gives `f₂(η) = 2/(1−2η)²`; Theorem 2.4 is better for
`η < 0.0858` and tends to HP's `1` as `η → 0`. For `e' = 1` it gives
`m·N_1 − r ≤ (1+O(m^{−1/2}))d`, against `2(m/(m−1))²d` in the abstract model of §4.

### 2.5 The constant bias saving (PROVED, unconditional)

**Theorem 2.5.** Let `p ≥ 11`, `0 < κ ≤ 3/2`, `u = (1+2κ)^{−1/2}`, and `A,B ⊆ F_p` with
`|A||B| ≥ (1/2+κ)p`. Then

    |S(A,B)| ≤ [ 1 − (1−u)² + (√((1/2+κ)p) + 1)/(2(p−1)) ]·|A||B|.

For `κ > 3/2` the bound with `κ = 3/2` applies (saving `1/4 − O(p^{−1/2})`). For small `κ`,
`(1−u)² = κ² − 3κ³ + O(κ⁴)`.

*Proof.* Let `m ≤ n` and `m₀ = ⌈(1/2+κ)p/n⌉ ≤ m`; `m₀ ≤ √((1/2+κ)p)+1 ≤ (p+1)/2` for `p ≥ 11`.
By averaging over `m₀`-subsets `A' ⊆ A` (as in Prop. 1.4) and the non-residue dilation, it
suffices to show `S(A',B) ≤ [1 − (1−u)² + u(1−u)m₀/d]·m₀n` for `|A'| = m₀`
(then use `u(1−u) ≤ 1/4`, `m₀/d ≤ 2(√((1/2+κ)p)+1)/(p−1)`). Put `θ = (1−u)/2 ∈ (0,1/4]` and
`e = ⌊θm₀⌋ ≤ (m₀−1)/2`, so `e+1 ≥ θm₀` and `m₀−2e ≥ um₀`. Corollary 2.3 with `r ≤ m₀`:
`N_− ≥ (e+1)[n − (d+m₀)/(um₀)]`. If the bracket is negative then `um₀n < d+m₀`, which with
`m₀n > d/u²` gives `(1−u)/u < m₀/d` and the claimed factor is `≥ 1`: nothing to prove. Otherwise
`N_− ≥ θ[m₀n − (d+m₀)/u]` and `S ≤ m₀n − 2N_− ≤ m₀n[1 − (1−u)(1 − (d+m₀)/(u m₀ n))]`; since
`(d+m₀)/(um₀n) < u(1+m₀/d)`, this is `≤ m₀n[1 − (1−u)² + u(1−u)m₀/d]`. ∎

**Remark 2.6 (prime sensitivity).** `D ≤ p−1` is used for: Taylor expansion and the order
criterion (`j!` units), the normalisation `u_s = F^{(s)}/(D)_s`, and `p ∤ Λ` (all factorials in
Lemma 2.2 have arguments `≤ D−e ≤ p−1`). Over `F_{p²}` with `A = B = F_p` every `b ∈ F_p` has
`e_b = 0` and `mn − r = p²−p > d`, so (★) is false there; correspondingly `(D)_j ≡ 0 (mod p)`
for some `j ≤ m−1` (verifier §F).

**Remark 2.7 (sharpness of the threshold 1/2).** For `A = {a_1,a_2}` the level set
`B = {b : e_b = 0}` has `|B| = (p+3)/4` for suitable `A` (CITED in-workspace:
`research/sigma-biclique-2026-09-05.md`, (b), `M_2(p) = (p+3)/4`), so `|A||B| = (p+3)/2` and
`S = |A||B| − r` with `r ≤ 2`: `|S|/|A||B| → 1` at `|A||B| = p/2 + O(1)`. Also for `m = 2`
Corollary 2.3 (with `e = 0`) gives `S ≤ d`, and adjoining arbitrary further points to that `B`
shows that no `η(κ) > 2κ/(1+2κ)+o(1)` is possible in general. So `η(κ)` lies between
`(1−(1+2κ)^{−1/2})² ≈ κ²` and `2κ/(1+2κ) ≈ 2κ`; closing this gap is OPEN.

**Proposition 2.9 (index-`k` subgroups; PROVED).** Let `k | p−1`, `d_k = (p−1)/k`,
`D = d_k+m−1 ≤ p−1`, `F = −1+Σ_l c_l(x+a_l)^D`, and redefine
`e_b = #{a∈A : a+b ≠ 0, (a+b)^{d_k} ≠ 1}`. Then (★) holds with `d` replaced by `d_k`, for every
`0 ≤ e` with `2e ≤ min(m−1, d_k)`. *Proof.* Only Lemma 1.1(d) changes. With `y_l = b+a_l`,
`ω_l = y_l^{d_k}` and `λ_l = c_l(1−ω_l^{−1})` for `l∈E(b)`, put `R = Σ_{l∈E(b)} λ_l(x+a_l)^D`.
For `1 ≤ j ≤ m−1`, `F^{(j)}(b) − R^{(j)}(b) = (D)_j Σ_{l: y_l≠0} c_l y_l^{m−1−j}[ω_l − (ω_l−1)[l∈E]]
= (D)_j Σ_{l: y_l≠0} c_l y_l^{m−1−j}`, which is `0` by Fact 0 when `b ∉ −A`, and `0` for
`j ≤ m−2` when `b ∈ −A`; similarly at `j = 0`. So `F − R` has order `≥ m−δ_b` at `b`, and `R` is a
combination of `e_b` pure `D`-th powers, which is all Steps 1–3 use. Step 4: `deg u_s = d_k−s ≥ 0`
for `s ≤ 2e ≤ d_k`, and Lemma 2.2 applies with `N'−c = d_k−e ≥ 0` and all factorial arguments
`≤ D−e ≤ p−1`. ∎ (Verifier §C5: `k = 3, 4`, 44 instances, no failure. We have not derived the
corresponding bias statement for characters of order `k`.)

**Remark 2.8 (what is and is not new here).** The only new input relative to HP is Steps 2–3:
the Taylor data of `F` at `b` form the syndrome of a Reed–Solomon-type code with error locators
`(b+a_k)^{−1}`, `k∈E(b)`; the Hankel minor vanishes to order `≈ (e+1−e_b)m` at `b` (higher when
there are *fewer* errors) while its degree is only `(e+1)d`. The weights `(e+1−e_b)` are exactly
what the bias application needs, because `Σ_b(e+1−e_b) = (e+1)n − N_−` is linear in the number
of negative sums. The abstract no-go of §4 shows that this algebraic step cannot be replaced by
HP applied to subsets.

## 3. Literature (CITED; searched 2026-09-26)

- Hanson–Petridis, arXiv:1905.09134 (PLMS 2021), Theorem 1.2: "`A+B ⊆ Z_d ∪ {0}` … Then
  `|A||B| ≤ d + |B ∩ (−A)|`"; §2 proof as re-derived in §1.1. The method "only works in prime
  fields" (text after Cor. 1.5). No robust form.
- Kalmynin, arXiv:2504.10202, *On additive irreducibility of multiplicative subgroups*
  (text: `sources/stepanov-robust-kalmynin-2504.10202.txt`). Definition 2 (d-critical pair:
  `A+B ⊆ μ_d∪{0}` and `|A||B| = d + |(−A)∩B|`); Lemma 4 (in the critical case
  `HP(x;A,d) = C Π_{b∈B}(x−b)^{α−ε(b)}`); Theorem 3 (Sárközy's conjecture). Works only with exact
  containment and equality; contains a differential operator "designed to annihilate" the
  root structure but no statement allowing exceptional sums.
- Rudnev–Tyrrell, arXiv:2607.24270 (v2 12 Aug 2026), Theorem 1.1: if `H = A+B` for a proper
  subgroup then a summand is a singleton or `|A|=|B|=2`, `|H|=4`. Exact decompositions only; its
  survey paragraph lists Yip [Yip24, Yip25], Kim–Yip–Yoo [KYY23, KYY26], Yip–Yoo [YY26a,b],
  Yoo, Cochrane — all exact-containment statements (abstracts of arXiv:2608.02568 and
  2607.25711 fetched; neither mentions exceptions or character sums).
- Yip, arXiv:2501.16620 (Combinatorica 46, 2026): "small perturbations" of shifted `k`-th powers
  over `ℕ` (product sets, Hajdu–Sárközy); not a Stepanov bound with exceptions.
- Yip's publication list (fetched from his research page, 72 items) contains no robust or
  character-sum version of HP; the workspace note `research/sigma-lit2026-2026-09-05.md`
  (§A) records the same for 2024–2026 and states that HP-type results "give nothing at all for
  sums with cancellation"; Theorem 2.5 above shows this is not so for a small proportion of
  negative sums.
- Heath-Brown–Konyagin (Q. J. Math. 2000) and Shkredov (Acta Arith. 164, 2014) are cited in
  HP's bibliography; we did not fetch them (UNVERIFIED CITATION for their content); from their
  titles and HP's description they concern subgroups/Gauss sums and exact sumset containment.
- Hankel-determinant input: Krattenthaler, *Advanced determinant calculus*, Thm 26 (3.12)
  (fetched, local copy).

Verdict: we found no prior statement of RHP, of (★), or of a constant saving for
`p/2 < |A||B| < p`. This is a negative search result, not a proof of novelty.

## 4. No-go: subset-HP alone does not imply RHP (PROVED)

HP applied to every `T ⊆ A` says: for `R_k = {x ∉ −A : χ(x+a_k) = 1}` (`|R_k| ≤ d`,
`|F_p∖(−A)| ≤ 2d+1`), every `t`-fold intersection has size `≤ d/t`. RHP concerns the number of
`x` lying in `≥ (1−η)m` of the `R_k`.

**Theorem 4.1.** Let `m ≥ 5` be odd, `j = (m−1)/2`, `k ≥ 1`, `w₀ = k·C(m,j)`,
`μ = m(m−3)k`, `d = w₀(m−1)²/2`. Let `U` consist of: for each `i ∈ [m]`, `w₀` points lying
exactly in the sets `R_l`, `l ≠ i`; for each `j`-subset `T ⊆ [m]`, `μ` points lying exactly in
`R_l`, `l∈T`; and `w₀+1` points in no set. Then `|U| = 2d+1`, `|R_l| = d` for all `l`, every
`t`-fold intersection has size `≤ d/t` (`t = 2` with equality), `Σ_{x∈U}(2cov(x)−m)² = m(2d+m)`
(the Paley second moment is `m(p−m)`), and exactly `m w₀ = 2md/(m−1)²` points have coverage
`m−1`. Hence `m·#{x : cov(x) ≥ m−1} = 2(m/(m−1))²·d > 2d`, while RHP with `η = 1/m` demands
`≤ f(1/m)d + g(m)` with `f(1/m) → 1`.

*Proof.* By symmetry a `t`-fold intersection has size
`I(t) = w₀(m−t) + μ·C(m−t, j−t)`. Using `μ·C(m,j) = m(m−3)w₀`:
`I(1) = w₀(m−1) + m(m−3)w₀·j/m = w₀(m−1)²/2 = d`;
`I(2) = w₀(m−2) + m(m−3)w₀·j(j−1)/(m(m−1)) = w₀(m−1)²/4 = d/2`;
`|U| = m w₀ + m(m−3)w₀ + w₀ + 1 = (m−1)²w₀ + 1 = 2d+1`.
For `t ≥ 3`: `μC(m−t,j−t) = μC(m,j)·C(j,t)/C(m,t) ≤ 2d·(j/m)^t < 2d·2^{−t}` and
`w₀(m−t) = 2d(m−t)/(m−1)²`, so `t·I(t)/d < t2^{1−t} + 2t(m−t)/(m−1)²`. For `t = 3` this is
`≤ 3/4 + 6/(m−1) < 1` (`m ≥ 26`); for `4 ≤ t ≤ 10`, `≤ 1/2 + 20m/(m−1)² < 1` (`m ≥ 43`); for
`t ≥ 11`, `≤ 11/1024 + m²/(2(m−1)²) < 1`. This settles `m ≥ 43`. For `5 ≤ m ≤ 41` the `m`
inequalities `t·I(t) ≤ d` are homogeneous in `k`, so their exact integer evaluation at `k = 1`
(verifier §E; all hold, for every odd `m` from `5` to `59`) is a proof for every `k`.
The second moment: `Σcov = md`,
`Σcov(cov−1) = 2C(m,2)I(2) = m(m−1)d/2`, so `Σ(2cov−m)² = 4Σcov² − 4mΣcov + m²|U| = m(2d+m)`. ∎

Values of `m·#{cov = m−1}/d = 2m²/(m−1)²`: `3.125` (`m=5`), `2.42` (`m=11`), `2.10` (`m=41`),
`2.07` (`m=59`). The verifier also cross-checks the symmetric intersection formula by explicit
enumeration of the materialised `m = 7` system (all `127` subsets `T`), and exhibits `m = 43` with
`2d+1` prime. Conclusion: any proof of RHP must use more than subset-HP
plus sizes and the second moment; Theorem 2.1 uses the Hankel (rank) structure of the Taylor
data, which the abstract system does not have.

## 5. Other attempts and where the method stops

(a) *Hankel route* (task 4a): succeeded — §2. The choice of the leading `(e+1)`-minor and of the
row expansion gives order `(e+1−e_b)(m − (3e+e_b)/2 − δ_b)`; other minors (rows `I`, columns `J`)
give `Σ_{i∈S}(m − i − max J)` and did not improve the constants in the cases we tried by hand.
(b) *Several auxiliary polynomials / polynomial coefficients* (task 4b): a naive parameter count
for `Ψ = g_0 + Σ_k g_k(x)(x+a_k)^d` with `deg g_k < L` needs `L ≳ m²/θ` to reach order
`(1−θ)m`, i.e. degree `≫ p` when `m ≍ √p`; HP's choice `g_k = c_k(x+a_k)^{m−1}` beats the count
by an exact cancellation, and the Hankel minors are the natural nonlinear continuation of that
choice (they are the Wronskian-type invariants of the `e`-dimensional span
`{(x+a_k)^D : k∈E}`). We did not pursue (b) further. HEURISTIC.
(c) *Limits* (OPEN): the method gives `f(η) = 1 + O(√η)` and `η(κ) ≍ κ²`; it gives nothing when
`|A||B| ≤ p/2` (necessarily, Remark 2.7) and no power saving. The loss `2e` in `m−2e` comes
from using entries `u_{i+j}` with `i+j` up to `2e`. For any `(e+1)`-minor (rows `I`, columns
`J`) the Step-3 bound in row `i` is `m−i−max J ≤ m−i−e`, so the Step-3 bound always loses at
least `e` per row, and optimising `e` then forces `f = 1+Θ(√η)` within this method
(HEURISTIC: the true orders could exceed the Step-3 bound). The extension
to index-`k` subgroups is Proposition 2.9.

## 6. Open obligations

1. Independent referee of Theorem 2.1 (Steps 1–5) and Lemma 2.2; the verifier checks every step
   numerically but is not a proof.
2. Improve `f(η)` to `1+O(η)` and `η(κ)` to `≍ κ` (Remark 2.7 gap), e.g. by optimising the
   minor or combining several `e`.
3. The bias statement for characters of order `k` from Proposition 2.9, and the `A − A`
   (clique) version: Theorem 2.5 with `B = −A` gives, for `|A|² ≥ (1/2+κ)p`, at least
   `≈ κ²|A|²/2` pairs with `χ(a−a') = −1` — to be written out with the diagonal handled.
4. Whether the Hankel structure can be pushed to growing `e` with a power saving (any such bound
   would have to fail over `F_{p²}`, which (★) does).
