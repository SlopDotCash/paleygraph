# Stepanov wave, worker `hybrid`: Hanson–Petridis counts plus semidefinite constraints (2026-09-29)

**Status:** PROVED: (i) validity and the standalone value of the natural positive-semidefinite (PSD)
constraints for a complete biclique `A + B ⊆ Q ∪ {0}` (§1): the spectral Gram family of the Paley matrix
(P1) gives exactly the Chung form `mn − r ≤ √(mn(p−m)(p−n)/p)` (so `mn ≲ p`) and, for cliques, the
ϑ-bound `m ≤ √p`; principal submatrices (P2) are a special case; the Gram matrix of the shifted characters
`x ↦ χ(x+a)` (P4, `k = 1`) gives `mn ≤ p − m + O(r)`; the bipartite ϑ-program is `≤ √p + 1`;
(ii) a **decoupling theorem** (Theorem 2.1): for every profile satisfying the exact first and second moment
identities (and, two-sided, the joint moment identity), an explicit symmetric operator `S̃` with
`S̃1 = 0`, `S̃² = pI − J`, `S̃1_A = f_A`, `S̃1_B = f_B` exists; hence every level-set compression and
spectral-Gram constraint (P1, P3) is automatically satisfiable and adds nothing to the count LP of the
balanced note; (iii) the `S_A`-symmetrised moment matrices of the lifted characters (P4, any degree) are
automatically PSD and contribute exactly the Krawtchouk–Weil linear inequalities (Proposition 2.3).
By exact rational certificates (§3): the hybrid program (all theorem rows of the balanced note + all
Krawtchouk–Weil rows + P1 + P3) has optimum **exactly at the Hanson–Petridis value** at
`p = 1009, 10009, 40009, 100049, 1000033`, for balanced bicliques (`mn/p = 0.4797 … 0.4998`, certified
`n_0 = ⌊d/m⌋`, LP optimum `d/m`) and for cliques (`m = ⌊(1+√(2p−1))/2⌋`); and the genuinely two-point,
Delsarte-type PSD program P5 (Gram matrices of the lifted characters indexed by pairs of points, not
covered by (ii)–(iii)) is feasible on those profiles, with exact `LDLᵀ` certificates at
`p = 1009, 10009, 40009, 100049` (bicliques and cliques). REFUTED: "adding spectral / PSD constraints of types
P1–P5 to the count constraints excludes Hanson–Petridis-tight balanced bicliques or cliques" (witnesses:
Tables 3.1, 3.2). OPEN: the target `|A||B| ≤ (1/2 − c)p + o(p)` (equivalently for cliques
`ω(G_p) ≤ √((1/2 − c)p)`). The obstruction (§4) is precise: the PSD families P1–P5, built from
`S² = pI − J` and Weil bounds, are field-agnostic and see `B` only through degree-≤2 identities, where
Hanson–Petridis-tightness is invisible; the only PSD routes with numerical evidence of beating
Hanson–Petridis (localised ϑ, degree-4 SOS; §5) use the Paley adjacency *inside* a localisation and are not
expressible in the aggregated profile variables. No combination of Stepanov/Hanson–Petridis with SDP/ϑ
constraints was found in the literature (§5; completeness not claimed).

Verifier: `experiments/stepanov_hybrid_2026_09_29.py` → `results/stepanov_hybrid_2026_09_29.json`
(run with `/opt/miniconda3/bin/python3`; about 6 min, of which 5 min are the two-point SDPs at
`p = 10009, 40009`; 22,148 checks, 0 failures;
needs `cvxpy` with the Clarabel/SCS solvers, which are installed). Local sources:
`sources/stepanov-hybrid-randomstrasse101-2603.29571-entry12.txt`,
`sources/stepanov-hybrid-mmp-1907.05971.txt`, `sources/stepanov-hybrid-sdp-abstracts.txt`, and the
pre-existing `sources/kunisky-2303.16475v1.html`.

---

## 0. Setting

`p` odd prime, `χ` the Legendre symbol (`χ(0) = 0`), `d = (p−1)/2`, `Q` the nonzero squares.
`S = [χ(x+y)]_{x,y∈F_p}` (sum form); for `p ≡ 1 mod 4` also `M = [χ(x−y)]` (the ±1 Paley matrix).
`A + B ⊆ Q ∪ {0}`, `m = |A|`, `n = |B|`, `r = |A ∩ (−B)|`, `s = |A ∩ B|`,
`f_A = S1_A`, i.e. `f_A(x) = Σ_{a∈A} χ(x+a) = m − δ_x − 2e_x` with `e_x = #{a : χ(x+a) = −1}` and
`δ_x = [x ∈ −A]`. The **profile** of `A` is `n_j = #{x ∉ −A : e_x = j}` (`0 ≤ j ≤ m`),
`n'_j = #{x ∈ −A : e_x = j}` (`0 ≤ j ≤ m−1`). Points of `B ∖ (−A)` have `e_x = 0`.

The balanced note (`research/stepanov-balanced-2026-09-29.md`, Theorem 3.5) showed that its list of count
constraints — the exact identities `Σ(n_j + n'_j) = p`, `Σn'_j = m`, `Σ_x f_A = 0`,
`Σ_x f_A² = m(p−m)`; the Hankel-minor inequalities (★)_e for `A` and `νA`; averaged subset
Hanson–Petridis; the pencil and derivative inequalities; Weil moments of orders 4, 6, 8 and Weil counts —
is satisfied by rational profiles with `n_0 = ⌊d/m⌋`, `m ≈ √(p/2)`. We call these the **theorem rows**. The
question of this note is whether PSD constraints, added to the theorem rows, exclude such profiles.

Throughout, "PSD" means positive semidefinite, and for a finite family `v_1,…,v_k ∈ R^{F_p}` and a symmetric
matrix `X`, `Gram_X(v) := [v_iᵀXv_j]`.

## 1. The PSD families and what each gives alone (task 1)

### 1.1 The Paley matrix

**Lemma 1.1 (PROVED).** (a) `S1 = 0`, `S² = pI − J`, `tr S = 0`. (b) `Π_± := ½(I − J/p ± S/√p)` are
orthogonal projections, `Π_+ + Π_− = I − J/p`, `Π_+Π_− = 0`, `S = √p(Π_+ − Π_−)`, `rank Π_± = d`.
(c) For `p ≡ 1 mod 4` the same holds for `M`.

*Proof.* (a) `(S²)_{xz} = Σ_y χ(x+y)χ(z+y) = Σ_u χ(u)χ(u+z−x)`, which is `p−1` for `z = x` and `−1`
otherwise (background fact `Σ_u χ(u)χ(u+c) = −1`, `c ≠ 0`); `Σ_y χ(x+y) = 0`; `tr S = Σ_x χ(2x) = 0`.
(b) Expand using `SJ = JS = 0`, `J² = pJ`, `S² = pI − J`; `tr Π_± = ½(p − 1 ± 0) = d`. (c) Same computation
with `u = x − y`; `M` is symmetric because `χ(−1) = 1`. ∎ (Verifier §A: all primes `5 ≤ p < 140`.)

### 1.2 P1: the spectral Gram family

**Definition.** For any family `v`, `G_±(v) := Gram_{2Π_±}(v) = Gram_{I−J/p}(v) ± Gram_S(v)/√p ⪰ 0`.
For `v = (1_A, 1_B)`, with `σ_A = 1_AᵀS1_A = Σ_{a,a'} χ(a+a')`, `σ_B` likewise and
`τ = 1_AᵀS1_B = mn − r`:

    G_± = [[ m − m²/p ± σ_A/√p ,  s − mn/p ± τ/√p ],
           [ s − mn/p ± τ/√p   ,  n − n²/p ± σ_B/√p ]]  ⪰ 0.

**Proposition 1.2 (PROVED).** (a) If `A + B ⊆ Q ∪ {0}` then P1 for `(1_A, 1_B)` implies
`mn − r ≤ √(mn(p−m)(p−n)/p)`. (b) This is the exact value of P1 as a program in the unknowns
`(σ_A, σ_B)`: for `m = n ≤ p/3`, `s = 0`, `τ ≥ 0` the matrices `G_±` are PSD with
`σ_A = σ_B = −m²/√p` **iff** `τ ≤ √p·m(1 − m/p)`. So P1 alone gives only `m² ≲ p` for balanced bicliques (Table 1.1), not `p/2`.
(c) For a clique `C` (`p ≡ 1 mod 4`), P1 for `1_C` with `M` reads `m − m²/p ± m(m−1)/√p ≥ 0`, i.e.
`m ≤ √p`: the Delsarte–Hoffman bound, equal to the Lovász number `ϑ(G_p) = √p` (CITED: Lovász 1979,
Thm 8, as quoted in MMP, arXiv:1907.05971, Proposition 1(ii), and Kunisky, arXiv:2303.16475, after
Theorem 1.2: "√p is the value of both the spectral “Hoffman bound” and the Lovász ϑ function bounds").

*Proof.* (a) Write `m' = m − m²/p`, `n' = n − n²/p`, `s' = s − mn/p`, `α = σ_A/√p`, `β = σ_B/√p`. PSD
`2×2` matrices satisfy `|G_{12}| ≤ √(G_{11}G_{22})`, so
`2τ/√p ≤ |s' + τ/√p| + |s' − τ/√p| ≤ √((m'+α)(n'+β)) + √((m'−α)(n'−β)) ≤ 2√(m'n')`
(the first step is the triangle inequality, the last is Cauchy–Schwarz). Square. (b) With `α = β = s'`
(`s' = −m²/p`), `G_+ = [[m'+s', s'+τ/√p],[s'+τ/√p, m'+s']]` is PSD iff `|s' + τ/√p| ≤ m' + s'`, and
`G_−` iff `τ/√p − s' ≤ m' − s'`; both iff `τ ≤ √p m'` (for `τ ≥ 0`). (c) The `1×1` Gram. ∎

**Table 1.1** (verifier §A; largest balanced `m = n` with `r = 0` allowed by P1 vs by Hanson–Petridis):

| p | P1: max m | P1: m²/p | HP: max m | HP: m²/p |
|---|---|---|---|---|
| 1009 | 30 | 0.892 | 22 | 0.480 |
| 10009 | 99 | 0.979 | 71 | 0.504 |
| 100049 | 315 | 0.992 | 224 | 0.502 |
| 1000033 | 999 | 0.998 | 707 | 0.500 |

(The HP column uses `m² ≤ d + m`, which allows `r ≤ m`.) P1 was checked exactly on 978 actual bicliques
(`p < 140`, random `A`, `B = B(A)` and random subsets), in `Q(√p)` arithmetic, together with the bound (a);
the ϑ-bound (c) on 300 greedy maximal cliques.

### 1.3 P2: principal submatrices (interlacing)

For `U ⊆ F_p`, P1 applied to the standard basis vectors `(e_u)_{u∈U}` gives
`(I − J/p ± M/√p)_{UU} ⪰ 0`, which contains the interlacing statement `‖M_U‖ ≤ √p` for the Paley matrix
restricted to `U = A ∪ (−B)`; and `M² = pI − J` gives `M_U² ⪯ (pI − J)_{UU}`. For a clique
`M_C = J − I` and these give `m ≤ √p` and `m(m−1) ≤ p − 1`. For a biclique the off-diagonal block of
`M_{A∪(−B)}` is `J` minus `r` entries (`M_{a,−b} = χ(a+b)`), and `‖M‖ = √p` gives
`mn − r = 1_AᵀM1_{−B} ≤ √p·√(mn)` (PROVED), i.e. `mn ≲ p` again; the internal blocks `M_A`, `M_{−B}` are
unconstrained by the biclique. P2 is a sub-family of P1 and is not used further.

### 1.4 P3: level-set compressions of the Paley matrix

**Lemma 1.3 (PROVED).** Let `F_p = ⊔_ℓ L_ℓ` be a partition on whose cells `f_A` and `f_B` are constant
(values `f^A_ℓ`, `f^B_ℓ`), `ν_ℓ = |L_ℓ|`, `a_ℓ = |A ∩ L_ℓ|`, `b_ℓ = |B ∩ L_ℓ|`,
`T_{ℓk} = 1_{L_ℓ}ᵀS1_{L_k}`. Then
(T1) `T1 = 0`; (T2) `T f^A = p a − m ν` and `T f^B = p b − n ν`;
(T3) `G_±` of the family `{1_{L_ℓ}} ∪ {1_A, 1_B}` is PSD, where the `S`-entries are `T`,
`1_{L_ℓ}ᵀS1_A = ν_ℓ f^A_ℓ`, `1_AᵀS1_A = Σ_ℓ a_ℓ f^A_ℓ`, `1_AᵀS1_B = Σ_ℓ b_ℓ f^A_ℓ = Σ_ℓ a_ℓ f^B_ℓ`;
(T4) `Σ ν_ℓ f^A_ℓ = 0`, `Σ ν_ℓ (f^A_ℓ)² = m(p−m)`, `Σ ν_ℓ f^A_ℓ f^B_ℓ = ps − mn`.

*Proof.* `S1 = 0`; `S f_A = S²1_A = p1_A − m1`; P1; `⟨S1_A, S1_B⟩ = 1_Aᵀ(pI − J)1_B`. ∎
(Verifier §B: exact on 137 actual configurations with the joint `(A,B)`-level partition, `p < 115`; `G_±`
PSD numerically.)

### 1.5 P4: Gram / moment matrices of the shifted characters

For `T ⊆ A` put `S_T := Σ_x Π_{a∈T} χ(x+a)` and `Z_{T,T'} := Σ_x Π_{a∈T}χ(x+a)·Π_{a∈T'}χ(x+a)`.

**Lemma 1.4 (PROVED).** (a) `Z_{T,T'} = S_{T△T'} − Σ_{c∈T∩T'} Π_{a∈T△T'} χ(a−c)`.
(b) `S_∅ = p`, `S_{a} = 0`, `S_{a,b} = −1`, and `|S_U| ≤ (|U|−1)√p` for `|U| ≥ 3` (Weil, for the
squarefree polynomial `Π_{a∈U}(x+a)`; background fact of the brief). (c) `Z_k = [Z_{T,T'}]_{|T|,|T'|≤k}`
is a Gram matrix, hence PSD. (d) (`k = 1`, the Gram of the vectors `x ↦ χ(x+a)`.) The `A×A` block is
`pI − J`, and splitting off the rows `x ∈ B` gives `pI − J − Σ_{b∈B} v_bv_bᵀ ⪰ 0`,
`v_b = (χ(a+b))_a`; at the vector `1`: `(n − r)m² + r(m−1)² ≤ m(p−m)`.

*Proof.* (a) `χ(x+c)² = 1` unless `x = −c`. (b) Background facts. (c), (d) Gram matrices of real vectors;
(d) at `1` is `Σ_{x∈B} f_A(x)² ≤ Σ_x f_A(x)²`. ∎ (Verifier §C: 4,760 entries of `Z_2` for random `A`,
`m ≤ 4`, 14 primes.)

Note that the joint Gram of `{χ(·+a)}_{a∈A} ∪ {χ(·+b)}_{b∈B}` is `pI − J` (plus `p` on `A ∩ B`) and does
not see the biclique relation at all: that relation lives in the coordinates `x ∈ B` of the vectors
`χ(·+a)`, i.e. in the *rows* of `V = [χ(x+a)]`, which is how (d) uses it.

### 1.6 P5: the two-point (Delsarte-type) Gram family

For `1 ≤ t ≤ m` let `ψ_t(x) = (Π_{a∈T}χ(x+a))_{|T|=t}`. The `p × p` matrix
`Y_t = [⟨ψ_t(x), ψ_t(y)⟩]_{x,y}` is a Gram matrix, and `⟨ψ_t(x), ψ_t(y)⟩ = K_t(h; ν)`, where `ν` is the
number of coordinates `a` with `χ(x+a)χ(y+a) ≠ 0`, `h` the number with `χ(x+a)χ(y+a) = −1`, and
`K_t(h; ν) = Σ_i (−1)^i C(h,i) C(ν−h, t−i)` is the Krawtchouk polynomial. Aggregated over the levels
`ℓ ∈ {N_j} ∪ {A_j}` of the profile, with `P` = the counts of ordered pairs `(x,y)` by
`(level x, level y, ν, h, ε)` (`ε` = the sign of `x` at the zero coordinate of `y ∈ −A`), the following hold
for every actual configuration (**Lemma 1.5, PROVED**):

- marginals: `Σ P = ν_ℓν_k` (minus the diagonal for `−A`-levels); for `N`-levels at least `ν_ℓ` pairs
  at `h = 0` (the pairs `x = y`); supports `|j − k| ≤ h ≤ min(j+k, 2ν − j − k)`, `h ≡ j + k`;
- `ε`-splits: a point with `e_x = j` has exactly `j` partners `y = −c ∈ −A` with `χ(x+c) = −1`;
  for `x, y ∈ −A` the two signs satisfy `ε' = χ(−1)ε`;
- `Γ^{(t)}_{ℓk} := Σ_{x∈L_ℓ, y∈L_k} ⟨ψ_t(x),ψ_t(y)⟩` (a linear function of `P`) and
  `g^{(t)}_ℓ = Σ_{x∈L_ℓ} e_t(s(x)) = ν_ℓ K_t(j; m)` or `ν_ℓ K_t(j; m−1)` satisfy
  `[[Γ^{(t)}, g^{(t)}],[g^{(t)ᵀ}, C(m,t)]] ⪰ 0` for every `t` (Gram of the level sums and the all-ones
  vector);
- exact rows: `Γ^{(1)}1 = 0` (because `Σ_x ψ_1(x) = 0`) and `Γ^{(2)}1 = −g^{(2)}` (because
  `Σ_x ψ_2(x) = −1`);
- Weil rows: `1ᵀΓ^{(t)}1 = Σ_{|T|=t} S_T² ≤ C(m,t)(t−1)²p` for `t ≥ 3`.

P5 constrains the distribution of *pairs* of points; it is not of the form treated in §2, and it implies
the Krawtchouk–Weil rows of Proposition 2.3 (Schur complement at the all-ones direction).

### 1.7 P6: ϑ-type programs and localisation

**Bipartite ϑ.** `ϑ_bip(p) := max Σ_{x,y} Y_{x,y'}` over PSD `Y` on `F_p ⊔ F_p'` with
`tr Y_{11} = tr Y_{22} = 1` and `Y_{x,y'} = 0` whenever `χ(x+y) = −1`. A biclique gives the feasible point
`vvᵀ`, `v = (1_A/√m, 1_B/√n)`, of value `√(mn)`. **Proposition 1.6 (PROVED).** `ϑ_bip(p) ≤ √p + 1`, so
this program gives `mn ≤ (√p+1)²`. *Proof.* On the support of `Y_{12}`, `J = S + D` with
`D_{xy} = [x+y = 0]`; `2⟨Y_{12}, S⟩ = ⟨Y, [[0,S],[S,0]]⟩ ≤ √p·tr Y = 2√p` since `‖S‖ = √p`, and
`⟨Y_{12}, D⟩ = Σ_x Y_{x,(−x)'} ≤ Σ_x √(Y_{xx}Y_{(−x)'(−x)'}) ≤ 1`. ∎ Numerically (verifier §H, after
symmetrising over `x ↦ qx + t`, `y ↦ qy − t`, `q ∈ Q`, which reduces the program to eight variables and
three `2×2` conditions) `ϑ_bip(p)²/p = 1.0000000` at `p = 101, 1009, 10009, 100049`, with or without
Schrijver's nonnegativity: the ϑ-type program for bicliques is exactly at the Chung level.

**Localisation** (CITED, §5): for cliques, `ω(G_p) ≤ 1 + ϑ(local graph on Q)` (MMP) and
`ω(G_p) ≤ 2 + ϑ(G_{p,2})` (Kunisky). These act on the Paley adjacency *among the common neighbours*, not on
the profile, and are discussed in §4–5.

## 2. Decoupling: what P1, P3, P4 add to the count constraints (task 2, theory)

**Theorem 2.1 (PROVED).** Let `ν_1,…,ν_N > 0` with `Σν_ℓ = p`. Let `𝒳` be a set of one or two symbols
(`A`, or `A` and `B`), and for `X ∈ 𝒳` let `m_X ≥ 1`, placements `0 ≤ a^X_ℓ ≤ ν_ℓ` with
`Σ_ℓ a^X_ℓ = m_X`, and values `f^X_ℓ ∈ R`. If `𝒳 = {A, B}` let
`s ∈ [Σ_ℓ max(0, a^A_ℓ + a^B_ℓ − ν_ℓ), Σ_ℓ min(a^A_ℓ, a^B_ℓ)]`, and put `s_{XX} = m_X`,
`s_{AB} = s_{BA} = s`. Assume

1. `Σ_ℓ ν_ℓ f^X_ℓ = 0` for each `X`;
2. `Σ_ℓ ν_ℓ f^X_ℓ f^Y_ℓ = p s_{XY} − m_X m_Y` for all `X, Y ∈ 𝒳`;
3. `Σ_ℓ a^B_ℓ f^A_ℓ = Σ_ℓ a^A_ℓ f^B_ℓ` (if `𝒳 = {A, B}`).

Then there are a finite-dimensional real inner-product space `H`, vectors `e_ℓ, 1_X ∈ H` with
`⟨e_ℓ, e_k⟩ = ν_ℓδ_{ℓk}`, `⟨e_ℓ, 1_X⟩ = a^X_ℓ`, `⟨1_X, 1_Y⟩ = s_{XY}`, and a symmetric operator `S̃` on
`H` with `S̃1 = 0` (`1 := Σ_ℓ e_ℓ`), `S̃² = pI − J` (`J = 1 1ᵀ`, so `⟨1,1⟩ = p`) and
`S̃1_X = Σ_ℓ f^X_ℓ e_ℓ` for each `X`. Consequently (T1)–(T4) of Lemma 1.3 hold for the compression
`T̃_{ℓk} = ⟨e_ℓ, S̃e_k⟩` in place of `T`, and `Gram_{I−J/p±S̃/√p}(w) ⪰ 0` for every family `w` in
`span{e_ℓ, 1_X}` (P1). `T̃` is an explicit rational matrix (formula in Step 4).

*Proof.* **Step 1 (model).** In `L²(Ω)`, `Ω = ⊔_ℓ [0, ν_ℓ)` with Lebesgue measure, let `e_ℓ` be the
indicator of the `ℓ`-th interval and `1_X` the indicator of a subset `E^X_ℓ` of measure `a^X_ℓ` in each
interval, chosen for `𝒳 = {A,B}` with `Σ_ℓ |E^A_ℓ ∩ E^B_ℓ| = s` (possible by the range condition, cell by
cell). Let `W = span{e_ℓ, 1_X}`.

**Step 2 (the invariant block).** Put `u_X = 1_X − (m_X/p)1` and `g_X = Σ_ℓ f^X_ℓ e_ℓ`. Then
`⟨u_X, 1⟩ = 0` and, by 1, `⟨g_X, 1⟩ = 0`. Let `V_0 = span{u_X, g_X : X ∈ 𝒳}` and define `S_0` on the
generators by `u_X ↦ g_X`, `g_X ↦ p u_X`. For all generators `w, w'`:

- `⟨S_0w, S_0w'⟩ = p⟨w, w'⟩`: for `(u_X, u_Y)` this is `⟨g_X, g_Y⟩ = p s_{XY} − m_Xm_Y = p⟨u_X, u_Y⟩`
  (hypothesis 2); for `(u_X, g_Y)` it is `⟨g_X, p u_Y⟩ = p⟨u_X, g_Y⟩`, i.e. `⟨g_X, u_Y⟩ = ⟨u_X, g_Y⟩`,
  which is trivial for `X = Y` and is hypothesis 3 for `X ≠ Y` (`⟨g_A, u_B⟩ = Σ_ℓ a^B_ℓ f^A_ℓ` because
  `⟨g_A, 1⟩ = 0`); for `(g_X, g_Y)` it is `p²⟨u_X,u_Y⟩ = p⟨g_X,g_Y⟩`, hypothesis 2 again;
- `⟨S_0w, w'⟩ = ⟨w, S_0w'⟩`: the same identities.

The first shows that `S_0` is well defined (if `Σc_iw_i = 0` then `‖Σc_iS_0w_i‖² = p‖Σc_iw_i‖² = 0`) and
maps `V_0` onto `V_0`; the second that it is symmetric; and `S_0² = p` on the generators.

**Step 3 (extension).** Let `W_1 = W ⊖ (V_0 ⊕ R1)` and `H = W ⊕ W_1'` with `ι: W_1 → W_1'` an isometric
copy. Put `S̃ = S_0` on `V_0`, `S̃1 = 0`, `S̃w = √p ι(w)` and `S̃ι(w) = √p w` for `w ∈ W_1`. Then `S̃` is
symmetric, `S̃² = p` on `V_0 ⊕ W_1 ⊕ W_1' = 1^⊥` and `S̃1 = 0`, i.e. `S̃² = pI − J`; and
`S̃1_X = S̃u_X = g_X`.

**Step 4 (compression).** `S̃` maps `W_1` into `W_1' ⊥ W`, so for `w, w' ∈ W`,
`⟨w, S̃w'⟩ = ⟨Pw, S_0Pw'⟩` with `P` the orthogonal projection onto `V_0`. With a basis `(w_i)` of `V_0`
taken from the generators, `Pw = Σ_i c_i(w)w_i` where `c(w) = G_0^{−1}(⟨w_i, w⟩)_i`, `G_0 = [⟨w_i,w_j⟩]`,
and `⟨w, S̃w'⟩ = c(w)ᵀ[⟨w_i, S_0w_j⟩]c(w')`: all entries are rational in the data. `(I ± S̃/√p)/2` restricted
to `1^⊥` are orthogonal projections, which gives P1; (T1) is `⟨e_ℓ, S̃1⟩ = 0`, (T2) is
`Σ_k T̃_{ℓk}f^X_k = ⟨e_ℓ, S̃g_X⟩ = p⟨e_ℓ, u_X⟩ = p a^X_ℓ − m_Xν_ℓ`, (T3)–(T4) are the hypotheses. ∎

**Corollary 2.2 (PROVED).** (a) *One-sided.* For every nonnegative real profile `(n_j, n'_j)` satisfying
`total`, `moment1`, `moment2`, and every placement of `A` in it, the level-set constraints (T1)–(T4) and P1
for `1_A` are satisfiable (take the cells to be the levels, `𝒳 = {A}`). Hence adding P1 + P3 to the theorem
rows does not change the feasible set of profiles. For cliques the ϑ-bound `m ≤ √p` is already implied by
`moment1`, `moment2` and nonnegativity: Cauchy–Schwarz on the complement of `C` gives
`m(p−m) = Σ f² ≥ m(m−1)² + (m(m−1))²/(p−m)`, i.e. `(p−m)² ≥ p(m−1)²`. (b) *Two-sided.* Given profiles for both sides, P1 + P3 are satisfiable as soon as a
joint profile satisfying 1–3 exists; a joint profile is a coupling of the two profiles with the cells of
`B` and of `A` prescribed and one linear condition `Σ ν f^A f^B = ps − mn`, which holds for a convex
combination of the product and the anti-monotone couplings whenever the target lies between their values.
This is verified exactly in every instance of §3 (general existence not claimed).

*Remarks.* (1) The operator `S̃` uses nothing about primes: the theorem holds verbatim with `p` replaced by
any `q`, and over `F_q`, `q = p²`, the actual Paley matrix satisfies `S² = qI − J`, while the subfield
`F_p` is a clique of size `√q` that attains `ϑ = √q`. So P1–P3 are field-agnostic in the sense of barrier
B4 of `research/sigma-barriers-2026-09-05.md`. (2) Verifier §D applies the construction to the level data
of 32 real configurations and checks (T1), (T2), the `S̃1_X` column, symmetry and the isometry identities
exactly, and `G_±` PSD numerically.

**Proposition 2.3 (P4 aggregated; PROVED).** Let `Z̄_k = (1/m!)Σ_π P_πZ_kP_πᵀ` (`π` permuting the labels of
`A`). Then `Z̄_k = Σ_j n_j E_{U_j}[ψψᵀ] + Σ_j n'_j E_{U'_j}[ψψᵀ]`, where `U_j` is uniform on the sign
vectors of weight `j` and `U'_j` uniform on the vectors with one zero coordinate (uniform position) and `j`
minus signs among the others. In particular `Z̄_k ⪰ 0` for **every** nonnegative profile. Its exact entries
(`|T△T'| ≤ 2`) are equivalent to `total`, `negA`, `moment1`, `moment2`, and the Weil bounds on its entries
(`|T△T'| = t ≥ 3`) are equivalent to the **Krawtchouk–Weil rows**

    KW_t:   | Σ_j n_j K_t(j; m) + Σ_j n'_j K_t(j; m−1) |  ≤  C(m,t)(t−1)√p,      3 ≤ t ≤ m.

*Proof.* `Z_k = Σ_x ψ_k(s(x))ψ_k(s(x))ᵀ` with `s(x) = (χ(x+a))_a`, so
`Z̄_k = Σ_x E_π[ψ(πs(x))ψ(πs(x))ᵀ]`, and `πs(x)` is uniform on the `S_m`-orbit of `s(x)`, which is
`U_j` or `U'_j` according to the level of `x`. By Lemma 1.4(a) each entry is `S_{T△T'}` minus a
correction that involves only the points of `−A`, whose orbit average is a linear function of `(n'_j)`; so a
bound on the orbit average of the entries is the same as a bound on the average of `S_U` over `|U| = t`, and
`Σ_{|U|=t} S_U = Σ_x e_t(s(x))` equals the left side of `KW_t` (`e_t` of a `±1` vector with `j` minus
signs is `K_t(j; ·)`). For `t = 0, 1, 2` the averages are `p`, `0`, `−1`, which are `total`, `moment1` and
(with `total`, `negA`, since `e_2(s) = (f² − |s|²)/2` and `|s(x)|² = m − δ_x`) `moment2`. ∎ (Verifier §C: the symmetrisation identity exactly in 56 cases, `m ≤ 4`, with
exact `LDLᵀ` of `Z̄_2`; the `KW_t` identity and bound on 341 actual `(A, t)`, max ratio `0.60`.)

So, in the aggregated profile variables, **P1 + P3 + P4 = theorem rows + KW rows**, for all `p` (one-sided).

## 3. The hybrid program and its optimum (task 2, numbers)

**The program.** `HYB(p, m)`: maximise `n_0` over profiles satisfying the theorem rows, the `KW_t` rows
(`3 ≤ t ≤ m`), `r = 0` (`n'_0 = 0`), plus P1 and P3 (placement variables for `A` and `B`, the compression
`T`, the Gram PSD conditions). Its clique version fixes `n'_0 = m` (the clique `C = −A'` at level 0),
`n_0 = 0`, and asks for feasibility at `m`.

By Corollary 2.2(a) and Proposition 2.3 the one-sided `HYB` optimum equals the LP optimum of theorem rows +
KW rows, and that is at most `d/m` (the `e = 0` row of (★) is Hanson–Petridis). Independently, the
level-set SDP of §1.4 solved numerically at `p = 1009`, `m = 22` (verifier §F; 45 levels, free symmetric
`T`, Schur-form Gram conditions; SCS) has optimum `22.909090`, equal to the LP optimum `504/22`.

**Theorem 3.1 (PROVED by exact certificates).** For each row of Table 3.1 there are nonnegative rationals
`(n_j, n'_j)` with the listed `n_0` satisfying every theorem row of the balanced note and every `KW_t` row
exactly, with all mass off `B` in the band `0.4m ≤ e_x ≤ 0.6m`; for bicliques a second certificate for the
`B`-side (roles of `m, n` exchanged), a joint profile satisfying hypotheses 1–3 of Theorem 2.1 (`s = 0`),
placements of `A` and `B` at band levels, an exact check of P1 in `Q(√p)`, and the rational compression
`T̃` with (T1), (T2) and the `S̃1_X` columns verified exactly. Hence the optimum of the one-sided `HYB` lies
in `[⌊d/m⌋, d/m]`, and the two-sided P1+P3 constraints do not exclude the certificate.

**Table 3.1** (verifier §E; `n_0` = certified number of complete points, `KW` = number of
Krawtchouk–Weil rows; for bicliques `n = n_0`, `r = 0`).

| p | kind | m | n_0 | mn/p | LP+KW optimum of n_0 (float) | d/m | rows (theorem + KW + fixings) | support off B (levels j of e_x) |
|---|---|---|---|---|---|---|---|---|
| 1009 | biclique | 22 | 22 | 0.4797 | 22.9091 | 22.9091 | 197 | 9, 10, 13; `n'_13 = 22` |
| 10009 | biclique | 71 | 70 | 0.4966 | 70.4789 | 70.4789 | 621 | 34, 35, 42; `n'_42 = 71` |
| 40009 | biclique | 141 | 141 | 0.4969 | 141.8723 | 141.8723 | 1223 | 69, 70, 84; `n'_84 = 141` |
| 100049 | biclique | 224 | 223 | 0.4993 | 223.3214 | 223.3214 | 1935 | 110–112, 134; `n'_134 = 224` |
| 1000033 | biclique | 707 | 707 | 0.4998 | 707.2362 | 707.2362 | 6089 | 352, 353, 424; `n'_424 = 707` |
| 1009 | clique | 22 | — | m²/p = 0.4797 | — | m(m−1) = 462 ≤ 504 | 202 | 9, 10, 13 |
| 10009 | clique | 71 | — | 0.5036 | — | 4970 ≤ 5004 | 635 | 34, 35, 42 |
| 40009 | clique | 141 | — | 0.4969 | — | 19740 ≤ 20004 | 1251 | 69, 70, 84 |
| 100049 | clique | 224 | — | 0.5015 | — | 49952 ≤ 50024 | 1980 | 110–112, 134 |
| 1000033 | clique | 707 | — | 0.4998 | — | 499142 ≤ 500016 | 6231 | 352, 353, 424 |

(For cliques `m = ⌊(1+√(2p−1))/2⌋`, the Hanson–Petridis maximum; `m + 1` violates the `e = 0` row.) The
`KW_t` rows never bind: over the ten certificates the largest ratio of a `KW_t` left side to its bound
`C(m,t)(t−1)⌊√p⌋` is `0.75` (at `p = 100049`); the certificates maximise a common relative slack of the
theorem rows, which is `0.234, 0.245, 0.248, 0.249, 0.250`, and at `p = 1009, 10009` the LP solution is the
same with and without the KW rows. The joint couplings mix the product and the
anti-monotone couplings with weights `λ = 0.29, 0.10, 0.14, 0.09, 0.14`; `dim V_0 = 4` in every case.

**Answer to task 2.** No. The hybrid optimum is within `m` of `(p−1)/2 + r` (it equals `d/m` up to the
integer part), not below it by a constant fraction of `p`.

**The two-point program P5.** Fixing the profile of Table 3.1, P5 (§1.6) becomes an SDP in the pair counts
`P`. Solved with Clarabel (maximising the smallest eigenvalue of the augmented Gram matrices on the
complement of the directions forced to be null by the exact rows, and imposing the Weil rows with a 0.1 %
margin), then rounded to denominator `10⁶`, repaired to satisfy every linear equality exactly (minimum-norm
rational correction), and checked in exact arithmetic: `P ≥ 0`, all marginals, `ε`-splits, both exact row
families, all Weil rows, and `[[Γ^{(t)}, g^{(t)}],[g^{(t)ᵀ}, C(m,t)]] ⪰ 0` for **every** `1 ≤ t ≤ m` by an
exact rational `LDLᵀ` (zero pivots allowed only with zero rows).

**Table 3.2** (verifier §G).

| p | m | kind | levels | pair variables | SDP margin | exact certificate | max Φ_t / Weil bound |
|---|---|---|---|---|---|---|---|
| 1009 | 22 | biclique | 5 | 137 | 3.7e-06 | yes (ranks 1–5) | 0.869 |
| 1009 | 22 | clique | 4 | 68 | 9.7e-06 | yes (ranks 1–4) | 0.999 |
| 10009 | 71 | biclique | 5 | 432 | 5.1e-07 | yes (ranks 1–5) | 0.999 |
| 10009 | 71 | clique | 4 | 203 | 1.0e-07 | yes (ranks 1–4) | 0.993 |
| 40009 | 141 | biclique | 5 | 845 | 1.1e-07 | yes (ranks 1–5) | 0.999 |
| 40009 | 141 | clique | 4 | 392 | 6.2e-09 | yes (ranks 1–4) | 0.994 |
| 100049 | 224 | biclique | 6 | 1939 | 3.8e-10 | yes (ranks 1–6) | 0.997 |
| 100049 | 224 | clique | 5 | 1043 | 1.5e-09 | yes (ranks 1–5) | 0.997 |

("SDP margin" is the smallest eigenvalue, after scaling by `C(m,t)p²`, on the complement of the forced null
directions; "ranks" are the exact ranks of the augmented matrices over `1 ≤ t ≤ m`; the last column is
`max_{t≥3} Φ_t/(C(m,t)(t−1)²p)` for the exact certificate. The margin-maximising solutions push some Weil
row up to the imposed cap `0.999`, i.e. the eigenvalue margin is traded against the Weil rows; the
certificates satisfy every Weil row exactly.)

So P5 does not exclude the tight profiles either. At `p = 1009` P5 is also feasible (numerically, pure
feasibility solve, Clarabel status `optimal`; verifier §G2) for two extreme one-point profiles with exact
theorem-row certificates: the LP profile maximising `Σ j n'_j` (the "`A` nearly
independent" profile of the balanced note's remark) and the one maximising the fourth moment. We found no
profile satisfying the theorem rows that P5 rejects; whether P5 is implied by the one-point constraints is
OPEN.

## 4. The obstruction (task 3)

The hybrid does not give a strict improvement at any tested `p`, so the task is to state exactly what
satisfies everything and why.

**The configuration.** For each row of Table 3.1: the rational profile (mass `n_0` at `e = 0`, the rest on
three or four band levels near `e = m/2`, all of `−A` on one band level); `A` and `B` placed on a band level
(`σ_A, σ_B ∈ {0, m, 2m}`, far inside the P1 range `|σ| ≲ 0.29 m√p`); the operator
`S̃ = S_0 ⊕ (swap dilation)` of Theorem 2.1, where `S_0` acts on the 4-dimensional span of
`u_A, f_A, u_B, f_B` by `u ↦ f ↦ pu`; and, for `p ≤ 100049`, the pair counts of Table 3.2. It satisfies
every theorem row, every Krawtchouk–Weil row, P1, P3 (T1–T4), P4 (all degrees, after symmetrisation) and,
where computed, P5.

**Why.**

1. *Exact information is degree ≤ 2 (description of the inputs).* The only character sums that the
   families P1–P5 evaluate exactly are those whose squarefree part has degree `≤ 2` (Lemma 1.4(b)); for
   degree `≥ 3` they use only Weil's bound (such sums genuinely vary with the polynomial, e.g. point counts
   of elliptic curves). In the aggregated variables the exact identities are `total`, `negA`, `moment1`,
   `moment2`, the joint moment and the P5 rows `t = 1, 2` (Lemma 1.4, Proposition 2.3, Lemma 1.5), and
   Hanson–Petridis-tight profiles satisfy them (they are rows of the LP).
2. *Weil information cannot see `B` (PROVED arithmetic).* In `KW_t` the points of `B` contribute `n·C(m,t)` against the
   allowance `C(m,t)(t−1)√p`; with `n ≤ √(p/2) + 1` this is a fraction `≤ 1/((t−1)√2) + o(1)` of the
   allowance for every `t ≥ 2`. In the P5 Weil rows `B` contributes `n²C(m,t) ≈ (p/2)C(m,t)` against
   `(t−1)²p·C(m,t)`. This is the square-root barrier (B1–B2 of the barrier map) in the form it takes here;
   the balanced note's Weil-moment ratio is the special case `t = 2k`.
3. *PSD structure built from `S² = pI − J` is abstract (PROVED).* Theorem 2.1: any profile with the right first two
   moments is realised by an abstract symmetric `S̃` with the Paley spectrum; P1–P3 cannot distinguish `S̃`
   from `S`, and (Remark 1) they hold verbatim over `F_{p²}`, where `ϑ = √q` is attained by the subfield
   clique. Proposition 2.3: the moment matrices of P4 are automatically PSD after symmetrisation.
4. *The only prime-field input is one-sided and not tight at the tight profiles (PROVED by the
   certificates).* Hanson–Petridis and its
   refinements (★), subset-HP, pencil, derivative-HP are the rows that fail over `F_{p²}`; the balanced note
   and Table 3.1 show they are consistent with the tight profiles, with common relative slack `≈ 1/4` on all
   non-identity rows.

**Precise obstruction.** Let `𝒞` be the class of constraints on `(profile, placements, T, Gram data, pair
counts)` that follow from (a) the degree-≤2 character identities, (b) Weil bounds for squarefree degree `≥ 3`,
(c) positive semidefiniteness of Gram matrices of vectors in the span of level indicators, `1_A`, `1_B` under
`I`, `S`, `Π_±`, or of the lifted characters indexed by points or by pairs of points, and (d) the Stepanov
rows of the balanced note. Every constraint in `𝒞` that we formulated is satisfied by the configuration
above, and for (a)–(c) at the level of single points and level sets this is a theorem (Theorem 2.1,
Proposition 2.3) valid for all `p`. **A hybrid that beats Hanson–Petridis must therefore contain an input
outside `𝒞`:** either a PSD constraint on the Paley adjacency among points of the complement that is not a
consequence of `S² = pI − J` (this is what localisation does: MMP's ϑ on the neighbourhood of `0`, which
reduces to an LP over the multiplicative group and involves Jacobi sums; Kunisky's `ϑ(G_{p,2})`), or a
prime-specific inequality that couples pairs of points or both sides of the biclique (e.g. a two-sided
auxiliary polynomial, obligation 2 of the balanced note). In the language of the balanced note: the target
is an anti-concentration statement for `f_A` at the square-root scale, and nothing in `𝒞` supplies it.

## 5. Literature (task 4)

Fetched 2026-09-29 (curl of the arXiv HTML/PDF/abstract pages; local copies listed at the top).

- **Randomstrasse101** (A. S. Bandeira, D. Dmitriev, K. Lucca, P. Nizić-Nikolac, A. Rödder, *Open
  Problems of 2025*, arXiv:2603.29571), Entry 12. Conjecture 25: `ω(G_p) = O(polylog p)`. After recalling
  `ω(G_p) ≤ √p` (spectral or `ϑ(Ḡ_p)`) and Hanson–Petridis `(1+o(1))√(p/2)`: Conjecture 26,
  `ϑ(Ḡ_{p,1}) ∼ √(p/2)` for the 1-localisation (observed empirically in MMP; "would recover the
  Hanson-Petridis bound"); Conjecture 27, `ϑ(Ḡ_{p,2}) ≤ (2/3)√p` for large `p`; the text before it
  records that Kunisky "observed empirically that ω(G_{p,2})∼(√(1/2)−ϵ)√p" (the symbol `ω` is as printed;
  from Kunisky's paper the quantity is the ϑ-value); Conjecture 28: degree-4 SOS gives `O(p^{1/2−ε})`.
- **Magsino–Mixon–Parshall**, *Linear programming bounds for cliques in Paley graphs*, arXiv:1907.05971
  (SPIE Wavelets and Sparsity XVIII, 2019): bound (3) `ω(G_p) ≤ ϑ_LS(L_p) + 1` for the local graph `L_p` on
  `Q_p`; Propositions 2–3 reduce `ϑ_L`, `ϑ_LS` of the circulant `L_p` to LPs; Proposition 5 is the dual
  certificate ("Then ω(G_p) ≤ f(0) + 1"); for `p < 3000`, `⌊LS(p)⌋ = ⌊HP(p)⌋` for most primes and
  `⌊LS(p)⌋ = ⌊HP(p)⌋ − 1` for 17 primes; Conjecture 6: "For infinitely many primes p ≡ 1 (mod 4), it holds
  that LS(p) < ⌊HP(p)⌋." They also note that the Gvozdenović–Laurent–Vallentin block-diagonal SDP hierarchy
  gives numerical bounds sharper than `HP(p)` for `p ≤ 809`.
- **Kunisky**, *Spectral pseudorandomness and the road to improved clique number bounds for Paley graphs*,
  arXiv:2303.16475 (local copy). Theorem 1.19 (conditional on his Conjecture 1.9, convergence of the
  minimum eigenvalue of the degree-`a` localisations to the Kesten–McKay edge):
  `ω(G_p) ≤ (√(2^a − 1)/2^{a−1})√p + o(√p)`; at `a = 3` this is `≈ 0.661√p` (below Hanson–Petridis);
  Theorem 1.18 proves the conjectures at `a = 1`. Appendix A, Conjecture A.1: for some `ε > 0` and all
  large `p`, `ϑ(G_{p,{0,g}}) ≤ (1/√2 − ε)√p` (2-localisation; numerics for `p ≤ 800`).
- **Kunisky–Yu**, arXiv:2211.02713 (CCC 2023): degree-4 SOS value `≥ Ω(p^{1/3})`; numerics suggest
  `O(p^{1/2−ε})`. **Kobzar–Mody**, arXiv:2304.08615: block-diagonal `L2` relaxations computed for
  `821 ≤ p ≤ 997` (Gvozdenović et al. up to 809), numerical evidence of scaling below `√p`, not conclusive.
  **Gaar–Pucher**, arXiv:2412.12958: the exact subgraph hierarchy stays at `ϑ` up to a threshold level for
  Paley graphs; a vertex-transitive variant does better numerically.
- Searches (2026-09-29): "Hanson-Petridis bound combined semidefinite programming Lovász theta Paley graph
  clique improvement"; "Stepanov method polynomial method combined with sum-of-squares or semidefinite
  relaxation clique number Paley graph 2025 2026"; a quoted variant with "Hanson" "Petridis" and
  "theta"/"semidefinite"/"linear programming". They returned only the items above, Yip
  (arXiv:2004.01175, prime powers) and blog/survey pages. **No work combining Stepanov/Hanson–Petridis
  polynomial constraints with SDP or ϑ constraints was found; completeness is not claimed.** The only
  convex-relaxation results near the Hanson–Petridis constant are the localised ϑ/LP bounds above, which are
  numerical (MMP: tie or beat by 1 for `p < 3000`; Kunisky: below `1/√2` for `p ≤ 800`) and unproved.

## 6. Open obligations

1. **The target (OPEN):** `mn ≤ (1/2 − c)p + o(p)` for balanced complete bicliques; `ω(G_p) ≤ √((1/2−c)p)`.
2. **An input outside `𝒞` (§4).** Candidates: (a) a localised ϑ with Stepanov cuts: the (★)-type rows need
   level counts, i.e. a lift of degree `≥ 3` of the ϑ matrix (Lasserre-type) in which "`e_x ≤ e`" is
   expressible — not attempted; (b) a proof of MMP's Conjecture 6 / Kunisky's Conjecture A.1 via their dual
   programs (MMP Proposition 5), which would be a pure-ϑ route; (c) a two-sided prime-specific inequality.
3. **P5 in general.** Exact certificates exist at `p = 1009, 10009, 40009, 100049`;
   a proof that P5 never excludes profiles satisfying the one-point rows (or an example where it does) is
   OPEN. We found no profile that P5 rejects.
4. **Corollary 2.2(b) in general:** existence of the joint coupling for all pairs of LP profiles (verified
   exactly only in the instances).
5. Independent review of Theorem 2.1, Proposition 2.3 and the verifier (§§A–H).
