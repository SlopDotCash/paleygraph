# Sigma pass, direction `dual`: the multiplicative-Fourier dual of the shifted-subgroup sum and cancellation across cosets (2026-09-05)

**Status:** PROVED: (i) the coset expansion `G·T_H(c) = Σ_{cosets ξH} χ(ξ)Ĥ(ξ)Ĥ_χ(cξ)` for every subgroup `H ≤ F_p^*` and every `c ∈ F_p` (no boundary term, also for `c ∈ −H` and `c = 0`), verified exactly in `Z[ζ_p]` for all 353 subgroups of all primes `p ≤ 200` and all `c` (Theorem 1.3); (ii) the exact `ℓ²` identities `Σ_{cosets}|Ĥ|² = p − |H|`, `Σ_{cosets}|Ĥ_χ(cξ)|² = p − ι_H|H|` with `ι_H = [H ⊆ QR]`, and that no bound using only `ℓ²`/`ℓ^∞` data of the two period families can fall below `√p·(1 − o(1))` for `H ⊆ QR`, nor below `√(p/|H|)(1−o(1))` in general (Props 2.1–2.4); in particular a Bourgain–Glibichuk–Konyagin saving `p^{−β}` on `max|Ĥ|` has `β ≤ ε/2` and the pointwise route is dead for `|H| ≤ √p`; (iii) the exact formula `D(H,c)² = (p·R_H(c) − ι_H|H|⁴)/|H|` for the coset `ℓ²`-weight in terms of a signed dilated additive energy `R_H(c)`, hence the exact relation `ρ(H,c) = (|T_H(c)|/√|H|)·(R_H(c)/|H|² − ι_H|H|²/p)^{−1/2}` (Theorem 3.1) and the sandwich `|T_H(c)|√(|H|/E(H,cH)) ≤ ρ(H,c) ≤ (|T_H(c)|/√|H|)(1 − |H|²/p)^{−1/2}` (Cor. 3.2, upper bound for `H ⊆ QR`); (iv) `sup ρ = ∞`: for every `m ≥ 2` and every prime `p ≡ 1 (mod m)`, `p ≥ 4^{m+1}m⁴`, the subgroup of order `m` has a shift with `ρ ≥ √m` (Theorem 3.3, via Weil); (v) the Gauss-period product relation, the identity `T_H(g^a) = Σ_k (−1)^k (a,k)_d` with cyclotomic numbers, the Jacobi-sum expansion, and the fact that the coset expansion and the Jacobi expansion are Fourier duals on `Z/d` (Props 4.1–4.4). REFUTED (Theorem 3.3 and witnesses): "square-root cancellation across cosets, `ρ(H,c) = O(1)` uniformly"; largest observed `ρ = 3.8006` (`ρ² = 879757/60906`) at `p = 2437, |H| = 29, c = 1226`. CONDITIONAL: `ρ ≤ R` for all `c` implies `M_H ≤ R·√(E(H)/|H|) ≪ R|H|^{3/4}` (Heath-Brown–Konyagin, CITED), so a bound `ρ ≤ p^{o(1)}` for `|H| ≥ p^ε` would settle Bourgain's Problem 5 with `δ = ε/4 − o(1)`. VERDICT: "coset cancellation" is a relabelling — in the multiplicative-dual (Jacobi) basis the cancellation ratio equals `|T_H(c)|/√|H|·(1+O(1/d))` exactly, and in the coset basis it equals `|T_H(c)|/√|H|` times the energy factor `(ε_H(c) − ι|H|²/p)^{−1/2}`, which for `H ⊆ QR` lies in `[C|H|^{−1/4}, (1−|H|²/p)^{−1/2}]` (observed range `[0.46, 1.41]` at the extremal shift, `|H| < √p`); the cyclotomic-number reduction is a tautology (`Σ_k(−1)^k(a,k)_d` is the definition of `T_H(g^a)` regrouped). OPEN: any nontrivial bound on `ρ`, equivalently on `M_H/√|H|`, below `√p`; the needed statement is cancellation in the Mellin-type sum `Σ_{ψ ∈ H^⊥} ψ(u) J(ψ,χ)` over the annihilator subgroup of characters (Statement CN below), for which no source examined gives anything. HEURISTIC (data, `p ≤ 3000`): `M_H/√(2|H| log d)` has median `0.84`, max `1.35`, consistent with a Gaussian model `T_H(c)/√|H| ~ N(0,1)` across cosets.

Verifier: `experiments/sigma_dual_2026_09_05.py` (310,324 checks, 0 failures, 16 s) → `results/sigma_dual_2026_09_05.json`.

## 0. Notation

`p` an odd prime, `χ` the Legendre symbol with `χ(0) = 0`, `e(x) = exp(2πi x/p)`, `ζ = e(1)`. `g` a primitive root. `H ≤ F_p^*` a subgroup, `f := |H|`, `d := (p−1)/f` the index, `H = ⟨g^d⟩`. `ι = ι_H := 1` if `H ⊆ QR` (equivalently `d` even, equivalently `f | (p−1)/2`), else `0`. `QR` = nonzero squares.

    T_H(c) := Σ_{h∈H} χ(h+c),   M_H := max_{c≠0} |T_H(c)|   (as in the `structured` note),
    G := Σ_{x∈F_p} χ(x) e(x),   |G| = √p,
    Ĥ(ξ) := Σ_{h∈H} e(ξh),      Ĥ_χ(η) := Σ_{h∈H} χ(h) e(ηh)   (the "twisted period"),
    η_j := Ĥ(g^j)  (j ∈ Z/d),  the Gauss periods of order `d`.

Cosets of `H` are `ξ_j H = g^j H`, `j ∈ Z/d`. "Σ_{cosets}" means a sum over one representative `ξ` of each coset; every summand below is shown to be independent of the representative.

Bourgain's Problem 5 (Chang 2010, quoted in the `structured` note, §1): bound `M_H` nontrivially for `|H| ∼ √p`; by Theorem 2.2 of the `structured` note, the subgroup case of the two-set conjecture is equivalent to `M_H ≤ p^{−δ}|H|` for all `|H| > p^ε`.

Additive-energy notation: `r(t) := #{(h,h') ∈ H²: h − h' = t}`, `r^χ(t) := Σ_{h−h'=t} χ(h)χ(h')` (so `|r^χ| ≤ r`, and `r^χ = r` when `H ⊆ QR`), `E(H) := Σ_t r(t)²`, `E(H,cH) := Σ_t r(−ct) r(t) = #{(h_1,h_2,h_3,h_4) ∈ H⁴ : h_1 + ch_3 = h_2 + ch_4}`.

## 1. The coset identity (Task 1) — PROVED

**Lemma 1.1.** For every `x ∈ F_p` (including `x = 0`): `χ(x) = (1/G) Σ_{ξ∈F_p} χ(ξ) e(ξx)`.

*Proof.* For `x ≠ 0` substitute `ξ = η/x`: `Σ_ξ χ(ξ)e(ξx) = Σ_η χ(η/x) e(η) = χ(x) Σ_η χ(η)e(η) = χ(x)G` (as `χ(1/x) = χ(x)`). For `x = 0` the left side is `Σ_ξ χ(ξ) = 0 = χ(0)G`. `G ≠ 0` since `|G|² = p` (standard). ∎

**Lemma 1.2 (coset covariance).** For `h ∈ H`: `Ĥ(hξ) = Ĥ(ξ)`, `Ĥ_χ(hη) = χ(h) Ĥ_χ(η)`, and `χ(hξ) = χ(h)χ(ξ)`. Hence `χ(ξ)Ĥ(ξ)Ĥ_χ(cξ)` depends only on the coset `ξH` (replacing `ξ` by `hξ` multiplies it by `χ(h)² = 1`), and `|Ĥ|`, `|Ĥ_χ|` are constant on cosets.

*Proof.* `h' ↦ hh'` permutes `H`; `Ĥ_χ(hη) = Σ_{h'} χ(h') e(ηhh') = Σ_{h''} χ(h''/h) e(ηh'') = χ(h)Ĥ_χ(η)`. ∎

**Theorem 1.3 (coset expansion).** For every `H ≤ F_p^*` and every `c ∈ F_p`,

    G · T_H(c) = Σ_{ξ ≠ 0} χ(ξ) e(ξc) Ĥ(ξ) = Σ_{cosets ξH} χ(ξ) Ĥ(ξ) Ĥ_χ(cξ).

No boundary term arises: the case `c ∈ −H` (where the summand `h = −c` of `T_H(c)` is `χ(0) = 0`) and the case `c = 0` (where `T_H(0) = ι f`) are covered by the same formula.

*Proof.* By Lemma 1.1 applied at `x = h + c` for each `h ∈ H` — valid also when `h + c = 0`, which is exactly why no correction term appears —
`G·T_H(c) = Σ_h Σ_ξ χ(ξ) e(ξ(h+c)) = Σ_ξ χ(ξ) e(ξc) Σ_h e(ξh) = Σ_{ξ≠0} χ(ξ)e(ξc)Ĥ(ξ)` (the term `ξ = 0` vanishes because `χ(0) = 0`). Group `ξ ≠ 0` by cosets: for `ξ = ξ_0 h`, `h ∈ H`, `Ĥ(ξ_0h) = Ĥ(ξ_0)` (Lemma 1.2), so
`Σ_{h∈H} χ(ξ_0h) e(ξ_0hc) Ĥ(ξ_0h) = χ(ξ_0)Ĥ(ξ_0) Σ_h χ(h) e((cξ_0)h) = χ(ξ_0)Ĥ(ξ_0)Ĥ_χ(cξ_0)`. By Lemma 1.2 the last expression does not depend on the representative `ξ_0`. ∎

*Remarks.* (a) For `H ⊆ QR`, `Ĥ_χ = Ĥ` and the identity reads `G·T_H(g^a) = Σ_{j∈Z/d} (−1)^j η_j η_{j+a}`, `G = Σ_j (−1)^j η_j`: a sign-twisted autocorrelation of the Gauss periods. (b) For `H ⊄ QR` (`d` odd), `f` is even, `K := H ∩ QR = ⟨g^{2d}⟩` has index 2 in `H`, and `Ĥ_χ(η) = K̂(η) − K̂(g^dη)`; also `T_H(c) = T_K(c) − T_K(c g^{−d})` (from `T_{sK}(c) = χ(s)T_K(c/s)`), so `M_H ≤ 2M_K` and the case `H ⊆ QR` suffices for the conjecture up to a factor 2. (c) (Exact inversion symmetry, `H ⊆ QR`.) `T_H(c) = χ(c) T_H(1/c)` for `c ≠ 0`: `T_H(c) = Σ_h χ(h)χ(1 + c/h) = Σ_{h'} χ(1 + ch') = χ(c) Σ_{h'} χ(1/c + h')`. Thus the exceptional shifts of a subgroup in `QR` are closed under `c ↦ 1/c`.

**Verification (exact).** The verifier represents elements of `Z[ζ]` as integer vectors on `ζ^0,…,ζ^{p−1}` reduced to the basis `ζ^1,…,ζ^{p−1}` by `1 + ζ + … + ζ^{p−1} = 0`; the product `Ĥ(ξ)Ĥ_χ(cξ)` is accumulated from the `f²` exponents `ξh + cξh'` with weights `χ(h')`. This arithmetic is cross-validated against sympy's remainder modulo the cyclotomic polynomial `Φ_p` (12 random products at `p ∈ {7,11,13}`). Theorem 1.3 is then checked for all 353 subgroups of all 45 primes `3 ≤ p ≤ 199` and all `c ∈ F_p`: 37,853 exact equalities in `Z[ζ]`, of which 8,965 have `c ∈ −H` and 353 have `c = 0`. The same identity is re-checked in floating point (`|·| < 10^{−7}p`) at 2,104 `(p,H,c)` triples in §2, and the coset invariance of `|Ĥ|`, `|Ĥ_χ|` at all 5,536 subgroups of all primes `≤ 3000`. The inversion symmetry (c) is checked exactly at all 23,298 pairs `(H ⊆ QR, c ≠ 0)`, `p ≤ 199`.

## 2. The `ℓ²` obstruction (Task 2) — PROVED

**Proposition 2.1 (exact `ℓ²` identities).**
(a) `Σ_{ξ∈F_p}|Ĥ(ξ)|² = pf`, `Σ_{ξ≠0}|Ĥ(ξ)|² = f(p−f)`, `Σ_{cosets}|Ĥ(ξ)|² = p − f`.
(b) `Σ_{η∈F_p}|Ĥ_χ(η)|² = pf`, `Ĥ_χ(0) = Σ_h χ(h) = ιf`, `Σ_{η≠0}|Ĥ_χ(η)|² = pf − ιf²`, and for every `c ≠ 0`: `Σ_{cosets}|Ĥ_χ(cξ)|² = p − ιf`.
(c) `Σ_{c≠0} T_H(c)² = pf − f² − ιf²`.

*Proof.* (a), (b): Parseval `Σ_ξ|Σ_x a(x)e(ξx)|² = pΣ_x|a(x)|²` for `a = 1_H` and `a = χ·1_H`; subtract the `0`-term (`Ĥ(0) = f`, `Ĥ_χ(0) = ιf` since `Σ_h χ(h) = f` if `H ⊆ QR` and `0` otherwise, `H ∩ QR` having index 2); each coset carries `f` equal values, and `ξ ↦ cξ` permutes the cosets. (c): `Σ_{c∈F_p} T_H(c)² = Σ_{h,h'} Σ_c χ(c+h)χ(c+h') = f(p−1) − f(f−1) = pf − f²` (background fact `Σ_x χ(x)χ(x+t) = −1`, `t ≠ 0`), and `T_H(0)² = ι f²`. ∎

The multisets `{|Ĥ(ξ)|}_{cosets}` and `{|Ĥ_χ(cξ)|}_{cosets}` do not depend on `c`, and for `H ⊆ QR` they coincide (verified at all 5,536 subgroups, `p ≤ 3000`; exact `Z[ζ]` check of (a),(b) at all 353 subgroups, `p ≤ 199`; (c) exact at all 5,536).

**Corollary 2.2 (Cauchy–Schwarz over cosets).** `|T_H(c)| ≤ (1/√p)·√(p−f)·√(p−ιf) < √p`; for `H ⊆ QR`, `|T_H(c)| ≤ (p−f)/√p`. This is the Weil-range bound (Gong's Theorem 2 / the elementary `√(p−f)` of the `structured` note, Prop. 3.3) and nothing more.

**Proposition 2.3 (norm-only information cannot beat `√p`).** Let `Φ: [0,∞)⁴ → [0,∞)` be any function such that `|Σ_{j=1}^d x_j y_j| ≤ Φ(‖x‖₂, ‖y‖₂, ‖x‖_∞, ‖y‖_∞)` for all `x, y ∈ C^d`. Then for all admissible data `(α,β,A,B)` (i.e. `A ≤ α ≤ A√d`, `B ≤ β ≤ B√d`):

    Φ(α,β,A,B) ≥ min(α²B/A, β²A/B) − AB.

Applied to `a_j := χ(ξ_j)Ĥ(ξ_j)`, `b_j := Ĥ_χ(cξ_j)` (`α² = p−f`, `β² = p−ιf`, `A = A_H := max_{ξ≠0}|Ĥ(ξ)|`, `B = B_H := max_{η≠0}|Ĥ_χ(η)|`), every bound on `|T_H(c)| = |Σ_j a_jb_j|/√p` that uses only these four numbers is at least
`[min(A/B, B/A)(p−f) − f²]/√p`, which for `H ⊆ QR` (`A = B`) equals `(p − f − f²)/√p ≥ √p(1 − 2f²/p)`, and in general is `≥ √(p/f)(1 − o(1)) − f²/√p` because `A, B ∈ [√(f(p−f)/(p−1)), f]`.

*Proof.* Put `m := ⌊α²/A²⌋ ≥ 1`, `m' := ⌊β²/B²⌋ ≥ 1`. Let `x` have `m` coordinates equal to `A`, one coordinate `√(α² − mA²) ∈ [0, A)` (omitted if `m = d`), zeros elsewhere; then `‖x‖₂ = α`, `‖x‖_∞ = A`. Define `y` likewise with `B, β`, on the same leading coordinates. Then `Σ x_jy_j ≥ min(m,m')AB ≥ min(α²/A² − 1, β²/B² − 1)AB`. The last claims use `A ≤ f`, `B ≤ f`, and the pigeonhole bounds of Prop. 2.4(b). ∎

Numerically (all 5,536 subgroups, `p ≤ 3000`): `min_H [min(α²B/A, β²A/B)/p] = 0.4707` (at `p = 1009`, `|H| = 144`, `d = 7`, `A = 14.56`, `B = 30.94`), i.e. the norm-only bound is `≥ 0.47·√p` in every case examined, and exactly `(p−f)/√p` whenever `H ⊆ QR`.

**Proposition 2.4 (the pointwise route and its ceiling).**
(a) For all `H`, `c ≠ 0`: `|T_H(c)| ≤ A_H·√(d(p−ιf)/p) < A_H√d`, and `|T_H(c)| ≤ B_H·√(d(p−f)/p)`.
(b) (Pigeonhole.) `A_H ≥ √(f(p−f)/(p−1))` and `B_H ≥ √((pf − ιf²)/(p−1))`; both are `≥ √f·(1 − f/p)^{1/2}`.
(c) Hence the bound in (a) is never below `√((p−f)(p−ιf)/p) ≥ √p(1 − f/p)`, whatever estimate for `A_H` is inserted. Quantitatively: if `|H| = p^ε` and `A_H ≤ p^{−β}|H|` (Bourgain–Glibichuk–Konyagin; CITED: Kurlberg 2007, Theorem 1.1, `sources/kurlberg-small-subgroups-2007.txt` lines 74–77: "Given α > 0, there exists β = β(α) > 0 such that if |H| > p^α … |Σ_{x∈H} ψ(x)| ≪ p^{−β}|H|"), then necessarily `β ≤ ε/2 + o(1)` by (b), and (a) gives `|T_H(c)| ≤ p^{−β}√(f(p−1))`, which is `≤ p^{−δ}f` only if `f ≥ p^{1+2δ−2β} ≥ p^{1+2δ−ε}`, i.e. only if `ε ≥ 1/2 + δ`. At `|H| = √p` one would need `β > 1/4 = ε/2`, which (b) forbids. The pointwise route is therefore dead for every `|H| ≤ √p`, and even the best conceivable pointwise input `A_H ≈ √f` reproduces exactly the Weil bound `√p`.

*Proof.* (a) From Theorem 1.3, `√p|T_H(c)| = |Σ_j a_jb_j| ≤ ‖a‖_∞‖b‖_1 ≤ A_H√d‖b‖₂ = A_H√(d(p−ιf))` by Cauchy–Schwarz and Prop. 2.1(b); symmetrically with the roles exchanged. (b) is Prop. 2.1(a),(b) divided by `p−1`. (c) Insert (b) into (a): `A_H√(d(p−ιf)/p) ≥ √(f(p−f)/(p−1))·√((p−1)(p−ιf)/(fp)) = √((p−f)(p−ιf)/p)`. ∎

(Both pigeonhole bounds and the mixed bound of (a) are checked at all 5,536 subgroups.)

## 3. The coset-cancellation ratio (Task 3)

**Definition.** For `c ≠ 0`,

    D(H,c)² := Σ_{cosets} |Ĥ(ξ)|² |Ĥ_χ(cξ)|²,      ρ(H,c) := |Σ_{cosets} χ(ξ)Ĥ(ξ)Ĥ_χ(cξ)| / D(H,c) = √p·|T_H(c)| / D(H,c).

(`D > 0` always: a Gauss period `Σ_{h} ζ^{ξh}` or twisted period `Σ_h χ(h)ζ^{ηh}` is a nontrivial `Z`-combination of distinct elements of the `Q`-basis `ζ^1,…,ζ^{p−1}` of `Q(ζ)`, hence nonzero.) "Square-root cancellation across cosets" is the hypothesis `ρ(H,c) = O(1)`.

**Theorem 3.1 (exact energy formula).** For `c ≠ 0` let `R_H(c) := Σ_{t∈F_p} r(−ct) r^χ(t) = f² + Σ_{t≠0} r(−ct)r^χ(t)`. Then

    D(H,c)² = (p·R_H(c) − ι f⁴)/f,   hence   ρ(H,c)² = f·T_H(c)² / (R_H(c) − ιf⁴/p),

i.e. `ρ(H,c) = (|T_H(c)|/√f)·(ε_H(c) − ιf²/p)^{−1/2}` with `ε_H(c) := R_H(c)/f²`. Moreover `|R_H(c)| ≤ E(H,cH) ≤ E(H)`, and `Σ_{c≠0}(R_H(c) − f²) = (ιf² − f)(f² − f)`, so `avg_{c≠0} ε_H(c) = 1 + (ι f − 1)(f−1)/(p−1)`.

*Proof.* `|Ĥ(ξ)|² = Σ_{h,h'} e(ξ(h−h')) = Σ_s r(s)e(ξs)` and `|Ĥ_χ(cξ)|² = Σ_t r^χ(t) e(cξt)`. Summing over all `ξ ∈ F_p`: `Σ_ξ |Ĥ(ξ)|²|Ĥ_χ(cξ)|² = Σ_{s,t} r(s)r^χ(t) Σ_ξ e(ξ(s+ct)) = p Σ_{s+ct=0} r(s)r^χ(t) = p R_H(c)`. The term `ξ = 0` equals `f²·(ιf)² = ιf⁴`, and each coset contributes `f` equal terms; so `f·D² = pR_H(c) − ιf⁴`. The `t = 0` term of `R_H` is `r(0)r^χ(0) = f·Σ_hχ(h)² = f²`. `|R_H(c)| ≤ Σ_t r(−ct)r(t) = E(H,cH)` since `|r^χ| ≤ r`; `E(H,cH) ≤ √(E(H)E(cH)) = E(H)` by Cauchy–Schwarz in `t` and `r_{cH−cH}(t) = r(t/c)`. Finally `Σ_{c≠0} Σ_{t≠0} r(−ct)r^χ(t) = Σ_{t≠0} r^χ(t)·Σ_{s≠0} r(s) = (Σ_{t≠0}r^χ(t))(f²−f)` and `Σ_{t≠0} r^χ(t) = (Σ_hχ(h))² − f = ιf² − f`. ∎

(Checked: `D²` from the FFT periods against the exact rational `(pR − ιf⁴)/f` at three shifts — the maximiser of `|T_H|`, the maximiser of `ρ`, and `c = 1` — for every one of the 5,536 subgroups, `p ≤ 3000`; the average identity `Σ_{cosets c} D² = (p−f)(pf−ιf²)/f` at every subgroup; `|R| ≤ E(H,cH) ≤ E(H)` at all 16,608 witness shifts.)

**Corollary 3.2 (sandwich).** For all `H`, `c ≠ 0`:

    ρ(H,c) ≥ |T_H(c)|·√(f/E(H,cH)) ≥ |T_H(c)|·√(f/E(H)),

and if `H ⊆ QR` and `f² < p` (then `R_H(c) ≥ f²` because `r^χ = r ≥ 0`):

    ρ(H,c) ≤ (|T_H(c)|/√f)·(1 − f²/p)^{−1/2}.

Consequently, for `H ⊆ QR`, `f ≤ p^{1/2−κ}`:  `(M_H/√f)·(f²/E(H))^{1/2} ≤ max_c ρ(H,c) ≤ (M_H/√f)(1 + O(p^{−2κ}))`, and by `E(H) ≪ f^{5/2}` for `f ≤ p^{2/3}` (CITED: Shkredov, arXiv:1311.5726, Corollary 5, `sources/sigma-subgroup-2026-09-05/shkredov-1311.5726.txt` lines 196–198, attributed there to Heath-Brown–Konyagin; the sharper `E(Γ) ≪ |Γ|^{49/20} log^{1/5}|Γ|` for `|Γ| ≤ √p` is Corollary 12 of Murphy–Rudnev–Shkredov–Shteinikov, arXiv:1712.00410, `mrss-1712.00410.txt` §6) the left factor is `≫ f^{−1/4}`:

    M_H/√f  ≍  max_c ρ(H,c)   up to a factor between `C f^{−1/4}` and `1 + o(1)`.

(Data, `p ≤ 3000`, `f < √p`: `max_H E(H)/f^{5/2} = 1.125`.) Both inequalities of the sandwich are checked at every subgroup (lower bound at the three witness shifts; upper bound at all cosets whenever `H ⊆ QR`, `f² < p`).

**Theorem 3.3 (`ρ` is unbounded; REFUTES `ρ = O(1)`).** Let `m ≥ 2` and let `p ≡ 1 (mod m)` be a prime with `p ≥ 4^{m+1} m⁴`; let `H` be the subgroup of order `m`. Then there is `c ∈ F_p^*` with `c + H ⊆ QR` and `R_H(c) = m²`, and for it

    T_H(c) = m,   D(H,c)² = pm − ιm³,   ρ(H,c) = m√p/√(pm − ιm³) ≥ √m.

Since for each `m` there are infinitely many such `p` (Dirichlet), `sup_{p,H,c} ρ(H,c) = ∞`.

*Proof.* Step 1 (Weil count). For `c ∈ F_p` put `Π(c) := Π_{h∈H}(1 + χ(c+h)) = Σ_{S⊆H} χ(f_S(c))`, `f_S(x) := Π_{h∈S}(x+h)`. `Σ_c χ(f_S(c))` equals `p` for `S = ∅`, `0` for `|S| = 1`, and has modulus `≤ (|S|−1)√p` for `|S| ≥ 2` (Weil, background fact of the brief: `f_S` has `|S|` distinct roots, so it is squarefree of degree `≥ 2` and not a constant times a square). Since `Σ_{|S|≥2}(|S|−1) = m2^{m−1} − 2^m + 1 = 2^{m−1}(m−2) + 1`, we get `Σ_c Π(c) ≥ p − (2^{m−1}(m−2)+1)√p`. On the other hand `Π(c) = 2^m` if `c + H ⊆ QR`, `Π(c) = 0` if some `χ(c+h) = −1`, and `0 ≤ Π(c) ≤ 2^{m−1}` for the `m` values `c ∈ −H`. Hence the number `N₊` of `c` with `c + H ⊆ QR` satisfies `2^m N₊ ≥ p − (2^{m−1}(m−2)+1)√p − m2^{m−1}`, i.e. `N₊ ≥ 2^{−m}p − ((m−2)/2 + 2^{−m})√p − m/2 ≥ 2^{−m}p − (m/2)√p − m/2`.
Step 2 (bad shifts). If `t ≠ 0` and `r(−ct)r^χ(t) ≠ 0` then `−ct ∈ (H−H)∖{0}` and `t ∈ (H−H)∖{0}`, so `c` lies in the set `−((H−H)∖{0})/((H−H)∖{0})` of size `≤ m²(m−1)²`. For every other `c ≠ 0`, `R_H(c) = f² = m²` (Theorem 3.1).
Step 3. If `√p ≥ 2^{m+1}m²` then `2^{−m}p ≥ 2m²√p`, so `N₊ ≥ (3/2)m²√p − m/2 ≥ 3·2^m m⁴ − m/2 > m²(m−1)²`, and some `c` with `c + H ⊆ QR` is not bad. For it `T_H(c) = m`, `D² = (pm² − ιm⁴)/m`, `ρ² = pm²/D² = pm/(p − ιm²) ≥ m`. ∎

The hypothesis `p ≥ 4^{m+1}m⁴` is far from necessary: in the scan, the first prime at which the recorded extremal shift of the order-`m` subgroup is monochromatic with `R_H = m²` is `p = 7, 31, 37, 131, 97, 449, 433, 883, 491, 2707, 1669, 2081` for `m = 2, …, 13` (witnesses `(p, c, ρ²)` in the JSON key `theorem33_witnesses`; e.g. `m = 10`: `p = 491`, `c = 90`, `ρ² = 10` exactly; `m = 13`: `p = 2081`, `c = 1070`, `ρ² = 27053/1912 = 14.15`, `ρ = 3.7615 ≥ √13 = 3.606`).

### 3.4 Numerical results, all 5,536 subgroups of all 429 primes `p ≤ 3000`, all cosets of `c`

(`|T_H|`, `D`, `ρ` are constant on cosets `cH` (Lemma 1.2 and Prop. 2.1(a) of the `structured` note), so `c` runs over the `d` coset representatives `g^a`; `T_H` is computed exactly, `D²` in double precision from FFT periods and verified against the exact rational of Theorem 3.1 at the witness shifts; `ρ²` at every witness is recorded as an exact rational.)

| band of `|H|` | n | max `ρ` | mean `max_c ρ` | max `ρ_max/√(2 log d)` | mean of it | max `M_H/√|H|` | mean `median_c ρ` | mean `RMS_c ρ` |
|---|---|---|---|---|---|---|---|---|
| 2–3 | 635 | 2.056 | 1.523 | 0.958 | 0.450 | 1.732 | 0.629 | 0.999 |
| 4–7 | 589 | 3.089 | 2.227 | 1.312 | 0.700 | 2.646 | 0.782 | 1.000 |
| 8–15 | 569 | 3.762 | 2.445 | 1.231 | 0.818 | 3.606 | 0.690 | 0.998 |
| 16–31 | 544 | 3.801 | 2.326 | 1.395 | 0.841 | 3.651 | 0.689 | 0.994 |
| 32–63 | 523 | 3.430 | 2.088 | 1.268 | 0.829 | 3.889 | 0.704 | 0.989 |
| 64–127 | 513 | 3.108 | 1.830 | 1.281 | 0.801 | 3.442 | 0.695 | 0.963 |
| 128–255 | 470 | 2.747 | 1.559 | 1.263 | 0.767 | 3.266 | 0.701 | 0.916 |
| 256–511 | 413 | 2.501 | 1.249 | 1.226 | 0.696 | 2.481 | 0.676 | 0.832 |
| 512–1023 | 301 | 1.795 | 0.776 | 1.078 | 0.502 | 1.778 | 0.534 | 0.572 |
| 1024–1500 | 121 | 0.062 | 0.057 | 0.053 | 0.048 | 0.031 | 0.028 | 0.040 |

(Since `√p < 55` throughout the scan, the rows with `|H| ≥ 64` lie entirely outside the regime `|H| < √p` of the conjecture, and the row 32–63 only partly inside it.)

*Largest observed `ρ` (REFUTATION witnesses for any constant `≤ 3.8`):*

| `p` | `|H|` | `d` | `c` | `T_H(c)` | `ε_H(c)` | `ρ` | `ρ²` exact | `M_H/√|H|` | `√(2 log d)` |
|---|---|---|---|---|---|---|---|---|---|
| 2437 | 29 | 84 | 1226 | 19 | 1.207 | 3.8006 | 879757/60906 | 3.528 | 2.977 |
| 2081 | 13 | 160 | 1070 | −13 | 1.000 | 3.7615 | 27053/1912 | 3.606 | 3.186 |
| 937 | 26 | 36 | 544 | −18 | 1.615 | 3.7336 | 151794/10889 | 3.530 | 2.677 |
| 2273 | 16 | 142 | 2073 | −14 | 1.000 | 3.7155 | 111377/8068 | 3.500 | 3.148 |
| 2551 | 30 | 85 | 1877 | 20 | 1.000 | 3.6515 | 40/3 | 3.651 | 2.981 |
| 1777 | 12 | 148 | 302 | −12 | 1.000 | 3.6136 | 21324/1633 | 3.464 | 3.161 |

*Answers to the two tests.* (i) `max_c ρ(H,c)` is **not** bounded independently of `p` and `|H|`: Theorem 3.3 proves `sup ρ = ∞`, and the data show `max ρ` rising from `2.06` (`|H| ≤ 3`) to `3.80` (`|H| ∈ [16,31]`). Within the range of the scan the growth is indistinguishable from `≈ 0.84·√(2 log d)`: over the 591 subgroups with `16 ≤ |H| < √p`, `M_H/√(2|H| log d)` has minimum `0.557`, median `0.839`, maximum `1.346` (at `p = 2081`, `|H| = 32`, `M_H = 22`), and `ρ_max/√(2 log d)` has mean between `0.80` and `0.84` and maximum at most `1.395` in each of the bands 8–15, 16–31, 32–63, 64–127. (ii) Yes: when `H` is small enough to have a monochromatic shift (`M_H = |H|`; 1,271 such subgroups in the scan, largest `|H| = 13`), `ρ(H,c) ≥ √|H|` at that shift whenever `R_H(c) = |H|²` (which is the typical case: 1,103 of the 1,271 monochromatic subgroups, and all 59 of order 8 except those at `p ∈ {457, 2689}`, where `R_H(c*) = 96`), exactly as Theorem 3.3 predicts; the table of all monochromatic subgroups with `|H| ≥ 8` is in the JSON (`witnesses.monochromatic`).

*Typical `ρ`.* Over subgroups with `|H| ≥ 8`, the median over `H` of `median_c ρ(H,c)` is `0.7071` and of `RMS_c ρ(H,c)` is `0.998`. The RMS is forced: `Σ_{cosets c} T_H(c)² = (pf − f² − ιf²)/f ≈ p` while `D² ≈ pf` on average (Theorem 3.1), so `RMS_c ρ ≈ 1` whenever `D²` is nearly constant in `c`. HEURISTIC reading: `T_H(c)/√|H|` behaves like a standard Gaussian across the `d` cosets, and `M_H ≈ √(2|H| log d)`; this is far stronger than the conjecture needs and is not proved.

*The energy factor at the extremal shift.* Over the 1,917 subgroups with `2 < |H| < √p`, the ratio `ρ(H,c*)/(M_H/√|H|) = (ε_H(c*) − ιf²/p)^{−1/2}` at the maximiser `c*` of `|T_H|` lies in `[0.459, 1.414]` with mean `0.964`; the minimum is at `p = 2113`, `|H| = 44`, where `c* = 1 ∈ H` and `ε_H(1) = E(H)/|H|² = 5.66` (the autocorrelation `a = 0` of the periods, which carries the full additive energy of `H`); the maximum `√2` occurs for `H ⊄ QR` with `R_H(c*) = f²/2` (signed energy below the trivial value). So in the regime of the conjecture the coset ratio and `M_H/√|H|` agree within a factor `≤ 2.2` in all cases examined.

### 3.5 Verdict on the formulation

In the coset basis, `ρ(H,c)` is `|T_H(c)|/√|H|` divided by the square root of the normalised signed energy `ε_H(c) − ιf²/p` (Theorem 3.1); the latter is `1 + (ιf−1)(f−1)/(p−1)` on average over `c`, is `≥ 1 − f²/p` for `H ⊆ QR`, and is `≤ E(H)/f² ≪ f^{1/2}` always. In the dual basis of `Z/d` (Prop. 4.4 below) the corresponding ratio is `|T_H(c)|/√|H|·(1 + O(1/d))` exactly, with a `c`-independent denominator. Hence "square-root cancellation across cosets" is not a new formulation: `ρ = O(1)` is false (Theorem 3.3), `ρ ≤ p^{o(1)}` is equivalent — up to the energy slack `f^{1/4}` — to `M_H ≤ √|H|·p^{o(1)}`, which is the Gaussian-model prediction and strictly stronger than Problem 5 (`M_H ≤ p^{−δ}|H|`). The only genuine content added by the coset viewpoint is the conditional implication of Cor. 3.2: `ρ ≤ R ⟹ M_H ≤ R√(E(H)/|H|) ≪ R|H|^{3/4}`, so a proof of `ρ ≤ p^{o(1)}` for `|H| ≥ p^ε` would give Problem 5 with `δ = ε/4 − o(1)`; but no route to such a bound is visible, since (§2) the norms of the two period families cannot supply it and (§4) the cyclotomic-number expression is the definition of `T_H` in disguise.

## 4. Gauss-period algebra, cyclotomic numbers, and the Jacobi dual (Task 4)

Throughout `H = ⟨g^d⟩`, `η_j = Ĥ(g^j)`, `ψ` the character of `F_p^*` of order `d` with `ψ(g) = ω := e^{2πi/d}`; the characters trivial on `H` are exactly `ψ^i`, `i ∈ Z/d` (the annihilator `H^⊥`). Cyclotomic numbers of order `d`: `(a,k)_d := #{h ∈ H : 1 + g^a h ∈ g^k H}` (`a, k ∈ Z/d`); `δ_a := [−1 ∈ g^aH]`. Note `Σ_k (a,k)_d = f − δ_a`.

**Proposition 4.1 (product relation).** `η_i η_{i+a} = f·δ_a + Σ_{k∈Z/d} (a,k)_d η_{i+k}`.

*Proof.* `η_iη_{i+a} = Σ_{h,h'} ζ^{g^i h + g^{i+a}h'} = Σ_{h''∈H} Σ_{h∈H} ζ^{g^i h (1 + g^a h'')}` (`h'' = h'/h`). For the `δ_a` values of `h''` with `1 + g^ah'' = 0` the inner sum is `f`; otherwise `1 + g^ah'' ∈ g^kH` for a unique `k` and the inner sum is `η_{i+k}`. ∎ (Exact in `Z[ζ]` for all `d`, all `(i,a) ∈ (Z/d)²`, all `p ≤ 100`: 87,948 checks.)

**Proposition 4.2 (the coset sum in cyclotomic numbers is the definition of `T_H`).** Let `e := d` if `d` is even and `e := 2d` if `d` is odd, and `K := ⟨g^e⟩ = H ∩ QR`. Then for `c = g^a`:

    (d even)   T_H(g^a) = Σ_{k∈Z/d} (−1)^k (a,k)_d ;
    (d odd)    T_H(g^a) = χ(g^a) Σ_{k∈Z/e} (−1)^k [ (−a,k)_e + (d−a,k)_e ].

*Proof.* `d` even: insert Prop. 4.1 into Theorem 1.3 (Remark 1.3(a)): `G·T_H(g^a) = Σ_j (−1)^j η_jη_{j+a} = fδ_aΣ_j(−1)^j + Σ_k (a,k)_d Σ_j (−1)^j η_{j+k} = 0 + Σ_k (a,k)_d (−1)^k G`, using `Σ_j(−1)^jη_{j+k} = (−1)^kG`. But this is a tautology: `Σ_k (−1)^k (a,k)_d = Σ_{h∈H} χ(1 + g^ah) = χ(g^a) T_H(g^{−a}) = T_H(g^a)` by the inversion symmetry of Remark 1.3(c) (equivalently the cyclotomic symmetry `(a,k) = (−a, k−a)`). `d` odd: `T_H(c) = χ(c)Σ_{x∈c^{−1}H} χ(1+x)` and `c^{−1}H = g^{−a}K ∪ g^{d−a}K`, on each of which `χ(g^kK) = (−1)^k`. ∎ (Exact for all `a`, all subgroups, `p ≤ 199`: 10,884 checks.)

**Proposition 4.3 (Jacobi expansion — the multiplicative-Fourier dual).** With `J(α,β) := Σ_{y≠0,1} α(y)β(1−y)` (so `J(1,χ) = −1`, `J(χ,χ) = −χ(−1)`, `|J(ψ^i,χ)| = √p` otherwise):

    T_H(c) = (χ(c)/d) Σ_{i∈Z/d} ψ^i(−c) J(ψ^i, χ),      Σ_{i∈Z/d}|J(ψ^i,χ)|² = (d−2)p + 2  (d even),  (d−1)p + 1  (d odd).

*Proof.* `1_H(x) = (1/d)Σ_i ψ^i(x)` for `x ≠ 0`; `Σ_{x} ψ^i(x)χ(x+c) = ψ^i(−c)χ(c)J(ψ^i,χ)` by `x = −cy` (for `i = 0` this reads `Σ_{x≠0}χ(x+c) = −χ(c)`). The moduli are standard. ∎ Consequently the **Jacobi-basis cancellation ratio** `ρ_J(H,c) := |Σ_i ψ^i(−c)J(ψ^i,χ)| / (Σ_i|J(ψ^i,χ)|²)^{1/2} = d|T_H(c)|/√((d−2)p+2) = (|T_H(c)|/√f)·(1 + O(1/d))` has a `c`-independent denominator and is exactly the normalised shifted sum. (Float checks, `p ≤ 199`: 2,104 expansions, 8,965 moduli, 353 `ℓ²` sums, all to `10^{−8}`.)

**Proposition 4.4 (the two expansions are Fourier duals on `Z/d`).** With Gauss sums `G(θ) := Σ_{x≠0}θ(x)e(x)` (`G(1) = −1`): `η_j = (1/d)Σ_i ω^{−ij}G(ψ^i)` and `Ĥ_χ(cg^j) = χ(g^j)(1/d)Σ_i ω^{−ij}ψ^{−i}(c)χ(c)G(ψ^iχ)`, hence

    Σ_{cosets} χ(ξ)Ĥ(ξ)Ĥ_χ(cξ) = (χ(c)/d) Σ_{i∈Z/d} ψ^i(c) G(ψ^i) G(ψ^{−i}χ) = (χ(c)G/d) Σ_i ψ^i(c) J(ψ^i, ψ^{−i}χ),

and `J(ψ^i,ψ^{−i}χ) = ψ^i(−1)J(ψ^i,χ)` recovers Prop. 4.3. *Proof.* `1_{g^jH}(x) = (1/d)Σ_iψ^i(x)ω^{−ij}`; `Σ_xψ^iχ(x)e(cx) = ψ^{−i}(c)χ(c)G(ψ^iχ)`; the `j`-sum of `ω^{−(i+i')j}` is `d[i' = −i]`; `G(α)G(β) = J(α,β)G(αβ)` when `αβ ≠ 1` (here `αβ = χ`), valid with the stated conventions also when `α = 1` or `β = 1`. The reflection formula is `J(α,ᾱβ̄) = α(−1)J(α,β)` (substitute `z = y/(1−y)`). ∎ (Float checks at 2,104 triples for the Gauss-product form and 8,965 reflection identities.) Thus the `ℓ²`-weight `D(H,c)²` of the coset basis and the flat weight `(d−2)p+2` of the Jacobi basis are the two Parseval normalisations of one and the same sum; the coset basis merely redistributes the weight according to `|η_j|²|η̃_{j+a}|²`, which is what produces the energy factor of Theorem 3.1.

### 4.5 The cyclotomic-number statement that would settle Problem 5, and why it is not a reduction

**Statement CN(ε,δ).** For every even `e | p−1` with `e ≤ p^{1−ε}` and every `a ∈ Z/e`:

    | Σ_{k∈Z/e} (−1)^k (a,k)_e |  ≤  p^{−δ} (p−1)/e.

By Prop. 4.2 and Remark 1.3(b), `CN(ε,δ)` is *equivalent* to `M_K ≤ p^{−δ}|K|` for all subgroups `K ⊆ QR` with `|K| ≥ p^ε`, and implies `M_H ≤ 2p^{−δ}|H|` for all `|H| ≥ 2p^ε` — that is, to Bourgain's Problem 5 in the `p^ε` range. The equivalence is a tautology: `(a,k)_e` counts the `h ∈ K` with `1 + g^ah` in the `k`-th coset, `(−1)^k = χ` on that coset, and the alternating sum is `Σ_{h∈K}χ(1+g^ah)`. Nothing about the individual cyclotomic numbers helps: the classical evaluation `(a,k)_e = (1/e²)Σ_{i,j∈Z/e} ω^{−ai−kj} ψ^i(−1) J(ψ^i,ψ^j)` (obtained exactly as in Prop. 4.3; all but three of the `e²` Jacobi sums have modulus `√p`) gives `|(a,k)_e − (p−2)/e²| < √p`, which is vacuous as soon as `e > p^{1/4}` (for `e > √p` most `(a,k)_e` are `0` or `1`); inserting it into `CN` gives `e√p`, worse than trivial. Inserting the Jacobi evaluation into the alternating sum instead picks out `j = e/2` (`ψ^{e/2} = χ`) and returns Prop. 4.3 verbatim:

    Σ_k (−1)^k (a,k)_e = (1/e) Σ_{i∈Z/e} ω^{−ai} ψ^i(−1) J(ψ^i,χ).

So the precise **open** input is:

**Statement J(ε,δ).** For every subgroup `Ψ ≤ \widehat{F_p^*}` of order `e ≤ p^{1−ε}` containing `χ`, and every `u ∈ F_p^*`: `|Σ_{ψ∈Ψ} ψ(u) J(ψ,χ)| ≤ p^{1−δ}`.

The trivial bound is `e√p`; `J(ε,δ)` asks for a cancellation factor `p^{1/2−δ}/e = |K|p^{−δ}/√p`, which for `δ < ε/2` is *less* than square-root cancellation (`e^{−1/2}`) in the `e`-term sum, and equals it at `δ = ε/2`. Over the *full* character group the sum is trivial to evaluate (`Σ_ψ ψ(u)J(ψ,χ) = (p−1)χ(1−1/u)`, `u ≠ 1`), and undoing the Mellin transform of a subgroup of characters gives back `Σ_{h∈K}` — the dual sum has `e = p^{1−ε}` terms of known modulus and unknown phase, i.e. more terms than the original. Katz's Sato–Tate theorems for finite-field Mellin transforms (CITED: `sources/katz-finite-field-mellin.txt`, introduction lines 205–232 and Theorem 1.1 at line 483) concern the distribution of such sums as the character ranges over *all* multiplicative characters of `k^×` with `#k → ∞`; they contain no statement about a proper subgroup `Ψ` of the characters of a fixed `F_p^*`, and no source examined in this workspace does. I am not aware of any result on `J(ε,δ)` for `e > p^{1/2}`; this is the honest endpoint of the "multiplicative dual" direction. Numerical evidence (`p ≤ 3000`) is the same data as §3.4: `|Σ_{ψ∈Ψ}ψ(u)J(ψ,χ)|/√(ep) = ρ_J(K,u)(1+O(1/e)) = |T_K(u)|/√|K|·(1+O(1/e))`, whose maximum over `u` is `≤ 1.35·√(2 log e)` in every case with `16 ≤ |K| < √p`, i.e. square-root cancellation up to the Gaussian extreme-value factor, and never a violation of it.

### 4.6 What a bound on `ρ` would have to be, and why none is available

By Cor. 3.2, for `H ⊆ QR`, `f ≤ p^{1/2−κ}`: any bound `max_c ρ(H,c) ≤ f^{1/2−θ}` with `θ > 1/4` yields `M_H ≪ f^{1−θ+1/4}`, a power saving, hence Problem 5 — so a *nontrivial* bound on `ρ` is at least as hard as the open problem; conversely the only unconditional bounds are the trivial ones `ρ ≤ √f·(1−f²/p)^{−1/2}` (from `M_H ≤ f`) and `ρ ≤ √(p/f)·(1+o(1))` (from Weil, `M_H < √p`), the second being better only for `f > √p`. The `ℓ²` weight of Theorem 3.1 is provably `D² ≤ (p/f)E(H) ≪ pf^{3/2}` (CITED energy bound) and `D² ≥ pf − f³` (`H ⊆ QR`), so the "denominator" side of Task 4 is closed: `Σ_{cosets}|ĤĤ_χ|² = p·f·(ε_H(c) − ιf²/p)` with `1 − f²/p ≤ ε_H(c) − ιf²/p ≤ E(H)/f² ≪ f^{1/2}`, and any further gain must come from the numerator, i.e. from `M_H` itself.

## 5. Verifier record

`experiments/sigma_dual_2026_09_05.py` (numpy, sympy for the cross-check; 16 s; 310,324 checks, 0 failures):

- B: sympy cross-validation of the exact `Z[ζ_p]` arithmetic (12 random products, `p ∈ {7,11,13}`).
- §1 (`p ≤ 199`, 45 primes, 353 subgroups): 37,853 exact coset identities (all `c ∈ F_p`; 8,965 with `c ∈ −H`, 353 with `c = 0`); 706 exact `ℓ²` identities (periods and twisted periods, over cosets); 87,948 exact period product relations (`p ≤ 100`); 10,884 exact cyclotomic alternating identities (both parities of `d`); 23,298 exact inversion symmetries.
- §2 (`p ≤ 199`): 2,104 Jacobi expansions, 4,208 Gauss-product dual / float coset identities, 8,965 Jacobi moduli, 353 Jacobi `ℓ²` sums, 8,965 reflection identities (all to `10^{−8}`), Gauss-sum modulus at each prime.
- §3 (`p ≤ 3000`, 429 primes, 5,536 subgroups, all cosets of `c`): `T_H` by FFT checked against direct sums; coset invariance of `|Ĥ|`, `|Ĥ_χ|`; equality of the two multisets when `H ⊆ QR`; float `ℓ²` identities; `D²` by cyclic correlation checked against direct sums and against the exact rational `(pR − ιf⁴)/f` at three shifts per subgroup; `D²`-average identity; second moment of `T_H`; `|R| ≤ E(H,cH) ≤ E(H)`; both sandwich inequalities; pigeonhole bounds; `Φ` exactness for `H ⊆ QR`; `Φ ≥ 0.4√p`; `M_H < √p`; Theorem 3.3 witnesses `ρ ≥ √m`.
- Output JSON: all 5,536 subgroup rows (`p, |H|, d, ι, M_H, c*, T_H(c*), ρ(c*), ε(c*), ρ² exact, ρ_max, c_ρ, …, E(H), Φ/√p, A, B`), the witness tables of §3.4, the band table, the monochromatic list, the Theorem 3.3 witnesses, and `Φ_min`.

## 6. Remaining obligations (OPEN)

1. Statement `J(ε,δ)` / `CN(ε,δ)`: any cancellation in `Σ_{ψ∈Ψ}ψ(u)J(ψ,χ)` over a subgroup `Ψ` of characters of order `e ∈ (p^{1/2}, p^{1−ε}]` — equivalently Problem 5. No tool in this note bears on it: the `ℓ²`/`ℓ^∞` data of the period families are exhausted (§2), and the coset/Jacobi bases are Parseval-equivalent (§4.4).
2. The Gaussian heuristic `M_H ≈ √(2|H| log d)` (HEURISTIC, §3.4) is consistent with all 5,536 subgroups at `p ≤ 3000`; a proof of even `M_H ≤ |H|^{1/2+o(1)}·p^{o(1)}` is far beyond reach and would imply Problem 5 with `δ = ε/2 − o(1)`.
3. The energy factor `ε_H(c*)` at the extremal shift: Cor. 3.2 bounds it by `E(H)/|H|² ≪ |H|^{1/2}`; the data suggest `ε_H(c*) = O(1)` for `|H| < √p` (max `5.66`, attained with `c* ∈ H`), i.e. the extremal shift is not additively special. Proving `ε_H(c*) = O(1)` would sharpen Cor. 3.2 to `max_c ρ ≍ M_H/√|H|` but would not, by itself, bound either side.
4. Not examined: Bourgain–Chang's multilinear Burgess note and Katz's "Gauss sums, Kloosterman sums, and monodromy groups" for a possible "horizontal" statement about Jacobi sums along a subgroup of characters (the local copy `sources/katz-gauss-kloosterman-monodromy.txt` was not read for this note).

## Sources

- Chang, M.-C., *Character sums in finite fields*, Contemp. Math. 518 (2010), Problem 5 — as quoted in `research/sigma-structured-2026-09-05.md` §1 from `sources/sigma-chang-character-sums-survey-DubProc.txt`.
- Kurlberg, P., *Bounds on exponential sums over small multiplicative subgroups* (2007), Theorem 1.1 — `sources/kurlberg-small-subgroups-2007.txt`, lines 74–77 (the Bourgain–Glibichuk–Konyagin bound).
- Shkredov, I. D., *On exponential sums over multiplicative subgroups of medium size*, arXiv:1311.5726, Corollary 5 (`E(Γ) ≪ |Γ|^{5/2}`, `|Γ| ≤ p^{2/3}`, attributed to Heath-Brown–Konyagin, Q. J. Math. 51 (2000)) — `sources/sigma-subgroup-2026-09-05/shkredov-1311.5726.txt`, lines 196–198, 244–252.
- Murphy, B., Rudnev, M., Shkredov, I. D., Shteinikov, Y. N., *On the few products, many sums problem*, arXiv:1712.00410, §6, Corollary 12 (`E(Γ) ≪ |Γ|^{49/20}log^{1/5}|Γ|`, `|Γ| ≤ √p`) — `sources/sigma-subgroup-2026-09-05/mrss-1712.00410.txt`.
- Katz, N. M., *Convolution and equidistribution: Sato–Tate theorems for finite-field Mellin transforms* — `sources/katz-finite-field-mellin.txt`, introduction (lines 205–232) and Theorem 1.1 (line 483); used only to record that its equidistribution is over all characters of `k^×`, `#k → ∞`.
- Weil's bound for `Σ_x χ(f(x))` and the second-moment identity: background facts of `research/sigma-brief-2026-09-05.md`.
- The equivalence Theorem 2.2 and Prop. 2.1 (coset covariance of `T_H`) of `research/sigma-structured-2026-09-05.md`; the dilation-invariance obstruction of `research/sigma-crux-2026-09-05.md` §4 (Prop. 4.6), of which Prop. 2.3 above is the coset-basis analogue.
