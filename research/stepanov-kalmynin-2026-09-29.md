# Stepanov wave, worker `kalmynin`: the exact-decomposition identities under containment (2026-09-29)

**Status:** PROVED — (i) the Taylor expansion of the Hanson–Petridis polynomial `F_A` to all
orders at every complete, anti-complete and `−A` point (Theorem 2.1). From it, Kalmynin's
Relations X and Y, the Yip–Yoo second-order identity and the top-coefficient identities hold for
**every** complete biclique `A + B ⊆ Q ∪ {0}`, each with an explicit additive defect term built
from the cofactor `Γ_A` of degree `g = d − mn + r` (Theorem 2.2, Corollary 2.3, Proposition 2.4).
(ii) The moment identities `Σ_{a,b}(a+b)^j = 0`, `1 ≤ j < d` (Kalmynin Lemma 3, Yip–Yoo (2.6),
Rudnev–Tyrrell Lemma 4.2) hold for a containment pair iff `A + B` covers the subgroup uniformly,
that is iff `A ⊕ B = Q` (Proposition 2.5). This is the **only** point at which the three proofs of
Sárközy's conjecture use `A + B = Q` rather than containment together with Hanson–Petridis
tightness (§1.4). (iii) A reach statement (Proposition 2.7): an identity that involves `Γ_A` only
through its Taylor data of order `≤ s` on `B` and its top `J` coefficients says nothing beyond
containment once `g ≥ (s+1)n + J`. Near HP-tightness in the balanced window (`g ≥ 2c·d`,
`n ≍ √p`), every such identity of order `s ≲ 2c·d/n ≍ c√p` is therefore void. Kalmynin, Yip–Yoo and
Rudnev–Tyrrell use `s ≤ 2`, `J ≤ 3`. (iv) Two transfers survive containment exactly: Kalmynin's
Möbius map for cliques (his Lemma 7) and the Rudnev–Tyrrell reciprocal map for bicliques. Each acts
on bad-partner profiles by an explicit linear map (Propositions 3.1–3.2). The "twisted"
Hanson–Petridis bound they produce is the Hanson–Petridis bound again (Proposition 3.3).
Exact rational LP certificates (Theorem 3.4) show that imposing every count constraint of the
balanced note on all transferred sets still does not exclude an HP-tight balanced profile. This
holds for cliques up to `m(m−1)/p = 0.4993` and for bicliques (joint profiles, with the Chung
identity) up to `mn/p = 0.4998`, for primes `p ≤ 1000033`. (v) At complete points, every Hankel
minor of the robust note carries only the Relation-X datum (Proposition 4.1).
REFUTED (explicit witnesses) — the containment-plus-tightness versions of Kalmynin's `α = β`
theorem and of the moment identities. For every `p ≡ 1 (mod 4)` and every non-residue `t`, the set
`A = {0,t}` with its complete set `B` (`|B| = (p−1)/4`, `r = 0`) is HP-critical, but `|A| ≠ |B|`
and `Σ(a+b)² ≠ 0` (e.g. `p = 13`: Proposition 2.6). Also refuted: the defect-free Relations X, Y and
Yip–Yoo (2.9) for `g > 0` (e.g. `p = 13`, `A = {0,1}`, `B = {3,9}`). OPEN — the target
`|A||B| ≤ (1/2 − c)p + o(p)` for balanced complete bicliques, and the prime Paley clique constant.
The precise obstruction is Theorem 5.1. Every ingredient these papers add beyond Hanson–Petridis is
one of: (a) the moment identities, which need `A ⊕ B = Q`; (b) a `g = 0` rigidity statement, whose
containment version carries an unconstrained defect; (c) a transfer, which no known count
constraint sees. The Hankel-minor robustness trick transports (b) and adds nothing.

Verifier: `experiments/stepanov_kalmynin_2026_09_29.py` → `results/stepanov_kalmynin_2026_09_29.json`
(run with `/opt/miniconda3/bin/python3`; about 60 s, one process; 1,228,382 exact checks, 0 failures).
Sources: Kalmynin arXiv:2504.10202v2 (`sources/sigma-kalmynin-2504.10202v2.txt`), Yip–Yoo
arXiv:2608.02568v2 of 9 Sep 2026 (fetched 2026-09-29 and saved as
`sources/stepanov-kalmynin-yipyoo-2608.02568.{pdf,txt}`), Rudnev–Tyrrell arXiv:2607.24270
(`sources/stepanov-robust-rudnev-tyrrell-2607.24270.txt`), Hanson–Petridis arXiv:1905.09134
(`sources/sigma-hanson-petridis-1905.09134.txt`), and in-workspace
`research/stepanov-balanced-2026-09-29.md` (cited as **[Bal]**) and
`research/stepanov-robust-2026-09-26.md` (**[Rob]**). Statements from the papers are restated in
our notation with their numbers; the wording is ours.

---

## 0. Notation

`p` odd prime, `χ` Legendre, `χ(0) = 0`, `d = (p−1)/2`, `Q` the nonzero squares, `ν` a
non-residue. `A = {a_1,…,a_m}`, `2 ≤ m ≤ (p+1)/2`, `D = d + m − 1 ≤ p − 1`,
`c_k = Π_{l≠k}(a_k − a_l)^{−1}`, and

    F_A(x) = −1 + Σ_k c_k (x + a_k)^D ,   deg F_A = d,  leading coefficient λ = C(D, m−1)

(Kalmynin's `HP(x; A, d)` for `μ_d = Q`). `h(X) = Π_{a∈A}(X + a)` (Yip–Yoo's `G`).
`C₁(A) = {x ∉ −A : x + A ⊆ Q}` (complete points), `C₀(A) = {x ∈ −A : x + A ⊆ Q ∪ {0}}`,
`AC(A) = {x ∉ −A : x + A ⊆ νQ}` (anti-complete). For `B ⊆ C₁(A) ∪ C₀(A)`: `n = |B|`,
`B₁ = B ∖ (−A)`, `B₀ = B ∩ (−A)`, `r = |B₀|`, `P₁, P₀, P = P₁P₀` the monic root polynomials and

    F_A = λ · P₁^m · P₀^{m−1} · Γ_A,   Γ_A monic,  deg Γ_A = g := d − mn + r,  Γ_A(b) ≠ 0 on B

([Bal] Theorem 1.1(3); `Γ_A` is [Bal]'s `G_A`). The pair is **critical** (Kalmynin's Definition 2,
"d-critical") iff `g = 0`. For `x ∈ F_p` let `σ_s(x) = Σ_k c_k (x + a_k)^{−1−s}` (`x ∉ −A`),
`S₁(b) = Σ_a 1/(a+b)`, `S₂(b) = Σ_a 1/(a+b)²`, `T₁(b) = Σ_{b'∈B∖b} 1/(b−b')`,
`T₂(b) = Σ_{b'∈B∖b} 1/(b−b')²`. `M_j = Σ_{a∈A,b∈B}(a+b)^j`. All binomials `C(D, j)`, `j ≤ D`, and
all `(D)_j` are units mod `p`.

## 1. What the exact-decomposition papers use (task 1; CITED)

### 1.1 Kalmynin, arXiv:2504.10202v2 (28 May 2025)

Setting: `d | p − 1`, `1 < d < p − 1`, `μ_d` the `d`-th roots of unity, `α = |A|`, `β = |B|`.

- **Lemmas 1–2** (pure algebra). `Σ_a c_a a^j = [j = α−1]` for `0 ≤ j ≤ α−1`, `c_a = 1/Π(a − a')`,
  shift invariance, and `Σ_a c_a a^{j+α−1} = h_j(A)` (complete homogeneous symmetric polynomial).
- **Definition 2.** `(A, B)` is `d`-critical if `A + B ⊆ μ_d ∪ {0}` and `|A||B| = d + |(−A) ∩ B|`.
  If `A + B` equals `μ_d` or `μ_d ∪ {0}`, the pair is critical: this follows from HP plus
  counting, and then every element has a unique representation.
- **Lemma 3** (moment identity). If `A + B = μ_d` or `μ_d ∪ {0}`, then `Σ_{a,b}(a+b)^k = 0` for
  `1 ≤ k < d`. *Uses:* unique representation of every element of `μ_d` (equality), and
  `Σ_{x∈μ_d} x^k = 0`.
- **Lemma 4** (critical factorisation). For a critical pair,
  `HP(x;A,d) = C Π_{b∈B}(x−b)^{α−ε(b)}` with `ε(b) = [b ∈ −A]`. *Uses:* containment (vanishing
  orders), `deg HP = d`, criticality (degree count). This is [Bal] Theorem 1.1(3) at `g = 0`.
- **Lemmas 5–6** (cliques, `B = −A`). If `(A,−A)` is critical and `d ∉ {2,6}`, then
  `p₂ = p₁²/α` and `p₃ = p₁³/α²`. *Uses:* comparison of the coefficients of `x^d, x^{d−2}, x^{d−3}`
  in Lemma 4, i.e. `g = 0`.
- **Lemma 7** (Möbius transfer). If `(A, −A)` is critical, then for every `a ∈ A` so is
  `(A^a, −A^a)`, where `A^a = {0} ∪ {1/(a−a') : a' ∈ A ∖ {a}}`. *Uses:* only
  `A − A ⊆ μ_d ∪ {0}` and `|A^a| = |A|`. **Containment-only.**
- **Corollary 2 and Theorem 5** (pure algebra). Lemmas 6–7 give, at every `a ∈ A`, relation (1)
  `Σ_{a'≠a}(a−a')^{−2} = α^{−1}(Σ_{a'≠a}(a−a')^{−1})²` and relation (2), the cubic analogue with
  `α^{−2}`. Theorem 5: over an algebraically closed field (with `α < √(char K)` in positive
  characteristic), sets satisfying (1) and (2) have `α = 5` and `A = {C, C ± D, C ± iD}`. The
  proof uses the differential operator
  `𝒟g = 4α(α−2)g′g‴ − 3(α−1)(α−2)g″² − α(α+1)g g⁗`, which kills `F_a(x) = (x−a)(1+s(x−a))^α`, and
  the factorisation `D(α,l) = −α(α−l+1)(α−l)(α−l−1)(α² − (l+5)α − l + 6)` (both re-verified
  symbolically, verifier §A10).
- **Theorem 4** (Lev–Sonn, containment form). If `d ∉ {2, 6}`, no `(A,−A)` is `d`-critical,
  except `p = 41`, `d = 20`. *Uses:* Lemmas 4–7 and Theorem 5 only, so **containment plus
  criticality, not equality.** For `d = (p−1)/2` it says: a Paley clique with `m(m−1) = d` exists
  only for `p ∈ {5, 13, 41}`. So it improves the Hanson–Petridis clique bound by at most one unit.
- **Lemma 8** (support of power sums) uses Lemma 3.
- **Lemma 9** (substitution `x ↦ 1/x + b` in Lemma 4) and **Lemma 10**. For `A + B = μ_d`, every
  `b ∈ B` satisfies

      (Relation X)  S₁(b) = α(α+1)/(d−1) · T₁(b),
      (Relation Y)  S₁(b)² + S₂(b) = α(α+1)(α+2)/((d−1)(d−2)) · (α T₁(b)² − T₂(b)).

  *Uses:* the critical factorisation with `r = 0` (`ε ≡ 0`) and `αβ = d`. The proof compares the
  coefficients of `x^{d−1}, x^{d−2}, x^{d−3}` of the transformed identity. These are the Taylor
  coefficients of order `α, α+1, α+2` of `HP` at `b` (§2).
- **Lemmas 11–12, 14–16** (residue theorem on `P¹`) are pure algebra, given the power-sum
  support that Lemma 8 extracts from Lemma 3 and given Relation X.
- **Lemma 13, Lemma 17, proof of Theorem 3.** These combine everything into congruences such as
  `2(3k−2)(k−1)α + (k+2)(k+3) ≡ 0 (mod p)` and then lift them to `Z`: the two sides are integers of
  size `< p`, so the congruence is an equality. *Uses:* `α = β`, `p = 2α² + 1`, and Lemmas 3/8.

### 1.2 Yip–Yoo, arXiv:2608.02568v2 (9 Sep 2026)

Setting: `H` a proper subgroup of `F_p^*` with `H = A + B`, `|A|, |B| ≥ 2`. `F(X) = Π_{b∈B}(X−b)`,
`G(X) = Π_{a∈A}(X+a)`.

- **Proposition 2.7.** `d = αβ`; `Σ_a c_a(X+a)^{d+α−1} − 1 = C(d+α−1, α−1)F(X)^α` (2.5); and
  `Σ_{i≤j} C(j,i)p_i(A)p_{j−i}(B) = 0` for `1 ≤ j < d` (2.6). *Uses:* for (2.5), containment and
  `αβ ≤ d = |A+B| ≤ αβ`; for (2.6), unique representation of `H` and orthogonality (their
  Lemma 2.2). **(2.6) is Kalmynin's Lemma 3.**
- **Proposition 2.8.** At every root `b` of `F`: `2(d−1)G′(b)/G(b) = α(α+1)F″(b)/F′(b)` (2.7). The
  analogue (2.8) holds at the roots of `G`. Globally,
  `α(α+1)GF″ + β(β+1)FG″ = 2(d−1)F′G′` (2.9). *Uses:* (2.5) for `A` and for `B`; `F` and `G`
  coprime (their Lemma 2.3, from `0 ∉ H`); and the gluing step: the left-minus-right polynomial has
  degree `≤ α+β−2` and `α+β` roots.
- **Proposition 2.9** (`α = β`). Take the least `s` with `p_s(A) ≠ 0`. Then (2.6) gives
  `p_s(B) = −βp_s(A)/α`, and the coefficient of `X^{α+β−s−2}` in (2.9) gives (2.14). *Uses (2.6),
  i.e. equality.*
- **Proposition 2.10** (third-order identity (2.19)) assumes `α = β` and `|H| = α²`. Lemma 3.1,
  Proposition 3.2 (`k` even, `Φ_α(k) ≡ 0`), Corollary 3.3 and §4 use (2.6) through Lemma 3.1.
- Their Remark 2.11 identifies (2.7) and (2.21) with Kalmynin's Relations X and Y.

### 1.3 Rudnev–Tyrrell, arXiv:2607.24270

- **Lemma 3.4** (modified HP). Let `S ⊆ μ_d`, `|S| = α`, and `T ⊆ F_p^*` with `S + t ⊆ tμ_d` for all
  `t ∈ T`. If `d = α(|T| + 1)`, then
  `Σ_s k_S(s)(x+s)^{d+α−1} − C_S x^{d+α−1} = C(α+d−1, α) x^{α−1} Π_{t∈T}(x−t)^α`,
  with `k_S(s) = c_S(s)/s` and `C_S = Σ k_S`.
- **Relation (1) in the proof of Lemma 3.5.** For distinct `b, b' ∈ B`:
  `1/(a+b) + 1/(b'−b) = (a+b′)/((a+b)(b′−b))`. So `S = 1/(A+b)` and `T = {1/(b′−b)}` satisfy
  `S + t ⊆ tμ_d`. *Uses:* only `A + B ⊆ μ_d`. **Containment-only.**
- **Lemmas 3.5, 3.8–3.9, Corollary 3.10.** These are Kalmynin's Lemma 9 and Relations X, Y
  (reciprocal transfer identities of first and second order). *Use* `d = α(|T|+1)`, i.e.
  criticality, and Kalmynin's `α = β` as input.
- **Lemma 4.2** (moment identity). `p_j(B) = −p_j(A)` for `j ≤ k`, proved from the uniqueness of
  representations in `A + B = μ_d`. **Equality.**

### 1.4 Where equality is used (summary; PROVED where it refers to §2)

| ingredient | Kalmynin | Yip–Yoo | Rudnev–Tyrrell | needs | containment version |
|---|---|---|---|---|---|
| critical factorisation | L4, L9 | (2.5) | L3.4–3.5 | containment + `g = 0` | [Bal] Thm 1.1: cofactor `Γ_A`, `deg = g` |
| top coefficients (`M₁ = 0`) | L5–6 | — | — | `g = 0` | Cor. 2.3 (defect `[x^{g−1}]Γ_A`) |
| Relations X, Y, (2.7), (2.21) | L10 | P2.8, P2.10 | L3.8–3.9 | `g = 0`, `r = 0` | Thm 2.2 (defect `Γ_A′/Γ_A`, …) |
| global identity (2.9) | — | P2.8 | — | `g = 0` both sides, `r = 0` | Prop. 2.4 (`Z = PhW`, `deg W ≤ 2g−2`) |
| moment identities | L3, L8 | (2.6), L3.1 | L4.2 | **`A ⊕ B = H`** | none: Prop. 2.5 (equivalent to `A ⊕ B = H`) |
| `α = β` | Thm 2 | P2.9 | (input) | moment identities | false for critical pairs: Prop. 2.6 |
| Möbius / reciprocal transfer | L7 | — | (1) | containment only | exact: Props. 3.1–3.3 |
| residues, `𝒟`, Newton | L11–16, Thm 5 | L2.1–2.6 | App. | nothing | valid for all sets |
| lifting congruences to `Z` | L13, L17 | (2.14), C3.3 | §7 | defect-free identities | unavailable when `g ≥ 3n+3` (Prop. 2.7) |

So, for the quadratic residues, Kalmynin's clique theorem (Theorem 4) is already a
containment-plus-criticality statement. The Sárközy/`α = β` part needs the moment identities, and
nothing weaker: the critical containment pairs of Proposition 2.6 satisfy every other identity in
the table (verified: Relations X and Y at 10,646 points, (2.9) for 123 pairs) and still have
`α ≠ β`.

## 2. The identities under containment (task 2)

### 2.1 Taylor expansion at complete points

**Theorem 2.1 (PROVED).** Let `x ∈ C₁(A)` (respectively `x ∈ AC(A)`) and put `ε = +1` (resp. `−1`).
Then `[Y⁰]F_A(x+Y) = ε − 1`, `[Y^j]F_A(x+Y) = 0` for `1 ≤ j ≤ m−1`, and for `0 ≤ s ≤ d − m`

    [Y^{m+s}] F_A(x+Y) = ε · C(D, m+s) · σ_s(x),    σ_s(x) = (−1)^{m−1} [Y^s] 1/h(x − Y).

In particular, for `x ∈ C₁(A)`:
`F_A(x+Y) = (−1)^{m−1} Y^m · Σ_s C(D,m+s) [Y^s](1/h(x−Y)) Y^s`. So `F_A(x+Y)/Y^m` is the
Hadamard product of `(−1)^{m−1}/h(x − Y)` with the weights `C(D, m+s)`. At `x = −a_{k₀} ∈ C₀(A)`:
`[Y^j] = 0` for `j ≤ m−2`, `[Y^{m−1}] = −C(D,m−1)c_{k₀}`, and
`[Y^{m−1+s}] = C(D,m−1+s)Σ_{k≠k₀}c_k(x+a_k)^{−s}` for `s ≥ 1`.

*Proof.* `[Y^j]F_A(x+Y) = C(D,j)Σ_k c_k y_k^{D−j} − [j=0]` with `y_k = x + a_k`. For `y ≠ 0` and any
integer `j`, `y^{D−j} = y^d · y^{m−1−j} = χ(y) y^{m−1−j}` (Euler). So on `C₁ ∪ AC` the coefficient is
`εC(D,j)Σ_k c_k y_k^{m−1−j} − [j=0]`. For `0 ≤ m−1−j ≤ m−2` this vanishes by
`Σ_k c_k y_k^i = [i = m−1]` (Kalmynin Lemma 1, translated). At `j = 0` it is `ε − 1`, and at
`j = m+s` the exponent is `−1−s`. Partial fractions,
`Σ_k c_k/(t + a_k) = (−1)^{m−1}/h(t)` ([Bal] proof of Theorem 1.1(1); Yip–Yoo Lemma 2.6), with
`t = x − Y` expanded in `Y`, give `Σ_s σ_s(x)Y^s = (−1)^{m−1}/h(x−Y)`. At `x = −a_{k₀}` the `k₀` term
is `0^{D−j} = 0`, and `Σ_{k≠k₀}c_k y_k^i = −[i=0]c_{k₀}` for `0 ≤ i ≤ m−2`. ∎
(Verifier §A1: 158,213 coefficients at 7,186 complete/anti-complete points and 434 `C₀` points, `p ≤ 61`.)

### 2.2 Relations X and Y with defect

**Theorem 2.2 (PROVED).** Let `A + B ⊆ Q ∪ {0}` and `b ∈ B₁`. Then

    (X_Γ)  (d−1)/(m+1) · S₁(b) = m Σ_{b'∈B₁∖b} 1/(b−b') + (m−1) Σ_{b'∈B₀} 1/(b−b') + Γ_A′(b)/Γ_A(b),

    (Y_Γ)  (d−1)(d−2)/(2(m+1)(m+2)) · (S₁² + S₂) − ½((d−1)/(m+1))² S₁²
             = −(m/2) Σ_{B₁∖b}(b−b')^{−2} − ((m−1)/2) Σ_{B₀}(b−b')^{−2} + ½[Γ_A″/Γ_A − (Γ_A′/Γ_A)²](b).

For `g = 0` and `r = 0`, (X_Γ) and (Y_Γ) are exactly Kalmynin's Relations X and Y (with `α = m`).
Eliminating `S₁²` via X turns Y_Γ into his form.

*Proof.* Let `Φ(Y) = F_A(b+Y)/Y^m`. By Theorem 2.1, `Φ = t₀ + t₁Y + t₂Y² + O(Y³)` with
`t_s = C(D,m+s)σ_s(b)`. Also `Σ_s σ_sY^s = σ₀·h(b)/h(b−Y) = σ₀ Π_a (1 − Y/(a+b))^{−1}`, so
`σ₁/σ₀ = S₁` and `σ₂/σ₀ = (S₁² + S₂)/2`. Since `C(D,m+1)/C(D,m) = (d−1)/(m+1)` and
`C(D,m+2)/C(D,m) = (d−1)(d−2)/((m+1)(m+2))`, the left sides are `ℓ₁(Φ) = t₁/t₀` and
`ℓ₂(Φ) = t₂/t₀ − t₁²/(2t₀²)`, where `ℓ₁, ℓ₂` are the `Y` and `Y²` coefficients of `log`. For `ψ`
with `ψ(0) ≠ 0` these are `ψ′/ψ(0)` and `½(ψ″/ψ − (ψ′/ψ)²)(0)`, and both are additive in products.
By the factorisation,
`Φ = λ·Π_{b'∈B₁∖b}(b−b'+Y)^m · Π_{b'∈B₀}(b−b'+Y)^{m−1} · Γ_A(b+Y)`, and every factor has nonzero
constant term. Finally `ℓ₁(c+Y) = 1/c` and `ℓ₂(c+Y) = −1/(2c²)`. ∎

*Checks (§A2).* Both identities hold at all 35,777 points of 3,149 complete bicliques
(`p ≤ 101`, `|A| ≤ 6`, `B` = all or part of `C₁ ∪ C₀`). They also hold on the exact decompositions
`μ₄ = {0, −1−i} + {1, i}` (13 primes). The defect-free Relation X **fails** at 12,192 of the 12,412
points with `g > 0` and `r = 0`. Witness: `p = 13`, `A = {0,1}`, `B = {3,9}`, `g = 2`.

**Corollary 2.3 (PROVED; top coefficient and summed X).** For every complete biclique,

    (d/m)·p₁(A) = −m·p₁(B₁) − (m−1)·p₁(B₀) + [x^{g−1}]Γ_A,

and if `r = 0`, `(d−1)/(m+1)·Σ_{a,b} 1/(a+b) = Σ_{b∈B} Γ_A′(b)/Γ_A(b)`. Hence critical pairs with
`r = 0` satisfy `M₁ = Σ(a+b) = 0` and `M_{d−1} = Σ(a+b)^{−1} = 0`.

*Proof.* `[x^{d−1}]F_A = C(D,m)Σ_k c_k a_k^m = C(D,m)p₁(A)` (Kalmynin Lemma 2 with `h₁ = p₁`), and
`C(D,m)/λ = d/m`. The right side is the `x^{d−1}` coefficient of `P₁^mP₀^{m−1}Γ_A`. For `g = r = 0`
and `mn = d`: `m(n p₁(A) + m p₁(B)) = 0`. Summing (X_Γ) over `B`, the double sum of `1/(b−b')`
cancels, and `(a+b)^{d−1} = χ(a+b)/(a+b)`. ∎ (§A2: 3,149 and 1,857 cases; 606 critical cases.)

### 2.3 The Yip–Yoo identity with defect

**Proposition 2.4 (PROVED).** Let `A + B ⊆ Q` (so `r = 0`), `P = Π_b(X−b)`, `h = Π_a(X+a)`.
Let `Γ_A`, `Γ_B` be the cofactors of `F_A` and of `F_B = λ_B Π_a(x−a)^n Γ_B` ([Bal] Theorem 1.1(4),
same `g`), and `Γ̃_B(X) = Γ_B(−X)`. Put
`𝒬 = m(m+1)hP″ + n(n+1)Ph″ − 2(d−1)P′h′` and

    Z = 𝒬·Γ_A·Γ̃_B + 2(m+1)·hP′·Γ_A′·Γ̃_B + 2(n+1)·Ph′·Γ_A·Γ̃_B′.

Then `Ph | Z` and `W := Z/(Ph)` has degree `≤ 2g − 2`. For `g = 0`, `W = 0` and `Z = 𝒬 = 0`
is Yip–Yoo's (2.9).

*Proof.* At `b ∈ B`, (X_Γ) with `T₁ = P″/(2P′)` gives
`𝒬(b) = −2(m+1)h(b)P′(b)Γ_A′(b)/Γ_A(b)`, so `Z(b) = 0`. At `z = −a`, (X_Γ) for the pair `(B, A)` reads
`(d−1)/(n+1)·Σ_b 1/(a+b) = nΣ_{a'≠a}1/(a−a') + Γ_B′(a)/Γ_B(a)`. With
`Σ_b 1/(a+b) = −P′(z)/P(z)`, `Σ_{a'≠a}1/(a−a') = −h″(z)/(2h′(z))` and
`Γ_B′(a)/Γ_B(a) = −Γ̃_B′(z)/Γ̃_B(z)`, this is `𝒬(z) = −2(n+1)P(z)h′(z)Γ̃_B′(z)/Γ̃_B(z)`, so
`Z(z) = 0`. The roots of `Ph` are simple and distinct (`r = 0`), and
`deg Z ≤ m + n − 2 + 2g`. ∎
(§A3: 932 pairs. `Ph | Z` always; `deg W ≤ 2g−2` always, with equality attained. `𝒬 = 0` for all
123 critical pairs, and `𝒬 ≠ 0` for all 809 pairs with `g > 0`.)

**What Yip–Yoo's gluing needs.** The gluing works because `deg 𝒬 < #roots` when `g = 0`. With a
defect, `W` is an unknown polynomial of degree `≤ 2g − 2`, and coefficient comparison in `Z = PhW`
at the top (Yip–Yoo §2.4) or at the bottom (§4) meets the coefficients of `W`, `Γ_A` and `Γ̃_B`
first.

### 2.4 The moment identities are exact decomposition

**Proposition 2.5 (PROVED).** Let `H ≤ F_p^*` have order `d`, and `μ : H → Z_{≥0}` with all
`μ(x) < p`. Then `Σ_{x∈H}μ(x)x^j ≡ 0` for all `1 ≤ j ≤ d−1` iff `μ` is constant. Consequently, for
`A + B ⊆ H`, the moment identities `M_j = 0` (`1 ≤ j < d`) hold iff every element of `H` has the same
number `c` of representations. Under Hanson–Petridis (`mn ≤ d` when `0 ∉ A + B`) this means
`c = 1`, i.e. `A ⊕ B = H`.

*Proof.* The vectors `(x^j)_{x∈H}`, `0 ≤ j ≤ d−1`, form a basis of `F_p^H` (Vandermonde, distinct
nodes). The form `⟨μ,ν⟩ = Σ_x μ(x)ν(x)` is nondegenerate, so the orthogonal complement of
`span{x^j : 1 ≤ j ≤ d−1}` is one-dimensional. The constant vector lies in it, since
`Σ_{x∈H}x^j = 0` for `d ∤ j`. Hence `μ ≡ c (mod p)`, and `0 ≤ μ < p` forces `μ = c`. The converse is
the same orthogonality. For `A + B`, `μ(x) ≤ min(m,n) < p` and `M_j = Σ_x μ(x)x^j`. ∎
(§A4: 240 random multiplicity functions and the `μ₄` decompositions.)

**Proposition 2.6 (PROVED; critical containment pairs refute the `α = β` extension).** Let
`p ≡ 1 (mod 4)` and `t ∈ νQ`. Then `A = {0,t}` has `C₀(A) = ∅` and `|C₁(A)| = (p−1)/4`. So
`(A, C₁(A))` is critical with `r = 0`, `|A| = 2` and `|B| = (p−1)/4`. For `p ≥ 13` this gives
`|A| ≠ |B|` and `A + B ≠ Q`, and the moment identities fail (Proposition 2.5). The pair
nevertheless satisfies Relations X, Y and (2.9) (Theorem 2.2, Proposition 2.4 with `Γ = 1`).

*Proof.* `x = 0` and `x = −t` fail (`χ(t) = χ(−t) = −1`). Substituting `x = ty` turns
`χ(x) = χ(x+t) = 1` into `χ(y) = χ(y+1) = −1`. Therefore
`4|C₁| = Σ_{y≠0,−1}(1−χ(y))(1−χ(y+1)) = (p−2) + 1 + 1 − 1 = p − 1`, using
`Σ_{y≠0,−1}χ(y) = −χ(−1) = −1`, `Σ_{y≠0,−1}χ(y+1) = −1` and `Σ_yχ(y)χ(y+1) = −1`.
Then `mn − r = d`. ∎

Witness: `p = 13`, `A = {0,2}`, `B = {1,10,12}`, `M₂ = 9 ≢ 0`. We enumerated exhaustively over
`A ∋ 0`, `|A| ≤ 4`, `p ∈ {13,17,19,29,37,41,53,61,73}`, which gives 192 critical `r = 0` pairs
(§A4b). In every one `M₁ = M_{d−1} = 0` (Corollary 2.3), and the first non-vanishing moment is `M₂`,
except for the nine `3 × 3` pairs at `p = 19`, where it is `M₃`. The pairs are symmetric, which
kills all odd moments when `m = 2`. The exhaustive search (§6) found only one family of critical
`r = 0` pairs with `min(m,n) ≥ 3`, at `p = 19` (nine `3 × 3` pairs). For them `M₃ = M₆ ≠ 0` and the
sums are not distinct. This is exactly the prime that Kalmynin treats separately.

### 2.5 The reach of coefficient comparison under a defect

The identities of §§2.2–2.3 prescribe the Taylor data of `Γ_A` on `B` as explicit functions of
`(A, B)`. By Theorem 2.1,
`Γ_A(b+Y) = [Σ_s C(D,m+s)σ_s(b)Y^s] / [λ(P₁(b+Y)/Y)^m P₀(b+Y)^{m−1}]` for `b ∈ B₁`. The top
coefficients of `Γ_A` are fixed by `[x^{d−j}]F_A = C(D,m−1+j)h_j(A)`. Kalmynin, Yip–Yoo and
Rudnev–Tyrrell turn such prescriptions into relations between `A` and `B` by **eliminating the
cofactor**. At `g = 0` it is known (`Γ = 1`), and every prescribed value is a relation.

**Proposition 2.7 (PROVED; reach).** Let `B ⊆ F_p` with `|B| = n`, `s, J ≥ 0`, `J ≤ g`. The
affine map from monic `Γ` of degree `g` to (the Taylor coefficients of order `≤ s` of `Γ` at each
`b ∈ B`; the coefficients of `x^{g−1}, …, x^{g−J}`) has linear part of rank `min(g, (s+1)n + J)`.
It is onto `F_p^{(s+1)n+J}` iff `g ≥ (s+1)n + J`. Consequently, if `g ≥ (s+1)n + J`, every
prescription of the order-`(s, J)` data of `Γ_A` by `(A, B)` is attained by some monic `Γ` of degree
`g`. Eliminating `Γ` then yields no relation in `(A, B)` beyond divisibility, i.e. beyond
`A + B ⊆ Q ∪ {0}` ([Bal] Theorem 1.1(3)).

*Proof.* Onto for `g ≥ (s+1)n + J`: fix the top `J` coefficients. The remaining part `R`
(`deg R < g − J`) must lie in a prescribed class modulo `Π_b(x−b)^{s+1}`, which has degree
`(s+1)n ≤ g − J` (Chinese remainder theorem). If `g < (s+1)n + J`, the domain has dimension
`g <` codomain. The rank formula: the kernel consists of polynomials of degree `< g − J`
divisible by `Π(x−b)^{s+1}`. ∎ (§A5: 120 random cases, rank formula exact.)

**Reading for the target.** Kalmynin's Relations X, Y and Lemma 5, and Yip–Yoo's (2.7), (2.21),
(2.9) and (2.19), use `s ≤ 2`, `J ≤ 3` on `B` (and on `−A` for the global identities). They
constrain `(A,B)` only while `g < 3n + 3` (resp. `3·max(m,n) + 3`), i.e. within `O(√p)` of the HP
bound in the balanced window. For `mn = (1/2 − c)p` one has `g ≈ 2c·d`. Any Kalmynin-type identity
carrying information there must involve Taylor data at `B` of order at least
`m + g/n − 1 ≈ (1 + 2c)·m`, that is, data of the same order as the vanishing order `m` that
Hanson–Petridis counts. This extends [Bal] Remark 1.2′ (`s = 0`, `J = 0`) to all orders and to
the top coefficients.

**Task 2, the constraint each surviving identity puts on `F_A = λP^mΓ_A`:**
(X_Γ) fixes `Γ_A′/Γ_A` on `B₁`. (Y_Γ) fixes the second log-coefficient. [Bal] Lemma 1.2 fixes
`Γ_A(b)`. Together they fix the order-2 Hermite data of `Γ_A` on `B₁`: `3|B₁|` affine conditions on
`g` free coefficients. Corollary 2.3 fixes `[x^{g−1}]Γ_A`, and Proposition 2.4 is the combination of
the order-1 data on `B` and on `−A`. These are identities in `F_p`, valid for every biclique. They
imply a count inequality (the only kind the LP certificates of [Bal] Theorem 3.5 can violate) only
through elimination of `Γ_A`, which Proposition 2.7 rules out once `g ≥ 3n + 3`. So they are not
testable against the certificates, and they add no count information in the near-tight regime.
The transfers of §3 are the surviving ingredients that *do* produce count constraints; §3.4 tests
them.

## 3. The containment-only transfers

### 3.1 Möbius transfer for cliques (Kalmynin Lemma 7)

**Proposition 3.1 (PROVED).** Let `p ≡ 1 (mod 4)` and `A` be a clique (`A − A ⊆ Q ∪ {0}`), `|A| = m`,
`a₀ ∈ A`, `A^{a₀} = {0} ∪ {1/(a₀−a) : a ∈ A∖{a₀}}`, and in difference form
`e_x(A) = #{a : χ(x−a) = −1}`. Then:
1. `A^{a₀}` is a clique of size `m` (so the same `g = d − m(m−1)`).
2. `x ↦ 1/(a₀ − x)` is a bijection `F_p ∖ A → F_p ∖ A^{a₀}`, and
   `e_{1/(a₀−x)}(A^{a₀}) = e_x(A)` if `χ(x−a₀) = 1`, and `m + 1 − e_x(A)` if `χ(x−a₀) = −1`.
3. Averaging over `a₀ ∈ A`, the profile `n_j = #{x ∉ A : e_x = j}` becomes
   `(Tn)_j = [(m−j)n_j + (m+1−j)n_{m+1−j}]/m` (`n_{m+1} := 0`), and the clique points stay fixed.
4. `T` fixes every profile with `j·n_j = (m+1−j)·n_{m+1−j}` for `1 ≤ j ≤ m`. In particular it fixes
   `n_j ∝ C(m,j)` (`j ≥ 1`), the profile of a random-like set, and any `n₀`.

*Proof.* (1) Kalmynin Lemma 7 (reproved: `1/(a₀−a′) − 1/(a₀−a″) = (a′−a″)/((a₀−a′)(a₀−a″)) ∈ Q`;
`χ(−1) = 1`). (2) For `a ≠ a₀`,
`1/(a₀−x) − 1/(a₀−a) = (x−a)/((a₀−x)(a₀−a))`, so `χ = χ(x−a)χ(a₀−x)`. Also
`χ(1/(a₀−x) − 0) = χ(a₀−x) = χ(x−a₀)`. So the `m−1` signs for `a ≠ a₀` are multiplied by `χ(x−a₀)`,
and the sign at the new element `0` equals the old sign at `a₀`. Bijectivity: `1/(a₀−x) = 0` never
occurs, and `1/(a₀−x) ∈ A^{a₀}∖{0}` iff `x ∈ A`. (3) `#{a₀ : χ(x−a₀) = 1} = m − e_x`,
`#{a₀ : χ(x−a₀) = −1} = e_x`. (4) Substitute. ∎ (§A6: 128 cliques, 647 transfers, 54,172 points;
Kalmynin's relations (1), (2) hold at `p = 41` and fail for all 45 near-critical cliques tested.)

### 3.2 Reciprocal transfer for bicliques (Rudnev–Tyrrell (1))

**Proposition 3.2 (PROVED).** Let `A + B ⊆ Q ∪ {0}`, `b ∈ B₁`, `S_b = {1/(a+b) : a ∈ A} ⊆ Q`, and
`ι_b(x) = 1/(x − b)` for `x ≠ b`. Then `ι_b(x) + 1/(a+b) = (x+a)/((x−b)(a+b))`, so

    χ(ι_b(x) + s_a) = χ(x+a)·χ(x−b),    e_{ι_b(x)}(S_b) = e_x(A) or m − δ_x − e_x(A)   (χ(x−b) = ±1),

and `δ_{ι_b(x)}(S_b) = δ_x`. The remaining point `0 = ι_b(∞)` has `e = 0` and `δ = 0`. In particular
`C₁(S_b) = {0} ∪ ι_b({x ∈ C₁(A) : χ(x−b)=1} ∪ {x ∈ AC(A) : χ(x−b) = −1})`: the transfer
exchanges complete and anti-complete points according to adjacency to `b`. For `p ≡ 1 (mod 4)` and
`r = 0`, let `u_x = e_x(A)` and `v_x = #{b : χ(x−b) = −1}`. The average over `b ∈ B` of the profile
of `S_b` is

    (TA)_j = Σ_{x∉B∪−A} [ (n−v_x)/n·[u_x=j] + v_x/n·[u_x=m−j] ] + Σ_{x∈B} [ (n−v_x)/n·[j=0] + v_x/n·[j=m] ],
    (TA)′_j = #{x ∈ −A : u_x = j}.

Here the extra point `0` of each `S_b` is booked in the `x ∈ B` term: `(n−1−v_b)/n + 1/n`.
Symmetrically, for `T_a = {1/(b+a)}` (profile in sum form, `x ↦ −x`), the average over `a ∈ A`
involves the same joint counts `(u_x, v_x)`. Summing the first-moment identity of `S_b` over `b`
gives the **Chung identity** `Σ_x f_A(x)g_B(x) = pr − mn`, with `f_A = Σ_aχ(·+a)` and
`g_B = Σ_bχ(·−b)`.

*Proof.* Direct computation, as in Proposition 3.1. `x ∈ −A` gives `v_x = 0` because
`χ(−a−b) = χ(a+b)`. The Chung identity is `Σ_x χ(x+a)χ(x−b) = −1 + p[a = −b]`. ∎
(§A7: 540,524 sign checks, 869 averaged profiles for `A` and for `B` compared with the direct
computation, 1,410 Chung identities.)

**Proposition 3.3 (PROVED; twisted HP = HP).** Let `S ⊆ Q`, `|S| = m`, and `T ⊆ F_p^*` with
`χ(s+t) = χ(t)` for all `s ∈ S`, `t ∈ T`. Then `m(|T| + 1) ≤ d`. For `S = S_b` the largest such `T` is
`{ι_b(x) : x ∈ C₁(A) ∖ {b}}`, so the inequality reads `m|C₁(A)| ≤ d`. This is the HP inequality
for `(A, C₁(A))`, and weaker than HP for `(A, C₁ ∪ C₀)`.

*Proof.* Rudnev–Tyrrell's computation for Lemma 3.4, without the hypothesis `d = α(|T|+1)`. Put
`F_tw = Σ_s k(s)(x+s)^D − C_S x^D`. The coefficient of `x^{D−i}` is `C(D,i)Σ_s c(s)s^{i−1}`, which is
`0` for `0 ≤ i ≤ m−1` and `C(D,m)` for `i = m`, so `deg F_tw = d − 1`. At `t ∈ T`, `(s+t)^d = t^d`
gives `[Y^j]F_tw(t+Y) = C(D,j)t^d[Σ_s k(s)(s+t)^{m−1−j} − C_S t^{m−1−j}] = 0` for `j ≤ m−1`. At `0`,
`s^d = 1` gives order `≥ m−1`. Counting roots: `m|T| + m − 1 ≤ d − 1`. For the maximal `T`: by
Proposition 3.2, `χ(ι_b(x)+s) = χ(ι_b(x))` for all `s` iff `χ(x+a) = 1` for all `a`. Points `x ∈ −A`
give `ι_b(x) + s = 0` for one `s` and are excluded. ∎ (§A7: `deg F_tw = d−1`, the orders, and
`T_max` for 1,410 pairs.)

### 3.3 What the transfers change: nothing at the level of the known count constraints

**Theorem 3.4 (PROVED by exact rational certificates).** Let `U(m)` be the family of count
constraints of [Bal] Theorem 3.5, valid for every `m`-set in `F_p`. It comprises the exact
identities, the Hankel inequalities (★)_e for `A` and `νA`, the averaged subset-HP inequalities, the
pencil inequalities, [Bal] Proposition 1.4, the Weil moment bounds of orders 4, 6, 8 and the Weil
counts.

*(Cliques.)* For each row of the first table, with `m` the largest integer satisfying
`m(m−1) ≤ d`, there is a nonnegative rational profile `(n_j)` on `F_p ∖ A` with the following
properties. It is supported on the band `⌈0.4m⌉ ≤ j ≤ m + 1 − ⌈0.4m⌉`, has `m` clique points, and
the profile together with `Tn`, `T²n` and `T³n` satisfies `U(m)`. Every inequality holds with the
stated relative slack.

| `p` | `m` | `m(m−1)/p` | `g` | rows | min slack | certificate |
|---|---|---|---|---|---|---|
| 1009 | 22 | 0.4579 | 42 | 557 | 0.0833 | `n₉ = 231, n₁₁ = 525, n₁₄ = 231` (`Tn` = `283.5, 262.5, 262.5, 178.5` on `9, 11, 12, 14`) |
| 10009 | 71 | 0.4966 | 34 | 1688 | 0.0068 | `T`-fixed: `n₂₉ ≈ 1089.9, n₃₆ ≈ 8113.0, n₄₃ ≈ 735.1` |
| 40009 | 141 | 0.4934 | 264 | 3298 | 0.0132 | results file |
| 100049 | 224 | 0.4993 | 72 | 5203 | 0.0014 | results file |
| 1000033 | 707 | 0.4991 | 874 | 16316 | 0.0017 | results file |

Sanity check: at `m + 1` (above the HP clique bound) the same LP is infeasible (`p = 1009, 10009`).

*(Bicliques, joint profile.)* For each row of the second table (`n = ⌊d/m⌋`, `r = 0`) there are
nonnegative rationals with the following properties. `X_{jk}` counts points `x ∉ B ∪ (−A)` with
`(u_x, v_x) = (j,k)`, `Y_k` counts the points of `B` and `Z_j` those of `−A`. All are supported on
the 40–60% bands. The four profiles (`A`; `−B`; the averaged `S_b`; the averaged `T_a`) satisfy
`U(m)`, `U(n)`, `U(m)` and `U(n)` respectively, and the Chung identity holds.

| `p` | `m` | `n` | `mn/p` | `g` | rows | min slack | `Y`, `Z` |
|---|---|---|---|---|---|---|---|
| 1009 | 22 | 22 | 0.4797 | 20 | 670 | 0.0397 | `Y₁₃ = 22`, `Z₁₁ = 22`; `X` on 6 cells |
| 10009 | 71 | 70 | 0.4966 | 34 | 2030 | 0.0068 | `Y₂₈ = 70`, `Z₄₂ = 71` |
| 40009 | 141 | 141 | 0.4969 | 123 | 4006 | 0.0061 | `Y₇₂`, `Z₇₉` |
| 100049 | 224 | 223 | 0.4993 | 72 | 6314 | 0.0014 | `Y₁₃₃`, `Z₁₃₄` |
| 1000033 | 707 | 707 | 0.4998 | 167 | 19854 | 0.00033 | `Y₂₈₃`, `Z₄₂₄` |

By weak LP duality, no nonnegative combination of these constraints, transfers included, proves
that an `m`-clique (resp. an `m × n` biclique with these parameters) does not exist. In particular,
none implies `m(m−1) ≤ (1/2 − c)p` (resp. `mn ≤ (1/2 − c)p`) for any `c` with
`(1/2 − c)p` below the tabulated `m(m−1)` (resp. `mn`).

*Method.* A HiGHS LP maximises the common slack. It is rounded to denominator 1000, the equalities
are repaired exactly (three pivots for cliques; six `X`-pivots for the joint LP, with `Y` and `Z`
integral), all auxiliary profiles are recomputed exactly from the definitions, and every row is
checked in rational arithmetic (27,062 clique rows and 32,874 joint rows). The band is a restriction,
not a theorem, and it makes the certificate stronger. Structural reason: `T` fixes binomial-type
profiles (Proposition 3.1(4)). The reciprocal map is a `B`-weighted reflection `j ↦ m − j`, and the
balanced constraints already hold for `νA`, i.e. for the reflected profile. A band-symmetric profile
is therefore invisible to both transfers. This settles [Bal] open obligation 4 for the joint
constraints coming from the critical-case papers: they do not make the LP infeasible.

## 4. The Hankel-minor robustness trick and the extra identities

**Proposition 4.1 (PROVED).** Let `x ∈ C₁(A)`, `1 ≤ e ≤ (m−1)/2`, `u_s = F_A^{(s)}/(D)_s`,
`H_{e+1} = det[u_{i+j}]_{0≤i,j≤e}` ([Rob] §2), `φ = t₀` and `κ = t₁/t₀ = (d−1)/(m+1)·S₁(x)`. Then

    H_{e+1}(x+Y) = φ^{e+1} Y^{(e+1)(m−e)} [ Δ_e + E_e κ Y + O(Y²) ],

where `Δ_e = det N₀`, `N₀ = [(m)_{i+j}/(D)_{i+j}]`, and `E_e = Σ_i det N^{(i)}`, with `N^{(i)}` equal
to `N₀` with row `i` replaced by `[(m+1)_{i+j}/(D)_{i+j}]_j`. These are universal constants.

*Proof.* `u_s(x+Y) = (φ/(D)_s)[(m)_sY^{m−s} + κ(m+1)_sY^{m+1−s} + O(Y^{m+2−s})]` because
`F_A(x+Y) = φY^m(1 + κY + …)`. Take `Y^{m−i}` out of row `i` and `Y^{−j}` out of column `j`, then
expand to first order by multilinearity. ∎ (§A8: 157 cases, `e = 1, 2`, exact.)

So at complete points the Hankel minors carry exactly the Relation-X datum `κ`. With
`H_{e+1} = P₁^{(e+1)(m−e)}P₀^{(e+1)(m−1−e)}K_e` and `deg K_e = (e+1)(g + e(n−1))` ([Bal]
Proposition 1.3(2)), the Hankel analogue of (X_Γ) reads, whenever `Δ_e ≢ 0 (mod p)` (so that
`K_e(b) ≠ 0`),
`(E_e/Δ_e)κ(b) = (e+1)(m−e)Σ_{B₁∖b}1/(b−b') + (e+1)(m−1−e)Σ_{B₀}1/(b−b') + K_e′(b)/K_e(b)`
(same proof as Theorem 2.2; checked for `e = 1`, `B = C₁(A)`, 110 points, §A8). The defect is
`K_e′/K_e`, with `deg K_e = (e+1)(g + e(n−1)) > g`, so Proposition 2.7 applies with an even larger
defect. The trick in [Rob] is a statement
about **vanishing orders** at points with `e_x ≤ e`: the rank-`e_x` perturbation `2ρ` is killed only in
orders. The leading coefficients at such points involve the bad partners `E(x)` (Step 3 of [Rob]
Theorem 2.1) and do not satisfy a `Γ`-free Relation X. **Answer to task 3's last question:** the
robustness trick applies formally, transporting Relation X with a larger defect. It creates no
defect-free identity, and it cannot, since even the non-robust identities (`e = 0`) are void at
`g ≍ p` (Proposition 2.7).

## 5. The obstruction

**Theorem 5.1 (status by item: (a) is a reading of the cited papers, CITED; (b), (c), (e) are
PROVED in §§2–3; (d) is descriptive, HEURISTIC).**
1. **(a)** Every identity that Kalmynin, Rudnev–Tyrrell and Yip–Yoo use beyond Hanson–Petridis is
   of one of four kinds: pure algebra (valid for all sets); a consequence of the critical
   factorisation `g = 0`; the moment identities; or the Möbius/reciprocal transfer (§1.4).
2. **(b)** The moment identities hold for a containment pair iff `A ⊕ B = Q` (Proposition 2.5).
   They fail for every pair with `g > 0` and for the critical families of Proposition 2.6. The
   `α = β` theorem does not extend to critical containment pairs.
3. **(c)** Every `g = 0` consequence they use is built from Taylor coefficients of order `≤ m + 2`
   on `B` (and on `−A`), glued into global identities by a degree count (Yip–Yoo), or from the top
   `≤ 3` coefficients. Its containment version holds with explicit defect
   terms (Theorem 2.2, Corollary 2.3, Proposition 2.4), and the defect is unconstrained as soon as
   `g ≥ 3n + 3` (Proposition 2.7).
4. **(d)** Their final step turns congruences with coefficients polynomial in `(α, k)` into integer
   equations by a size argument. This needs defect-free identities, and with `g ≥ 3n + 3` none of
   the order-`≤ 2` identities is defect-free.
5. **(e)** The transfers preserve `(m, n, g)` and act linearly on profiles, and the LP with all
   transferred constraints stays feasible at HP-tightness (Theorem 3.4).

Hence the target `|A||B| ≤ (1/2 − c)p + o(p)` cannot be reached by containment versions of these
ingredients as they stand. A proof along these lines would need at least one of the following:
1. identities of order `≳ (1 + 2c)m` at the points of `B`, i.e. comparable to the Hanson–Petridis
   vanishing order itself, together with a new mechanism that turns them into integer
   *inequalities* (Proposition 2.7 and §2.5);
2. a count constraint violated by the certificates of Theorem 3.4 and [Bal] Theorem 3.5, which by
   [Bal] §3.3 is a square-root-scale anti-concentration statement for `f_A`;
3. a genuinely joint input beyond the Chung identity and the transfers, which Theorem 3.4 shows is
   insufficient.

What survives, in summary: under containment, near-tightness gives no new identity. At exact
tightness (`g = 0`) the containment-only content is Kalmynin's Theorem 4 for cliques (at most one
unit off HP) and Corollary 2.3 / Proposition 2.4 for bicliques. None of these is an inequality.

## 6. Data: Hanson–Petridis-tight containment pairs (§A9)

We searched exhaustively over `A ∋ 0` with `2 ≤ |A| ≤ 5` (within a budget of `3·10⁵` sets per
`(p, m)`) and all primes `11 ≤ p ≤ 53`, for pairs `(A, C₁ ∪ C₀)` with `m|C₁| + (m−1)|C₀| = d`. A
tight pair must have `B = C(A)`, since removing a point lowers `mn − r`. Results:
- `r = 0`: `m = 2` for `p ≡ 1 (mod 4)` (Proposition 2.6), plus `p = 13` (`3 × 2`), `p = 17`
  (`4 × 2`) and `p = 19` (`3 × 3`). None has `A + B = Q`, and `(2.9)` holds for the `3 × 3` pairs
  at `p = 19`.
- `r > 0`, `min(m,n) ≥ 3`: `p = 13` (`3×3`, `r = 3`), `17` (`3×3`), `23` (`4×3`, `3×4`), `29`
  (`5×3`, `3×5`), `31` (`4×4`), `37` (`3×7`, `r = 3`) and `41` (the 5-cliques, `r = 5`).
- `p = 43, 47, 53`: only `m = 2`. This is consistent with [Bal]/`sigma-biclique`: `P(p) = (p+3)/2`
  for `43 ≤ p ≤ 1000`.

## 7. Literature check (2026-09-29)

The three papers above were read from the local copies. Yip–Yoo v2 was fetched from arXiv on
2026-09-29; the arXiv abstract page lists an earlier title, *…: a self-contained approach*, while the
v2 PDF is titled *…via differential identities*. None of the three states a containment version
of the Relations X/Y or of the moment identities, or a bound below Hanson–Petridis for
non-decompositions. The only containment result is Kalmynin's Theorem 4, at exact criticality.
Two web searches ("Hanson–Petridis bound improvement sumset contained in quadratic residues
Kalmynin critical pairs 2026"; "Hanson Petridis clique number Paley graph improved constant
sqrt(p/2) 2026") returned only these papers, Randomstrasse101 (arXiv:2603.29571), Yip's
prime-power papers and older material. **No improvement of the prime-field constant was found.
Completeness is not claimed.**

## 8. Open obligations

1. **The target (OPEN).** `|A||B| ≤ (1/2 − c)p + o(p)` for balanced complete bicliques, and
   `ω(G_p) ≤ (1/√2 − c)√p`.
2. **High-order Kalmynin relations.** Proposition 2.7 leaves open whether Taylor data of order
   `≈ (1+2c)m` on `B` (codimension `(s+1)n − g > 0`) yield anything beyond containment. They are
   polynomial identities of degree `≍ √p` in the elements of `A`, `B`. No mechanism converting them
   to inequalities is known.
3. **Iterated / higher transfers.** Theorem 3.4 uses one level of the reciprocal transfer and three
   iterates of the Möbius map. Transfers of the transferred twisted pairs `(S_b, T)` need triple
   joint profiles (`u_x`, `v_x`, adjacency to `b`); this was not attempted.
4. **Classify critical containment pairs.** Beyond the `m = 2` family, do `r = 0` critical pairs with
   `min(m,n) ≥ 3` exist for `p > 19`? Does (2.9) plus Kalmynin's coefficient comparisons exclude
   them without the moment identities? That would be the containment extension of Sárközy's
   conjecture at `g = 0`. It is not attempted here: it is an additive (`g = 0` vs `g ≥ 1`) statement
   and irrelevant to the constant.
5. **Robust points.** An explicit formula for the leading Taylor coefficient of `H_{e+1}` at points
   with `1 ≤ e_x ≤ e` (§4) was not derived.
6. Independent review of Theorem 2.2, Propositions 2.4, 2.5, 2.7, 3.1–3.3, 4.1 and of the two LP
   certificate routines (`clique_lp`, `joint_lp`).
