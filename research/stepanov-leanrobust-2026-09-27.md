# Stepanov wave, worker `leanrobust`: Lean 4 formalization of the Hankel-minor Stepanov inequality (2026-09-27)

**Status:** PROVED and machine-checked **locally** in Lean 4 (compile results verified by the orchestrator, §4) (no `sorry`; every `#print axioms` line
is exactly `propext, Classical.choice, Quot.sound`), see §4 for the per-environment table:
(1) **Theorem 2.1** of `research/stepanov-robust-2026-09-26.md` — inequality (★), the exact degree
`(e+1)(d−e)` of the Hankel determinant and its order at every `b` — for **every** `e` with
`2e+1 ≤ |A|` and every `A ⊆ F_p` with `|A| ≤ (p+1)/2`, **unconditionally**; the leading coefficient
`Λ = det[C(D−i−j, m−1)]` is shown `≢ 0 (mod p)` by a new elementary argument (§3.4) that does not
use Krattenthaler's determinant; (2) **Corollary 2.3** (both displays); (3) **Theorem 2.5** (the
constant-bias bound for `|A||B| ≥ (1/2+κ)p`), exactly as stated in the robust note. Nothing is
conditional and there is no `_partial` file. REFUTED: only a prose remark of the robust note
("(★) is tight only at `e = 0`"; witness §5). OPEN: server-side (prove2me) verification — no
workspace/key on this machine; compilation in Mathlib `c5ea003`; see §7.

Verifier: `experiments/stepanov_leanrobust_2026_09_27.py` → `results/stepanov_leanrobust_2026_09_27.json`
(standard library only; recompiles the Lean file in both environments, read-only, plus ~2·10⁵ exact
modular checks). Lean source: `experiments/stepanov_leanrobust_lean/StepanovRobust.lean`.

---

## 0. Environment (as for `leanhp`)

- No prove2me workspace, no Mathlib `c5ea003`, no API key on this machine; nothing was submitted and no
  network request was made.
- Every compile calls the toolchain's `lean` binary directly with `LEAN_PATH` = the prebuilt package
  directories of an existing project (exactly what `lake env lean` does), running no `lake` command
  in, and writing nothing to, those projects (`lean <file>` without `-o` writes no `.olean`):
  - **E1**: Lean v4.30.0-rc2, Mathlib `5450b53` (`~/proximityprize/.lake/packages`);
  - **E2**: Lean v4.29.1, Mathlib `5e932f9` (`~/TheLeaningOfEverything/.lake/packages`).
- The machine was under very heavy load from other sessions (load average 150–360 on 16 cores)
  during this work; Lean used 1–10 % CPU, so wall-clock compile times (§4) are 10–30× the CPU time.
  The verifier therefore may exceed ten minutes wall-clock on a loaded machine (§8).

## 1. Files

| path | content |
|---|---|
| `experiments/stepanov_leanrobust_lean/StepanovRobust.lean` | the whole development (≈1380 lines, `import Mathlib`), namespace `StepanovRobust` |
| `experiments/stepanov_leanrobust_2026_09_27.py` | verifier (Lean re-compilation in E1/E2 + exact checks) |
| `results/stepanov_leanrobust_2026_09_27.json` | verifier output |

The development re-proves (copied from `leanhp`'s `StepanovHP.lean`, unchanged in substance) the
Lagrange-weight lemmas `sum_wt_eval`, `sum_wt_shift_pow`; everything else is new.

## 2. Definitions and statements (verbatim Lean)

Preamble: `import Mathlib`, `set_option autoImplicit false`, `open Polynomial Finset`,
`namespace StepanovRobust`, `variable {p : ℕ} [Fact p.Prime]`.

```lean
noncomputable def wt (A : Finset (ZMod p)) (a : ZMod p) : ZMod p :=
  (Lagrange.basis A id a).coeff (A.card - 1)            -- c_a = Π_{a'≠a} (a - a')⁻¹
noncomputable def uPoly (A : Finset (ZMod p)) (d s : ℕ) : (ZMod p)[X] :=
  C (if s = 0 then -1 else 0) + ∑ a ∈ A, C (wt A a) * (X + C a) ^ (d + A.card - 1 - s)
noncomputable def hankel (A : Finset (ZMod p)) (d e : ℕ) :
    Matrix (Fin (e + 1)) (Fin (e + 1)) (ZMod p)[X] :=
  Matrix.of fun i j => uPoly A d ((i : ℕ) + j)
noncomputable def hankelLead (p : ℕ) (d m e : ℕ) : ZMod p :=
  (Matrix.of fun i j : Fin (e + 1) =>
    (((d + m - 1 - ((i : ℕ) + j)).choose (m - 1) : ℕ) : ZMod p)).det
def charSum (A B : Finset (ZMod p)) : ℤ :=
  ∑ a ∈ A, ∑ b ∈ B, quadraticChar (ZMod p) (a + b)
```
`uPoly A d s` is the robust note's `u_s = F^{(s)}/(D)_s` (`F = −1 + Σ_a c_a(x+a)^D`, `D = d+m−1`)
written without derivatives (for `s ≥ 1`, `F^{(s)} = (D)_s Σ_a c_a (x+a)^{D−s}`); the identification
with derivatives is not formalized and is not needed, since (★) does not mention `u_s`.
`quadraticChar (ZMod p)` is the Legendre symbol with `χ(0) = 0`.

**Theorem 2.1 (★) — PROVED (local), unconditional.** With `e_b = #{a ∈ A : a + b non-square}`,
`δ_b = [−b ∈ A]`, `d = (p−1)/2`:
```lean
theorem hankel_stepanov_paley (hp : p ≠ 2) (A : Finset (ZMod p))
    (hA : A.card ≤ (p + 1) / 2) (e : ℕ) (he : 2 * e + 1 ≤ A.card) :
    ∑ b : ZMod p, ∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) ≤ (e + 1) * ((p - 1) / 2 - e)
```
The inner sum is `Σ_{i=e_b}^{e}(m − e − i − δ_b) = (e+1−e_b)(m − (3e+e_b)/2 − δ_b)`; it is empty when
`e_b > e`, and the `ℕ`-subtractions never truncate because `m − δ_b − i − e ≥ m − 1 − 2e ≥ 0`. So this
is exactly (★) of the robust note.

**Theorem 2.1 (degree and orders) — PROVED (local).**
```lean
theorem hankel_det_natDegree_eq (hp : p ≠ 2) (A : Finset (ZMod p))
    (hA : A.card ≤ (p + 1) / 2) (e : ℕ) (he : 2 * e + 1 ≤ A.card) :
    (hankel A ((p - 1) / 2) e).det.natDegree = (e + 1) * ((p - 1) / 2 - e)
theorem hankel_order_paley (hp : p ≠ 2) (A : Finset (ZMod p)) (hA1 : 1 ≤ A.card) (e : ℕ)
    (hed : 2 * e ≤ (p - 1) / 2) (b : ZMod p) :
    (X - C b) ^ (∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e))) ∣ (hankel A ((p - 1) / 2) e).det
```

**Lemma 2.2, non-vanishing part — PROVED (local).**
```lean
theorem hankelLead_ne_zero (d m e : ℕ) (hed : 2 * e ≤ d) (hem : 2 * e + 1 ≤ m)
    (hDp : d + m - 1 < p) : hankelLead p d m e ≠ 0
```
(`d` arbitrary here, not only `(p−1)/2`.) The hypothesis `d + m − 1 < p` cannot be dropped: for
`p = 3, d = 2, m = 3, e = 1` (`D = 4`) one has `Λ ≡ 0 (mod 3)` (verifier witness).

**General-exponent core (for reuse) — PROVED (local).** For any `d ≥ 1` such that every
`x ∈ F_p` has `x = 0 ∨ x^d = 1 ∨ x^d = −1`, with `badSet A d b = {a ∈ A : (b+a)^d = −1}`:
```lean
theorem hankel_stepanov_core (A : Finset (ZMod p)) (d e : ℕ) (hd1 : 1 ≤ d) (hed : 2 * e ≤ d)
    (hdx : ∀ x : ZMod p, x = 0 ∨ x ^ d = 1 ∨ x ^ d = -1)
    (hΛ : hankelLead p d A.card e ≠ 0) :
    ∑ b : ZMod p, ∑ i ∈ Finset.Icc (badSet A d b).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) ≤ (e + 1) * (d - e)
```
(`hankel_stepanov_paley` = this with `d = (p−1)/2`, `badSet = {¬IsSquare}` (Euler's criterion) and
`hΛ` discharged by `hankelLead_ne_zero`.)

**Abstract determinant lemma (Steps 2–3) — PROVED (local).** Over any integral domain `R`:
```lean
theorem pow_dvd_det_add_mul {n : ℕ} {κ : Type*} [Fintype κ]
    (V : Matrix (Fin n) κ R) (W : Matrix κ (Fin n) R) (Q : Matrix (Fin n) (Fin n) R)
    (π : R) (t : Fin n → ℕ) (ht : Antitone t) (hQ : ∀ i j, π ^ t i ∣ Q i j) :
    π ^ (∑ i ∈ Finset.univ.filter (fun i : Fin n => Fintype.card κ ≤ (i : ℕ)), t i) ∣
      (Q + V * W).det
```

**Corollary 2.3 — PROVED (local).** With `n = |B|`, `N_− = Σ_{b∈B} e_b`, `r = #{b ∈ B : −b ∈ A}`:
```lean
theorem bias_inequality (hp : p ≠ 2) (A : Finset (ZMod p)) (hA : A.card ≤ (p + 1) / 2) (e : ℕ)
    (he : 2 * e + 1 ≤ A.card) (B : Finset (ZMod p)) :
    ((A.card : ℤ) - 2 * e) * ((e + 1) * B.card -
        ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ)) ≤
      (e + 1) * ((((p - 1) / 2 : ℕ) : ℤ) - e + (B.filter fun b => -b ∈ A).card)
theorem bias_bound (hp : p ≠ 2) (A : Finset (ZMod p)) (hA : A.card ≤ (p + 1) / 2) (e : ℕ)
    (he : 2 * e + 1 ≤ A.card) (B : Finset (ZMod p)) :
    (charSum A B : ℝ) ≤ (A.card : ℝ) * B.card - (B.filter fun b => -b ∈ A).card -
      2 * (e + 1) * (B.card - ((((p - 1) / 2 : ℕ) : ℝ) - e + (B.filter fun b => -b ∈ A).card) /
        ((A.card : ℝ) - 2 * e))
```
together with `charSum_eq : charSum A B = |A||B| − r − 2N_−`.

**Theorem 2.5 — PROVED (local).**
```lean
theorem constant_bias (hp11 : 11 ≤ p) (κ : ℝ) (hκ0 : 0 < κ) (hκ1 : κ ≤ 3 / 2)
    (A B : Finset (ZMod p)) (hAB : (1 / 2 + κ) * p ≤ (A.card : ℝ) * B.card) :
    |(charSum A B : ℝ)| ≤
      (1 - (1 - 1 / Real.sqrt (1 + 2 * κ)) ^ 2 +
        (Real.sqrt ((1 / 2 + κ) * p) + 1) / (2 * ((p : ℝ) - 1))) * ((A.card : ℝ) * B.card)
```
This is the robust note's Theorem 2.5 verbatim (`u = (1+2κ)^{−1/2} = 1/√(1+2κ)`). The intermediate
per-subset bound (`|A| ≤ (p+1)/2`, `S ≤ [1 − (1−u)² + u(1−u)|A|/d]|A||B|`) is `subset_bound`, and
the `|A| ≤ |B|` upper bound after sub-sampling is `upper_bound`.

## 3. The proof as formalized

Notation: `m = |A|`, `K = m−1`, `D = d+K`, `c_a = wt A a`, `y = b + a`, `E = E(b) = badSet`,
`δ = δ_b`, `ρ_s = Σ_{a∈E} c_a (X+a)^{D−s}`, `ε_s = u_s − 2ρ_s` (`epsPoly`).

### 3.1 Step 1 (`epsPoly_coeff_comp`, `epsPoly_dvd`)
Taylor coefficients by the binomial theorem only (`coeff_X_add_C_pow`), no derivatives:
`[x^j] ε_s(x+b) = −[s=j=0] + C(D−s, j)·(Σ_a c_a y^{D−s−j} − 2Σ_{a∈E} c_a y^{D−s−j})`.
For `s+j+δ < m`, termwise `c_a y^{D−s−j} − 2[a∈E]c_a y^{D−s−j} = c_a y^{m−1−s−j}` (if `y^d = −1`:
`y^{D−s−j} = −y^{m−1−s−j}`; if `y^d = 1`: equal; if `y = 0`: then `δ = 1`, both exponents `≥ 1`,
both sides `0`), and `Σ_a c_a y^{m−1−s−j} = [s+j = 0]` (Lagrange, `sum_wt_shift_pow`). Hence the
coefficient is `0`, and `(X−b)^{m−δ−s} ∣ ε_s` (via `X^n ∣ P(x+b)` and composing with `x−b`).
No hypothesis `D < p` is used here.

### 3.2 Steps 2–3 (`hankel_eq_add_mul`, `pow_dvd_det_add_mul`, `hankel_order`)
`[u_{i+j}] = [ε_{i+j}] + V·W` with `V_{i,a} = 2c_a(X+a)^{D−2e}(X+a)^{e−i}`, `W_{a,j} = (X+a)^{e−j}`,
`a ∈ E` (no fractions: `(D−2e)+(e−i)+(e−j) = D−i−j`). Row `i` of `[ε_{i+j}]` is divisible by
`(X−b)^{t_i}`, `t_i = m−δ−(i+e)`, antitone in `i`. For the abstract lemma: multilinearity of `det`
in rows (`AlternatingMap.map_add_univ`) gives `det(Q+VW) = Σ_S det M_S` (rows in `S` from `Q`).
If `|S| + |κ| < n`, then `M_S = X·Y` with inner index `S ⊕ κ` of size `< n`; mapping to the fraction
field, `rank(XY) ≤ |S⊕κ| < n`, so `det M_S = 0` (`Matrix.rank_mul_le_right`, `rank_of_isUnit`). If
`|S| ≥ n − |κ|`, pull `π^{t_i}` out of each row in `S` (`Matrix.det_mul_column`) and use
`Σ_{i∈S} t_i ≥ Σ_{i ≥ |κ|} t_i` for antitone `t` (`sum_filter_le_sum_of_antitone`).
This gives `(X−b)^{Σ_{i=e_b}^{e}(m−δ−i−e)} ∣ det`.

### 3.3 Steps 4–5 (`hankel_det_natDegree_le`, `hankel_det_coeff`, `sum_le_natDegree_of_dvd`)
`deg u_s ≤ d−s` with `[x^{d−s}]u_s = C(D−s, m−1)` (Lagrange again). Leibniz expansion of `det` and
`Σ_i (d − σ(i) − i) = (e+1)(d−e)` give `natDegree ≤ (e+1)(d−e)` and top coefficient `Λ`.
Then `Σ_b ord_b ≤ #roots ≤ natDegree` (`le_rootMultiplicity_iff`, `count_roots`, `card_roots'`).

### 3.4 `Λ ≢ 0 (mod p)` without Krattenthaler (`hankelLead_ne_zero`) — PROVED
Assume `2e ≤ d`, `e ≤ K` (implied by `2e+1 ≤ m`), `D < p`. For `0 ≤ i, j ≤ e`, with falling
factorials `x^{\underline k}`:

(a) `C(D−i−j, K)·K!·(d−i)! = (D−i−e)!·(d−i)^{\underline j}·(D−i−j)^{\underline{e−j}}` (integers;
`choose_factor`), from `C(n,K)K!(n−K)! = n!`, `(d−i−j)!(d−i)^{\underline j} = (d−i)!` and
`(D−i−e)!(D−i−j)^{\underline{e−j}} = (D−i−j)!`.

(b) Hence in `F_p`: `C(D−i−j,K) = α_i·P_j(i)` with `α_i = (D−i−e)!/(K!(d−i)!) ≠ 0` (all
factorials of numbers `< p`) and `P_j(x) = (d−x)^{\underline j}(D−j−x)^{\underline{e−j}}`
(`leadPoly`, via `descPochhammer`), a polynomial of degree `≤ e`. So
`Λ = (Π_i α_i)·det[P_j(i)]_{i,j}`.

(c) If `det[P_j(i)] = 0`, some `c ≠ 0` has `Σ_j c_j P_j(i) = 0` for `i = 0..e`; `Q = Σ c_j P_j` has
degree `≤ e` and the `e+1` distinct roots `0, …, e` of `F_p` (`e < p`), so `Q = 0`. Evaluate at
`x = d−k`: `P_j(d−k) = k^{\underline j}(K+k−j)^{\underline{e−j}}`, which is `0` for `j > k` and equals
`k!·K^{\underline{e−k}} ≢ 0` for `j = k`. Induction on `k` gives `c_k = 0` for all `k`: contradiction.

This replaces the citation of Krattenthaler (3.12) in the robust note's Lemma 2.2 for the purpose of
Theorem 2.1 (it proves `p ∤ Λ`, not the product formula). The same argument needs only
`e ≤ m−1` rather than `2e+1 ≤ m`; only the latter is machine-checked.

### 3.5 Corollary 2.3 and Theorem 2.5
`bias_inequality`: pointwise (`pt_bound`), for each `b ∈ B`,
`(m−2e)(e+1−e_b) − (e+1)δ_b ≤ Σ_{i=e_b}^{e}(m−δ_b−i−e)` (both cases `e_b ≤ e`, `e_b > e`), summed over
`B ⊆ F_p` and combined with (★). `charSum_eq` from `χ(x) = 1 − [x=0] − 2[x non-square]`.
`constant_bias`: symmetry `S(A,B) = S(B,A)` (WLOG `|A| ≤ |B|`); dilation by a non-square `ν`
(`FiniteField.exists_nonsquare`) gives `S(νA,νB) = −S(A,B)`; sub-sampling `m₀ = ⌈(1/2+κ)p/|B|⌉₊`
with the double-counting identity `Σ_{|T|=k+1, T⊆A} Σ_{a∈T} g(a) = C(|A|−1,k)·Σ_{a∈A} g(a)`
(`sum_powersetCard_sum`); `m₀ ≤ √((1/2+κ)p)+1 ≤ (p+1)/2` for `p ≥ 11`; per subset, `e = ⌊θm₀⌋₊`,
`θ = (1−u)/2`, and the two cases of the robust note (`um₀n ≥ d+m₀`: the algebraic identity
`… = m₀(1−u)(d+m₀)(u²m₀n − d) ≥ 0`; otherwise the constant is `≥ 1` and `S ≤ m₀n`).

## 4. Compile and axiom results

Filled in by the orchestrator on 2026-09-29 (the worker was stopped before its final run).
The first independent compile found that in **E2** the proof of `subset_bound` exceeded the
default `maxHeartbeats 200000`, so `subset_bound`, `upper_bound` and `constant_bias` carried
`sorryAx` there, while E1 was clean. The fix is `set_option maxHeartbeats 1000000 in` on
`subset_bound` alone (a resource limit; no logic changed). The verifier's source scan also
flagged the English word "admit" in a doc comment; it now strips comments before scanning.
After both fixes, `experiments/stepanov_leanrobust_2026_09_27.py`:

| environment | exit | errors | `#print axioms` lines | all exactly `propext, Classical.choice, Quot.sound` |
|---|---|---|---|---|
| E1: Lean v4.30.0-rc2 / Mathlib `5450b53` | 0 | 0 | 17 / 17 | yes |
| E2: Lean v4.29.1 / Mathlib `5e932f9` | 0 | 0 | 17 / 17 | yes |

Together with 202,148 exact checks (0 failures). Machine-checked, in both environments:
`hankel_stepanov_paley` (Theorem 2.1, unconditional, with the elementary proof of
`Λ ≢ 0 (mod p)`), `hankel_det_natDegree_eq`, `hankel_order_paley`, `bias_inequality` and
`bias_bound` (Corollary 2.3), and `constant_bias` (Theorem 2.5, exactly as stated in the
robust note). Not done: compilation in the prove2me environment (Mathlib `c5ea003`) and a
server verdict (no workspace or key on this machine).

## 5. Exact checks (`experiments/stepanov_leanrobust_2026_09_27.py`)

PENDING.

## 6. Fidelity notes

- (★) is stated as a sum over **all** `b ∈ F_p` of `Σ_{i=e_b}^{e}`; for `e_b > e` the inner sum is empty,
  so this is the robust note's sum over `{b : e_b ≤ e}`.
- Hypotheses: exactly those of the robust note (`p` odd prime, `|A| ≤ (p+1)/2`, `0 ≤ e ≤ (|A|−1)/2`;
  Theorem 2.5: `p ≥ 11`, `0 < κ ≤ 3/2`). `hp : p ≠ 2` is necessary in the Paley statements
  (Euler's criterion); Theorem 2.5 derives it from `p ≥ 11`.
- The formal proof of Theorem 2.1 differs from the note's in two technical places: derivatives are
  replaced by Taylor coefficients (so `D < p` is used **only** in `Λ ≢ 0`), and the rank argument is
  done over the fraction field of `F_p[x]` for the matrices `M_S` of the row expansion.

## 7. Remaining obligations

1. Compile in Mathlib `c5ea003` and obtain server verdicts (needs the prove2me workspace and key).
2. Formalize the exact value of `Λ` (Krattenthaler's product) if it is ever needed; the proofs here
   need only `p ∤ Λ`.
3. Theorem 2.4 (RHP) and Proposition 2.9 (index-`k` subgroups) of the robust note are **not**
   formalized (`hankel_stepanov_core` is stated for `d` with `x^d ∈ {0, ±1}` only).

## 8. Reproduction

`/opt/miniconda3/bin/python3 experiments/stepanov_leanrobust_2026_09_27.py` (flags `--no-lean`,
`--e1-only`). No network; writes only `results/stepanov_leanrobust_2026_09_27.json`. The arithmetic
part takes about one minute; each Lean run needs ≈ 1–2 CPU-minutes (E1 and E2 run in parallel).
