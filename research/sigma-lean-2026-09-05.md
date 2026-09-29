# Sigma pass, direction `lean`: Lean 4 / Mathlib formalization of the Paley two-set programme (2026-09-05)

**Status:** PROVED (machine-checked in Lean 4.30.0 / Mathlib `c5ea003`, axioms exactly `propext, Classical.choice, Quot.sound`): (i) shift orthogonality `∑_x χ(x)χ(x+c) = −1` for `c ≠ 0`, (ii) the exact second moment `∑_x (∑_{b∈B} χ(x−b))² = |B|(p−|B|)`, (iii) the Chung/Vinogradov bound `(∑_{a∈A}∑_{b∈B} χ(a−b))² ≤ |A||B|(p−|B|)`, (iv) the interval-to-least-non-residue implication. REFUTED: nothing (no claim tested false; every exact check passed). OPEN: the two-set Paley conjecture (goal), and, on the platform, the three literature statements (Karatsuba's amplification bound, both sentences of Hanson–Petridis Corollary 1.5), which are formalized faithfully as `sorry` statements and labelled CITED. **Nothing was uploaded to prove2me**: this session had no API key and was instructed not to authenticate; the complete proposal package is prepared for the user to upload (`~/prove2me_workspace/proposals/README.md`). Independent review of the faithfulness of the three literature statements is still an obligation.

## 1. Environment and method

- Workspace `/Users/shawwalters/prove2me_workspace`, pinned by `lean-toolchain` / `lakefile.lean` to `leanprover/lean4:v4.30.0` and Mathlib `c5ea00351c28e24afc9f0f84379aa41082b1188f`. This is the platform environment "Mathlib c5ea003 (Lean v4.30.0)" (non-default; default is v4.33.1 / `0df444a…`), as recorded in the workspace's `environments.json` snapshot.
- Mathlib cache present: `lake env lean Solutions/SmokeTest.lean` exits 0 in 9.2 s wall. No `lake exe cache get` and no download of any kind was performed.
- Every Lean run below is `lake env lean <file>` from the workspace root (`lake env lean` elaborates one file against the cached oleans; it never builds Mathlib). All runs are reproduced by `experiments/sigma_lean_2026_09_05.py` (stdlib only), which wrote `results/sigma_lean_2026_09_05.json`: **36820 checks, 0 failures, 79.3 s wall**.
- A `credentials.json` (210 bytes, mode 600) exists in the workspace, created by the concurrent Proximity-Prize session. It was not opened or used; no authenticated endpoint was called. Files `Thm_ProximityPrize_*` / `Sol_ProximityPrize_*` were not touched.

## 2. The eight statements (verbatim Lean, status, axioms, paths)

Common preamble of every file: `import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic` (plus `Mathlib.Order.Interval.Finset.Nat`, `Mathlib.Analysis.SpecialFunctions.Pow.Real` or `Mathlib.Data.Real.Sqrt` where used), `set_option autoImplicit false`, `open Finset`. Throughout, `quadraticChar (ZMod p)` is Mathlib's quadratic character; it is *definitionally* the Legendre symbol (`legendreSym p a = quadraticChar (ZMod p) ↑a` holds by `rfl`, checked), with `χ(0) = 0` and integer values.

### 2.1 PROVED — `Theorems/Thm_paley_shift_orthogonality.lean`, proof `Solutions/Sol_paley_shift_orthogonality.lean`
```lean
theorem paley_shift_orthogonality (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (c : ZMod p) (hc : c ≠ 0) :
    ∑ x : ZMod p, quadraticChar (ZMod p) x * quadraticChar (ZMod p) (x + c) = -1 := by sorry
```
Compile: exit 0, empty output (3.29 s). `#print axioms solution`: `['Classical.choice', 'Quot.sound', 'propext']`. Statement match with `theorem solution`: normalized text identical; `example : @paley_shift_orthogonality = @solution := rfl` elaborates (defeq of the full Π-types, using proof irrelevance).

### 2.2 PROVED — `Theorems/Thm_paley_second_moment.lean`, proof `Solutions/Sol_paley_second_moment.lean`
```lean
theorem paley_second_moment (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (B : Finset (ZMod p)) :
    ∑ x : ZMod p, (∑ b ∈ B, quadraticChar (ZMod p) (x - b)) ^ 2
      = (B.card : ℤ) * ((p : ℤ) - B.card) := by sorry
```
Compile: exit 0, empty output (3.25 s). Axioms: `['Classical.choice', 'Quot.sound', 'propext']`. Match: text + `rfl` both pass.

### 2.3 PROVED — `Theorems/Thm_paley_chung_bound.lean`, proof `Solutions/Sol_paley_chung_bound.lean`
```lean
theorem paley_chung_bound (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (A B : Finset (ZMod p)) :
    (∑ a ∈ A, ∑ b ∈ B, quadraticChar (ZMod p) (a - b)) ^ 2
      ≤ (A.card : ℤ) * B.card * ((p : ℤ) - B.card) := by sorry
```
Compile: exit 0, empty output (3.3 s). Axioms: `['Classical.choice', 'Quot.sound', 'propext']`. Match: text + `rfl` both pass.

### 2.4 PROVED — `Theorems/Thm_paley_interval_nonresidue.lean`, proof `Solutions/Sol_paley_interval_nonresidue.lean`
```lean
theorem paley_interval_nonresidue (p : ℕ) [Fact p.Prime] (N : ℕ) (hN : 2 * N < p)
    (h : ∑ a ∈ Icc 1 N, ∑ b ∈ Icc 1 N,
        quadraticChar (ZMod p) ((a : ZMod p) + (b : ZMod p)) < (N : ℤ) ^ 2) :
    ∃ n : ℕ, 2 ≤ n ∧ n ≤ 2 * N ∧ quadraticChar (ZMod p) (n : ZMod p) = -1 := by sorry
```
Compile: exit 0, empty output (4.02 s). Axioms: `['Classical.choice', 'Quot.sound', 'propext']`. Match: text + `rfl` both pass. (No oddness hypothesis is needed.)

### 2.5 OPEN (goal) — `Theorems/Thm_paley_two_set_conjecture.lean`
```lean
theorem paley_two_set_conjecture :
    ∀ ε : ℝ, 0 < ε → ε < 1 →
      ∃ δ : ℝ, 0 < δ ∧ ∃ p₀ : ℕ, ∀ (p : ℕ) [Fact p.Prime], p₀ < p →
        ∀ A B : Finset (ZMod p), (p : ℝ) ^ ε < A.card → (p : ℝ) ^ ε < B.card →
          ((|∑ a ∈ A, ∑ b ∈ B, quadraticChar (ZMod p) (a + b)| : ℤ) : ℝ)
            ≤ (p : ℝ) ^ (-δ) * A.card * B.card := by sorry
```
Elaborates with exactly one `declaration uses sorry` warning. Source: Satake, arXiv:2011.02907v2, Conjecture 7 (quoted verbatim in the docstring; local copy `sources/satake-2011.02907.txt`, lines 204–212). Encoding decisions: `ε = α`, `δ = β`, `p₀ = p(α)`; `a + b` instead of `s − t` (equivalent: `B ↦ −B` ranges over all subsets); `0 < ε < 1` omits only `α = 1`, for which `P(1, β)` is vacuous because `|S| > p` is impossible.

### 2.6 CITED (literature; open on platform) — `Theorems/Thm_paley_karatsuba_amplification.lean`
```lean
theorem paley_karatsuba_amplification :
    ∀ δ : ℝ, 0 < δ → ∃ C : ℝ, 0 < C ∧ ∀ (p : ℕ) [Fact p.Prime], p ≠ 2 →
      ∀ A B : Finset (ZMod p), (p : ℝ) ^ (1 / 2 + δ) < A.card → (p : ℝ) ^ δ < B.card →
        ((|∑ a ∈ A, ∑ b ∈ B, quadraticChar (ZMod p) (a + b)| : ℤ) : ℝ)
          ≤ C * (p : ℝ) ^ (-(0.05 * δ ^ 2)) * A.card * B.card := by sorry
```
Elaborates (one sorry warning). Source: M.-C. Chang, *Character sums in finite fields* (survey), §4, Theorem 4.4 (Karacuba), read on page 9 of the local PDF `sources/sigma-chang-character-sums-survey-DubProc.pdf` (rendered image checked: exponent is `p^{-0.05δ²}`, hypotheses `|A| > p^{1/2+δ}`, `|B| > p^{δ}`, χ non-trivial); Chang cites [Kar3] = A. A. Karacuba, *A certain arithmetic sum*, Soviet Math. Dokl. 12 (1971), no. 4, 1172–1174 (bibliography line 895 of the local text). Encoding decisions: specialized to the quadratic character; `p ≠ 2` added because for `p = 2` the quadratic character is trivial (`quadraticChar_eq_one_of_char_two`) and the source assumes non-triviality; the Vinogradov `≪` is read as a constant `C = C(δ) > 0` uniform in `p`, `A`, `B` (for fixed `δ` this is equivalent to the asymptotic statement, since each fixed prime contributes finitely many constraints).

### 2.7 CITED — `Theorems/Thm_paley_hanson_petridis_difference_bound.lean`
```lean
theorem paley_hanson_petridis_difference_bound (p : ℕ) [Fact p.Prime] (d : ℕ)
    (hd : d ∣ p - 1) (hd' : d < p - 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a - a' = 0 ∨ (a - a') ^ d = 1) :
    A.card * (A.card - 1) ≤ d := by sorry
```
Elaborates (one sorry warning). Source: Hanson–Petridis, arXiv:1905.09134v3 = Proc. Lond. Math. Soc. (3) 121 (2020) 287–292 (journal ref and DOI 10.1112/plms.12322 verified on the arXiv abstract page), Corollary 1.5, first sentence, local text lines 200–201: "Let p be a prime, d properly dividing p − 1 and suppose A ⊆ F_p is such that A − A ⊆ Z_d ∪ {0}. Then |A|(|A| − 1) ≤ d." Encoding: `Z_d = {z : z^d = 1}` (their §1 definition), "properly dividing" = `d ∣ p − 1 ∧ d < p − 1`; the product is in `ℕ` (`|A| = 0` gives `0 ≤ d`). Corner check: `d = 1` (if counted as proper) is consistent — for odd `p`, `A − A ⊆ {0, 1}` forces `|A| ≤ 1`.

### 2.8 CITED — `Theorems/Thm_paley_hanson_petridis_clique_number.lean`
```lean
theorem paley_hanson_petridis_clique_number (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1)
    (A : Finset (ZMod p)) (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    (A.card : ℝ) ≤ (Real.sqrt (2 * (p : ℝ) - 1) + 1) / 2 := by sorry
```
Elaborates (one sorry warning). Source: ibid., Corollary 1.5, second sentence: "In particular, for p ≡ 1 (mod 4), we have ω(G_p) ≤ (√(2p − 1) + 1)/2." Encoding: a clique is a finset all of whose differences of distinct elements satisfy `IsSquare` (a nonzero square, since `a ≠ a'`), i.e. `χ(a − a') = 1` (`quadraticChar_one_iff_isSquare`); the clique-number bound is stated for every clique. Consistency: `|A|(|A|−1) ≤ (p−1)/2` ⇔ `|A| ≤ (1 + √(2p−1))/2`.

Not drafted: Chang's small-doubling theorem (Duke Math. J. 145 (2008); abstract result (2) verified on Project Euclid: `|A|,|B| > p^{4/9+ε}`, `|B+B| < K|B|` ⇒ `|Σχ(x+y)| < p^{-τ}|A||B|`). Only the abstract was checked, not the theorem index or the dependence of `τ` on `(ε, K)` inside the paper, so no Lean statement was written; it is cited in the mission description's timeline only.

## 3. Proofs of the four PROVED statements (mathematical form; the Lean files are the certificates)

Let `p` be an odd prime, `χ` the quadratic character with `χ(0)=0`, so `χ(x)² = 1` for `x ≠ 0` and `∑_x χ(x) = 0` (`quadraticChar_sum_zero`, which needs `ringChar (ZMod p) ≠ 2`, i.e. `p ≠ 2`).

**(i) Shift orthogonality.** For `x ≠ 0`: `1 + c x⁻¹ = x⁻¹(x + c)`, so `χ(1 + c x⁻¹) = χ(x⁻¹)χ(x+c) = χ(x)χ(x+c)` (as `χ(x)χ(x⁻¹) = 1` and `χ(x)² = 1`). For `x = 0` both `χ(0)χ(c)` and `χ(1) − 1` vanish. Hence pointwise `χ(x)χ(x+c) = χ(1 + c x⁻¹) − [x = 0]`. The map `x ↦ 1 + c x⁻¹` is a bijection of `𝔽_p` (inverse `y ↦ c (y − 1)⁻¹`, with `0⁻¹ = 0`; `Function.bijective_iff_has_inverse`), so summing gives `∑_y χ(y) − 1 = −1`.

**(ii) Second moment.** Pair correlation: `∑_x χ(x−b)χ(x−b') = p − 1` if `b = b'` (all terms 1 except `x = b`), and `= −1` if `b ≠ b'` (substitute `x = y + b'` and apply (i) with `c = b' − b`). Expand the square and exchange sums: `∑_x F_B(x)² = ∑_{b,b'∈B} ([b=b'] p − 1) = |B|(p − |B|)`.

**(iii) Chung bound.** Cauchy–Schwarz on `A` with weights 1 (`Finset.sum_mul_sq_le_sq_mul_sq`): `(∑_{a∈A} F_B(a))² ≤ |A| ∑_{a∈A} F_B(a)²`; the terms are squares, so the sum extends to all of `𝔽_p` (`sum_le_sum_of_subset_of_nonneg`); then (ii).

**(iv) Interval link.** Suppose no `n ∈ [2, 2N]` has `χ(n) = −1`. For `a, b ∈ [1, N]`, `0 < a + b ≤ 2N < p`, so `a + b ≢ 0 (mod p)` (`ZMod.val_cast_of_lt`), hence `χ(a+b) ∈ {1, −1}` (`quadraticChar_dichotomy`) and by assumption `= 1`. So the double sum equals `N²`, contradicting `< N²`.

## 4. Exact-computation checks (all passed; from `results/sigma_lean_2026_09_05.json`)

Primes tested: 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53 (seed 20260905).

- Shift orthogonality: 364 checks (every `c ≠ 0`, every `p`).
- Second moment: 13408 checks (all `2^p` subsets for `p ≤ 13`, 300 random subsets for larger `p`).
- Chung bound: 22272 checks (all pairs `(A, B)` for `p ≤ 7`, 400 random pairs otherwise). Largest observed ratio `S²/(|A||B|(p−|B|))`: p=3: 1/2 on (A,B,S)=[[0], [1], -1]; p=5: 4/6 on (A,B,S)=[[0], [2, 3], -2]; p=7: 9/12 on (A,B,S)=[[0], [1, 2, 4], -3]; p=11: 49/90 on (A,B,S)=[[5, 6, 8], [0, 6, 8, 9, 10], -7].
- Interval ⇒ non-residue, and its converse identity (`no non-residue in [2,2N]` ⇒ sum `= N²`): 364 checks; least non-residues p=3: 2, p=5: 2, p=7: 3, p=11: 2, p=13: 2, p=17: 3, p=19: 2, p=23: 5, p=29: 2, p=31: 3, p=37: 2, p=41: 3, p=43: 2, p=47: 5, p=53: 2.
- Hanson–Petridis difference bound: 273 checks, exhaustive over all `A ⊆ 𝔽_p` with `A − A ⊆ Z_d ∪ {0}` for `p ≤ 13` and every proper divisor `d` of `p − 1`; extremal sizes (p=3, d=1: max|A|=1), (p=5, d=1: max|A|=1), (p=5, d=2: max|A|=2), (p=7, d=1: max|A|=1), (p=7, d=2: max|A|=2), (p=7, d=3: max|A|=1), (p=11, d=1: max|A|=1), (p=11, d=2: max|A|=2), (p=11, d=5: max|A|=1), (p=13, d=1: max|A|=1), (p=13, d=2: max|A|=2), (p=13, d=3: max|A|=1), (p=13, d=4: max|A|=2), (p=13, d=6: max|A|=3).
- Hanson–Petridis clique bound: exact clique numbers by Bron–Kerbosch: p=5: ω=2 ≤ 2.0, p=13: ω=3 ≤ 3.0, p=17: ω=3 ≤ 3.372, p=29: ω=4 ≤ 4.275, p=37: ω=4 ≤ 4.772, p=41: ω=5 ≤ 5.0, p=53: ω=5 ≤ 5.623.
- Not numerically testable: the goal (asymptotic, open) and Karatsuba's bound (asymptotic with an unspecified constant).

## 5. Proposal package (prepared, not uploaded)

All in `/Users/shawwalters/prove2me_workspace/proposals/`:
- `paley-mission-proposal.json` — mission metadata (`name`, 1,482-word `description`, `mission_type: OpenProblem`, `env: c5ea003…`, `field_ids: []` to be filled), 8 draft items in the platform's `/submit-problem` shape (`kind, theorem_name, theorem_title, formal_statement, natural_language_statement, preamble, source, tags`; `formal_statement` is read from the Lean file by `paley_build_proposal.py` and checked equal by the verifier), `main_item_theorem_name`, `item_order`, and 7 milestones (`milestone_title` starting with the source index, `milestone_description` verbatim from the source).
- `paley-mission-description.md` — the seven-section introduction per `references/mission_description.md`, with every citation re-verified this session (Volostnov's arXiv:1712.09355 was misattributed to "Volostnov–Shkredov" in the previous draft; fixed).
- `README.md` — exact `curl`/`jq` commands for: token exchange (`POST /agent/refresh`), field lookup, `POST /mission-proposals`, `POST …/items` for all eight, `PATCH` for `main_item_id`/`item_order`, `POST …/milestones` ×7, read-backs, hand-off, and post-launch `POST /verify` of the four solutions with the explanations `Solutions/Sol_paley_*.md`.

Why nothing was uploaded: no API key was available to this session and the orchestrator forbade authenticated calls. Two documented deviations from the captain playbook are left to the user: flat `paley_` names instead of a `namespace` (file names were fixed by the orchestrator), and no read-backs (sub-agents were forbidden by the brief).

## 6. Remaining obligations

1. Independent read-backs of items 2.5–2.8 (blind auditor per `references/mission_auditor.md`) before the human confirms the proposal.
2. Confirm on the platform that no mission/theorem named `paley_*` exists in the chosen environment (Step 2 of the README), and choose field ids.
3. If the default environment (Lean v4.33.1 / Mathlib `0df444a`) is preferred, re-verify all eight files there; only `c5ea003` was tested.
4. Formalize the literature layer: Weil's bound for `χ(f(x))` and Hölder amplification (Karatsuba), Stepanov's method (Hanson–Petridis). Both are open on the platform and untouched here.
5. The goal remains open; nothing in this direction bears on its truth.

## 7. Citation ledger

VERIFIED from local copies / arXiv abstract pages: Satake Conjecture 7, Remarks 8, 9, 16, 17, Theorem 10, Theorem 14 (`sources/satake-2011.02907.txt`); Hanson–Petridis §1 "Theorem (Vinogradov)", "Conjecture", Theorem 1.2, Corollary 1.5, journal ref + DOI (`sources/sigma-hanson-petridis-1905.09134.txt`, arXiv abs); Chang survey Problems 4.2–4.3, Theorem 4.4, Remark after 4.4, Corollary 1.3, [Kar3] (`sources/sigma-chang-character-sums-survey-DubProc.{pdf,txt}`); Chang Duke 2008 abstract result (2) (Project Euclid); Volostnov arXiv:1712.09355 title/author/abstract (arXiv abs). UNVERIFIED (labelled as such in the package): the theorem number in F. Chung, *Several generalizations of Weil's sums*, J. Number Theory 49 (1994) (cited via Satake [8]; the `(p−|B|)` form is proved here anyway); the Contemporary Mathematics volume number of Chang's survey (search result only: AMS Contemp. Math., 2010, pp. 83–98).
