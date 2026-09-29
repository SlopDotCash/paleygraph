# Stepanov wave, worker `leanhp`: Lean 4 formalization of Hanson–Petridis for prime Paley graphs (2026-09-26)

**Status:** PROVED and machine-checked **locally** in Lean 4 (no `sorry`; `#print axioms` = exactly
`propext, Classical.choice, Quot.sound` for every theorem), in two Mathlib versions, Lean v4.30.0-rc2 /
Mathlib `5450b53` and Lean v4.29.1 / Mathlib `5e932f9`: (1) `paley_hp_sharp_example` (for `p ≡ 1 mod 4`,
`B = {b : b, b+1 squares}` has `|B| = (p+3)/4`, and `A = {0,1}`, `B` attain equality in HP);
(2) `hanson_petridis_paley` (HP Theorem 1.2 for `d = (p−1)/2`); (2′) `hanson_petridis_theorem_1_2`
(HP Theorem 1.2 for every proper divisor `d` of `p − 1`); (3) `paley_clique_number_bound`
(`|A|(|A|−1) ≤ (p−1)/2`); and, as corollaries, full proofs of the two HP statements that the sigma pass
had uploaded to prove2me as CITED/Open (`paley_hanson_petridis_difference_bound`,
`paley_hanson_petridis_clique_number`). **NOT DONE: nothing was compiled in Mathlib `c5ea003` and nothing
was submitted to prove2me**, because the prove2me workspace `~/prove2me_workspace` (its Mathlib `c5ea003`
checkout and its `credentials.json`) no longer exists on this machine (§0); zero HTTP requests were made.
REFUTED: nothing. OPEN: server-side verification (needs the workspace and key restored by the user).

## 0. Environment blocker (read first)

- `/Users/shawwalters/prove2me_workspace` does not exist (`ls` → "No such file or directory"). Searches:
  `find` over every top-level directory of `$HOME` (depth ≤ 6) for `mathlib` package checkouts, and
  Spotlight (`mdfind`) for `Sol_paley_chung_bound.lean`, `lean-toolchain`, `lake-locked.sh`, found no
  prove2me workspace and no Mathlib at `c5ea003`. The only Lean projects with prebuilt Mathlib are
  `~/proximityprize` (Lean v4.30.0-rc2, Mathlib `5450b53e5ddc75d46418fabb605edbf36bd0beb6`) and
  `~/TheLeaningOfEverything` (Lean v4.29.1, Mathlib `5e932f97dd25535344f80f9dd8da3aab83df0fe6`).
  `~/.Trash` is not readable from this session (macOS permission), so whether the workspace was
  moved there is unknown.
- No `credentials.json` for prove2me was found (searched: every non-hidden top-level directory of
  `$HOME` and `~/.config`, depth ≤ 4; the three files of that name found belong to unrelated services
  and were not opened). Registering or re-authenticating was not
  attempted (out of scope and not permitted).
- Obtaining Mathlib `c5ea003` would need `git` (excluded by the brief) and a network download (which
  needs the user's permission); not attempted.
- What was done instead: every Lean file was elaborated with the toolchain's `lean` binary directly,
  with `LEAN_PATH` set to the prebuilt package directories of the two projects above. This is exactly
  what `lake env lean <file>` does, but it runs no `lake` command in those projects and writes nothing
  there (they belong to other sessions). Wrappers: `scratchpad/leanhp/runlean.sh` (5450b53) and
  `runlean2.sh` (5e932f9); the verifier reproduces both without the wrappers.
- Because the proofs compile unchanged in two Mathlib revisions of the Lean v4.29.1 and v4.30.0-rc2
  generations, it is *likely* but **not verified** that they compile in `c5ea003` (Lean v4.30.0).

## 1. Files

All under `experiments/stepanov_leanhp_lean/`:

| file | content |
|---|---|
| `StepanovHP.lean` | development file: Lagrange weights, auxiliary polynomial, HP core, Theorem 1.2, Paley case, clique bounds (317 lines) |
| `StepanovSharp.lean` | shift orthogonality `Σχ(x)χ(x+c) = −1` and the sharp example (135 lines) |
| `make_solutions.py` | generates the prove2me-format files below from the two development files |
| `Theorems/Thm_<name>.lean` | statement + `:= by sorry` (what `/submit-problem` would receive) |
| `Solutions/Sol_<name>.lean` | self-contained proof: imports + development namespace + `theorem solution` with the identical signature, closed by the development theorem, then `#print axioms solution` |

`<name>` ∈ {`hanson_petridis_paley`, `hanson_petridis_theorem_1_2`, `paley_clique_number_bound`,
`paley_hp_sharp_example`, `paley_hanson_petridis_difference_bound`, `paley_hanson_petridis_clique_number`}.

## 2. The statements (verbatim Lean; preamble `import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic`,
`import Mathlib.NumberTheory.LegendreSymbol.Basic` (or `Mathlib.Data.Real.Sqrt` for the last), `set_option autoImplicit false`, `open Finset`)

**Target 1 — PROVED (local).**
```lean
theorem paley_hp_sharp_example (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) :
    let A : Finset (ZMod p) := {0, 1}
    let B : Finset (ZMod p) := Finset.univ.filter (fun b => IsSquare b ∧ IsSquare (b + 1))
    B.card = (p + 3) / 4 ∧ (∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) ∧
      A.card * B.card = (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card
```
**Target 2 — PROVED (local).**
```lean
theorem hanson_petridis_paley (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (A B : Finset (ZMod p))
    (hAB : ∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) :
    A.card * B.card ≤ (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card
```
(`hp : p ≠ 2` is necessary: for `p = 2`, `A = B = {0,1}` gives `4 > 0 + 2`.)

**Target 2′ (general divisor, HP Theorem 1.2 exactly) — PROVED (local).**
```lean
theorem hanson_petridis_theorem_1_2 (p : ℕ) [Fact p.Prime] (d : ℕ) (hd : d ∣ p - 1) (hd' : d < p - 1)
    (A B : Finset (ZMod p)) (hAB : ∀ a ∈ A, ∀ b ∈ B, a + b = 0 ∨ (a + b) ^ d = 1) :
    A.card * B.card ≤ d + (B.filter (fun b => -b ∈ A)).card
```
**Target 3 — PROVED (local).**
```lean
theorem paley_clique_number_bound (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    A.card * (A.card - 1) ≤ (p - 1) / 2
```
(`hp` is kept for fidelity to the Paley graph; the proof uses only `p ≠ 2`, which it implies.)

**Existing private prove2me statements** (verbatim from `research/sigma-lean-2026-09-05.md` §2.7–2.8,
theorem ids `58b3552d…`, `cb279d48…`, status Open/CITED there) — **PROVED (local)**:
```lean
theorem paley_hanson_petridis_difference_bound (p : ℕ) [Fact p.Prime] (d : ℕ)
    (hd : d ∣ p - 1) (hd' : d < p - 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a - a' = 0 ∨ (a - a') ^ d = 1) :
    A.card * (A.card - 1) ≤ d
theorem paley_hanson_petridis_clique_number (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1)
    (A : Finset (ZMod p)) (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    (A.card : ℝ) ≤ (Real.sqrt (2 * (p : ℝ) - 1) + 1) / 2
```
Caveat: the original `Theorems/Thm_paley_hanson_petridis_*.lean` files were in the deleted workspace;
the statements above are transcribed from the sigma note, and the server's stored `formal_statement`
was not re-read (no key). A `WA` verdict would indicate a transcription mismatch, not a false proof.

## 3. The proof (as formalized)

Notation: `p` prime, `1 ≤ d`, `A, B ⊆ F_p`, `M = |A|`, `D = d + M − 1`, hypothesis
`H(a,b)`: `a + b = 0 ∨ (a + b)^d = 1`, `r = |{b ∈ B : −b ∈ A}|`.

**Lemma 1 (Lagrange weights; `sum_wt_eval`, `sum_wt_shift_pow`).** Put `c_a = [x^{M−1}] L_a(x)`, where
`L_a` is the Lagrange basis polynomial of the nodes `A` (Mathlib `Lagrange.basis A id a`). For every
`f` with `deg f < M`, `f = Σ_a f(a) L_a` (`Lagrange.eq_interpolate`), so
`Σ_a c_a f(a) = [x^{M−1}] f`. With `f = (x + t)^n`, `n < M`: `Σ_a c_a (t + a)^n = [n = M − 1]`.
(Non-vanishing of `c_a` is never needed.)

**Lemma 2 (coefficients; `coeff_auxF_comp`).** `F(x) = −1 + Σ_a c_a (x + a)^D`. For every `b` and `j`,
`[x^j] F(x + b) = −[j = 0] + binom(D, j) Σ_a c_a (b + a)^{D−j}` (`Polynomial.coeff_X_add_C_pow`).

**Lemma 3 (`F ≠ 0`, `deg F ≤ d`; `auxF_ne_zero`, `auxF_natDegree_le`).** Take `b = 0`. For `N > d`:
if `N ≤ D` then `D − N < M − 1`, so the sum vanishes by Lemma 1; if `N > D`, `binom(D, N) = 0`. For
`N = d`: `D − d = M − 1`, the sum is 1 and the coefficient is `binom(D, d)`, which is nonzero in `F_p`
when `D < p`, because `binom(D,d)·d!·(D−d)! = D!` and `p ∤ D!` (`Nat.Prime.dvd_factorial`). This is the
only place where primality of the field (as opposed to `F_{p²}`) enters, exactly as in HP.

**Lemma 4 (root multiplicity; `le_rootMultiplicity_auxF`).** If `d + M ≤ p` and `H(a, b)` holds for all
`a ∈ A`, then `mult_b(F) ≥ M − [−b ∈ A]`. Proof: `mult_b F = natTrailingDegree F(x + b)`
(`rootMultiplicity_eq_natTrailingDegree`), so it suffices that `[x^j]F(x+b) = 0` for `j < M − [−b∈A]`.
For each `a`, `(b+a)^{D−j} = (b+a)^{M−1−j}`: if `(b+a)^d = 1` split the exponent `d + (M−1−j)`; if
`b + a = 0` then `−b = a ∈ A`, so `j < M − 1` and both powers are `0`. Lemma 1 with `t = b`,
`n = M−1−j` gives `Σ_a c_a (b+a)^{M−1−j} = [j = 0]`, and `−[j=0] + binom(D,j)[j=0] = 0`.
This replaces HP's derivative computation `F^{(j)}(b) = (D)_j/((M−1)⋯(M−j)) · G^{(j)}(b)` by the
equivalent Taylor-coefficient computation; no factorial other than `D!` needs to be invertible.

**Core (`hp_core`).** If `d + M ≤ p` and `H` holds on `A × B`: the multiset `Σ_{b∈B}` of the lower
bounds is dominated by `F.roots` (`count_roots`), so
`Σ_{b∈B} (M − [−b∈A]) ≤ |roots F| ≤ deg F ≤ d` (`card_roots'`), i.e. `M|B| − r ≤ d`.

**Removing `d + M ≤ p` (`hp_of_two_mul_lt`, needs `2d < p`).** HP use `M ≤ |A + B| ≤ d + 1`. The
formalization instead argues: if `d + M > p` and `B ∋ b`, then `|A ∖ {−b}| ≥ M − 1 ≥ d + 1`; pick
`A' ⊆ A ∖ {−b}` with `|A'| = d + 1` (so `d + |A'| ≤ p`) and apply the core to `(A', {b})`:
`d + 1 ≤ d + 0`, contradiction. (HP's WLOG `|A| ≤ |B|` is not needed.)

**Theorem 1.2 (`hanson_petridis`).** `d ∣ p−1`, `d < p−1` give `p − 1 = dk` with `k ≥ 2`, so `1 ≤ d`
and `2d < p`. **Paley case**: `IsSquare x ⇒ x = 0 ∨ x^{(p−1)/2} = 1` (`x = r²`,
`r^{p−1} = 1`, `ZMod.pow_card_sub_one_eq_one`), and `d = (p−1)/2` satisfies `1 ≤ d`, `2d < p` for odd `p`.

**Clique corollaries.** `B = −A`: `a + (−a') = a − a'` is a square (0 when `a = a'`), and every
`b ∈ −A` has `−b ∈ A`, so `r = |A|`; HP gives `|A|² ≤ d + |A|`. The real form follows from
`2|A|(|A|−1) ≤ p − 1 ⇒ (2|A|−1)² ≤ 2p − 1` and `Real.abs_le_sqrt`.

**Sharp example.** For `x ≠ 0`: `χ(x)χ(x+c) = χ(1 + c x^{-1})`; with the `x = 0` correction and the
bijection `x ↦ 1 + c x^{-1}`, `Σχ(x)χ(x+c) = −1`. Pointwise
`(1+χ(b))(1+χ(b+1)) = 4·[b, b+1 squares] − 2[b = 0] − 2[b = −1]` (uses `χ(−1) = 1`,
`ZMod.exists_sq_eq_neg_one_iff`); summing, `p − 1 = 4|B| − 4`. Then `B ∩ (−{0,1}) = {0, −1}`,
`r = 2`, and `2·(p+3)/4 = (p−1)/2 + 2`.

## 4. Compile and axiom results

From the last run of `experiments/stepanov_leanhp_2026_09_26.py --with-dev --second-env`
(`results/stepanov_leanhp_2026_09_26.json`, 33.8 s wall, 0 failures). "E1" = Lean v4.30.0-rc2 /
Mathlib `5450b53`; "E2" = Lean v4.29.1 / Mathlib `5e932f9`. Every Solutions file prints exactly one
line, `'solution' depends on axioms: [propext, Classical.choice, Quot.sound]`, and nothing else
(no error, no warning, no `sorry`).

| file | E1 exit | E2 exit | axioms of `solution` |
|---|---|---|---|
| `Solutions/Sol_paley_hp_sharp_example.lean` | 0 | 0 | propext, Classical.choice, Quot.sound |
| `Solutions/Sol_hanson_petridis_paley.lean` | 0 | 0 | propext, Classical.choice, Quot.sound |
| `Solutions/Sol_hanson_petridis_theorem_1_2.lean` | 0 | 0 | propext, Classical.choice, Quot.sound |
| `Solutions/Sol_paley_clique_number_bound.lean` | 0 | 0 | propext, Classical.choice, Quot.sound |
| `Solutions/Sol_paley_hanson_petridis_difference_bound.lean` | 0 | 0 | propext, Classical.choice, Quot.sound |
| `Solutions/Sol_paley_hanson_petridis_clique_number.lean` | 0 | 0 | propext, Classical.choice, Quot.sound |
| `StepanovHP.lean` (5 `#print axioms`) | 0 | 0 (earlier run) | same three, each |
| `StepanovSharp.lean` (2 `#print axioms`) | 0 | 0 (earlier run) | same three, each |
| `Theorems/Thm_*.lean` (6 files) | 0, one `declaration uses sorry` warning each | — | (statements only) |

Statement match: for each of the six names, the text after `theorem <name>` in `Thm_<name>.lean` and
after `theorem solution` in `Sol_<name>.lean` is identical up to whitespace (checked by the verifier),
and no Solutions file imports `Theorems.*` or contains `sorry`.
**Mathlib `c5ea003`: not run** (§0). Server verdicts: **none** (§5).

## 5. prove2me: what was (not) done

- Requests made: **none** (no method, no path, no id). No theorem was created, no proof verified.
- Prepared for the user, in `experiments/stepanov_leanhp_lean/`: the six `Theorems/Thm_*.lean` (the
  four new ones would go to `POST /submit-problem` with `{"problems": [...], "env":
  "c5ea00351c28e24afc9f0f84379aa41082b1188f", "private": true}`), and the six `Solutions/Sol_*.lean`
  for `POST /verify` (`proof_type=prove`); the two existing private theorems `58b3552d…`
  (`paley_hanson_petridis_difference_bound`) and `cb279d48…` (`paley_hanson_petridis_clique_number`)
  need only `/verify`. Before any submission: restore `~/prove2me_workspace`, run
  `lake env lean` on each file there (the verifier does this automatically when the workspace exists),
  and check for name collisions with `GET /theorems?theorem_name=…&env=c5ea003…`.

## 6. Exact checks (`experiments/stepanov_leanhp_2026_09_26.py`, integers mod p only)

102532 checks, 0 failures (20 of them are the Lean runs and 6 the statement matches):

- HP Theorem 1.2 for every prime `3 ≤ p ≤ 47` and every proper divisor `d` of `p − 1`: all `A ⊆ F_p`
  for `p ≤ 13` (random `A`, `|A| ≤ 6`, 400 per `(p, d)` above), each with its maximal `B` (the binding
  case: each `b` adds `|A| − [−b∈A] ≥ 0` to `|A||B| − r`) plus random sub-`B`s: 86552 checks.
  Equality with `|A| ≥ 2` was observed for 50 of the 63 pairs `(p, d)` (search exhaustive only for
  `p ≤ 13`; sampled above, so the other 13 are not claimed to lack equality), e.g. `p = 7, d = 3`:
  `A = B = {0,1}`, `|A||B| = 4 = 3 + 1`.
- The squares-with-0 set equals `Z_{(p−1)/2} ∪ {0}` (14 primes): this is the Paley reduction.
- Internal claims of the proof on random `(A, d)`: Lagrange identity `Σ_a c_a (t+a)^n = [n = M−1]`
  for all `t`, `n < M` (13618), `deg F = d` exactly with leading coefficient `binom(D,d)` (230), and
  `mult_b F ≥ M − [−b ∈ A]` whenever `A + b ⊆ Z_d ∪ {0}` (395).
- The brief's bad-partner formula (Paley `d`; `b ∉ −A`; `E(b) = {a : χ(b+a) = −1}`):
  `[x^j] F(x+b) = −2·binom(D,j)·Σ_{a∈E(b)} c_a (b+a)^{M−1−j}` for `0 ≤ j ≤ M−1` (1662 checks), i.e.
  `F^{(j)}(b) = −2 (D)_j Σ_{E(b)} c_a (b+a)^{M−1−j}` as stated in the brief. PROVED (by hand): by
  Lemma 2, `[x^j]F(x+b) = −[j=0] + binom(D,j) Σ_a c_a χ(b+a)(b+a)^{M−1−j}`; write
  `χ = 1 − 2·1_E` and use Lemma 1 (`Σ_a c_a (b+a)^{M−1−j} = [j=0]`, `binom(D,0) = 1`).
- Sharp example `|B| = (p+3)/4`, `A + B ⊆ Q ∪ {0}`, `|A||B| = d + 2`: all 21 primes `p ≡ 1 (4)`, `p ≤ 197`.
- Clique numbers by exact Bron–Kerbosch for the 13 primes `p ≡ 1 (4)`, `p ≤ 109`: `ω(G_p)` =
  2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 6, 5, 6 for p = 5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109;
  `ω(ω−1) ≤ (p−1)/2` and `(2ω−1)² ≤ 2p−1` in every case (equality in the first at p = 5, 13, 41).
- `p = 2` witness that `hp : p ≠ 2` cannot be dropped: `A = B = {0,1}` ⊆ `F_2`, `4 > 0 + 2`.

## 7. Remaining obligations

1. Compile the six Solutions in Mathlib `c5ea003` (restore the workspace; rerun the verifier).
2. Private upload and server verdicts (needs the user's API key).
3. Priority: two web searches (2026-09-27) found no earlier Lean/Mathlib formalization of
   Hanson–Petridis; this is not an exhaustive check, so "first machine-checked proof" is UNVERIFIED.
4. The HP theorem remains a statement about **complete** bicliques; nothing here addresses the robust
   version RHP of the brief. The formal Lemmas 1–4 (in particular `coeff_auxF_comp`, the exact Taylor
   coefficients of `F` at any point `b`) are, however, stated without the biclique hypothesis and can
   be reused by a formal attempt at RHP.

## 8. Reproduction

- `python3 experiments/stepanov_leanhp_2026_09_26.py` — arithmetic checks, textual statement match, and
  `lean` on the six `Solutions/` and six `Theorems/` files in E1 (or `lake env lean` in
  `~/prove2me_workspace`, i.e. Mathlib `c5ea003`, automatically if that workspace exists again).
  Flags: `--with-dev` adds the two development files, `--second-env` re-runs the Solutions in E2,
  `--no-lean` skips Lean. No network access; writes only `results/stepanov_leanhp_2026_09_26.json`.
- `python3 experiments/stepanov_leanhp_lean/make_solutions.py` regenerates `Theorems/` and `Solutions/`
  from `StepanovHP.lean` / `StepanovSharp.lean`.
