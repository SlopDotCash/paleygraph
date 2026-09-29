import Mathlib

/-!
# The Hankel-minor Stepanov inequality over prime fields (robust Hanson–Petridis)

Worker `leanrobust`, Stepanov wave, 2026-09-27.  Formalizes Theorem 2.1, Lemma 2.2 (non-vanishing
part), Corollary 2.3 and Theorem 2.5 of `research/stepanov-robust-2026-09-26.md`.

* `pow_dvd_det_add_mul`: `π ^ (Σ_{i ≥ |κ|} t i) ∣ det (Q + V * W)` when row `i` of `Q` is divisible
  by `π ^ t i` (`t` antitone) and `V * W` factors through `κ` (Steps 2–3 of Theorem 2.1).
* `epsPoly_dvd`: `(X - b)^{|A| - δ_b - s} ∣ u_s - 2 ρ_s` (Step 1), with
  `u_s = -[s=0] + Σ_a c_a (X + a)^{D - s}` (`= F^{(s)}/(D)_s`, written without derivatives).
* `hankel_det_coeff`, `hankel_det_natDegree_le`: degree `≤ (e+1)(d-e)`, top coefficient
  `Λ = det[C(D-i-j, m-1)]` (Step 4); `hankelLead_ne_zero`: `Λ ≠ 0` in `ZMod p` (Lemma 2.2,
  non-vanishing part, proved by a factorisation + interpolation argument, not by Krattenthaler).
* `hankel_stepanov_paley`: inequality (★) (Theorem 2.1), unconditional.
* `bias_inequality`, `bias_bound`: Corollary 2.3.  `constant_bias`: Theorem 2.5.
-/

set_option autoImplicit false

open Polynomial Finset

namespace StepanovRobust

section DetLemma

variable {R : Type*} [CommRing R] [IsDomain R]

/-- A product `X * Y` through fewer than `n` intermediate indices is singular. -/
theorem det_mul_eq_zero_of_card_lt {n : ℕ} {ι : Type*} [Fintype ι]
    (X : Matrix (Fin n) ι R) (Y : Matrix ι (Fin n) R) (h : Fintype.card ι < n) :
    (X * Y).det = 0 := by
  classical
  apply IsFractionRing.injective R (FractionRing R)
  rw [map_zero, RingHom.map_det, RingHom.mapMatrix_apply, Matrix.map_mul]
  by_contra hne
  have hu : IsUnit ((X.map (algebraMap R (FractionRing R))) * (Y.map (algebraMap R (FractionRing R)))) :=
    (Matrix.isUnit_iff_isUnit_det _).2 (isUnit_iff_ne_zero.2 hne)
  have h1 := Matrix.rank_of_isUnit _ hu
  have h2 := Matrix.rank_mul_le_right (X.map (algebraMap R (FractionRing R)))
    (Y.map (algebraMap R (FractionRing R)))
  have h3 := Matrix.rank_le_card_height (Y.map (algebraMap R (FractionRing R)))
  rw [Fintype.card_fin] at h1
  omega

/-- Rows in `S` taken from `Q`, the others from `V * W`: if `|S| + |κ| < n` the determinant
vanishes. -/
theorem det_piecewise_eq_zero {n : ℕ} {κ : Type*} [Fintype κ]
    (V : Matrix (Fin n) κ R) (W : Matrix κ (Fin n) R) (Q : Matrix (Fin n) (Fin n) R)
    (S : Finset (Fin n)) (h : S.card + Fintype.card κ < n) :
    Matrix.det (S.piecewise Q (V * W) : Matrix (Fin n) (Fin n) R) = 0 := by
  classical
  let X : Matrix (Fin n) (S ⊕ κ) R := fun i c => match c with
    | Sum.inl s => if i = (s : Fin n) then 1 else 0
    | Sum.inr l => if i ∈ S then 0 else V i l
  let Y : Matrix (S ⊕ κ) (Fin n) R := fun c j => match c with
    | Sum.inl s => Q s j
    | Sum.inr l => W l j
  have hXY : (S.piecewise Q (V * W) : Matrix (Fin n) (Fin n) R) = X * Y := by
    ext i j
    rw [Matrix.mul_apply, Fintype.sum_sum_type]
    by_cases hi : i ∈ S
    · rw [Finset.piecewise_eq_of_mem _ _ _ hi]
      simp only [X, Y, if_pos hi, zero_mul, Finset.sum_const_zero, add_zero]
      rw [Finset.sum_eq_single ⟨i, hi⟩]
      · simp
      · intro s _ hs
        rw [if_neg, zero_mul]
        intro h'
        exact hs (Subtype.ext h'.symm)
      · simp
    · rw [Finset.piecewise_eq_of_notMem _ _ _ hi]
      simp only [X, Y, if_neg hi, Matrix.mul_apply]
      have : ∀ s : S, (if i = (s : Fin n) then (1 : R) else 0) * Q s j = 0 := by
        intro s
        rw [if_neg, zero_mul]
        intro h'
        exact hi (h' ▸ s.2)
      simp [this]
  rw [hXY]
  apply det_mul_eq_zero_of_card_lt
  simp only [Fintype.card_sum, Fintype.card_coe]
  omega

/-- For antitone `t` and `|S| ≥ n - k`, the sum of `t` over `S` dominates the sum over `i ≥ k`. -/
theorem sum_filter_le_sum_of_antitone {n k : ℕ} (t : Fin n → ℕ) (ht : Antitone t)
    (S : Finset (Fin n)) (hS : n ≤ S.card + k) :
    ∑ i ∈ Finset.univ.filter (fun i : Fin n => k ≤ (i : ℕ)), t i ≤ ∑ i ∈ S, t i := by
  classical
  set T := Finset.univ.filter (fun i : Fin n => k ≤ (i : ℕ)) with hT
  have hTcard : T.card + k ≤ n ∨ T = ∅ := by
    by_cases hkn : k ≤ n
    · left
      have h1 : (T.map Fin.valEmbedding) ⊆ Finset.Ico k n := by
        intro x hx
        simp only [Finset.mem_map, hT, Finset.mem_filter, Finset.mem_univ, true_and,
          Fin.valEmbedding_apply] at hx
        obtain ⟨i, hi, rfl⟩ := hx
        simp only [Finset.mem_Ico]
        exact ⟨hi, i.2⟩
      have h2 := Finset.card_le_card h1
      rw [Finset.card_map, Nat.card_Ico] at h2
      omega
    · right
      rw [Finset.eq_empty_iff_forall_notMem]
      intro i hi
      simp only [hT, Finset.mem_filter, Finset.mem_univ, true_and] at hi
      have := i.2
      omega
  rw [← Finset.sum_inter_add_sum_diff T S t, ← Finset.sum_inter_add_sum_diff S T t,
    Finset.inter_comm T S]
  apply Nat.add_le_add_left
  rcases (T \ S).eq_empty_or_nonempty with hE | ⟨i₀, hi₀⟩
  · rw [hE, Finset.sum_empty]; exact Nat.zero_le _
  · have hi₀' : k ≤ (i₀ : ℕ) := by
      simp only [Finset.mem_sdiff, hT, Finset.mem_filter, Finset.mem_univ, true_and] at hi₀
      exact hi₀.1
    have hkn : k < n := lt_of_le_of_lt hi₀' i₀.2
    set c := t ⟨k, hkn⟩
    have hup : ∀ i ∈ T \ S, t i ≤ c := by
      intro i hi
      simp only [Finset.mem_sdiff, hT, Finset.mem_filter, Finset.mem_univ, true_and] at hi
      exact ht (show (⟨k, hkn⟩ : Fin n) ≤ i from hi.1)
    have hlow : ∀ i ∈ S \ T, c ≤ t i := by
      intro i hi
      simp only [Finset.mem_sdiff, hT, Finset.mem_filter, Finset.mem_univ, true_and] at hi
      have hik : (i : ℕ) ≤ k := by have := hi.2; omega
      exact ht (show i ≤ (⟨k, hkn⟩ : Fin n) from hik)
    have hcard : (T \ S).card ≤ (S \ T).card := by
      have e1 := Finset.card_sdiff_add_card_inter T S
      have e2 := Finset.card_sdiff_add_card_inter S T
      rw [Finset.inter_comm] at e2
      rcases hTcard with h | h
      · omega
      · rw [h] at e1 ⊢; simp
    calc ∑ i ∈ T \ S, t i ≤ (T \ S).card • c := Finset.sum_le_card_nsmul _ _ _ hup
      _ ≤ (S \ T).card • c := by rw [smul_eq_mul, smul_eq_mul]; exact Nat.mul_le_mul_right _ hcard
      _ ≤ ∑ i ∈ S \ T, t i := Finset.card_nsmul_le_sum _ _ _ hlow

/-- **Row-filtered divisibility of a low-rank perturbation.** If `Q`'s row `i` is divisible by
`π ^ t i` with `t` antitone and `V * W` has inner dimension `k = |κ|`, then
`π ^ (∑_{i ≥ k} t i)` divides `det (Q + V * W)`. -/
theorem pow_dvd_det_add_mul {n : ℕ} {κ : Type*} [Fintype κ]
    (V : Matrix (Fin n) κ R) (W : Matrix κ (Fin n) R) (Q : Matrix (Fin n) (Fin n) R)
    (π : R) (t : Fin n → ℕ) (ht : Antitone t) (hQ : ∀ i j, π ^ t i ∣ Q i j) :
    π ^ (∑ i ∈ Finset.univ.filter (fun i : Fin n => Fintype.card κ ≤ (i : ℕ)), t i) ∣
      (Q + V * W).det := by
  classical
  have hexp : (Q + V * W).det =
      ∑ S : Finset (Fin n), Matrix.det (S.piecewise Q (V * W) : Matrix (Fin n) (Fin n) R) :=
    (Matrix.detRowAlternating : (Fin n → R) [⋀^Fin n]→ₗ[R] R).map_add_univ Q (V * W)
  rw [hexp]
  apply Finset.dvd_sum
  intro S _
  by_cases hS : S.card + Fintype.card κ < n
  · rw [det_piecewise_eq_zero V W Q S hS]
    exact dvd_zero _
  · rw [not_lt] at hS
    choose Q' hQ' using hQ
    let v : Fin n → R := fun i => if i ∈ S then π ^ t i else 1
    let N : Matrix (Fin n) (Fin n) R := fun i j => if i ∈ S then Q' i j else (V * W) i j
    have hM : (S.piecewise Q (V * W) : Matrix (Fin n) (Fin n) R) =
        Matrix.of fun i j => v i * N i j := by
      ext i j
      by_cases hi : i ∈ S
      · simp only [Finset.piecewise_eq_of_mem _ _ _ hi, Matrix.of_apply, v, N, if_pos hi]
        exact hQ' i j
      · simp only [Finset.piecewise_eq_of_notMem _ _ _ hi, Matrix.of_apply, v, N, if_neg hi,
          one_mul]
    rw [hM, Matrix.det_mul_column]
    apply Dvd.dvd.mul_right
    have hprod : ∏ i, v i = π ^ (∑ i ∈ S, t i) := by
      simp only [v]
      rw [Finset.prod_ite_mem, Finset.univ_inter, Finset.prod_pow_eq_pow_sum]
    rw [hprod]
    exact pow_dvd_pow π (sum_filter_le_sum_of_antitone t ht S (by omega))

end DetLemma

section Hankel

variable {p : ℕ} [Fact p.Prime]

/-- Lagrange leading-coefficient weights `c_a = [x^{|A|-1}] L_a(x)` (as in `StepanovHP`). -/
noncomputable def wt (A : Finset (ZMod p)) (a : ZMod p) : ZMod p :=
  (Lagrange.basis A id a).coeff (A.card - 1)

theorem sum_wt_eval (A : Finset (ZMod p)) (f : (ZMod p)[X]) (hf : f.degree < A.card) :
    ∑ a ∈ A, wt A a * f.eval a = f.coeff (A.card - 1) := by
  have h := Lagrange.eq_interpolate (s := A) (v := id) (Set.injOn_id _) hf
  conv_rhs => rw [h]
  rw [Lagrange.interpolate_apply, finset_sum_coeff]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [coeff_C_mul, wt, mul_comm]
  rfl

/-- `Σ_a c_a (t + a)^n = [n = |A| - 1]` for `n < |A|`. -/
theorem sum_wt_shift_pow (A : Finset (ZMod p)) (t : ZMod p) (n : ℕ) (hn : n < A.card) :
    ∑ a ∈ A, wt A a * (t + a) ^ n = if n = A.card - 1 then 1 else 0 := by
  have hdeg : ((X + C t) ^ n : (ZMod p)[X]).degree < A.card := by
    rw [degree_lt_iff_coeff_zero]
    intro m hm
    rw [coeff_X_add_C_pow, Nat.choose_eq_zero_of_lt (by omega), Nat.cast_zero, mul_zero]
  have h := sum_wt_eval A ((X + C t) ^ n) hdeg
  simp only [eval_pow, eval_add, eval_X, eval_C] at h
  have h' : ∑ a ∈ A, wt A a * (t + a) ^ n = ∑ a ∈ A, wt A a * (a + t) ^ n :=
    Finset.sum_congr rfl fun a _ => by rw [add_comm t a]
  rw [h', h, coeff_X_add_C_pow]
  split_ifs with hn'
  · subst hn'; simp
  · rw [Nat.choose_eq_zero_of_lt (by omega), Nat.cast_zero, mul_zero]

/-- `Σ_{a∈E} C(w a) (X + a)^N`: Taylor coefficients at `b`. -/
theorem coeff_shiftSum_comp (E : Finset (ZMod p)) (w : ZMod p → ZMod p) (N : ℕ) (b : ZMod p)
    (j : ℕ) :
    ((∑ a ∈ E, C (w a) * (X + C a) ^ N).comp (X + C b)).coeff j =
      (∑ a ∈ E, w a * (b + a) ^ (N - j)) * (N.choose j : ZMod p) := by
  rw [Polynomial.sum_comp, finset_sum_coeff, Finset.sum_mul]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [mul_comp, C_comp, pow_comp, add_comp, X_comp, C_comp, coeff_C_mul,
    show (X + C b + C a : (ZMod p)[X]) = X + C (b + a) by rw [C_add, add_assoc],
    coeff_X_add_C_pow, mul_assoc]

/-- `u_s = F^{(s)}/(D)_s` for `F = -1 + Σ_a c_a (X+a)^D`, `D = d + |A| - 1`, written without
derivatives: `u_s = -[s = 0] + Σ_a c_a (X + a)^{D - s}`. -/
noncomputable def uPoly (A : Finset (ZMod p)) (d s : ℕ) : (ZMod p)[X] :=
  C (if s = 0 then -1 else 0) + ∑ a ∈ A, C (wt A a) * (X + C a) ^ (d + A.card - 1 - s)

/-- `ρ_s = Σ_{a∈E} c_a (X + a)^{D - s}`. -/
noncomputable def rhoPoly (A : Finset (ZMod p)) (d : ℕ) (E : Finset (ZMod p)) (s : ℕ) :
    (ZMod p)[X] :=
  ∑ a ∈ E, C (wt A a) * (X + C a) ^ (d + A.card - 1 - s)

/-- The non-residue set of `b`: `E(b) = {a ∈ A : (b + a)^d = -1}`. -/
noncomputable def badSet (A : Finset (ZMod p)) (d : ℕ) (b : ZMod p) : Finset (ZMod p) :=
  A.filter fun a => (b + a) ^ d = -1

/-- `ε_s = u_s - 2 ρ_s`. -/
noncomputable def epsPoly (A : Finset (ZMod p)) (d : ℕ) (b : ZMod p) (s : ℕ) : (ZMod p)[X] :=
  uPoly A d s - C 2 * rhoPoly A d (badSet A d b) s

/-- **Step 1 (local structure).** The Taylor coefficients of `ε_s` at `b` vanish below
`|A| - δ_b - s`. -/
theorem epsPoly_coeff_comp (A : Finset (ZMod p)) (d : ℕ) (hd1 : 1 ≤ d) (b : ZMod p)
    (hb : ∀ a ∈ A, b + a = 0 ∨ (b + a) ^ d = 1 ∨ (b + a) ^ d = -1) (s j : ℕ)
    (hsj : s + j + (if -b ∈ A then 1 else 0) < A.card) :
    ((epsPoly A d b s).comp (X + C b)).coeff j = 0 := by
  classical
  have hsj' : s + j < A.card := by split_ifs at hsj <;> omega
  simp only [epsPoly, uPoly, rhoPoly, sub_comp, add_comp, mul_comp, C_comp]
  rw [coeff_sub, coeff_add, coeff_C_mul, coeff_C, coeff_shiftSum_comp, coeff_shiftSum_comp]
  set N := d + A.card - 1 - s with hN
  have hNj : N - j = d + (A.card - 1 - s - j) := by omega
  rw [hNj]
  have key : (∑ a ∈ A, wt A a * (b + a) ^ (d + (A.card - 1 - s - j))) -
      2 * ∑ a ∈ badSet A d b, wt A a * (b + a) ^ (d + (A.card - 1 - s - j)) =
      ∑ a ∈ A, wt A a * (b + a) ^ (A.card - 1 - s - j) := by
    rw [badSet, Finset.sum_filter, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun a ha => ?_
    by_cases hP : (b + a) ^ d = -1
    · rw [if_pos hP, pow_add, hP]; ring
    · rw [if_neg hP, mul_zero, sub_zero]
      rcases hb a ha with h0 | h1 | h2
      · have hmem : -b ∈ A := by
          have : a = -b := by linear_combination h0
          rw [← this]; exact ha
        rw [if_pos hmem] at hsj
        rw [h0, zero_pow (by omega), zero_pow (by omega)]
      · rw [pow_add, h1, one_mul]
      · exact absurd h2 hP
  have hsum := sum_wt_shift_pow A b (A.card - 1 - s - j) (by omega)
  have e1 : (∑ a ∈ A, wt A a * (b + a) ^ (d + (A.card - 1 - s - j))) * (N.choose j : ZMod p) -
      2 * ((∑ a ∈ badSet A d b, wt A a * (b + a) ^ (d + (A.card - 1 - s - j))) *
        (N.choose j : ZMod p)) =
      ((∑ a ∈ A, wt A a * (b + a) ^ (d + (A.card - 1 - s - j))) -
      2 * ∑ a ∈ badSet A d b, wt A a * (b + a) ^ (d + (A.card - 1 - s - j))) *
        (N.choose j : ZMod p) := by ring
  rw [add_sub_assoc, e1, key, hsum]
  by_cases hs0 : s = 0
  · by_cases hj0 : j = 0
    · subst hs0; subst hj0; simp
    · rw [if_neg hj0, if_neg (by omega)]; simp
  · rw [if_neg (show A.card - 1 - s - j ≠ A.card - 1 by omega)]
    simp [hs0]

theorem X_sub_C_pow_dvd_of_coeff_comp (P : (ZMod p)[X]) (b : ZMod p) (n : ℕ)
    (h : ∀ j < n, (P.comp (X + C b)).coeff j = 0) : (X - C b) ^ n ∣ P := by
  have h1 : X ^ n ∣ P.comp (X + C b) := by rw [X_pow_dvd_iff]; exact h
  obtain ⟨Rq, hR⟩ := h1
  refine ⟨Rq.comp (X - C b), ?_⟩
  have hP : P = (P.comp (X + C b)).comp (X - C b) := by
    rw [comp_assoc]; simp
  rw [hP, hR, mul_comp, pow_comp, X_comp]

theorem epsPoly_dvd (A : Finset (ZMod p)) (d : ℕ) (hd1 : 1 ≤ d) (b : ZMod p)
    (hb : ∀ a ∈ A, b + a = 0 ∨ (b + a) ^ d = 1 ∨ (b + a) ^ d = -1) (s : ℕ) :
    (X - C b) ^ (A.card - (if -b ∈ A then 1 else 0) - s) ∣ epsPoly A d b s :=
  X_sub_C_pow_dvd_of_coeff_comp _ b _ fun j hj =>
    epsPoly_coeff_comp A d hd1 b hb s j (by split_ifs at hj ⊢ <;> omega)

/-- The Hankel matrix `[u_{i+j}]_{0 ≤ i,j ≤ e}`. -/
noncomputable def hankel (A : Finset (ZMod p)) (d e : ℕ) :
    Matrix (Fin (e + 1)) (Fin (e + 1)) (ZMod p)[X] :=
  Matrix.of fun i j => uPoly A d ((i : ℕ) + j)

/-- **Step 2 (rank).** `[u_{i+j}] = [ε_{i+j}] + V * W` with inner dimension `|E(b)|`. -/
theorem hankel_eq_add_mul (A : Finset (ZMod p)) (d e : ℕ) (hed : 2 * e ≤ d) (hA1 : 1 ≤ A.card)
    (b : ZMod p) :
    hankel A d e =
      (Matrix.of fun i j : Fin (e + 1) => epsPoly A d b ((i : ℕ) + j)) +
      (Matrix.of fun (i : Fin (e + 1)) (a : badSet A d b) =>
          C (2 * wt A a) * (X + C (a : ZMod p)) ^ (d + A.card - 1 - 2 * e) *
            (X + C (a : ZMod p)) ^ (e - i)) *
      (Matrix.of fun (a : badSet A d b) (j : Fin (e + 1)) => (X + C (a : ZMod p)) ^ (e - j)) := by
  ext i j
  have hi := i.isLt
  have hj := j.isLt
  simp only [hankel, Matrix.add_apply, Matrix.of_apply, Matrix.mul_apply, epsPoly, rhoPoly]
  have hterm : ∀ a ∈ badSet A d b, C (2 * wt A a) * (X + C a) ^ (d + A.card - 1 - 2 * e) *
      (X + C a) ^ (e - (i : ℕ)) * (X + C a) ^ (e - (j : ℕ)) =
      C 2 * (C (wt A a) * (X + C a) ^ (d + A.card - 1 - ((i : ℕ) + j))) := by
    intro a _
    rw [mul_assoc (C (2 * wt A a) * _), ← pow_add, mul_assoc, ← pow_add, C_mul,
      show d + A.card - 1 - 2 * e + (e - (i : ℕ) + (e - (j : ℕ))) =
        d + A.card - 1 - ((i : ℕ) + j) by omega]
    ring
  rw [Finset.sum_coe_sort (badSet A d b) (fun a => C (2 * wt A a) *
      (X + C a) ^ (d + A.card - 1 - 2 * e) * (X + C a) ^ (e - (i : ℕ)) *
      (X + C a) ^ (e - (j : ℕ))), Finset.sum_congr rfl hterm, ← Finset.mul_sum, sub_add_cancel]

end Hankel

section Degree

variable {p : ℕ} [Fact p.Prime]

theorem coeff_uPoly (A : Finset (ZMod p)) (d s N : ℕ) :
    (uPoly A d s).coeff N = (if N = 0 then (if s = 0 then (-1 : ZMod p) else 0) else 0) +
      (∑ a ∈ A, wt A a * a ^ (d + A.card - 1 - s - N)) *
        ((d + A.card - 1 - s).choose N : ZMod p) := by
  rw [uPoly, coeff_add, coeff_C]
  congr 1
  have h := coeff_shiftSum_comp A (wt A) (d + A.card - 1 - s) 0 N
  simpa using h

theorem natDegree_uPoly_le (A : Finset (ZMod p)) (d s : ℕ) (hs : s ≤ d)
    (hA1 : 1 ≤ A.card) : (uPoly A d s).natDegree ≤ d - s := by
  rw [natDegree_le_iff_coeff_eq_zero]
  intro N hN
  rw [coeff_uPoly, if_neg (by omega), zero_add]
  by_cases hND : N ≤ d + A.card - 1 - s
  · have h := sum_wt_shift_pow A 0 (d + A.card - 1 - s - N) (by omega)
    simp only [zero_add] at h
    rw [h, if_neg (by omega), zero_mul]
  · rw [Nat.choose_eq_zero_of_lt (by omega), Nat.cast_zero, mul_zero]

theorem coeff_uPoly_top (A : Finset (ZMod p)) (d s : ℕ) (hd1 : 1 ≤ d) (hs : s ≤ d)
    (hA1 : 1 ≤ A.card) :
    (uPoly A d s).coeff (d - s) = ((d + A.card - 1 - s).choose (A.card - 1) : ZMod p) := by
  rw [coeff_uPoly]
  have h0 : (if d - s = 0 then (if s = 0 then (-1 : ZMod p) else 0) else 0) = 0 := by
    split_ifs with h1 h2
    · omega
    · rfl
    · rfl
  rw [h0, zero_add, show d + A.card - 1 - s - (d - s) = A.card - 1 by omega]
  have h := sum_wt_shift_pow A 0 (A.card - 1) (by omega)
  simp only [zero_add, if_true] at h
  rw [h, one_mul, show d + A.card - 1 - s = (d - s) + (A.card - 1) by omega, Nat.choose_symm_add]

theorem coeff_prod_of_natDegree_le' {ι : Type*} [DecidableEq ι] (s : Finset ι)
    (f : ι → (ZMod p)[X]) (n : ι → ℕ) (h : ∀ i ∈ s, (f i).natDegree ≤ n i) :
    (∏ i ∈ s, f i).coeff (∑ i ∈ s, n i) = ∏ i ∈ s, (f i).coeff (n i) := by
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.prod_insert ha, Finset.sum_insert ha, Finset.prod_insert ha,
      coeff_mul_add_eq_of_natDegree_le (h a (Finset.mem_insert_self a s))
        ((natDegree_prod_le _ _).trans
          (Finset.sum_le_sum fun i hi => h i (Finset.mem_insert_of_mem hi))),
      ih (fun i hi => h i (Finset.mem_insert_of_mem hi))]

theorem sum_perm_sub (e d : ℕ) (hed : 2 * e ≤ d) (σ : Equiv.Perm (Fin (e + 1))) :
    ∑ i : Fin (e + 1), (d - (((σ i : Fin (e + 1)) : ℕ) + i)) = (e + 1) * (d - e) := by
  have h1 : ∑ i : Fin (e + 1), ((σ i : Fin (e + 1)) : ℕ) = ∑ i : Fin (e + 1), (i : ℕ) :=
    Equiv.sum_comp σ (fun i => (i : ℕ))
  have h2 : (∑ i : Fin (e + 1), (i : ℕ)) * 2 = (e + 1) * e := by
    rw [Fin.sum_univ_eq_sum_range (fun i => i) (e + 1), Finset.sum_range_id_mul_two]
    simp
  have h3 : ∑ i : Fin (e + 1), (d - (((σ i : Fin (e + 1)) : ℕ) + i)) +
      (∑ i : Fin (e + 1), ((σ i : Fin (e + 1)) : ℕ) + ∑ i : Fin (e + 1), (i : ℕ)) =
      ∑ _i : Fin (e + 1), d := by
    rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    have := (σ i).isLt
    have := i.isLt
    omega
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, smul_eq_mul, h1] at h3
  have h4 : (e + 1) * d = (e + 1) * (d - e) + (e + 1) * e := by
    rw [← mul_add]; congr 1; omega
  generalize (e + 1) * d = P1 at h3 h4
  generalize (e + 1) * e = P2 at h2 h4
  generalize (e + 1) * (d - e) = P3 at h4 ⊢
  omega

theorem smul_units_coeff (u : ℤˣ) (P : (ZMod p)[X]) (N : ℕ) :
    (u • P).coeff N = u • P.coeff N := by
  rw [Units.smul_def, Units.smul_def, zsmul_eq_mul, zsmul_eq_mul, ← C_eq_intCast, coeff_C_mul]

/-- The leading coefficient `Λ = det[C(D - i - j, m - 1)]_{0 ≤ i,j ≤ e}`, `D = d + m - 1`. -/
noncomputable def hankelLead (p : ℕ) (d m e : ℕ) : ZMod p :=
  (Matrix.of fun i j : Fin (e + 1) =>
    (((d + m - 1 - ((i : ℕ) + j)).choose (m - 1) : ℕ) : ZMod p)).det

theorem hankel_det_natDegree_le (A : Finset (ZMod p)) (d e : ℕ) (hed : 2 * e ≤ d)
    (hA1 : 1 ≤ A.card) : (hankel A d e).det.natDegree ≤ (e + 1) * (d - e) := by
  rw [Matrix.det_apply]
  refine natDegree_sum_le_of_forall_le _ _ fun σ _ => ?_
  rw [Units.smul_def]
  refine (natDegree_smul_le _ _).trans ?_
  refine (natDegree_prod_le _ _).trans ?_
  rw [← sum_perm_sub e d hed σ]
  refine Finset.sum_le_sum fun i _ => ?_
  simp only [hankel, Matrix.of_apply]
  exact natDegree_uPoly_le A d _ (by have := (σ i).isLt; have := i.isLt; omega) hA1

/-- **Step 4 (top coefficient).** -/
theorem hankel_det_coeff (A : Finset (ZMod p)) (d e : ℕ) (hd1 : 1 ≤ d) (hed : 2 * e ≤ d)
    (hA1 : 1 ≤ A.card) :
    (hankel A d e).det.coeff ((e + 1) * (d - e)) = hankelLead p d A.card e := by
  rw [Matrix.det_apply, hankelLead, Matrix.det_apply, finset_sum_coeff]
  refine Finset.sum_congr rfl fun σ _ => ?_
  rw [smul_units_coeff]
  congr 1
  rw [← sum_perm_sub e d hed σ, coeff_prod_of_natDegree_le']
  · refine Finset.prod_congr rfl fun i _ => ?_
    simp only [hankel, Matrix.of_apply]
    rw [coeff_uPoly_top A d _ hd1 (by have := (σ i).isLt; have := i.isLt; omega) hA1]
  · intro i _
    simp only [hankel, Matrix.of_apply]
    exact natDegree_uPoly_le A d _ (by have := (σ i).isLt; have := i.isLt; omega) hA1

theorem sum_count_le_card {α : Type*} [DecidableEq α] (s : Multiset α) (B : Finset α) :
    ∑ b ∈ B, s.count b ≤ Multiset.card s := by
  have h1 : ∑ b ∈ B, s.count b ≤ ∑ b ∈ B ∪ s.toFinset, s.count b :=
    Finset.sum_le_sum_of_subset subset_union_left
  have h2 : ∑ b ∈ s.toFinset, s.count b = ∑ b ∈ B ∪ s.toFinset, s.count b :=
    Finset.sum_subset subset_union_right fun x _ hx =>
      Multiset.count_eq_zero.2 (by simpa using hx)
  rw [Multiset.toFinset_sum_count_eq] at h2
  omega

/-- **Step 5 (count).** -/
theorem sum_le_natDegree_of_dvd (H : (ZMod p)[X]) (hH : H ≠ 0) (f : ZMod p → ℕ)
    (hf : ∀ b, (X - C b) ^ f b ∣ H) : ∑ b : ZMod p, f b ≤ H.natDegree := by
  classical
  calc ∑ b, f b ≤ ∑ b, H.roots.count b := Finset.sum_le_sum fun b _ => by
          rw [count_roots]; exact (le_rootMultiplicity_iff hH).2 (hf b)
    _ ≤ Multiset.card H.roots := sum_count_le_card _ _
    _ ≤ H.natDegree := card_roots' H

/-- **Step 3 (order at `b`).** -/
theorem hankel_order (A : Finset (ZMod p)) (d e : ℕ) (hd1 : 1 ≤ d) (hed : 2 * e ≤ d)
    (hA1 : 1 ≤ A.card) (b : ZMod p)
    (hb : ∀ a ∈ A, b + a = 0 ∨ (b + a) ^ d = 1 ∨ (b + a) ^ d = -1) :
    (X - C b) ^ (∑ i ∈ Finset.univ.filter (fun i : Fin (e + 1) => (badSet A d b).card ≤ (i : ℕ)),
        (A.card - (if -b ∈ A then 1 else 0) - ((i : ℕ) + e))) ∣ (hankel A d e).det := by
  have h := pow_dvd_det_add_mul
    (Matrix.of fun (i : Fin (e + 1)) (a : badSet A d b) =>
          C (2 * wt A a) * (X + C (a : ZMod p)) ^ (d + A.card - 1 - 2 * e) *
            (X + C (a : ZMod p)) ^ (e - i))
    (Matrix.of fun (a : badSet A d b) (j : Fin (e + 1)) => (X + C (a : ZMod p)) ^ (e - j))
    (Matrix.of fun i j : Fin (e + 1) => epsPoly A d b ((i : ℕ) + j)) (X - C b)
    (fun i : Fin (e + 1) => A.card - (if -b ∈ A then 1 else 0) - ((i : ℕ) + e))
    (fun i i' hii' => by
      have : (i : ℕ) ≤ i' := hii'
      simp only
      omega)
    (fun i j => by
      simp only [Matrix.of_apply]
      exact (pow_dvd_pow _ (by have := j.isLt; omega)).trans (epsPoly_dvd A d hd1 b hb _))
  rw [Fintype.card_coe] at h
  rw [hankel_eq_add_mul A d e hed hA1 b]
  exact h

theorem sum_Icc_eq_sum_fin (f : ℕ → ℕ) (k e : ℕ) :
    ∑ i ∈ Finset.Icc k e, f i =
      ∑ i ∈ Finset.univ.filter (fun i : Fin (e + 1) => k ≤ (i : ℕ)), f i := by
  rw [Finset.sum_filter, Fin.sum_univ_eq_sum_range (fun i => if k ≤ i then f i else 0) (e + 1),
    ← Finset.sum_filter]
  congr 1
  ext i
  simp only [Finset.mem_Icc, Finset.mem_filter, Finset.mem_range]
  omega

/-- **Theorem 2.1, general exponent, conditional on `Λ ≠ 0`.** -/
theorem hankel_stepanov_core (A : Finset (ZMod p)) (d e : ℕ) (hd1 : 1 ≤ d) (hed : 2 * e ≤ d)
    (hdx : ∀ x : ZMod p, x = 0 ∨ x ^ d = 1 ∨ x ^ d = -1)
    (hΛ : hankelLead p d A.card e ≠ 0) :
    ∑ b : ZMod p, ∑ i ∈ Finset.Icc (badSet A d b).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) ≤ (e + 1) * (d - e) := by
  rcases Nat.eq_zero_or_pos A.card with hA0 | hA1
  · simp [hA0]
  have hH : (hankel A d e).det ≠ 0 := by
    intro h0
    apply hΛ
    rw [← hankel_det_coeff A d e hd1 hed hA1, h0, coeff_zero]
  refine le_trans (le_of_eq ?_) ((sum_le_natDegree_of_dvd _ hH _ fun b =>
    hankel_order A d e hd1 hed hA1 b fun a _ => hdx (b + a)).trans
      (hankel_det_natDegree_le A d e hed hA1))
  refine Finset.sum_congr rfl fun b _ => ?_
  exact sum_Icc_eq_sum_fin (fun i => A.card - (if -b ∈ A then 1 else 0) - (i + e)) _ e

end Degree

section Paley

variable {p : ℕ} [Fact p.Prime]

theorem isSquare_zero_or_pow_half (hp : p ≠ 2) {x : ZMod p} (hx : IsSquare x) :
    x = 0 ∨ x ^ ((p - 1) / 2) = 1 := by
  obtain ⟨r, rfl⟩ := hx
  by_cases hr : r = 0
  · left; simp [hr]
  · right
    have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
    rw [← sq, ← pow_mul, show 2 * ((p - 1) / 2) = p - 1 by omega]
    exact ZMod.pow_card_sub_one_eq_one hr

theorem one_ne_neg_one_zmod (hp : p ≠ 2) : (1 : ZMod p) ≠ -1 := by
  intro h
  have h2 : ((2 : ℕ) : ZMod p) = 0 := by
    push_cast
    linear_combination h
  rw [ZMod.natCast_eq_zero_iff] at h2
  have := Nat.le_of_dvd (by norm_num) h2
  have := (Fact.out : p.Prime).two_le
  omega

theorem paley_dichotomy (hp : p ≠ 2) (x : ZMod p) :
    x = 0 ∨ x ^ ((p - 1) / 2) = 1 ∨ x ^ ((p - 1) / 2) = -1 := by
  by_cases hx : x = 0
  · exact Or.inl hx
  · right
    have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
    rw [show (p - 1) / 2 = p / 2 by omega]
    exact ZMod.pow_div_two_eq_neg_one_or_one (p := p) hx

theorem badSet_paley (hp : p ≠ 2) (A : Finset (ZMod p)) (b : ZMod p) :
    badSet A ((p - 1) / 2) b = A.filter (fun a => ¬ IsSquare (a + b)) := by
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  have hp3 : 3 ≤ p := by have := (Fact.out : p.Prime).two_le; omega
  ext a
  simp only [badSet, Finset.mem_filter]
  constructor
  · rintro ⟨ha, h⟩
    refine ⟨ha, fun hsq => ?_⟩
    rw [add_comm] at h
    rcases isSquare_zero_or_pow_half hp hsq with h0 | h1
    · rw [h0, zero_pow (by omega)] at h
      exact one_ne_zero (neg_eq_zero.1 h.symm)
    · rw [h] at h1
      exact one_ne_neg_one_zmod hp h1.symm
  · rintro ⟨ha, h⟩
    refine ⟨ha, ?_⟩
    have hne : a + b ≠ 0 := fun h0 => h (h0 ▸ ⟨0, by ring⟩)
    rw [add_comm, show (p - 1) / 2 = p / 2 by omega]
    rcases ZMod.pow_div_two_eq_neg_one_or_one (p := p) hne with h1 | h1
    · exact absurd ((ZMod.euler_criterion (p := p) hne).2 h1) h
    · exact h1

/-- **Theorem 2.1 (Paley, conditional form).** -/
theorem hankel_stepanov_paley_of_lead (hp : p ≠ 2) (A : Finset (ZMod p)) (e : ℕ)
    (hed : 2 * e ≤ (p - 1) / 2) (hΛ : hankelLead p ((p - 1) / 2) A.card e ≠ 0) :
    ∑ b : ZMod p, ∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) ≤ (e + 1) * ((p - 1) / 2 - e) := by
  have hp3 : 3 ≤ p := by
    have := (Fact.out : p.Prime).two_le
    omega
  have h := hankel_stepanov_core A ((p - 1) / 2) e (by omega) hed (paley_dichotomy hp) hΛ
  simp only [badSet_paley hp] at h
  exact h

end Paley

section Lead

variable {p : ℕ} [Fact p.Prime]

theorem factorial_ne_zero_zmod {n : ℕ} (hn : n < p) : ((n.factorial : ℕ) : ZMod p) ≠ 0 := by
  rw [Ne, CharP.cast_eq_zero_iff (ZMod p) p]
  intro h
  have := (Nat.Prime.dvd_factorial (Fact.out : p.Prime)).1 h
  omega

theorem descFactorial_ne_zero_zmod {n k : ℕ} (hk : k ≤ n) (hn : n < p) :
    ((n.descFactorial k : ℕ) : ZMod p) ≠ 0 := by
  intro h
  apply factorial_ne_zero_zmod hn
  rw [← Nat.factorial_mul_descFactorial hk, Nat.cast_mul, h, mul_zero]

/-- Factorisation of the Hankel entries: `C(D-i-j, K)·K!·(d-i)! = (D-i-e)!·(d-i)^{(j)}·(D-i-j)^{(e-j)}`
with falling factorials, `D = d + K`. -/
theorem choose_factor (d K e i j : ℕ) (hi : i ≤ e) (hj : j ≤ e) (hed : 2 * e ≤ d) :
    (d + K - (i + j)).choose K * K.factorial * (d - i).factorial =
      (d + K - i - e).factorial * ((d - i).descFactorial j * (d + K - i - j).descFactorial (e - j)) := by
  have h1 := Nat.choose_mul_factorial_mul_factorial (show K ≤ d + K - (i + j) by omega)
  have h2 := Nat.factorial_mul_descFactorial (show j ≤ d - i by omega)
  have h3 := Nat.factorial_mul_descFactorial (show e - j ≤ d + K - i - j by omega)
  rw [show d + K - (i + j) - K = d - i - j by omega, show d + K - (i + j) = d + K - i - j by omega]
    at h1
  rw [show d + K - i - j - (e - j) = d + K - i - e by omega] at h3
  rw [show d + K - (i + j) = d + K - i - j by omega, ← h2]
  calc (d + K - i - j).choose K * K.factorial * ((d - i - j).factorial * (d - i).descFactorial j)
      = ((d + K - i - j).choose K * K.factorial * (d - i - j).factorial) *
          (d - i).descFactorial j := by ring
    _ = (d + K - i - j).factorial * (d - i).descFactorial j := by rw [h1]
    _ = (d + K - i - e).factorial *
          ((d - i).descFactorial j * (d + K - i - j).descFactorial (e - j)) := by
        rw [← h3]; ring

/-- The column polynomials `P_j(x) = (d-x)^{(j)} (D-j-x)^{(e-j)}` (falling factorials). -/
noncomputable def leadPoly (p d K e j : ℕ) : (ZMod p)[X] :=
  (descPochhammer (ZMod p) j).comp (C (d : ZMod p) - X) *
    (descPochhammer (ZMod p) (e - j)).comp (C ((d + K - j : ℕ) : ZMod p) - X)

theorem natDegree_C_sub_X_le (a : ZMod p) : (C a - X).natDegree ≤ 1 := by
  refine (natDegree_sub_le _ _).trans ?_
  simp

theorem natDegree_leadPoly_le (d K e j : ℕ) (hj : j ≤ e) : (leadPoly p d K e j).natDegree ≤ e := by
  unfold leadPoly
  refine (natDegree_mul_le).trans ?_
  have h1 : ((descPochhammer (ZMod p) j).comp (C (d : ZMod p) - X)).natDegree ≤ j := by
    refine (natDegree_comp_le).trans ?_
    rw [descPochhammer_natDegree]
    calc j * (C (d : ZMod p) - X).natDegree ≤ j * 1 :=
          Nat.mul_le_mul_left _ (natDegree_C_sub_X_le _)
      _ = j := mul_one j
  have h2 : ((descPochhammer (ZMod p) (e - j)).comp
      (C ((d + K - j : ℕ) : ZMod p) - X)).natDegree ≤ e - j := by
    refine (natDegree_comp_le).trans ?_
    rw [descPochhammer_natDegree]
    calc (e - j) * (C ((d + K - j : ℕ) : ZMod p) - X).natDegree ≤ (e - j) * 1 :=
          Nat.mul_le_mul_left _ (natDegree_C_sub_X_le _)
      _ = e - j := mul_one _
  omega

theorem eval_leadPoly (d K e j : ℕ) (x : ZMod p) :
    (leadPoly p d K e j).eval x =
      (descPochhammer (ZMod p) j).eval ((d : ZMod p) - x) *
        (descPochhammer (ZMod p) (e - j)).eval (((d + K - j : ℕ) : ZMod p) - x) := by
  simp [leadPoly, eval_comp]

theorem eval_leadPoly_node (d K e i j : ℕ) (hi : i ≤ e) (hj : j ≤ e) (hed : 2 * e ≤ d) :
    (leadPoly p d K e j).eval ((i : ℕ) : ZMod p) =
      (((d - i).descFactorial j : ℕ) : ZMod p) *
        (((d + K - i - j).descFactorial (e - j) : ℕ) : ZMod p) := by
  rw [eval_leadPoly]
  have e1 : (d : ZMod p) - (i : ZMod p) = ((d - i : ℕ) : ZMod p) := by
    rw [Nat.cast_sub (by omega)]
  have e2 : ((d + K - j : ℕ) : ZMod p) - (i : ZMod p) = ((d + K - i - j : ℕ) : ZMod p) := by
    rw [← Nat.cast_sub (by omega), show d + K - j - i = d + K - i - j by omega]
  rw [e1, e2, descPochhammer_eval_eq_descFactorial, descPochhammer_eval_eq_descFactorial]

theorem eval_leadPoly_witness (d K e k j : ℕ) (hk : k ≤ e) (hj : j ≤ e) (hed : 2 * e ≤ d)
    (heK : e ≤ K) :
    (leadPoly p d K e j).eval ((d - k : ℕ) : ZMod p) =
      ((k.descFactorial j : ℕ) : ZMod p) * (((K + k - j).descFactorial (e - j) : ℕ) : ZMod p) := by
  rw [eval_leadPoly]
  have e1 : (d : ZMod p) - ((d - k : ℕ) : ZMod p) = ((k : ℕ) : ZMod p) := by
    rw [Nat.cast_sub (by omega)]; ring
  have e2 : ((d + K - j : ℕ) : ZMod p) - ((d - k : ℕ) : ZMod p) = ((K + k - j : ℕ) : ZMod p) := by
    rw [← Nat.cast_sub (by omega), show d + K - j - (d - k) = K + k - j by omega]
  rw [e1, e2, descPochhammer_eval_eq_descFactorial, descPochhammer_eval_eq_descFactorial]

/-- A matrix `[P_j(x_i)]` at `e+1` distinct nodes is nonsingular when the `P_j` have degree `≤ e`
and admit a triangular family of witness points. -/
theorem det_eval_ne_zero (e : ℕ) (P : Fin (e + 1) → (ZMod p)[X])
    (hdeg : ∀ j, (P j).natDegree ≤ e) (hep : e < p) (z : Fin (e + 1) → ZMod p)
    (hz : ∀ k j : Fin (e + 1), (k : ℕ) < j → (P j).eval (z k) = 0)
    (hzk : ∀ k, (P k).eval (z k) ≠ 0) :
    (Matrix.of fun i j : Fin (e + 1) => (P j).eval ((i : ℕ) : ZMod p)).det ≠ 0 := by
  classical
  intro hdet
  obtain ⟨c, hc0, hc⟩ := Matrix.exists_mulVec_eq_zero_iff.2 hdet
  set Qp : (ZMod p)[X] := ∑ j, C (c j) * P j with hQp
  have hinj : Function.Injective (fun i : Fin (e + 1) => ((i : ℕ) : ZMod p)) := by
    intro i i' h
    simp only at h
    rw [ZMod.natCast_eq_natCast_iff', Nat.mod_eq_of_lt (by omega),
      Nat.mod_eq_of_lt (by omega)] at h
    exact Fin.ext h
  have hQ : Qp = 0 := by
    apply eq_zero_of_natDegree_lt_card_of_eval_eq_zero Qp hinj
    · intro i
      have h := congrFun hc i
      simp only [Matrix.mulVec, dotProduct, Matrix.of_apply, Pi.zero_apply] at h
      rw [hQp, eval_finset_sum]
      rw [← h]
      refine Finset.sum_congr rfl fun j _ => ?_
      rw [eval_mul, eval_C, mul_comm]
    · rw [Fintype.card_fin]
      refine Nat.lt_succ_of_le (natDegree_sum_le_of_forall_le _ _ fun j _ => ?_)
      exact (natDegree_C_mul_le _ _).trans (hdeg j)
  apply hc0
  have hall : ∀ n : ℕ, ∀ k : Fin (e + 1), (k : ℕ) = n → c k = 0 := by
    intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      intro k hk
      have hev : Qp.eval (z k) = 0 := by rw [hQ, eval_zero]
      rw [hQp, eval_finset_sum] at hev
      simp only [eval_mul, eval_C] at hev
      rw [Finset.sum_eq_single k] at hev
      · exact (mul_eq_zero.1 hev).resolve_right (hzk k)
      · intro j _ hjk
        have hne : (j : ℕ) ≠ k := fun h => hjk (Fin.ext h)
        rcases Nat.lt_or_gt_of_ne hne with hlt | hgt
        · rw [ih j (by omega) j rfl, zero_mul]
        · rw [hz k j hgt, mul_zero]
      · intro h; exact absurd (Finset.mem_univ k) h
  funext k
  exact hall k k rfl

/-- **Lemma 2.2 (non-vanishing of the leading coefficient), PROVED without Krattenthaler.**
For `1 ≤ d`, `2e ≤ d`, `2e + 1 ≤ m` and `D = d + m - 1 < p`,
`Λ = det[C(D - i - j, m - 1)]_{0 ≤ i,j ≤ e} ≠ 0` in `ZMod p`. -/
theorem hankelLead_ne_zero (d m e : ℕ) (hed : 2 * e ≤ d) (hem : 2 * e + 1 ≤ m)
    (hDp : d + m - 1 < p) : hankelLead p d m e ≠ 0 := by
  set K := m - 1 with hK
  have hm : d + m - 1 = d + K := by omega
  let α : Fin (e + 1) → ZMod p := fun i =>
    (((d + K - i - e).factorial : ℕ) : ZMod p) *
      ((((K.factorial : ℕ) : ZMod p) * (((d - i).factorial : ℕ) : ZMod p)))⁻¹
  have hα : ∀ i, α i ≠ 0 := by
    intro i
    have hi := i.isLt
    refine mul_ne_zero (factorial_ne_zero_zmod (by omega)) (inv_ne_zero (mul_ne_zero
      (factorial_ne_zero_zmod (by omega)) (factorial_ne_zero_zmod (by omega))))
  have hentry : ∀ i j : Fin (e + 1),
      (((d + m - 1 - ((i : ℕ) + j)).choose (m - 1) : ℕ) : ZMod p) =
        α i * (leadPoly p d K e j).eval ((i : ℕ) : ZMod p) := by
    intro i j
    have hi := i.isLt
    have hj := j.isLt
    rw [hm, ← hK, eval_leadPoly_node d K e i j (by omega) (by omega) hed]
    have hcast := congrArg (Nat.cast : ℕ → ZMod p)
      (choose_factor d K e i j (by omega) (by omega) hed)
    simp only [Nat.cast_mul] at hcast
    have hu : (((K.factorial : ℕ) : ZMod p) * (((d - i).factorial : ℕ) : ZMod p)) ≠ 0 :=
      mul_ne_zero (factorial_ne_zero_zmod (by omega)) (factorial_ne_zero_zmod (by omega))
    simp only [α]
    rw [show ∀ a b c : ZMod p, a * b⁻¹ * c = (a * c) * b⁻¹ from fun a b c => by ring,
      eq_mul_inv_iff_mul_eq₀ hu]
    linear_combination hcast
  have hmat : (Matrix.of fun i j : Fin (e + 1) =>
      (((d + m - 1 - ((i : ℕ) + j)).choose (m - 1) : ℕ) : ZMod p)) =
      Matrix.of fun i j : Fin (e + 1) =>
        α i * (leadPoly p d K e j).eval ((i : ℕ) : ZMod p) := by
    ext i j
    exact hentry i j
  unfold hankelLead
  rw [hmat, Matrix.det_mul_column]
  refine mul_ne_zero (Finset.prod_ne_zero_iff.2 fun i _ => hα i) ?_
  refine det_eval_ne_zero e (fun j => leadPoly p d K e j)
    (fun j => natDegree_leadPoly_le d K e j (by have := j.isLt; omega)) (by omega)
    (fun k => ((d - k : ℕ) : ZMod p)) ?_ ?_
  · intro k j hkj
    have hk := k.isLt
    have hj := j.isLt
    simp only
    rw [eval_leadPoly_witness d K e k j (by omega) (by omega) hed (by omega),
      (Nat.descFactorial_eq_zero_iff_lt).2 hkj, Nat.cast_zero, zero_mul]
  · intro k
    have hk := k.isLt
    simp only
    rw [eval_leadPoly_witness d K e k k (by omega) (by omega) hed (by omega),
      Nat.descFactorial_self, show K + (k : ℕ) - k = K by omega]
    exact mul_ne_zero (factorial_ne_zero_zmod (by omega))
      (descFactorial_ne_zero_zmod (by omega) (by omega))

end Lead

section Main

variable {p : ℕ} [Fact p.Prime]

/-- **Theorem 2.1 of `stepanov-robust-2026-09-26.md` (Paley case), unconditional.**
For `A ⊆ F_p` with `|A| ≤ (p+1)/2` and `2e + 1 ≤ |A|`, with `e_b = #{a ∈ A : a + b non-square}`
and `δ_b = [-b ∈ A]`:
`Σ_b Σ_{i = e_b}^{e} (|A| - δ_b - (i + e)) ≤ (e + 1)(d - e)`, `d = (p-1)/2`. -/
theorem hankel_stepanov_paley (hp : p ≠ 2) (A : Finset (ZMod p))
    (hA : A.card ≤ (p + 1) / 2) (e : ℕ) (he : 2 * e + 1 ≤ A.card) :
    ∑ b : ZMod p, ∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) ≤ (e + 1) * ((p - 1) / 2 - e) := by
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  have hp3 : 3 ≤ p := by
    have := (Fact.out : p.Prime).two_le
    omega
  exact hankel_stepanov_paley_of_lead hp A e (by omega)
    (hankelLead_ne_zero _ _ e (by omega) he (by omega))

end Main

section MainExtras

variable {p : ℕ} [Fact p.Prime]

/-- **Theorem 2.1, degree statement:** `deg H_{e+1} = (e+1)(d-e)` exactly (Paley `d`). -/
theorem hankel_det_natDegree_eq (hp : p ≠ 2) (A : Finset (ZMod p))
    (hA : A.card ≤ (p + 1) / 2) (e : ℕ) (he : 2 * e + 1 ≤ A.card) :
    (hankel A ((p - 1) / 2) e).det.natDegree = (e + 1) * ((p - 1) / 2 - e) := by
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  have hp3 : 3 ≤ p := by
    have := (Fact.out : p.Prime).two_le
    omega
  refine le_antisymm (hankel_det_natDegree_le A _ e (by omega) (by omega)) ?_
  refine le_natDegree_of_ne_zero ?_
  rw [hankel_det_coeff A _ e (by omega) (by omega) (by omega)]
  exact hankelLead_ne_zero _ _ e (by omega) he (by omega)

/-- **Theorem 2.1, order statement (Paley):** `(X - b)^{Σ_{i=e_b}^{e} (m - δ_b - (i+e))}` divides
`H_{e+1} = det[u_{i+j}]`. -/
theorem hankel_order_paley (hp : p ≠ 2) (A : Finset (ZMod p)) (hA1 : 1 ≤ A.card) (e : ℕ)
    (hed : 2 * e ≤ (p - 1) / 2) (b : ZMod p) :
    (X - C b) ^ (∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e))) ∣ (hankel A ((p - 1) / 2) e).det := by
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  have hp3 : 3 ≤ p := by
    have := (Fact.out : p.Prime).two_le
    omega
  rw [← badSet_paley hp,
    sum_Icc_eq_sum_fin (fun i => A.card - (if -b ∈ A then 1 else 0) - (i + e))]
  exact hankel_order A _ e (by omega) hed hA1 b (fun a _ => paley_dichotomy hp (b + a))

end MainExtras

section Bias

variable {p : ℕ} [Fact p.Prime]

/-- Pointwise lower bound for the inner sum of (★). -/
theorem pt_bound (m e k δ : ℕ) (hme : 2 * e + 1 ≤ m) :
    ((m : ℤ) - 2 * e) * ((e + 1 : ℤ) - k) - (e + 1) * δ ≤
      ((∑ i ∈ Finset.Icc k e, (m - δ - (i + e)) : ℕ) : ℤ) := by
  rw [Nat.cast_sum]
  by_cases hk : k ≤ e
  · have h1 : ∀ i ∈ Finset.Icc k e, ((m : ℤ) - 2 * e - δ) ≤ ((m - δ - (i + e) : ℕ) : ℤ) := by
      intro i hi
      rw [Finset.mem_Icc] at hi
      omega
    have h2 := Finset.card_nsmul_le_sum (Finset.Icc k e) (fun i => ((m - δ - (i + e) : ℕ) : ℤ)) _ h1
    rw [Nat.card_Icc, nsmul_eq_mul, Nat.cast_sub (by omega)] at h2
    push_cast at h2
    have hk0 : (0 : ℤ) ≤ k := by positivity
    have hδ0 : (0 : ℤ) ≤ δ := by positivity
    nlinarith [mul_nonneg hk0 hδ0]
  · rw [Finset.Icc_eq_empty hk, Finset.sum_empty]
    have h1 : (1 : ℤ) ≤ (m : ℤ) - 2 * e := by omega
    have h2 : (e + 1 : ℤ) - k ≤ 0 := by omega
    have hδ0 : (0 : ℤ) ≤ δ := by positivity
    have h3 := mul_le_mul_of_nonneg_left h2 (by linarith : (0 : ℤ) ≤ (m : ℤ) - 2 * e)
    have h4 := mul_nonneg (by positivity : (0 : ℤ) ≤ (e : ℤ) + 1) hδ0
    linarith

/-- **Corollary 2.3 (integer form).** With `n = |B|`, `N_- = Σ_{b∈B} e_b` and
`r = #{b ∈ B : -b ∈ A}`: `(m - 2e)((e+1)n - N_-) ≤ (e+1)(d - e + r)`. -/
theorem bias_inequality (hp : p ≠ 2) (A : Finset (ZMod p)) (hA : A.card ≤ (p + 1) / 2) (e : ℕ)
    (he : 2 * e + 1 ≤ A.card) (B : Finset (ZMod p)) :
    ((A.card : ℤ) - 2 * e) * ((e + 1) * B.card -
        ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ)) ≤
      (e + 1) * ((((p - 1) / 2 : ℕ) : ℤ) - e + (B.filter fun b => -b ∈ A).card) := by
  classical
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  have hed : e ≤ (p - 1) / 2 := by omega
  have hstar := hankel_stepanov_paley hp A hA e he
  have hB := (Finset.sum_le_sum_of_subset (f := fun b : ZMod p =>
      ∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e))) (Finset.subset_univ B)).trans hstar
  have hBz : ((∑ b ∈ B, ∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) : ℕ) : ℤ) ≤
      ((e + 1 : ℕ) : ℤ) * ((((p - 1) / 2 : ℕ) : ℤ) - e) := by
    rw [← Nat.cast_sub hed, ← Nat.cast_mul]
    exact_mod_cast hB
  rw [Nat.cast_sum] at hBz
  have hpt : ∀ b ∈ B, ((A.card : ℤ) - 2 * e) *
      ((e + 1 : ℤ) - (A.filter fun a => ¬ IsSquare (a + b)).card) -
        (e + 1) * ((if -b ∈ A then 1 else 0 : ℕ) : ℤ) ≤
      ((∑ i ∈ Finset.Icc (A.filter (fun a => ¬ IsSquare (a + b))).card e,
        (A.card - (if -b ∈ A then 1 else 0) - (i + e)) : ℕ) : ℤ) := by
    intro b _
    exact pt_bound _ _ _ _ he
  have hsum := Finset.sum_le_sum hpt
  have hL : ∑ b ∈ B, (((A.card : ℤ) - 2 * e) *
      ((e + 1 : ℤ) - (A.filter fun a => ¬ IsSquare (a + b)).card) -
        (e + 1) * ((if -b ∈ A then 1 else 0 : ℕ) : ℤ)) =
      ((A.card : ℤ) - 2 * e) * ((e + 1) * B.card -
        ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ)) -
      (e + 1) * ((B.filter fun b => -b ∈ A).card : ℤ) := by
    rw [Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_sub_distrib,
      Finset.sum_const, Finset.card_filter, Nat.cast_sum]
    simp only [nsmul_eq_mul]
    ring
  rw [hL] at hsum
  push_cast at hBz hsum ⊢
  linarith

/-- The two-set character sum `S(A,B) = Σ_{a∈A} Σ_{b∈B} χ(a + b)`, `χ` the quadratic character
of `F_p` (with `χ(0) = 0`). -/
def charSum (A B : Finset (ZMod p)) : ℤ :=
  ∑ a ∈ A, ∑ b ∈ B, quadraticChar (ZMod p) (a + b)

theorem quadraticChar_eq_ite (x : ZMod p) :
    quadraticChar (ZMod p) x =
      1 - (if x = 0 then 1 else 0) - 2 * (if IsSquare x then 0 else 1) := by
  by_cases hx : x = 0
  · subst hx
    simp
  · by_cases hs : IsSquare x
    · rw [(quadraticChar_one_iff_isSquare hx).2 hs, if_neg hx, if_pos hs]
      norm_num
    · rw [quadraticChar_neg_one_iff_not_isSquare.2 hs, if_neg hx, if_neg hs]
      norm_num

/-- `S(A,B) = mn - r - 2 N_-`. -/
theorem charSum_eq (A B : Finset (ZMod p)) :
    charSum A B = (A.card : ℤ) * B.card - (B.filter fun b => -b ∈ A).card -
      2 * ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ) := by
  classical
  unfold charSum
  rw [Finset.sum_comm]
  have hb : ∀ b : ZMod p, ∑ a ∈ A, quadraticChar (ZMod p) (a + b) =
      (A.card : ℤ) - (if -b ∈ A then 1 else 0) -
        2 * ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ) := by
    intro b
    simp_rw [quadraticChar_eq_ite]
    rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, Finset.sum_const, ← Finset.mul_sum]
    have h1 : ∑ a ∈ A, (if a + b = 0 then (1 : ℤ) else 0) = if -b ∈ A then 1 else 0 := by
      have : ∀ a, (if a + b = 0 then (1 : ℤ) else 0) = if -b = a then 1 else 0 := by
        intro a
        by_cases h : a + b = 0
        · rw [if_pos h, if_pos (by linear_combination -h)]
        · rw [if_neg h, if_neg (fun h' => h (by rw [← h']; ring))]
      rw [Finset.sum_congr rfl fun a _ => this a, Finset.sum_ite_eq]
    have h2 : ∑ a ∈ A, (if IsSquare (a + b) then (0 : ℤ) else 1) =
        ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ) := by
      rw [Finset.card_filter, Nat.cast_sum]
      refine Finset.sum_congr rfl fun a _ => ?_
      split_ifs <;> simp_all
    rw [h1, h2, nsmul_eq_mul, mul_one]
  rw [Finset.sum_congr rfl fun b _ => hb b, Finset.sum_sub_distrib, Finset.sum_sub_distrib,
    Finset.sum_const, ← Finset.mul_sum, Finset.card_filter, Nat.cast_sum, nsmul_eq_mul]
  push_cast
  ring

/-- **Corollary 2.3 (bias form).** -/
theorem bias_bound (hp : p ≠ 2) (A : Finset (ZMod p)) (hA : A.card ≤ (p + 1) / 2) (e : ℕ)
    (he : 2 * e + 1 ≤ A.card) (B : Finset (ZMod p)) :
    (charSum A B : ℝ) ≤ (A.card : ℝ) * B.card - (B.filter fun b => -b ∈ A).card -
      2 * (e + 1) * (B.card - ((((p - 1) / 2 : ℕ) : ℝ) - e + (B.filter fun b => -b ∈ A).card) /
        ((A.card : ℝ) - 2 * e)) := by
  have h1 := bias_inequality hp A hA e he B
  have h2 := charSum_eq A B
  have hpos : (0 : ℝ) < (A.card : ℝ) - 2 * e := by
    have : (2 * e + 1 : ℝ) ≤ A.card := by exact_mod_cast he
    linarith
  have h1' : ((A.card : ℝ) - 2 * e) * ((e + 1) * B.card -
        ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℝ)) ≤
      (e + 1) * ((((p - 1) / 2 : ℕ) : ℝ) - e + (B.filter fun b => -b ∈ A).card) := by
    exact_mod_cast h1
  have h2' : (charSum A B : ℝ) = (A.card : ℝ) * B.card - (B.filter fun b => -b ∈ A).card -
      2 * ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℝ) := by
    exact_mod_cast h2
  set N := ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℝ)
  set Rr := ((((p - 1) / 2 : ℕ) : ℝ) - e + (B.filter fun b => -b ∈ A).card)
  have h3 : (e + 1 : ℝ) * B.card - N ≤ (e + 1) * Rr / ((A.card : ℝ) - 2 * e) := by
    rw [le_div_iff₀ hpos]; linarith
  have h4 : (e + 1 : ℝ) * Rr / ((A.card : ℝ) - 2 * e) =
      (e + 1) * (Rr / ((A.card : ℝ) - 2 * e)) := mul_div_assoc _ _ _
  rw [h2']
  linarith

end Bias

section ConstantBias

variable {p : ℕ} [Fact p.Prime]

theorem card_powersetCard_mem {α : Type*} [DecidableEq α] (A : Finset α) (k : ℕ) (a : α)
    (ha : a ∈ A) :
    ((A.powersetCard (k + 1)).filter (fun T => a ∈ T)).card = (A.card - 1).choose k := by
  have hA : A = insert a (A.erase a) := (Finset.insert_erase ha).symm
  have hnot : a ∉ A.erase a := Finset.notMem_erase a A
  have h1 : ((A.erase a).powersetCard (k + 1)).filter (fun T => a ∈ T) = ∅ := by
    rw [Finset.filter_eq_empty_iff]
    intro T hT haT
    exact hnot ((Finset.mem_powersetCard.1 hT).1 haT)
  have h2 : (((A.erase a).powersetCard k).image (insert a)).filter (fun T => a ∈ T) =
      ((A.erase a).powersetCard k).image (insert a) := by
    rw [Finset.filter_eq_self]
    intro T hT
    obtain ⟨S, _, rfl⟩ := Finset.mem_image.1 hT
    exact Finset.mem_insert_self a S
  have hinj : Set.InjOn (insert a) (((A.erase a).powersetCard k : Finset (Finset α)) :
      Set (Finset α)) := by
    intro S hS S' hS' h
    have hS1 := (Finset.mem_powersetCard.1 hS).1
    have hS'1 := (Finset.mem_powersetCard.1 hS').1
    have haS : a ∉ S := fun h' => hnot (hS1 h')
    have haS' : a ∉ S' := fun h' => hnot (hS'1 h')
    rw [← Finset.erase_insert haS, ← Finset.erase_insert haS']
    exact congrArg (fun T => Finset.erase T a) h
  conv_lhs => rw [hA, Finset.powersetCard_succ_insert hnot, Finset.filter_union, h1, h2,
    Finset.empty_union]
  rw [Finset.card_image_of_injOn hinj, Finset.card_powersetCard, Finset.card_erase_of_mem ha]

/-- Double counting: `Σ_{T ⊆ A, |T| = k+1} Σ_{a∈T} g a = C(|A|-1, k) Σ_{a∈A} g a`. -/
theorem sum_powersetCard_sum {α : Type*} [DecidableEq α] (A : Finset α) (k : ℕ) (g : α → ℝ) :
    ∑ T ∈ A.powersetCard (k + 1), ∑ a ∈ T, g a =
      ((A.card - 1).choose k : ℝ) * ∑ a ∈ A, g a := by
  have h1 : ∀ T ∈ A.powersetCard (k + 1), ∑ a ∈ T, g a = ∑ a ∈ A, if a ∈ T then g a else 0 := by
    intro T hT
    rw [← Finset.sum_filter, Finset.filter_mem_eq_inter,
      Finset.inter_eq_right.2 (Finset.mem_powersetCard.1 hT).1]
  rw [Finset.sum_congr rfl h1, Finset.sum_comm, Finset.mul_sum]
  refine Finset.sum_congr rfl fun a ha => ?_
  rw [← Finset.sum_filter, Finset.sum_const, card_powersetCard_mem A k a ha, nsmul_eq_mul]

theorem card_filter_neg_mem_le (A B : Finset (ZMod p)) :
    (B.filter fun b => -b ∈ A).card ≤ A.card := by
  refine Finset.card_le_card_of_injOn (fun b => -b) ?_ ?_
  · intro b hb
    exact (Finset.mem_filter.1 hb).2
  · intro b _ b' _ h
    simpa using h

theorem charSum_le_card_mul (A B : Finset (ZMod p)) :
    (charSum A B : ℝ) ≤ (A.card : ℝ) * B.card := by
  have h := charSum_eq A B
  have h1 : (0 : ℤ) ≤ ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℤ) :=
    Finset.sum_nonneg fun b _ => by positivity
  have h2 : (0 : ℤ) ≤ ((B.filter fun b => -b ∈ A).card : ℤ) := by positivity
  have h3 : charSum A B ≤ (A.card : ℤ) * B.card := by linarith
  exact_mod_cast h3

-- Root fix 2026-09-29: this proof exceeds the default 200000 heartbeats under Mathlib 5e932f9
-- (Lean v4.29.1); raising the limit changes no logic.
set_option maxHeartbeats 1000000 in
/-- The per-subset bound behind Theorem 2.5: if `1 ≤ |A| ≤ (p+1)/2` and `|A||B| ≥ (1/2+κ)p`,
`0 < κ ≤ 3/2`, `u = (1+2κ)^{-1/2}`, then
`S(A,B) ≤ [1 - (1-u)² + u(1-u)|A|/d]·|A||B|`. -/
theorem subset_bound (hp : p ≠ 2) (κ : ℝ) (hκ0 : 0 < κ) (hκ1 : κ ≤ 3 / 2)
    (A B : Finset (ZMod p)) (hA1 : 1 ≤ A.card) (hA : A.card ≤ (p + 1) / 2)
    (hAB : (1 / 2 + κ) * p ≤ (A.card : ℝ) * B.card) :
    (charSum A B : ℝ) ≤
      (1 - (1 - 1 / Real.sqrt (1 + 2 * κ)) ^ 2 +
        (1 / Real.sqrt (1 + 2 * κ)) * (1 - 1 / Real.sqrt (1 + 2 * κ)) * A.card /
          (((p : ℝ) - 1) / 2)) * ((A.card : ℝ) * B.card) := by
  classical
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  have hp3 : 3 ≤ p := by have := (Fact.out : p.Prime).two_le; omega
  set s := Real.sqrt (1 + 2 * κ) with hs
  have hs2 : s ^ 2 = 1 + 2 * κ := Real.sq_sqrt (by linarith)
  have hs0 : 0 ≤ s := Real.sqrt_nonneg _
  have hs1 : 1 < s := by nlinarith
  have hs4 : s ≤ 2 := by nlinarith
  set u := 1 / s with hu
  have hu0 : 0 < u := by positivity
  have hus : u * s = 1 := by rw [hu, one_div, inv_mul_cancel₀ (by linarith : (0 : ℝ) < s).ne']
  have hu1 : u < 1 := by nlinarith
  have hu2 : 1 / 2 ≤ u := by nlinarith
  have hu2k : u ^ 2 * (1 + 2 * κ) = 1 := by rw [← hs2, ← mul_pow, hus, one_pow]
  set M : ℝ := (A.card : ℝ) with hM
  set n : ℝ := (B.card : ℝ) with hn
  have hM1 : (1 : ℝ) ≤ M := by rw [hM]; exact_mod_cast hA1
  have hn0 : (0 : ℝ) ≤ n := by positivity
  set dR : ℝ := ((p : ℝ) - 1) / 2 with hdR
  have hdnat : (((p - 1) / 2 : ℕ) : ℝ) = dR := by
    rw [hdR]
    have : 2 * ((p - 1) / 2) + 1 = p := by omega
    have h' : (2 : ℝ) * (((p - 1) / 2 : ℕ) : ℝ) + 1 = p := by exact_mod_cast this
    linarith
  have hp3r : (3 : ℝ) ≤ p := by exact_mod_cast hp3
  have hdR0 : 0 < dR := by rw [hdR]; linarith
  -- `u² M n ≥ d`
  have hkey : dR ≤ u ^ 2 * (M * n) := by
    have h1 : (1 + 2 * κ) * dR ≤ M * n := by rw [hdR]; linarith
    have h2 : u ^ 2 * ((1 + 2 * κ) * dR) ≤ u ^ 2 * (M * n) :=
      mul_le_mul_of_nonneg_left h1 (sq_nonneg u)
    have h3 : u ^ 2 * ((1 + 2 * κ) * dR) = dR := by rw [← mul_assoc, hu2k, one_mul]
    linarith
  -- the exponent `e = ⌊θ M⌋`, `θ = (1-u)/2`
  set θ : ℝ := (1 - u) / 2 with hθ
  have hθ0 : 0 ≤ θ := by rw [hθ]; linarith
  set e : ℕ := ⌊θ * M⌋₊ with he
  have he1 : (e : ℝ) ≤ θ * M := Nat.floor_le (mul_nonneg hθ0 (by linarith))
  have he2 : θ * M < e + 1 := Nat.lt_floor_add_one _
  have hθ4 : θ ≤ 1 / 4 := by rw [hθ]; linarith
  have hem : 2 * e + 1 ≤ A.card := by
    have h0 := mul_le_mul_of_nonneg_right hθ4 (by linarith : (0 : ℝ) ≤ M)
    have h1 : (2 * e : ℝ) < M := by linarith
    rw [hM] at h1
    have h2 : 2 * e < A.card := by exact_mod_cast h1
    omega
  have hmain := bias_inequality hp A hA e hem B
  have hS := charSum_eq A B
  have hr := card_filter_neg_mem_le A B
  set N : ℝ := ∑ b ∈ B, ((A.filter fun a => ¬ IsSquare (a + b)).card : ℝ) with hN
  set R : ℝ := ((B.filter fun b => -b ∈ A).card : ℝ) with hR
  have hmain' : (M - 2 * e) * ((e + 1) * n - N) ≤ (e + 1) * (dR - e + R) := by
    rw [← hdnat, hM, hn, hN, hR]
    exact_mod_cast hmain
  have hS' : (charSum A B : ℝ) = M * n - R - 2 * N := by
    rw [hM, hn, hN, hR]
    exact_mod_cast hS
  have hRM : R ≤ M := by rw [hR, hM]; exact_mod_cast hr
  have hR0 : 0 ≤ R := by positivity
  have hN0 : 0 ≤ N := Finset.sum_nonneg fun b _ => by positivity
  have he0 : (0 : ℝ) ≤ e := by positivity
  have huM : u * M ≤ M - 2 * e := by rw [hθ] at he1; linarith
  -- `u M X ≤ (e+1)(d+M)` with `X = (e+1)n - N`
  have hX : u * M * ((e + 1) * n - N) ≤ (e + 1) * (dR + M) := by
    by_cases hXs : 0 ≤ (e + 1) * n - N
    · calc u * M * ((e + 1) * n - N) ≤ (M - 2 * e) * ((e + 1) * n - N) :=
            mul_le_mul_of_nonneg_right huM hXs
        _ ≤ (e + 1) * (dR - e + R) := hmain'
        _ ≤ (e + 1) * (dR + M) :=
            mul_le_mul_of_nonneg_left (by linarith) (by positivity)
    · rw [not_le] at hXs
      have : u * M * ((e + 1) * n - N) ≤ 0 :=
        mul_nonpos_of_nonneg_of_nonpos (by positivity) hXs.le
      have : 0 ≤ (e + 1) * (dR + M) := by positivity
      linarith
  rw [hS']
  have hc : (1 - (1 - u) ^ 2 + u * (1 - u) * M / dR) * (M * n) =
      M * n - (1 - u) ^ 2 * (M * n) + u * (1 - u) * M * (M * n) / dR := by ring
  rw [hc]
  by_cases hcase : dR + M ≤ u * M * n
  · -- Case A: `N ≥ θ(Mn - (d+M)/u)`, i.e. `u M · 2N ≥ (1-u) M (uMn - d - M)`
    have hN2 : (1 - u) * M * (u * M * n - dR - M) ≤ u * M * (2 * N) := by
      have h1 : θ * M * (u * M * n - dR - M) ≤ (e + 1) * (u * M * n - dR - M) :=
        mul_le_mul_of_nonneg_right he2.le (by linarith)
      rw [hθ] at h1
      linarith
    -- goal: `Mn - R - 2N ≤ Mn - (1-u)²Mn + u(1-u)M·Mn/d`
    have hgoal : u * M * dR * (M * n - 2 * N) ≤
        u * M * dR * (M * n - (1 - u) ^ 2 * (M * n)) + u * M * (u * (1 - u) * M * (M * n)) := by
      have hid : u * M * dR * (M * n - (1 - u) ^ 2 * (M * n)) +
          u * M * (u * (1 - u) * M * (M * n)) - u * M * dR * (M * n) +
          dR * ((1 - u) * M * (u * M * n - dR - M)) =
          M * (1 - u) * (dR + M) * (u ^ 2 * (M * n) - dR) := by ring
      have hpos : 0 ≤ M * (1 - u) * (dR + M) * (u ^ 2 * (M * n) - dR) :=
        mul_nonneg (mul_nonneg (mul_nonneg (by linarith) (by linarith)) (by linarith))
          (by linarith)
      linarith [mul_le_mul_of_nonneg_left hN2 hdR0.le]
    have huMd : 0 < u * M * dR := mul_pos (mul_pos hu0 (by linarith)) hdR0
    have hsplit : u * M * dR * (M * n - (1 - u) ^ 2 * (M * n) + u * (1 - u) * M * (M * n) / dR) =
        u * M * dR * (M * n - (1 - u) ^ 2 * (M * n)) + u * M * (u * (1 - u) * M * (M * n)) := by
      rw [mul_add]
      congr 1
      rw [mul_div_assoc', div_eq_iff hdR0.ne']
      ring
    have h3 : M * n - 2 * N ≤ M * n - (1 - u) ^ 2 * (M * n) + u * (1 - u) * M * (M * n) / dR :=
      le_of_mul_le_mul_left (by rw [hsplit]; exact hgoal) huMd
    linarith
  · -- Case B: the constant is `≥ 1`
    rw [not_le] at hcase
    have h1 : (1 - u) * dR ≤ u * M := by
      have := mul_lt_mul_of_pos_left hcase hu0
      linarith
    have h2 : (1 - u) ^ 2 * (M * n) ≤ u * (1 - u) * M * (M * n) / dR := by
      rw [le_div_iff₀ hdR0]
      have h0 : 0 ≤ (1 - u) * (M * n) := mul_nonneg (by linarith) (mul_nonneg (by linarith) hn0)
      linarith [mul_le_mul_of_nonneg_left h1 h0]
    linarith

/-- Upper bound `S(A,B) ≤ c·|A||B|` when `|A| ≤ |B|`. -/
theorem upper_bound (hp : p ≠ 2) (hp11 : 11 ≤ p) (κ : ℝ) (hκ0 : 0 < κ) (hκ1 : κ ≤ 3 / 2)
    (A B : Finset (ZMod p)) (hmn : A.card ≤ B.card)
    (hAB : (1 / 2 + κ) * p ≤ (A.card : ℝ) * B.card) :
    (charSum A B : ℝ) ≤
      (1 - (1 - 1 / Real.sqrt (1 + 2 * κ)) ^ 2 +
        (Real.sqrt ((1 / 2 + κ) * p) + 1) / (2 * ((p : ℝ) - 1))) * ((A.card : ℝ) * B.card) := by
  classical
  have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
  set s := Real.sqrt (1 + 2 * κ) with hs
  have hs2 : s ^ 2 = 1 + 2 * κ := Real.sq_sqrt (by linarith)
  have hs0 : 0 ≤ s := Real.sqrt_nonneg _
  have hs1 : 1 < s := by nlinarith
  set u := 1 / s with hu
  have hu0 : 0 < u := by positivity
  have hus : u * s = 1 := by rw [hu, one_div, inv_mul_cancel₀ (by linarith : (0 : ℝ) < s).ne']
  have hu1 : u < 1 := by nlinarith
  set Q : ℝ := (1 / 2 + κ) * p with hQ
  have hp11r : (11 : ℝ) ≤ p := by exact_mod_cast hp11
  have hQ0 : 0 < Q := by rw [hQ]; positivity
  set sQ := Real.sqrt Q with hsQ
  have hsQ2 : sQ ^ 2 = Q := Real.sq_sqrt hQ0.le
  have hsQ0 : 0 < sQ := Real.sqrt_pos.2 hQ0
  set M : ℝ := (A.card : ℝ) with hM
  set n : ℝ := (B.card : ℝ) with hn
  have hMn : M ≤ n := by rw [hM, hn]; exact_mod_cast hmn
  have hMnQ : Q ≤ M * n := hAB
  have hn00 : (0 : ℝ) ≤ n := by rw [hn]; positivity
  have hM0 : 0 < M := by
    by_contra h
    rw [not_lt] at h
    have : M * n ≤ 0 := by linarith [mul_nonneg (neg_nonneg.2 h) hn00]
    linarith
  have hn0 : 0 < n := lt_of_lt_of_le hM0 hMn
  have hnsQ : sQ ≤ n := by
    rw [hsQ, Real.sqrt_le_left hn0.le]
    linarith [mul_le_mul_of_nonneg_right hMn hn0.le]
  set x : ℝ := Q / n with hx
  have hx0 : 0 < x := by positivity
  have hxsQ : x ≤ sQ := by
    rw [hx, div_le_iff₀ hn0, ← hsQ2]
    linarith [mul_le_mul_of_nonneg_left hnsQ hsQ0.le]
  have hxM : x ≤ M := by
    rw [hx, div_le_iff₀ hn0]
    linarith
  set m₀ : ℕ := ⌈x⌉₊ with hm₀
  have hm₀1 : 1 ≤ m₀ := Nat.ceil_pos.2 hx0
  have hm₀A : m₀ ≤ A.card := Nat.ceil_le.2 hxM
  have hxm₀ : x ≤ m₀ := Nat.le_ceil x
  have hm₀lt : (m₀ : ℝ) < x + 1 := Nat.ceil_lt_add_one hx0.le
  have hm₀sQ : (m₀ : ℝ) ≤ sQ + 1 := by linarith
  -- `m₀ ≤ (p+1)/2`
  have hsQp : sQ ≤ ((p : ℝ) - 1) / 2 := by
    rw [hsQ, Real.sqrt_le_left (by linarith)]
    rw [hQ]
    linarith [mul_le_mul_of_nonneg_right hκ1 (by linarith : (0 : ℝ) ≤ p),
      mul_nonneg (by linarith : (0 : ℝ) ≤ p - 11) (by linarith : (0 : ℝ) ≤ p + 1)]
  have hm₀p : m₀ ≤ (p + 1) / 2 := by
    have h1 : (m₀ : ℝ) < ((p : ℝ) + 1) / 2 := by linarith
    have h2 : (2 * m₀ : ℝ) < p + 1 := by linarith
    have h3 : 2 * m₀ < p + 1 := by exact_mod_cast h2
    omega
  have hm₀n : Q ≤ (m₀ : ℝ) * n := by
    rw [hx, div_le_iff₀ hn0] at hxm₀
    linarith
  -- the per-subset constant is at most the final constant
  set c : ℝ := 1 - (1 - u) ^ 2 + (sQ + 1) / (2 * ((p : ℝ) - 1)) with hc
  have hcc : ∀ A' : Finset (ZMod p), A'.card = m₀ →
      (charSum A' B : ℝ) ≤ c * ((m₀ : ℝ) * n) := by
    intro A' hA'
    have h := subset_bound hp κ hκ0 hκ1 A' B (by omega) (by omega)
      (by rw [hA']; exact hm₀n)
    rw [hA'] at h
    refine h.trans (mul_le_mul_of_nonneg_right ?_ (by positivity))
    rw [hc]
    have huu : u * (1 - u) ≤ 1 / 4 := by linarith [sq_nonneg (u - 1 / 2)]
    have hp1 : (0 : ℝ) < (p : ℝ) - 1 := by linarith
    have h1 : u * (1 - u) * m₀ / (((p : ℝ) - 1) / 2) = 2 * (u * (1 - u)) * m₀ / ((p : ℝ) - 1) := by
      field_simp
    have h2 : 2 * (u * (1 - u)) * m₀ / ((p : ℝ) - 1) ≤ (sQ + 1) / (2 * ((p : ℝ) - 1)) := by
      rw [div_le_div_iff₀ hp1 (by positivity)]
      have hm0 : (0 : ℝ) ≤ m₀ := by positivity
      have hq1 := mul_le_mul_of_nonneg_right huu
        (mul_nonneg hm0 hp1.le : (0 : ℝ) ≤ (m₀ : ℝ) * ((p : ℝ) - 1))
      have hq2 := mul_le_mul_of_nonneg_right hm₀sQ hp1.le
      linarith
    show 1 - (1 - u) ^ 2 + u * (1 - u) * (m₀ : ℝ) / (((p : ℝ) - 1) / 2) ≤ c
    rw [hc]
    linarith
  -- averaging over `m₀`-subsets
  obtain ⟨k, hk⟩ : ∃ k, m₀ = k + 1 := ⟨m₀ - 1, by omega⟩
  have hcs : ∀ T : Finset (ZMod p),
      (charSum T B : ℝ) = ∑ a ∈ T, ∑ b ∈ B, (quadraticChar (ZMod p) (a + b) : ℝ) := by
    intro T
    simp only [charSum, Int.cast_sum]
  have havg : ∑ T ∈ A.powersetCard (k + 1), ∑ a ∈ T, ∑ b ∈ B,
      (quadraticChar (ZMod p) (a + b) : ℝ) =
      ((A.card - 1).choose k : ℝ) * (charSum A B : ℝ) := by
    rw [hcs A]
    exact sum_powersetCard_sum A k (fun a => ∑ b ∈ B, (quadraticChar (ZMod p) (a + b) : ℝ))
  have hsum : ∑ T ∈ A.powersetCard (k + 1), ∑ a ∈ T, ∑ b ∈ B,
      (quadraticChar (ZMod p) (a + b) : ℝ) ≤
      ((A.powersetCard (k + 1)).card : ℝ) * (c * ((m₀ : ℝ) * n)) := by
    rw [← nsmul_eq_mul]
    apply Finset.sum_le_card_nsmul
    intro T hT
    rw [← hcs T]
    exact hcc T (by rw [(Finset.mem_powersetCard.1 hT).2, hk])
  rw [havg, Finset.card_powersetCard] at hsum
  -- `C(m, m₀)·m₀ = m·C(m-1, m₀-1)`
  have hchoose : ((A.card.choose (k + 1) : ℕ) : ℝ) * m₀ = M * ((A.card - 1).choose k : ℕ) := by
    have h : A.card * (A.card - 1).choose k = A.card.choose (k + 1) * (k + 1) := by
      have h0 := Nat.add_one_mul_choose_eq (A.card - 1) k
      rwa [Nat.sub_add_cancel (by omega : 1 ≤ A.card)] at h0
    rw [hk, hM]
    have h' : ((A.card * (A.card - 1).choose k : ℕ) : ℝ) =
        ((A.card.choose (k + 1) * (k + 1) : ℕ) : ℝ) := by rw [h]
    push_cast at h' ⊢
    linarith
  have hpos : (0 : ℝ) < ((A.card - 1).choose k : ℕ) := by
    exact_mod_cast Nat.choose_pos (by omega)
  have h4 : ((A.card - 1).choose k : ℝ) * (charSum A B : ℝ) ≤
      ((A.card - 1).choose k : ℝ) * (c * (M * n)) := by
    calc ((A.card - 1).choose k : ℝ) * (charSum A B : ℝ)
        ≤ ((A.card.choose (k + 1) : ℕ) : ℝ) * (c * ((m₀ : ℝ) * n)) := hsum
      _ = c * n * (((A.card.choose (k + 1) : ℕ) : ℝ) * m₀) := by ring
      _ = c * n * (M * ((A.card - 1).choose k : ℕ)) := by rw [hchoose]
      _ = ((A.card - 1).choose k : ℝ) * (c * (M * n)) := by ring
  exact le_of_mul_le_mul_left h4 hpos

theorem charSum_comm (A B : Finset (ZMod p)) : charSum A B = charSum B A := by
  unfold charSum
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun a _ => ?_
  rw [add_comm]

theorem charSum_dilate (ν : ZMod p) (hν : ¬ IsSquare ν) (A B : Finset (ZMod p)) :
    charSum (A.image (ν * ·)) (B.image (ν * ·)) = -charSum A B := by
  classical
  have hν0 : ν ≠ 0 := fun h => hν (h ▸ ⟨0, by ring⟩)
  have hinj : Function.Injective (ν * · : ZMod p → ZMod p) := mul_right_injective₀ hν0
  have hχ : quadraticChar (ZMod p) ν = -1 := quadraticChar_neg_one_iff_not_isSquare.2 hν
  unfold charSum
  rw [Finset.sum_image fun x _ y _ h => hinj h, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Finset.sum_image fun x _ y _ h => hinj h, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun b _ => ?_
  beta_reduce
  rw [show ν * a + ν * b = ν * (a + b) by ring, map_mul, hχ]
  ring

theorem card_image_dilate (ν : ZMod p) (hν0 : ν ≠ 0) (A : Finset (ZMod p)) :
    (A.image (ν * ·)).card = A.card :=
  Finset.card_image_of_injective _ (mul_right_injective₀ hν0)

/-- **Theorem 2.5 (constant bias saving above `p/2`).** For `p ≥ 11`, `0 < κ ≤ 3/2`,
`u = (1+2κ)^{-1/2}` and `A, B ⊆ F_p` with `|A||B| ≥ (1/2+κ)p`:
`|S(A,B)| ≤ [1 - (1-u)² + (√((1/2+κ)p) + 1)/(2(p-1))]·|A||B|`. -/
theorem constant_bias (hp11 : 11 ≤ p) (κ : ℝ) (hκ0 : 0 < κ) (hκ1 : κ ≤ 3 / 2)
    (A B : Finset (ZMod p)) (hAB : (1 / 2 + κ) * p ≤ (A.card : ℝ) * B.card) :
    |(charSum A B : ℝ)| ≤
      (1 - (1 - 1 / Real.sqrt (1 + 2 * κ)) ^ 2 +
        (Real.sqrt ((1 / 2 + κ) * p) + 1) / (2 * ((p : ℝ) - 1))) * ((A.card : ℝ) * B.card) := by
  classical
  have hp : p ≠ 2 := by omega
  set c : ℝ := 1 - (1 - 1 / Real.sqrt (1 + 2 * κ)) ^ 2 +
    (Real.sqrt ((1 / 2 + κ) * p) + 1) / (2 * ((p : ℝ) - 1)) with hc
  have key : ∀ A' B' : Finset (ZMod p), (1 / 2 + κ) * p ≤ (A'.card : ℝ) * B'.card →
      (charSum A' B' : ℝ) ≤ c * ((A'.card : ℝ) * B'.card) := by
    intro A' B' h
    rcases le_total A'.card B'.card with hle | hle
    · exact upper_bound hp hp11 κ hκ0 hκ1 A' B' hle h
    · have h' : (1 / 2 + κ) * p ≤ (B'.card : ℝ) * A'.card := by linarith
      have h2 : (charSum B' A' : ℝ) ≤ c * ((B'.card : ℝ) * A'.card) :=
        upper_bound hp hp11 κ hκ0 hκ1 B' A' hle h'
      rw [charSum_comm]
      linarith
  have hringChar : ringChar (ZMod p) ≠ 2 := by rw [ZMod.ringChar_zmod_n]; exact hp
  obtain ⟨ν, hν⟩ := FiniteField.exists_nonsquare hringChar
  have hν0 : ν ≠ 0 := fun h => hν (h ▸ ⟨0, by ring⟩)
  have hup := key A B hAB
  have hlow := key (A.image (ν * ·)) (B.image (ν * ·))
    (by rw [card_image_dilate ν hν0, card_image_dilate ν hν0]; exact hAB)
  rw [charSum_dilate ν hν, card_image_dilate ν hν0, card_image_dilate ν hν0] at hlow
  push_cast at hlow
  rw [abs_le]
  constructor <;> linarith

end ConstantBias

end StepanovRobust

#print axioms StepanovRobust.pow_dvd_det_add_mul
#print axioms StepanovRobust.epsPoly_dvd
#print axioms StepanovRobust.hankel_order
#print axioms StepanovRobust.hankel_det_coeff
#print axioms StepanovRobust.hankel_det_natDegree_le
#print axioms StepanovRobust.hankel_stepanov_core
#print axioms StepanovRobust.hankelLead_ne_zero
#print axioms StepanovRobust.hankel_stepanov_paley_of_lead
#print axioms StepanovRobust.hankel_stepanov_paley
#print axioms StepanovRobust.hankel_det_natDegree_eq
#print axioms StepanovRobust.hankel_order_paley
#print axioms StepanovRobust.charSum_eq
#print axioms StepanovRobust.bias_inequality
#print axioms StepanovRobust.bias_bound
#print axioms StepanovRobust.subset_bound
#print axioms StepanovRobust.upper_bound
#print axioms StepanovRobust.constant_bias
