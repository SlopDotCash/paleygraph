import Mathlib.LinearAlgebra.Lagrange
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.Algebra.Polynomial.RingDivision
import Mathlib.Data.Nat.Prime.Factorial
import Mathlib.Data.Real.Sqrt
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

section
open Polynomial Finset

namespace StepanovHP

variable {p : ℕ} [Fact p.Prime]

/-- Lagrange leading-coefficient weights `c_a = [x^{|A|-1}] L_a(x)`. -/
noncomputable def wt (A : Finset (ZMod p)) (a : ZMod p) : ZMod p :=
  (Lagrange.basis A id a).coeff (A.card - 1)

/-- `Σ_a c_a f(a) = [x^{|A|-1}] f` for every `f` of degree `< |A|`. -/
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

/-- The Stepanov auxiliary polynomial `F(x) = -1 + Σ_a c_a (x + a)^{d + |A| - 1}`. -/
noncomputable def auxF (A : Finset (ZMod p)) (d : ℕ) : (ZMod p)[X] :=
  -1 + ∑ a ∈ A, C (wt A a) * (X + C a) ^ (d + A.card - 1)

theorem coeff_auxF_comp (A : Finset (ZMod p)) (d : ℕ) (b : ZMod p) (j : ℕ) :
    ((auxF A d).comp (X + C b)).coeff j =
      -(if j = 0 then 1 else 0) +
        (∑ a ∈ A, wt A a * (b + a) ^ (d + A.card - 1 - j)) * ((d + A.card - 1).choose j : ZMod p) := by
  simp only [auxF, add_comp, neg_comp, one_comp, Polynomial.sum_comp, mul_comp, C_comp, pow_comp, X_comp]
  rw [coeff_add, coeff_neg, coeff_one, finset_sum_coeff, Finset.sum_mul]
  congr 1
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [coeff_C_mul, show (X + C b + C a : (ZMod p)[X]) = X + C (b + a) by rw [C_add, add_assoc],
    coeff_X_add_C_pow, mul_assoc]

theorem coeff_auxF (A : Finset (ZMod p)) (d : ℕ) (j : ℕ) :
    (auxF A d).coeff j =
      -(if j = 0 then 1 else 0) +
        (∑ a ∈ A, wt A a * a ^ (d + A.card - 1 - j)) * ((d + A.card - 1).choose j : ZMod p) := by
  have h := coeff_auxF_comp A d 0 j
  simpa using h

theorem choose_ne_zero_zmod {n k : ℕ} (hk : k ≤ n) (hn : n < p) :
    ((n.choose k : ℕ) : ZMod p) ≠ 0 := by
  rw [Ne, CharP.cast_eq_zero_iff (ZMod p) p]
  intro h
  have h2 : p ∣ n.factorial := by
    rw [← Nat.choose_mul_factorial_mul_factorial hk, mul_assoc]
    exact Dvd.dvd.mul_right h _
  have := (Nat.Prime.dvd_factorial (Fact.out : p.Prime)).1 h2
  omega

theorem auxF_coeff_d (A : Finset (ZMod p)) (d : ℕ) (hd1 : 1 ≤ d) (hA : 1 ≤ A.card) :
    (auxF A d).coeff d = ((d + A.card - 1).choose d : ZMod p) := by
  rw [coeff_auxF, if_neg (by omega), show d + A.card - 1 - d = A.card - 1 by omega]
  have h := sum_wt_shift_pow A 0 (A.card - 1) (by omega)
  simp only [zero_add, if_true] at h
  rw [h]
  ring

theorem auxF_ne_zero (A : Finset (ZMod p)) (d : ℕ) (hd1 : 1 ≤ d) (hA : 1 ≤ A.card)
    (hAd : d + A.card ≤ p) : auxF A d ≠ 0 := by
  intro h0
  have h := auxF_coeff_d A d hd1 hA
  rw [h0, coeff_zero] at h
  exact choose_ne_zero_zmod (p := p) (k := d) (n := d + A.card - 1) (by omega) (by omega) h.symm

theorem auxF_natDegree_le (A : Finset (ZMod p)) (d : ℕ) (hd1 : 1 ≤ d) :
    (auxF A d).natDegree ≤ d := by
  rw [natDegree_le_iff_coeff_eq_zero]
  intro N hN
  rw [coeff_auxF, if_neg (by omega)]
  by_cases hND : N ≤ d + A.card - 1
  · have h := sum_wt_shift_pow A 0 (d + A.card - 1 - N) (by omega)
    simp only [zero_add] at h
    rw [h, if_neg (by omega)]
    ring
  · rw [Nat.choose_eq_zero_of_lt (by omega)]
    ring

/-- Root multiplicity of the auxiliary polynomial at `b`. -/
theorem le_rootMultiplicity_auxF (A : Finset (ZMod p)) (d : ℕ) (hd1 : 1 ≤ d)
    (hAd : d + A.card ≤ p) (b : ZMod p)
    (hb : ∀ a ∈ A, b + a = 0 ∨ (b + a) ^ d = 1) :
    (if -b ∈ A then A.card - 1 else A.card) ≤ rootMultiplicity b (auxF A d) := by
  rcases Nat.eq_zero_or_pos A.card with hA0 | hApos
  · split_ifs <;> omega
  have hF := auxF_ne_zero A d hd1 hApos hAd
  rw [rootMultiplicity_eq_natTrailingDegree]
  apply le_natTrailingDegree
  · intro h0
    apply hF
    have : ((auxF A d).comp (X + C b)).comp (X - C b) = auxF A d := by
      rw [comp_assoc]; simp
    rw [← this, h0, zero_comp]
  · intro j hj
    have hjA : j < A.card := by split_ifs at hj <;> omega
    rw [coeff_auxF_comp]
    have key : ∀ a ∈ A, wt A a * (b + a) ^ (d + A.card - 1 - j)
        = wt A a * (b + a) ^ (A.card - 1 - j) := by
      intro a ha
      rcases hb a ha with h0 | h1
      · have hmem : -b ∈ A := by
          have : a = -b := by linear_combination h0
          rw [← this]; exact ha
        rw [if_pos hmem] at hj
        rw [h0, zero_pow (by omega), zero_pow (by omega)]
      · rw [show d + A.card - 1 - j = d + (A.card - 1 - j) by omega, pow_add, h1, one_mul]
    rw [Finset.sum_congr rfl key, sum_wt_shift_pow A b (A.card - 1 - j) (by omega)]
    by_cases hj0 : j = 0
    · subst hj0; simp
    · rw [if_neg hj0, if_neg (show A.card - 1 - j ≠ A.card - 1 by omega)]; ring

theorem sum_count_le_card {α : Type*} [DecidableEq α] (s : Multiset α) (B : Finset α) :
    ∑ b ∈ B, s.count b ≤ Multiset.card s := by
  have h1 : ∑ b ∈ B, s.count b ≤ ∑ b ∈ B ∪ s.toFinset, s.count b :=
    Finset.sum_le_sum_of_subset subset_union_left
  have h2 : ∑ b ∈ s.toFinset, s.count b = ∑ b ∈ B ∪ s.toFinset, s.count b :=
    Finset.sum_subset subset_union_right fun x _ hx =>
      Multiset.count_eq_zero.2 (by simpa using hx)
  rw [Multiset.toFinset_sum_count_eq] at h2
  omega

/-- Core of Hanson–Petridis under the size hypothesis `d + |A| ≤ p`. -/
theorem hp_core (d : ℕ) (hd1 : 1 ≤ d) (A B : Finset (ZMod p)) (hAd : d + A.card ≤ p)
    (hAB : ∀ a ∈ A, ∀ b ∈ B, a + b = 0 ∨ (a + b) ^ d = 1) :
    A.card * B.card ≤ d + (B.filter (fun b => -b ∈ A)).card := by
  classical
  rcases Nat.eq_zero_or_pos A.card with hA0 | hApos
  · rw [hA0, zero_mul]; omega
  set F := auxF A d with hFdef
  have hF : F ≠ 0 := auxF_ne_zero A d hd1 hApos hAd
  have hdeg : F.natDegree ≤ d := auxF_natDegree_le A d hd1
  have hm : ∀ b ∈ B, (if -b ∈ A then A.card - 1 else A.card) ≤ F.roots.count b := by
    intro b hb
    rw [count_roots]
    refine le_rootMultiplicity_auxF A d hd1 hAd b fun a ha => ?_
    rw [add_comm]; exact hAB a ha b hb
  have hsum : ∑ b ∈ B, (if -b ∈ A then A.card - 1 else A.card) ≤ d :=
    calc ∑ b ∈ B, (if -b ∈ A then A.card - 1 else A.card)
        ≤ ∑ b ∈ B, F.roots.count b := Finset.sum_le_sum hm
      _ ≤ Multiset.card F.roots := sum_count_le_card _ _
      _ ≤ F.natDegree := card_roots' F
      _ ≤ d := hdeg
  have hsplit : ∑ b ∈ B, (if -b ∈ A then A.card - 1 else A.card)
      + (B.filter (fun b => -b ∈ A)).card = A.card * B.card := by
    rw [Finset.card_filter, ← Finset.sum_add_distrib, mul_comm, ← smul_eq_mul, ← Finset.sum_const]
    refine Finset.sum_congr rfl fun b _ => ?_
    split_ifs <;> omega
  omega

/-- Hanson–Petridis with the size hypothesis removed, assuming `2d < p`. -/
theorem hp_of_two_mul_lt (d : ℕ) (hd1 : 1 ≤ d) (h2d : 2 * d < p) (A B : Finset (ZMod p))
    (hAB : ∀ a ∈ A, ∀ b ∈ B, a + b = 0 ∨ (a + b) ^ d = 1) :
    A.card * B.card ≤ d + (B.filter (fun b => -b ∈ A)).card := by
  classical
  by_cases hA : d + A.card ≤ p
  · exact hp_core d hd1 A B hA hAB
  · rcases B.eq_empty_or_nonempty with hB | ⟨b, hb⟩
    · simp [hB]
    · exfalso
      have hcard : d + 1 ≤ (A.erase (-b)).card := by
        have := Finset.pred_card_le_card_erase (s := A) (a := -b)
        omega
      obtain ⟨A', hA'sub, hA'card⟩ := Finset.exists_subset_card_eq hcard
      have h := hp_core d hd1 A' {b} (by omega) (fun a ha b' hb' => by
        rw [Finset.mem_singleton] at hb'
        subst hb'
        exact hAB a (Finset.mem_of_mem_erase (hA'sub ha)) _ hb)
      have hfilt : ({b} : Finset (ZMod p)).filter (fun b' => -b' ∈ A') = ∅ := by
        rw [Finset.filter_eq_empty_iff]
        intro x hx hmem
        rw [Finset.mem_singleton] at hx
        subst hx
        exact Finset.notMem_erase (-x) A (hA'sub hmem)
      rw [hfilt, Finset.card_empty, Finset.card_singleton, hA'card] at h
      omega

/-- **Hanson–Petridis, Theorem 1.2**: `A + B ⊆ Z_d ∪ {0}`, `d` a proper divisor of `p - 1`. -/
theorem hanson_petridis (p : ℕ) [Fact p.Prime] (d : ℕ) (hd : d ∣ p - 1) (hd' : d < p - 1)
    (A B : Finset (ZMod p)) (hAB : ∀ a ∈ A, ∀ b ∈ B, a + b = 0 ∨ (a + b) ^ d = 1) :
    A.card * B.card ≤ d + (B.filter (fun b => -b ∈ A)).card := by
  obtain ⟨k, hk⟩ := hd
  have hd1 : 1 ≤ d := by
    rcases Nat.eq_zero_or_pos d with h | h
    · subst h; omega
    · exact h
  have hk2 : 2 ≤ k := by
    rcases k with _ | _ | k
    · omega
    · omega
    · omega
  have h2d : 2 * d < p := by
    have hp2 := (Fact.out : p.Prime).two_le
    have : d * 2 ≤ d * k := Nat.mul_le_mul_left d hk2
    omega
  exact hp_of_two_mul_lt d hd1 h2d A B hAB

theorem isSquare_zero_or_pow_half (hp : p ≠ 2) {x : ZMod p} (hx : IsSquare x) :
    x = 0 ∨ x ^ ((p - 1) / 2) = 1 := by
  obtain ⟨r, rfl⟩ := hx
  by_cases hr : r = 0
  · left; simp [hr]
  · right
    have hodd : p % 2 = 1 := Nat.odd_iff.mp ((Fact.out : p.Prime).odd_of_ne_two hp)
    rw [← sq, ← pow_mul, show 2 * ((p - 1) / 2) = p - 1 by omega]
    exact ZMod.pow_card_sub_one_eq_one hr

/-- **Hanson–Petridis for the Paley graph** (`d = (p-1)/2`). -/
theorem hanson_petridis_paley (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (A B : Finset (ZMod p))
    (hAB : ∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) :
    A.card * B.card ≤ (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card := by
  have hp3 : 3 ≤ p := by
    have := (Fact.out : p.Prime).two_le
    omega
  exact hp_of_two_mul_lt ((p - 1) / 2) (by omega) (by omega) A B
    fun a ha b hb => isSquare_zero_or_pow_half hp (hAB a ha b hb)

/-- **Clique bound** (Hanson–Petridis, Corollary 1.5, Paley case). -/
theorem paley_clique_number_bound (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    A.card * (A.card - 1) ≤ (p - 1) / 2 := by
  classical
  have hp2 : p ≠ 2 := by omega
  have h := hanson_petridis_paley p hp2 A (A.image Neg.neg) (by
    intro a ha b hb
    obtain ⟨a', ha', rfl⟩ := Finset.mem_image.1 hb
    by_cases hEq : a = a'
    · subst hEq; exact ⟨0, by ring⟩
    · rw [← sub_eq_add_neg]; exact hA a ha a' ha' hEq)
  have hB : (A.image Neg.neg).card = A.card := Finset.card_image_of_injective _ neg_injective
  have hfilt : (A.image Neg.neg).filter (fun b => -b ∈ A) = A.image Neg.neg := by
    rw [Finset.filter_eq_self]
    intro b hb
    obtain ⟨a', ha', rfl⟩ := Finset.mem_image.1 hb
    simpa using ha'
  rw [hfilt, hB] at h
  rw [Nat.mul_sub_one]
  omega

/-- **Hanson–Petridis, Corollary 1.5** (first sentence): `A - A ⊆ Z_d ∪ {0}` gives `|A|(|A|-1) ≤ d`. -/
theorem difference_bound (p : ℕ) [Fact p.Prime] (d : ℕ) (hd : d ∣ p - 1) (hd' : d < p - 1)
    (A : Finset (ZMod p)) (hA : ∀ a ∈ A, ∀ a' ∈ A, a - a' = 0 ∨ (a - a') ^ d = 1) :
    A.card * (A.card - 1) ≤ d := by
  classical
  have h := hanson_petridis p d hd hd' A (A.image Neg.neg) (by
    intro a ha b hb
    obtain ⟨a', ha', rfl⟩ := Finset.mem_image.1 hb
    rw [← sub_eq_add_neg]; exact hA a ha a' ha')
  have hB : (A.image Neg.neg).card = A.card := Finset.card_image_of_injective _ neg_injective
  have hfilt : (A.image Neg.neg).filter (fun b => -b ∈ A) = A.image Neg.neg := by
    rw [Finset.filter_eq_self]
    intro b hb
    obtain ⟨a', ha', rfl⟩ := Finset.mem_image.1 hb
    simpa using ha'
  rw [hfilt, hB] at h
  rw [Nat.mul_sub_one]
  omega

/-- **Hanson–Petridis, Corollary 1.5** (second sentence): `ω(G_p) ≤ (√(2p-1) + 1)/2`. -/
theorem clique_number_real (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1)
    (A : Finset (ZMod p)) (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    (A.card : ℝ) ≤ (Real.sqrt (2 * (p : ℝ) - 1) + 1) / 2 := by
  have h := paley_clique_number_bound p hp A hA
  rcases Nat.eq_zero_or_pos A.card with h0 | hpos
  · rw [h0, Nat.cast_zero]; positivity
  · have h2 : 2 * (A.card * (A.card - 1)) + 1 ≤ p := by omega
    have h3 : 2 * ((A.card : ℝ) * ((A.card : ℝ) - 1)) + 1 ≤ (p : ℝ) := by
      have h2' : ((2 * (A.card * (A.card - 1)) + 1 : ℕ) : ℝ) ≤ (p : ℝ) := by exact_mod_cast h2
      have hc : ((A.card - 1 : ℕ) : ℝ) = (A.card : ℝ) - 1 := by
        rw [Nat.cast_sub hpos, Nat.cast_one]
      push_cast at h2'
      rw [hc] at h2'
      exact h2'
    have h4 : (2 * (A.card : ℝ) - 1) ^ 2 ≤ 2 * (p : ℝ) - 1 := by nlinarith
    have h5 := Real.abs_le_sqrt h4
    have h6 := le_abs_self (2 * (A.card : ℝ) - 1)
    linarith

end StepanovHP

end

open Finset

theorem solution (p : ℕ) [Fact p.Prime] (d : ℕ)
    (hd : d ∣ p - 1) (hd' : d < p - 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a - a' = 0 ∨ (a - a') ^ d = 1) :
    A.card * (A.card - 1) ≤ d :=
  StepanovHP.difference_bound p d hd hd' A hA

#print axioms solution
