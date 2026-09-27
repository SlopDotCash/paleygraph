import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

section
open Finset

namespace StepanovHPSharp

variable {p : ℕ} [Fact p.Prime]

/-- Shift orthogonality `Σ_x χ(x) χ(x + c) = -1` for `c ≠ 0`. -/
theorem shift_orth (hp : p ≠ 2) (c : ZMod p) (hc : c ≠ 0) :
    ∑ x : ZMod p, quadraticChar (ZMod p) x * quadraticChar (ZMod p) (x + c) = -1 := by
  have hchar : ringChar (ZMod p) ≠ 2 := by rwa [ZMod.ringChar_zmod_n]
  have hpt : ∀ x : ZMod p, quadraticChar (ZMod p) x * quadraticChar (ZMod p) (x + c) =
      quadraticChar (ZMod p) (1 + c * x⁻¹) - if x = 0 then 1 else 0 := by
    intro x
    by_cases hx : x = 0
    · simp [hx]
    · rw [if_neg hx, sub_zero, ← map_mul]
      have : x * (x + c) = x ^ 2 * (1 + c * x⁻¹) := by
        have hxx := mul_inv_cancel₀ hx
        linear_combination (-(c * x)) * hxx
      rw [this, map_mul, quadraticChar_sq_one' hx, one_mul]
  simp_rw [hpt]
  rw [Finset.sum_sub_distrib]
  have hinj : Function.Injective (fun x : ZMod p => 1 + c * x⁻¹) := by
    intro x y h
    simpa [hc] using h
  have hbij : ∑ x : ZMod p, quadraticChar (ZMod p) (1 + c * x⁻¹)
      = ∑ y : ZMod p, quadraticChar (ZMod p) y :=
    Fintype.sum_bijective _ (Finite.injective_iff_bijective.1 hinj) _ _ (fun _ => rfl)
  rw [hbij, quadraticChar_sum_zero hchar]
  simp

theorem sharp_card (hp : p % 4 = 1) :
    (univ.filter (fun b : ZMod p => IsSquare b ∧ IsSquare (b + 1))).card = (p + 3) / 4 := by
  have hp2 : p ≠ 2 := by omega
  have hchar : ringChar (ZMod p) ≠ 2 := by rwa [ZMod.ringChar_zmod_n]
  have hneg1 : IsSquare (-1 : ZMod p) := ZMod.exists_sq_eq_neg_one_iff.2 (by omega)
  have hm10 : (-1 : ZMod p) ≠ 0 := neg_ne_zero.2 one_ne_zero
  have hχm1 : quadraticChar (ZMod p) (-1) = 1 := (quadraticChar_one_iff_isSquare hm10).2 hneg1
  have hpt : ∀ b : ZMod p,
      (1 + quadraticChar (ZMod p) b) * (1 + quadraticChar (ZMod p) (b + 1)) =
        4 * (if IsSquare b ∧ IsSquare (b + 1) then 1 else 0)
          - 2 * (if b = 0 then 1 else 0) - 2 * (if b = -1 then 1 else 0) := by
    intro b
    by_cases hb0 : b = 0
    · subst hb0
      have h0m1 : (0 : ZMod p) ≠ -1 := fun h => hm10 h.symm
      simp [h0m1]
    · by_cases hbm : b = -1
      · subst hbm
        simp [hm10, hχm1, hneg1]
      · have hb1 : b + 1 ≠ 0 := fun h => hbm (eq_neg_of_add_eq_zero_left h)
        have e1 : IsSquare b ↔ quadraticChar (ZMod p) b = 1 :=
          (quadraticChar_one_iff_isSquare hb0).symm
        have e2 : IsSquare (b + 1) ↔ quadraticChar (ZMod p) (b + 1) = 1 :=
          (quadraticChar_one_iff_isSquare hb1).symm
        rw [if_neg hb0, if_neg hbm]
        rcases quadraticChar_dichotomy hb0 with h1 | h1 <;>
          rcases quadraticChar_dichotomy hb1 with h2 | h2 <;>
          simp only [e1, e2, h1, h2] <;> norm_num
  have hsumL : ∑ b : ZMod p,
      (1 + quadraticChar (ZMod p) b) * (1 + quadraticChar (ZMod p) (b + 1)) = (p : ℤ) - 1 := by
    have hexp : ∀ b : ZMod p,
        (1 + quadraticChar (ZMod p) b) * (1 + quadraticChar (ZMod p) (b + 1)) =
          1 + quadraticChar (ZMod p) b + quadraticChar (ZMod p) (b + 1)
            + quadraticChar (ZMod p) b * quadraticChar (ZMod p) (b + 1) := fun b => by ring
    simp_rw [hexp]
    rw [Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
      shift_orth hp2 1 one_ne_zero, quadraticChar_sum_zero hchar]
    have hshift : ∑ b : ZMod p, quadraticChar (ZMod p) (b + 1) = ∑ b : ZMod p, quadraticChar (ZMod p) b :=
      Fintype.sum_equiv (Equiv.addRight (1 : ZMod p)) _ _ (fun _ => rfl)
    rw [hshift, quadraticChar_sum_zero hchar]
    simp only [Finset.sum_const, Finset.card_univ, ZMod.card, nsmul_eq_mul, mul_one]
    ring
  have hsumR : ∑ b : ZMod p, (4 * (if IsSquare b ∧ IsSquare (b + 1) then (1 : ℤ) else 0)
      - 2 * (if b = 0 then 1 else 0) - 2 * (if b = -1 then 1 else 0))
        = 4 * ((univ.filter (fun b : ZMod p => IsSquare b ∧ IsSquare (b + 1))).card : ℤ) - 4 := by
    rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum,
      ← Finset.mul_sum, Finset.sum_boole]
    simp
    ring
  have key : 4 * ((univ.filter (fun b : ZMod p => IsSquare b ∧ IsSquare (b + 1))).card : ℤ) - 4
      = (p : ℤ) - 1 := by
    rw [← hsumR, ← hsumL]
    exact Finset.sum_congr rfl fun b _ => (hpt b).symm
  omega

/-- **Sharpness of Hanson–Petridis** at `|A| = 2`, `p ≡ 1 (mod 4)`. -/
theorem paley_hp_sharp_example (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) :
    let A : Finset (ZMod p) := {0, 1}
    let B : Finset (ZMod p) := Finset.univ.filter (fun b => IsSquare b ∧ IsSquare (b + 1))
    B.card = (p + 3) / 4 ∧ (∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) ∧
      A.card * B.card = (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card := by
  intro A B
  have hcard : B.card = (p + 3) / 4 := sharp_card hp
  have hneg1 : IsSquare (-1 : ZMod p) := ZMod.exists_sq_eq_neg_one_iff.2 (by omega)
  have h01 : (0 : ZMod p) ≠ 1 := zero_ne_one
  have h0m1 : (0 : ZMod p) ≠ -1 := fun h => (neg_ne_zero.2 (one_ne_zero (α := ZMod p))) h.symm
  refine ⟨hcard, ?_, ?_⟩
  · intro a ha b hb
    simp only [A, B, Finset.mem_insert, Finset.mem_singleton, Finset.mem_filter,
      Finset.mem_univ, true_and] at ha hb
    rcases ha with rfl | rfl
    · simpa using hb.1
    · rw [add_comm]; exact hb.2
  · have hA : A.card = 2 := Finset.card_pair h01
    have hfilt : B.filter (fun b => -b ∈ A) = {0, -1} := by
      ext x
      simp only [A, B, Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_insert,
        Finset.mem_singleton]
      constructor
      · rintro ⟨-, h | h⟩
        · left; exact neg_eq_zero.1 h
        · right; rw [← h, neg_neg]
      · rintro (rfl | rfl)
        · exact ⟨⟨⟨0, by ring⟩, ⟨1, by ring⟩⟩,
            Or.inl neg_zero⟩
        · exact ⟨⟨hneg1, ⟨0, by ring⟩⟩, Or.inr (neg_neg 1)⟩
    rw [hA, hcard, hfilt, Finset.card_pair h0m1]
    omega

end StepanovHPSharp

end

open Finset

theorem solution (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) :
    let A : Finset (ZMod p) := {0, 1}
    let B : Finset (ZMod p) := Finset.univ.filter (fun b => IsSquare b ∧ IsSquare (b + 1))
    B.card = (p + 3) / 4 ∧ (∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) ∧
      A.card * B.card = (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card :=
  StepanovHPSharp.paley_hp_sharp_example p hp

#print axioms solution
