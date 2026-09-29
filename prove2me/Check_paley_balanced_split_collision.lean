import Mathlib.Algebra.Ring.Commute
import Mathlib.Tactic.Ring

namespace PaleyBalancedSplit

theorem two_equations_force_collision {F : Type*} [Field F]
    (A B c d : F) (hA : A ≠ 0)
    (h₁ : A * c + B * d = 0) (h₂ : A * d + B * c = 0) :
    c = d ∨ c = -d := by
  have hz : A * (c * c - d * d) = 0 := by
    calc
      A * (c * c - d * d) = c * (A * c + B * d) - d * (A * d + B * c) := by ring
      _ = 0 := by rw [h₁, h₂]; ring
  have hs : c * c = d * d :=
    sub_eq_zero.mp ((mul_eq_zero.mp hz).resolve_left hA)
  exact mul_self_eq_mul_self_iff.mp hs

theorem two_balanced_splits_force_collision {F : Type*} [Field F]
    (a b c d e f : F) (ha : a ≠ 0) (hb : b ≠ 0)
    (h₁ : a * b * c + d * e * f = 0)
    (h₂ : a * b * d + c * e * f = 0) : c = d ∨ c = -d := by
  apply two_equations_force_collision (a * b) (e * f) c d (mul_ne_zero ha hb)
  · simpa only [mul_assoc, mul_left_comm, mul_comm] using h₁
  · simpa only [mul_assoc, mul_left_comm, mul_comm] using h₂

theorem no_two_balanced_splits {F : Type*} [Field F]
    (a b c d e f : F) (ha : a ≠ 0) (hb : b ≠ 0)
    (hne : c ≠ d) (hnopp : c ≠ -d) :
    ¬(a * b * c + d * e * f = 0 ∧ a * b * d + c * e * f = 0) := by
  rintro ⟨h₁, h₂⟩
  exact (two_balanced_splits_force_collision a b c d e f ha hb h₁ h₂).elim hne hnopp

end PaleyBalancedSplit

#print axioms PaleyBalancedSplit.two_equations_force_collision
#print axioms PaleyBalancedSplit.two_balanced_splits_force_collision
#print axioms PaleyBalancedSplit.no_two_balanced_splits
