import Mathlib.NumberTheory.LegendreSymbol.JacobiSymbol
import Mathlib.Tactic

/-!
An exact obstruction to a proposed Gaussian eighth-moment bound.
This file does NOT prove or disprove the Paley graph conjecture.
-/

set_option autoImplicit false
set_option maxRecDepth 4096
set_option maxHeartbeats 2000000

namespace PaleyResearch

instance : Fact (Nat.Prime 1009) := ⟨by norm_num⟩

def rowSum (B : Finset ℤ) (x : ℕ) : ℤ :=
  ∑ b ∈ B, legendreSym 1009 ((x : ℤ) - b)

def evenMoment (B : Finset ℤ) (r : ℕ) : ℤ :=
  ∑ x ∈ Finset.range 1009, rowSum B x ^ (2 * r)

-- Negatives of the first thirty nonzero squares; all representatives are distinct.
def witness : Finset ℤ :=
  {1008, 1005, 1000, 993, 984, 973, 960, 945, 928, 909,
   888, 865, 840, 813, 784, 753, 720, 685, 648, 609,
   568, 525, 480, 433, 384, 333, 280, 225, 168, 109}

theorem witness_card : witness.card = 30 := by
  decide

theorem witness_representatives : ∀ b ∈ witness, 0 ≤ b ∧ b < 1009 := by
  decide

theorem witness_small : witness.card ^ 2 < 1009 := by
  rw [witness_card]
  decide

theorem witness_zero_row : rowSum witness 0 = 30 := by
  norm_num [rowSum, witness]

theorem witness_moment_lower : (30 : ℤ) ^ 8 ≤ evenMoment witness 4 := by
  have h : rowSum witness 0 ^ (2 * 4) ≤ evenMoment witness 4 := by
    unfold evenMoment
    apply Finset.single_le_sum (f := fun x => rowSum witness x ^ (2 * 4))
    · intro x hx
      rw [pow_mul]
      exact pow_nonneg (sq_nonneg (rowSum witness x)) 4
    · decide
  rw [witness_zero_row] at h
  exact h

theorem gaussian_eighth_moment_bound_fails :
    ¬ evenMoment witness 4 ≤ 105 * 1009 * (witness.card : ℤ) ^ 4 := by
  have h := witness_moment_lower
  rw [witness_card]
  norm_num at h ⊢
  omega

#print axioms witness_zero_row
#print axioms witness_moment_lower
#print axioms gaussian_eighth_moment_bound_fails

end PaleyResearch
