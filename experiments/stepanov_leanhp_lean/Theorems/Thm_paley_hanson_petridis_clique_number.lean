import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.Data.Real.Sqrt

set_option autoImplicit false

open Finset

theorem paley_hanson_petridis_clique_number (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1)
    (A : Finset (ZMod p)) (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    (A.card : ℝ) ≤ (Real.sqrt (2 * (p : ℝ) - 1) + 1) / 2 := by sorry
