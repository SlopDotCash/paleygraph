import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

open Finset

theorem paley_clique_number_bound (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a ≠ a' → IsSquare (a - a')) :
    A.card * (A.card - 1) ≤ (p - 1) / 2 := by sorry
