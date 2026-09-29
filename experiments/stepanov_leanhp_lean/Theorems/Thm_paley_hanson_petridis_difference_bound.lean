import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

open Finset

theorem paley_hanson_petridis_difference_bound (p : ℕ) [Fact p.Prime] (d : ℕ)
    (hd : d ∣ p - 1) (hd' : d < p - 1) (A : Finset (ZMod p))
    (hA : ∀ a ∈ A, ∀ a' ∈ A, a - a' = 0 ∨ (a - a') ^ d = 1) :
    A.card * (A.card - 1) ≤ d := by sorry
