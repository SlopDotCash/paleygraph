import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

open Finset

theorem hanson_petridis_theorem_1_2 (p : ℕ) [Fact p.Prime] (d : ℕ) (hd : d ∣ p - 1) (hd' : d < p - 1)
    (A B : Finset (ZMod p)) (hAB : ∀ a ∈ A, ∀ b ∈ B, a + b = 0 ∨ (a + b) ^ d = 1) :
    A.card * B.card ≤ d + (B.filter (fun b => -b ∈ A)).card := by sorry
