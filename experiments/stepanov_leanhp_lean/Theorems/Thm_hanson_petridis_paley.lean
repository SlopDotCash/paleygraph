import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

open Finset

theorem hanson_petridis_paley (p : ℕ) [Fact p.Prime] (hp : p ≠ 2) (A B : Finset (ZMod p))
    (hAB : ∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) :
    A.card * B.card ≤ (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card := by sorry
