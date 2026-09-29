import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.NumberTheory.LegendreSymbol.Basic

set_option autoImplicit false

open Finset

theorem paley_hp_sharp_example (p : ℕ) [Fact p.Prime] (hp : p % 4 = 1) :
    let A : Finset (ZMod p) := {0, 1}
    let B : Finset (ZMod p) := Finset.univ.filter (fun b => IsSquare b ∧ IsSquare (b + 1))
    B.card = (p + 3) / 4 ∧ (∀ a ∈ A, ∀ b ∈ B, IsSquare (a + b)) ∧
      A.card * B.card = (p - 1) / 2 + (B.filter (fun b => -b ∈ A)).card := by sorry
