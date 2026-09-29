import Mathlib.Combinatorics.Enumerative.DoubleCounting
import Mathlib.Tactic.Linarith
import Lean.Elab.Tactic.Omega

open Finset

theorem solution {α β : Type*} (S : Finset α) (P : Finset β)
    (R : α → β → Prop) [∀ a b, Decidable (R a b)]
    (q t K : ℕ) (hq : 2 ≤ q) (ht : 1 ≤ t) (hK : K < q)
    (hP : P.card + 1 = q * t)
    (hS : ∀ a ∈ S, (q - 1) * t ≤ (P.bipartiteAbove R a).card)
    (hcap : ∀ b ∈ P, (S.bipartiteBelow R b).card ≤ K) :
    S.card ≤ K := by
  have hcount := Finset.card_mul_le_card_mul R hS hcap
  by_contra h
  have hlarge : K + 1 ≤ S.card := by omega
  have hmult := Nat.mul_le_mul_right ((q - 1) * t) hlarge
  have hbudget := congrArg (fun z : ℕ => z * K) hP
  have hqsub : q - 1 + 1 = q := by omega
  have htsub : t - 1 + 1 = t := by omega
  have hsmall : K ≤ q - 1 := by omega
  have hbound := Nat.mul_le_mul_right (t - 1) hsmall
  have hqmul := congrArg (fun z : ℕ => z * (t * K)) hqsub
  have htK := congrArg (fun z : ℕ => z * K) htsub
  have htq := congrArg (fun z : ℕ => z * (q - 1)) htsub
  nlinarith only [hmult, hcount, hbudget, hbound, hqmul, htK, htq, hqsub, hq]

#print axioms solution
