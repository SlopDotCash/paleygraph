import Mathlib.LinearAlgebra.Isomorphisms
import Mathlib.LinearAlgebra.Quotient.Card
import Mathlib.Tactic.FieldSimp
import Mathlib.Data.Fintype.Card
import Mathlib.Tactic.Linarith

import Mathlib.Combinatorics.Enumerative.DoubleCounting
import Mathlib.Tactic.Linarith
import Lean.Elab.Tactic.Omega



open Finset

private theorem projection_incidence_bound {α β : Type*} (S : Finset α) (P : Finset β)
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



theorem paley_functional_kernel_card {F V : Type*} [Field F]
    [AddCommGroup V] [Module F V] (f : V →ₗ[F] F) (hf : f ≠ 0) :
    Nat.card V = Nat.card f.ker * Nat.card F := by
  have hs : Function.Surjective f := by
    obtain ⟨v, hv⟩ := DFunLike.ne_iff.mp hf
    have hv' : f v ≠ 0 := by simpa using hv
    intro y
    refine ⟨(y / f v) • v, ?_⟩
    simp [map_smul, smul_eq_mul, div_mul_cancel₀ _ hv']
  rw [Submodule.card_eq_card_quotient_mul_card f.ker,
      Nat.card_congr (f.quotKerEquivOfSurjective hs).toEquiv]


theorem paley_functional_survival_card {F V : Type*} [Field F]
    [AddCommGroup V] [Module F V] [Fintype F] [Fintype V]
    (f : V →ₗ[F] F) (hf : f ≠ 0) :
    ∃ t : ℕ, 1 ≤ t ∧ Nat.card V = Nat.card F * t ∧
      Nat.card {v : V // f v ≠ 0} = (Nat.card F - 1) * t := by
  classical
  let t := Fintype.card f.ker
  have ht : 1 ≤ t := Fintype.card_pos_iff.mpr ⟨0⟩
  have hcard : Fintype.card V = Fintype.card F * t := by
    simpa [Nat.card_eq_fintype_card, t, Nat.mul_comm] using
      paley_functional_kernel_card f hf
  have hzero : (Finset.univ.filter (fun v : V => f v = 0)).card = t := by
    simp [t, Fintype.card_subtype, LinearMap.mem_ker]
  have hsplit := Finset.card_filter_add_card_filter_not
    (s := Finset.univ) (fun v : V => f v = 0)
  rw [hzero, Finset.card_univ] at hsplit
  have hq : 1 ≤ Fintype.card F := Fintype.card_pos_iff.mpr ⟨0⟩
  have hqsub : Fintype.card F - 1 + 1 = Fintype.card F := by omega
  have hqmul := congrArg (fun z : ℕ => z * t) hqsub
  refine ⟨t, ht, ?_, ?_⟩
  · simpa [Nat.card_eq_fintype_card] using hcard
  · rw [Nat.card_eq_fintype_card, Nat.card_eq_fintype_card, Fintype.card_subtype]
    nlinarith only [hsplit, hcard, hqmul]

theorem paley_nonzero_functional_count {F V α : Type*} [Field F]
    [AddCommGroup V] [Module F V] [Fintype F] [Fintype V] [DecidableEq F]
    (S : Finset α) (f : α → V →ₗ[F] F)
    (hf : ∀ b ∈ S, f b ≠ 0) (K : ℕ) (hK : K < Nat.card F)
    (hcap : ∀ a : V, a ≠ 0 →
      (S.filter (fun b => f b a ≠ 0)).card ≤ K) :
    S.card ≤ K := by
  classical
  by_cases hS : S.Nonempty
  · obtain ⟨b0, hb0⟩ := hS
    obtain ⟨t, ht, htotal, _⟩ := paley_functional_survival_card (f b0) (hf b0 hb0)
    let P : Finset V := Finset.univ.erase 0
    let R : α → V → Prop := fun b a => f b a ≠ 0
    have hq : 2 ≤ Nat.card F := by
      simpa [Nat.card_eq_fintype_card] using (Fintype.one_lt_card : 1 < Fintype.card F)
    have hP : P.card + 1 = Nat.card F * t := by
      have he := Finset.card_erase_add_one (s := (Finset.univ : Finset V))
        (Finset.mem_univ (0 : V))
      simpa [P, Finset.card_univ, ← Nat.card_eq_fintype_card, htotal] using he
    apply projection_incidence_bound S P R (Nat.card F) t K hq ht hK hP
    · intro b hb
      obtain ⟨tb, _, htotalb, hsurvive⟩ := paley_functional_survival_card (f b) (hf b hb)
      have heq : tb = t := Nat.eq_of_mul_eq_mul_left (by omega : 0 < Nat.card F)
        (htotalb.symm.trans htotal)
      subst tb
      have hset : P.bipartiteAbove R b = Finset.univ.filter (fun a : V => f b a ≠ 0) := by
        ext a
        simp only [Finset.mem_bipartiteAbove, P, R, Finset.mem_erase,
          Finset.mem_univ, and_true, Finset.mem_filter, true_and]
        constructor
        · exact fun h => h.2
        · intro ha
          refine ⟨?_, ha⟩
          intro hz
          subst a
          exact ha (map_zero (f b))
      rw [hset]
      simpa [Nat.card_eq_fintype_card, Fintype.card_subtype] using hsurvive.ge
    · intro a ha
      exact hcap a (Finset.mem_erase.mp ha).1
  · simp [Finset.not_nonempty_iff_eq_empty.mp hS]

#print axioms paley_nonzero_functional_count

theorem paley_common_nonzero_projection {F V α : Type*} [Field F]
    [AddCommGroup V] [Module F V] [Fintype F] [Fintype V] [DecidableEq F]
    (S : Finset α) (f : α → V →ₗ[F] F) (hS : S.Nonempty)
    (hf : ∀ b ∈ S, f b ≠ 0) (hsize : S.card ≤ Nat.card F) :
    ∃ a : V, a ≠ 0 ∧ ∀ b ∈ S, f b a ≠ 0 := by
  classical
  by_contra h
  have hpositive : 0 < S.card := Finset.card_pos.mpr hS
  have hcap : ∀ a : V, a ≠ 0 →
      (S.filter (fun b => f b a ≠ 0)).card ≤ S.card - 1 := by
    intro a ha
    have hbad : ∃ b ∈ S, f b a = 0 := by
      by_contra hn
      apply h
      refine ⟨a, ha, ?_⟩
      intro b hb hzero
      exact hn ⟨b, hb, hzero⟩
    have hstrict : (S.filter (fun b => f b a ≠ 0)).card < S.card := by
      apply Finset.card_lt_card
      apply Finset.filter_ssubset.mpr
      simpa using hbad
    omega
  have hbound := paley_nonzero_functional_count S f hf (S.card - 1)
    (by omega) hcap
  omega

#print axioms paley_common_nonzero_projection
