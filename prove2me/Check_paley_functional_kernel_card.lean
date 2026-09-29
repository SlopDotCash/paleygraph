import Mathlib.LinearAlgebra.Isomorphisms
import Mathlib.LinearAlgebra.Quotient.Card
import Mathlib.Tactic.FieldSimp
import Mathlib.Data.Fintype.Card
import Mathlib.Tactic.Linarith

-- Local next-input check, not a platform submission.
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

#print axioms paley_functional_kernel_card

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

#print axioms paley_functional_survival_card
