import Mathlib.LinearAlgebra.Dual.Lemmas
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

namespace PaleyMCAProjection

variable {F V I Y : Type*} [Field F] [AddCommGroup V] [Module F V]

def RowsBad [Fintype I] (U : Y → Submodule F V) (W Z : I → V) (γ : F) : Prop :=
  ∃ y : Y, (∀ i : I, W i + γ • Z i ∈ U y) ∧ ∃ i : I, Z i ∉ U y

def ScalarBad (U : Y → Submodule F V) (w z : V) (γ : F) : Prop :=
  ∃ y : Y, w + γ • z ∈ U y ∧ z ∉ U y

def project [Fintype I] (a : I → F) (W : I → V) : V := ∑ i, a i • W i

theorem projection_preserves_bad [Fintype F] [Fintype I] [Nonempty I]
    (U : Y → Submodule F V) (W Z : I → V) :
    ∃ a : I → F, a ≠ 0 ∧ ∀ γ : F,
      RowsBad U W Z γ → ScalarBad U (project a W) (project a Z) γ := by
  classical
  let S : Finset F := Finset.univ.filter (RowsBad U W Z)
  by_cases hS : S.Nonempty
  · letI : Nonempty S := hS.to_subtype
    have witnesses : ∀ γ : S, ∃ (y : Y) (i : I),
        (∀ j : I, W j + (γ : F) • Z j ∈ U y) ∧ Z i ∉ U y := by
      intro γ
      have hb : RowsBad U W Z γ := (Finset.mem_filter.mp γ.property).2
      obtain ⟨y, hfold, i, hi⟩ := hb
      exact ⟨y, i, hfold, hi⟩
    choose y i hfold hnot using witnesses
    have separators : ∀ γ : S, ∃ φ : V →ₗ[F] F,
        φ (Z (i γ)) ≠ 0 ∧ (U (y γ)).map φ = ⊥ := by
      intro γ
      exact Submodule.exists_dual_map_eq_bot_of_notMem (hnot γ) inferInstance
    choose φ hφ hzero using separators
    let forms : S → (I → F) →ₗ[F] F := fun γ =>
      { toFun := fun a => ∑ j, a j * φ γ (Z j)
        map_add' := by intro a b; simp [add_mul, Finset.sum_add_distrib]
        map_smul' := by intro c a; simp [Finset.mul_sum, mul_assoc] }
    have hforms : ∀ γ ∈ (Finset.univ : Finset S), forms γ ≠ 0 := by
      intro γ _ hz
      have hv := congrArg (fun f : (I → F) →ₗ[F] F => f (Pi.single (i γ) 1)) hz
      have hbad : φ γ (Z (i γ)) = 0 := by simpa [forms, Pi.single_apply] using hv
      exact hφ γ hbad
    have hsize : (Finset.univ : Finset S).card ≤ Nat.card F := by
      simpa [Nat.card_eq_fintype_card, Fintype.card_coe, S] using
        (Finset.card_filter_le (s := (Finset.univ : Finset F)) (RowsBad U W Z))
    obtain ⟨a, ha, hall⟩ := paley_common_nonzero_projection
      (Finset.univ : Finset S) forms Finset.univ_nonempty hforms hsize
    refine ⟨a, ha, ?_⟩
    intro γ hbad
    let g : S := ⟨γ, by simp [S, hbad]⟩
    refine ⟨y g, ?_, ?_⟩
    · have heq : project a W + γ • project a Z =
          ∑ j, a j • (W j + γ • Z j) := by
        simp [project, smul_add, Finset.sum_add_distrib, Finset.smul_sum,
          smul_smul, mul_comm]
      rw [heq]
      exact (U (y g)).sum_mem (fun j _ => (U (y g)).smul_mem (a j) (hfold g j))
    · intro hz
      have hmem : φ g (project a Z) ∈ (U (y g)).map (φ g) :=
        ⟨project a Z, hz, rfl⟩
      rw [hzero g] at hmem
      have hvalue : forms g a = 0 := by
        simpa [forms, project, map_sum, map_smul, smul_eq_mul] using hmem
      exact hall g (Finset.mem_univ g) hvalue
  · refine ⟨fun _ => 1, ?_, ?_⟩
    · intro hz
      have he := congrFun hz (Classical.arbitrary I)
      exact one_ne_zero he
    · intro γ hbad
      exact False.elim (hS ⟨γ, by simp [S, hbad]⟩)

theorem rowsBad_iff_input_failure [Fintype I]
    (U : Y → Submodule F V) (W Z : I → V) (γ : F) :
    RowsBad U W Z γ ↔ ∃ y : Y, (∀ i : I, W i + γ • Z i ∈ U y) ∧
      ∃ i : I, W i ∉ U y ∨ Z i ∉ U y := by
  constructor
  · rintro ⟨y, hfold, i, hi⟩
    exact ⟨y, hfold, i, Or.inr hi⟩
  · rintro ⟨y, hfold, i, hi⟩
    refine ⟨y, hfold, i, ?_⟩
    intro hZ
    rcases hi with hW | hZ'
    · apply hW
      simpa using (U y).sub_mem (hfold i) ((U y).smul_mem γ hZ)
    · exact hZ' hZ

theorem scalarBad_iff_input_failure
    (U : Y → Submodule F V) (w z : V) (γ : F) :
    ScalarBad U w z γ ↔ ∃ y : Y, w + γ • z ∈ U y ∧ (w ∉ U y ∨ z ∉ U y) := by
  constructor
  · rintro ⟨y, hfold, hz⟩
    exact ⟨y, hfold, Or.inr hz⟩
  · rintro ⟨y, hfold, hw | hz⟩
    · refine ⟨y, hfold, ?_⟩
      intro hZ
      apply hw
      simpa using (U y).sub_mem hfold ((U y).smul_mem γ hZ)
    · exact ⟨y, hfold, hz⟩

theorem uniform_bad_count_bound_iff [Fintype F] [Fintype I] [Nonempty I]
    (U : Y → Submodule F V) (K : ℕ) :
    (∀ W Z : I → V, Nat.card {γ : F // RowsBad U W Z γ} ≤ K) ↔
      (∀ w z : V, Nat.card {γ : F // ScalarBad U w z γ} ≤ K) := by
  constructor
  · intro h w z
    simpa [RowsBad, ScalarBad] using h (fun _ => w) (fun _ => z)
  · intro h W Z
    obtain ⟨a, _, ha⟩ := projection_preserves_bad U W Z
    let inclusion : {γ : F // RowsBad U W Z γ} →
        {γ : F // ScalarBad U (project a W) (project a Z) γ} :=
      fun γ => ⟨γ, ha γ γ.property⟩
    have hinj : Function.Injective inclusion := by
      intro x y hxy
      exact Subtype.ext (congrArg
        (fun t : {γ : F // ScalarBad U (project a W) (project a Z) γ} => t.val) hxy)
    exact (Nat.card_le_card_of_injective inclusion hinj).trans (h (project a W) (project a Z))

variable {D : Type*}

def agreementSpace (C : Submodule F (D → F)) (T : Set D) : Submodule F (D → F) where
  carrier := {w | ∃ c ∈ C, Set.EqOn w c T}
  zero_mem' := ⟨0, C.zero_mem, fun _ _ => rfl⟩
  add_mem' := by
    rintro w z ⟨c, hc, hw⟩ ⟨d, hd, hz⟩
    refine ⟨c + d, C.add_mem hc hd, ?_⟩
    intro x hx
    change w x + z x = c x + d x
    rw [hw hx, hz hx]
  smul_mem' := by
    rintro a w ⟨c, hc, hw⟩
    refine ⟨a • c, C.smul_mem a hc, ?_⟩
    intro x hx
    change a • w x = a • c x
    rw [hw hx]

def CodeRowsBad [Fintype I] (C : Submodule F (D → F)) (A : Set (Set D))
    (W Z : I → D → F) (γ : F) : Prop :=
  ∃ T ∈ A, (∀ i : I, ∃ c ∈ C, Set.EqOn (W i + γ • Z i) c T) ∧
    ∃ i : I, ¬∃ c ∈ C, Set.EqOn (Z i) c T

def CodeScalarBad (C : Submodule F (D → F)) (A : Set (Set D))
    (w z : D → F) (γ : F) : Prop :=
  ∃ T ∈ A, (∃ c ∈ C, Set.EqOn (w + γ • z) c T) ∧
    ¬∃ c ∈ C, Set.EqOn z c T

theorem codeRowsBad_iff [Fintype I] (C : Submodule F (D → F)) (A : Set (Set D))
    (W Z : I → D → F) (γ : F) :
    CodeRowsBad C A W Z γ ↔ RowsBad (fun T : A => agreementSpace C T) W Z γ := by
  simp [CodeRowsBad, RowsBad, agreementSpace, and_comm, and_left_comm]

theorem codeScalarBad_iff (C : Submodule F (D → F)) (A : Set (Set D))
    (w z : D → F) (γ : F) :
    CodeScalarBad C A w z γ ↔ ScalarBad (fun T : A => agreementSpace C T) w z γ := by
  simp [CodeScalarBad, ScalarBad, agreementSpace, and_comm, and_left_comm]

theorem code_projection_preserves_bad [Fintype F] [Fintype I] [Nonempty I]
    (C : Submodule F (D → F)) (A : Set (Set D)) (W Z : I → D → F) :
    ∃ a : I → F, a ≠ 0 ∧ ∀ γ : F, CodeRowsBad C A W Z γ →
      CodeScalarBad C A (project a W) (project a Z) γ := by
  simpa only [codeRowsBad_iff, codeScalarBad_iff] using
    projection_preserves_bad (fun T : A => agreementSpace C T) W Z

theorem code_uniform_bad_count_bound_iff [Fintype F] [Fintype I] [Nonempty I]
    (C : Submodule F (D → F)) (A : Set (Set D)) (K : ℕ) :
    (∀ W Z : I → D → F, Nat.card {γ : F // CodeRowsBad C A W Z γ} ≤ K) ↔
      (∀ w z : D → F, Nat.card {γ : F // CodeScalarBad C A w z γ} ≤ K) := by
  simpa only [codeRowsBad_iff, codeScalarBad_iff] using
    uniform_bad_count_bound_iff (I := I) (fun T : A => agreementSpace C T) K

theorem codeRowsBad_iff_input_failure [Fintype I]
    (C : Submodule F (D → F)) (A : Set (Set D)) (W Z : I → D → F) (γ : F) :
    CodeRowsBad C A W Z γ ↔
      ∃ T ∈ A, (∀ i : I, ∃ c ∈ C, Set.EqOn (W i + γ • Z i) c T) ∧
        ∃ i : I, (¬∃ c ∈ C, Set.EqOn (W i) c T) ∨
          (¬∃ c ∈ C, Set.EqOn (Z i) c T) := by
  rw [codeRowsBad_iff, rowsBad_iff_input_failure]
  simp [agreementSpace, and_comm, and_left_comm]

end PaleyMCAProjection

#print axioms PaleyMCAProjection.projection_preserves_bad
#print axioms PaleyMCAProjection.rowsBad_iff_input_failure
#print axioms PaleyMCAProjection.scalarBad_iff_input_failure
#print axioms PaleyMCAProjection.uniform_bad_count_bound_iff
#print axioms PaleyMCAProjection.code_projection_preserves_bad
#print axioms PaleyMCAProjection.code_uniform_bad_count_bound_iff
#print axioms PaleyMCAProjection.codeRowsBad_iff_input_failure
