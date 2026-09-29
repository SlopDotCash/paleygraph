import Mathlib

set_option autoImplicit false
set_option maxRecDepth 8192
set_option maxHeartbeats 8000000

noncomputable section

namespace PaleySharpConstant

def H : Finset ℕ :=
  {1, 2, 4, 8, 16, 32, 64, 128,
  256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
  52347, 65536, 104694, 131072, 209388, 262144, 418776, 524288,
  837552, 1048576, 1675104, 1688191, 2097152, 2506113, 3324035, 3350208,
  3350209, 3376382, 4194304, 4603265, 5012226, 5025313, 5651841, 5862865,
  6176129, 6281641, 6438273, 6491029, 6569345, 6595723, 6634881, 6648070,
  6667649, 6684033, 6692225, 6696321, 6698369, 6699393, 6699905, 6700161,
  6700289, 6700353, 6700385, 6700401, 6700409, 6700413, 6700415, 6700416}

theorem subgroup_data : Nat.Prime 6700417 ∧ H.card = 64 ∧
    (∀ h ∈ H, 0 < h ∧ h < 6700417) ∧ 1 ∈ H ∧ 2 ∈ H ∧
    (∀ a ∈ H, ∀ b ∈ H, a * b % 6700417 ∈ H) ∧
    H.card ^ 4 ≤ 4 * 6700417 ∧ 6700417 ≤ H.card ^ 4 := by
  constructor
  · norm_num
  · decide

def cosLower (u : ℝ) : ℝ :=
  let c := 1 - u ^ 2 / 32
  let d := 2 * c ^ 2 - 1
  if 0 ≤ d then 2 * d ^ 2 - 1 else -1

theorem cosLower_le_cos (x u : ℝ) (hu : 0 ≤ u) (hu4 : u ≤ 4)
    (hx : |x| ≤ u) : cosLower u ≤ Real.cos x := by
  have hu2 : u ^ 2 ≤ 16 := by nlinarith [sq_le_sq₀ hu (by norm_num : (0:ℝ) ≤ 4) |>.mpr hu4]
  have hx2 : x ^ 2 ≤ u ^ 2 := by
    exact (sq_le_sq).mpr (by simpa [abs_of_nonneg hu] using hx)
  have hc : 0 ≤ 1 - u ^ 2 / 32 := by linarith
  have hq : 1 - u ^ 2 / 32 ≤ Real.cos (x / 4) := by
    have h := Real.one_sub_sq_div_two_le_cos (x := x / 4)
    nlinarith
  have hq2 : (1 - u ^ 2 / 32) ^ 2 ≤ Real.cos (x / 4) ^ 2 :=
    (sq_le_sq₀ hc (hc.trans hq)).mpr hq
  have hh : 2 * (1 - u ^ 2 / 32) ^ 2 - 1 ≤ Real.cos (x / 2) := by
    rw [show x / 2 = 2 * (x / 4) by ring, Real.cos_two_mul]
    linarith
  unfold cosLower
  dsimp only
  split_ifs with hd
  · have hh2 := (sq_le_sq₀ hd (hd.trans hh)).mpr hh
    rw [show x = 2 * (x / 2) by ring, Real.cos_two_mul]
    linarith
  · exact Real.neg_one_le_cos x

def radius (h : ℕ) : ℝ := (63 : ℝ) * ((min h (6700417 - h) : ℕ) : ℝ) / (10 * 6700417)

def periodReal : ℝ := ∑ h ∈ H, Real.cos (2 * Real.pi * (h : ℝ) / 6700417)

def period : ℂ := ∑ h ∈ H,
  Complex.exp (((2 * Real.pi * (h : ℝ) / 6700417 : ℝ) : ℂ) * Complex.I)

theorem period_re : period.re = periodReal := by
  simp only [period, periodReal, Complex.re_sum, Complex.exp_ofReal_mul_I_re]

theorem representatives : ∀ h ∈ H, h ≤ 6700417 := by decide

theorem rational_lower_sum : (43 : ℝ) < ∑ h ∈ H, cosLower (radius h) := by
  norm_num [H, cosLower, radius]

theorem fold_cos (h : ℕ) (hh : h ≤ 6700417) :
    Real.cos (2 * Real.pi * (h : ℝ) / 6700417) =
    Real.cos (2 * Real.pi * ((min h (6700417 - h) : ℕ) : ℝ) / 6700417) := by
  by_cases hs : h ≤ 6700417 - h
  · rw [min_eq_left hs]
  · rw [min_eq_right (le_of_not_ge hs), Nat.cast_sub hh]
    norm_num only [Nat.cast_ofNat]
    have he : 2 * Real.pi * ((6700417 : ℝ) - h) / 6700417 =
        2 * Real.pi - (2 * Real.pi * (h : ℝ) / 6700417) := by ring
    rw [he, Real.cos_two_pi_sub]

theorem radius_bounds (h : ℕ) (hh : h ≤ 6700417) :
    0 ≤ radius h ∧ radius h ≤ 4 ∧
    |2 * Real.pi * ((min h (6700417 - h) : ℕ) : ℝ) / 6700417| ≤ radius h := by
  have hn : 2 * min h (6700417 - h) ≤ 6700417 := by
    have := min_le_left h (6700417 - h)
    have := min_le_right h (6700417 - h)
    omega
  have hn' : (2 : ℝ) * ((min h (6700417 - h) : ℕ) : ℝ) ≤ 6700417 := by exact_mod_cast hn
  have ht : (0 : ℝ) ≤ ((min h (6700417 - h) : ℕ) : ℝ) := Nat.cast_nonneg _
  have hpi := Real.pi_lt_d2.le
  have hx : (0 : ℝ) ≤ 2 * Real.pi * ((min h (6700417 - h) : ℕ) : ℝ) / 6700417 := by positivity
  unfold radius
  rw [abs_of_nonneg hx]
  constructor
  · positivity
  constructor
  · nlinarith
  · nlinarith [mul_le_mul_of_nonneg_right hpi ht]

theorem periodReal_gt_43 : (43 : ℝ) < periodReal := by
  apply rational_lower_sum.trans_le
  unfold periodReal
  apply Finset.sum_le_sum
  intro h hh
  have hr := radius_bounds h (representatives h hh)
  rw [fold_cos h (representatives h hh)]
  exact cosLower_le_cos _ _ hr.1 hr.2.1 hr.2.2

theorem log_ratio_lt_12 : Real.log ((6700417 : ℝ) / 64) < 12 := by
  have h := Real.log_lt_log (by norm_num : (0 : ℝ) < 6700417 / 64)
    (by norm_num : (6700417 : ℝ) / 64 < 2 ^ 17)
  rw [Real.log_pow] at h
  have ht := Real.log_two_lt_d9
  norm_num at h
  linarith

theorem sharp_constant_fails :
    ¬ periodReal ^ 2 ≤ 2 * 64 * Real.log ((6700417 : ℝ) / 64) := by
  have h := periodReal_gt_43
  have hl := log_ratio_lt_12
  nlinarith

theorem period_norm_gt_43 : (43 : ℝ) < ‖period‖ := by
  have h := Complex.re_le_norm period
  rw [period_re] at h
  exact periodReal_gt_43.trans_le h

theorem sharp_period_bound_fails :
    ¬ ‖period‖ ≤ Real.sqrt (2 * 64 * Real.log ((6700417 : ℝ) / 64)) := by
  have hn := period_norm_gt_43
  have hl := log_ratio_lt_12
  have hlog : 0 ≤ Real.log ((6700417 : ℝ) / 64) :=
    Real.log_nonneg (by norm_num)
  intro hb
  have hsq := (sq_le_sq₀ (norm_nonneg period) (Real.sqrt_nonneg _)).mpr hb
  rw [Real.sq_sqrt (by positivity)] at hsq
  nlinarith

#print axioms subgroup_data
#print axioms cosLower_le_cos
#print axioms periodReal_gt_43
#print axioms log_ratio_lt_12
#print axioms sharp_constant_fails
#print axioms sharp_period_bound_fails

end PaleySharpConstant
