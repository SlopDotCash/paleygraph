import Mathlib.Tactic

/-! Exact finite certificate for a failed unrestricted dyadic-subgroup Gaussian comparison.
The Fourier/normalization interpretation is proved in the accompanying mathematical note.
This is a finite quartic-window obstruction to an exact Gaussian coefficient, not to the target bound. -/

set_option autoImplicit false
set_option maxRecDepth 16384
set_option maxHeartbeats 4000000

namespace PaleyQuarticSubgroupCertificate

def H : Finset ℕ :=
  {1, 2, 4, 8, 16, 32, 64, 128,
  256, 512, 1024, 2048, 4096, 8192, 16384, 32768,
  52347, 65536, 104694, 131072, 209388, 262144, 418776, 524288,
  837552, 1048576, 1675104, 1688191, 2097152, 2506113, 3324035, 3350208,
  3350209, 3376382, 4194304, 4603265, 5012226, 5025313, 5651841, 5862865,
  6176129, 6281641, 6438273, 6491029, 6569345, 6595723, 6634881, 6648070,
  6667649, 6684033, 6692225, 6696321, 6698369, 6699393, 6699905, 6700161,
  6700289, 6700353, 6700385, 6700401, 6700409, 6700413, 6700415, 6700416}

def normalizedCollisions : ℕ :=
  ∑ a ∈ H, ∑ b ∈ H, if (1 + a + 6700417 - b) % 6700417 ∈ H then 1 else 0

theorem prime_modulus : Nat.Prime 6700417 := by norm_num

theorem card_H : H.card = 64 := by decide

theorem positive_representatives : ∀ a ∈ H, 0 < a ∧ a < 6700417 := by decide

theorem contains_one : 1 ∈ H := by decide

theorem product_closed : ∀ a ∈ H, ∀ b ∈ H, a * b % 6700417 ∈ H := by decide

theorem all_roots : ∀ a ∈ H, a ^ 64 % 6700417 = 1 := by decide

theorem normalized_count : normalizedCollisions = 201 := by decide

theorem quartic_window : H.card ^ 4 ≤ 4 * 6700417 ∧ 6700417 ≤ H.card ^ 4 := by
  rw [card_H]
  decide

theorem gaussian_comparison_failure_arithmetic :
    3 * (6700417 - 1) * H.card ^ 2 <
      6700417 * H.card * normalizedCollisions - H.card ^ 4 := by
  rw [card_H, normalized_count]
  decide

#print axioms product_closed
#print axioms all_roots
#print axioms normalized_count
#print axioms gaussian_comparison_failure_arithmetic

end PaleyQuarticSubgroupCertificate
