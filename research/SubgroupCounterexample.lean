import Mathlib.Tactic

/-! Exact finite certificate for a failed unrestricted dyadic-subgroup Gaussian comparison.
The Fourier/normalization interpretation is proved in the accompanying mathematical note.
This file does not claim a counterexample in the asymptotic prize regime. -/

set_option autoImplicit false
set_option maxRecDepth 16384
set_option maxHeartbeats 4000000

namespace PaleySubgroupCertificate

def H : Finset ℕ :=
  {1, 2, 4, 8, 16, 32, 64, 128,
  255, 256, 257, 510, 512, 514, 1020, 1024,
  1028, 2040, 2048, 2056, 4080, 4096, 4112, 8160,
  8192, 8224, 16320, 16384, 16448, 32640, 32641, 32768,
  32769, 32896, 32897, 49089, 49153, 49217, 57313, 57345,
  57377, 61425, 61441, 61457, 63481, 63489, 63497, 64509,
  64513, 64517, 65023, 65025, 65027, 65280, 65281, 65282,
  65409, 65473, 65505, 65521, 65529, 65533, 65535, 65536}

def normalizedCollisions : ℕ :=
  ∑ a ∈ H, ∑ b ∈ H, if (1 + a + 65537 - b) % 65537 ∈ H then 1 else 0

theorem prime_modulus : Nat.Prime 65537 := by norm_num

theorem card_H : H.card = 64 := by decide

theorem positive_representatives : ∀ a ∈ H, 0 < a ∧ a < 65537 := by decide

theorem contains_one : 1 ∈ H := by decide

theorem product_closed : ∀ a ∈ H, ∀ b ∈ H, a * b % 65537 ∈ H := by decide

theorem all_roots : ∀ a ∈ H, a ^ 64 % 65537 = 1 := by decide

theorem normalized_count : normalizedCollisions = 309 := by decide

theorem gaussian_comparison_failure_arithmetic :
    3 * (65537 - 1) * H.card ^ 2 <
      65537 * H.card * normalizedCollisions - H.card ^ 4 := by
  rw [card_H, normalized_count]
  decide

#print axioms product_closed
#print axioms all_roots
#print axioms normalized_count
#print axioms gaussian_comparison_failure_arithmetic

end PaleySubgroupCertificate
