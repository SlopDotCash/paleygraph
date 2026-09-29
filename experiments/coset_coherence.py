#!/usr/bin/env python3
"""Exact coset multiplication identities and a synthetic-spectrum obstruction.

No floating-point computations are used. Synthetic vectors are NOT actual
Gauss periods. Their coordinates lie in Q(sqrt(D)); the research note proves
their autocorrelation identities and explains the arithmetic they violate.
"""
from collections import Counter
from fractions import Fraction
from math import isqrt
from pathlib import Path
import hashlib
import json

from paley_exact import prime
from subgroup_moments import primitive_root


def legendre(x, p):
    residue = pow(x % p, (p - 1) // 2, p)
    assert residue in (0, 1, p - 1)
    return -1 if residue == p - 1 else residue


def coset_kernel(p, n, prescribed_generator=None):
    assert prime(p) and (p - 1) % n == 0 and n % 2 == 0
    m, g = (p - 1) // n, primitive_root(p)
    h = pow(g, m, p) if prescribed_generator is None else prescribed_generator
    H = sorted({pow(h, j, p) for j in range(n)})
    assert len(H) == n and p - 1 in H
    # x^n identifies xH. Multiplication by g shifts the quotient index by one.
    lookup, value, step = {}, 1, pow(g, n, p)
    for j in range(m):
        assert value not in lookup
        lookup[value] = j
        value = value * step % p
    assert value == 1
    kernel = Counter(lookup[pow((1 + x) % p, n, p)] for x in H if x != p - 1)
    assert sum(kernel.values()) == n - 1
    exceptional = lookup[pow(2, n, p)]
    assert all(count % 2 == (j == exceptional) for j, count in kernel.items())
    assert kernel[exceptional] % 2 == 1

    # Verify the entire coefficient identity in the group ring of F_p:
    # (sum_{x in H} [x])^2 = n[0] + sum_{t in H, t!=-1} sum_{x in H} [(1+t)x].
    lhs = Counter((x + y) % p for x in H for y in H)
    rhs = Counter({0: n})
    rhs.update((1 + t) * x % p for t in H if t != p - 1 for x in H)
    assert lhs == rhs
    energy = sum(count * count for count in lhs.values())
    kernel_square_sum = sum(count * count for count in kernel.values())
    assert energy == n * n + n * kernel_square_sum
    assert kernel_square_sum >= 2 * n - 3
    return {"p": p, "n": n, "m": m, "primitive_root": g,
            "subgroup_generator": h, "kernel": sorted(kernel.items()),
            "kernel_at_zero": kernel[0], "kernel_square_sum": kernel_square_sum,
            "nonzero_kernel_size": len(kernel),
            "kernel_multiplicity_histogram": dict(sorted(Counter(kernel.values()).items())),
            "energy2": energy, "full_group_ring_identity_checked": True}


def synthetic(p, n, verify_small_autocorrelations=False):
    m = (p - 1) // n
    assert prime(p) and prime(m) and m % 4 == 3 and m > n * n
    kernel_data = coset_kernel(p, n)
    a, D = Fraction(n + 1, m), Fraction(p - n * n, m)
    # Vector: v_0 = n-a; v_j = -a+sqrt(D)*chi_m(j), j != 0.
    assert 0 < D < n and 0 < a < 1
    assert n - a > n - 1
    assert (n - a) ** 2 > D
    assert (n - a) ** 2 < n * n
    assert n - 2 * a > 0 and D < (n - 2 * a) ** 2
    assert D.denominator == m and isqrt(m) ** 2 != m
    # Exact symbolic autocorrelation. The coefficient of sqrt(D) vanishes
    # because chi is odd. Its rational coefficients are as follows.
    off_diagonal = m * a * a - 2 * n * a - D
    diagonal = m * a * a - 2 * n * a + n * n + (m - 1) * D
    assert off_diagonal == -n and diagonal == p - n
    assert -m * a + n == -1

    checked = 0
    if verify_small_autocorrelations:
        chi = [legendre(j, m) for j in range(m)]
        for shift in range(m):
            actual = sum(chi[j] * chi[(j + shift) % m] for j in range(m))
            assert actual == (m - 1 if shift == 0 else -1)
        checked = m

    # At index zero the period equation would read v_0^2=n+sum_s k_s v_s.
    # Compute its residual A+B*sqrt(D) exactly; D is irrational, so any
    # nonzero coefficient proves this is not an actual period vector.
    k0 = kernel_data["kernel_at_zero"]
    A = (n - a) ** 2 - n + (n - 1) * a - n * k0
    B = -sum(count * legendre(j, m) for j, count in kernel_data["kernel"])
    assert A != 0 or B != 0
    assert A * A != B * B * D
    return {"p": p, "n": n, "m": m, "quartic_window": n**4 // 4 <= p <= n**4,
            "a": str(a), "D": str(D), "spike": str(n - a),
            "coordinate_formula": "v_0=n-a; v_j=-a+sqrt(D)*Legendre(j,m) for j!=0",
            "nonprincipal_coordinate_squared_scale": str(D),
            "all_coordinates_bounded_by_n": True,
            "sum": -1, "autocorrelation_zero_shift": p - n,
            "autocorrelation_other_shifts": -n,
            "direct_character_autocorrelations_checked": checked,
            "nonlinear_residual_at_zero": {"rational_coefficient": str(A),
                                           "sqrt_D_coefficient": B,
                                           "squared_nonzero_check": str(A*A-B*B*D)},
            "actual_period_kernel": kernel_data}


def main():
    source = Path(__file__).resolve()
    examples = [(1049, 8), (17393, 16), (263009, 32), (4201409, 64)]
    data = {"source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "claim": "Obstruction to using only coset autocorrelations; not counterexamples to Paley.",
            "exact_synthetic_examples": [synthetic(p, n, n <= 16) for p, n in examples],
            "original_resonant_example": coset_kernel(6700417, 64, 2)}
    original = data["original_resonant_example"]
    assert original["energy2"] == 12864
    assert original["kernel_at_zero"] == 3
    assert original["kernel_square_sum"] == 137
    assert original["kernel_multiplicity_histogram"] == {2: 28, 3: 1, 4: 1}
    path = source.parents[1] / "results" / "coset_coherence.json"
    path.write_text(json.dumps(data, indent=2) + "\n")
    print("Exact group-ring and kernel-energy identities passed in five prime fields.")
    print("Four synthetic vectors pass all symbolic autocorrelations and fail the nonlinear identity.")
    print(f"Original kernel histogram: {original['kernel_multiplicity_histogram']}")
    print(path)


if __name__ == "__main__":
    main()
