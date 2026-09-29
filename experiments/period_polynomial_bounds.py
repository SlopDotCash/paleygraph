#!/usr/bin/env python3
"""Exact signed moments and polynomial certificates for the resonant subgroup.

The final acceptance checks use integers and Fraction only. NumPy is used
for bounded int64 quotient convolutions; all moment products use Python
integers. No numerical eigenvalues, numerical integration, or float-valued
Fourier sums enter the certificate.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import json

import numpy as np

from quotient_moments import quotient


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(poly, x):
    value = F(0)
    for coefficient in reversed(poly):
        value = value * x + coefficient
    return value


def inner(left, right, moments):
    return sum((a * b * moments[i + j]
                for i, a in enumerate(left) for j, b in enumerate(right)), F(0))


def shifted(poly, origin):
    return [sum((poly[j] * comb(j, i) * origin ** (j - i)
                 for j in range(i, len(poly))), F(0)) for i in range(len(poly))]


def bernstein(poly, left, right):
    degree = len(poly) - 1
    power = [v * (right - left) ** j for j, v in enumerate(shifted(poly, left))]
    return [sum((power[j] * F(comb(k, j), comb(degree, j))
                 for j in range(k + 1)), F(0)) for k in range(degree + 1)]


def rational_list(values):
    return [str(x) for x in values]


def signed_moments(root):
    carry_path = root / "results/dyadic_carries.json"
    quotient_path = root / "results/quotient_moments.json"
    carry = json.loads(carry_path.read_text())
    old = json.loads(quotient_path.read_text())
    assert carry["quotient_moments_sha256"] == digest(quotient_path)
    p, n, generator, depth = 6700417, 64, 2, 12
    assert (carry["p"], carry["n"]) == (p, n)
    subgroup, representatives, neighbors = quotient(p, n, generator)
    assert np.array_equal(subgroup, p - subgroup[::-1])
    counts = np.zeros(len(representatives), dtype=np.int64)
    counts[0] = 1
    zero_counts = [1] + [None] * (2 * depth)
    steps = []
    for r in range(1, depth + 1):
        bound = n * int(counts.max())
        assert bound <= np.iinfo(np.int64).max
        previous = counts
        counts = counts[neighbors].sum(axis=1, dtype=np.int64)
        assert np.all(counts >= 0)
        assert int(counts[0]) + n * sum(map(int, counts[1:])) == n**r
        odd = int(previous[0]) * int(counts[0]) + n * sum(
            int(a) * int(b) for a, b in zip(previous[1:], counts[1:]))
        even = int(counts[0])**2 + n * sum(int(a)**2 for a in counts[1:])
        zero_counts[2*r-1], zero_counts[2*r] = odd, even
        assert even == old["target"]["moments"][r-1]["energy"]
        steps.append({"r": r, "int64_pre_step_upper_bound": bound,
                      "odd_zero_count": odd, "even_zero_count": even})
    carried = [a + b for a, b in zip(carry["intrinsic_periodic_zero_counts"],
                                    carry["prime_alias_counts"])]
    assert zero_counts == carried
    moments = [F(p * count - n**k, n) for k, count in enumerate(zero_counts)]
    assert all(x.denominator == 1 for x in moments)
    assert moments[:4] == [F((p-1)//n), F(-1), F(p-n), F(3*p-n*n)]
    return subgroup, moments, {"p": p, "n": n, "generator": generator,
                              "coset_count": (p-1)//n, "maximum_moment_order": 24,
                              "zero_counts": zero_counts,
                              "signed_coset_moments": [int(x) for x in moments],
                              "quotient_steps": steps,
                              "all_orders_match_independent_carry_counts": True,
                              "input_hashes": {str(carry_path.relative_to(root)): digest(carry_path),
                                               str(quotient_path.relative_to(root)): digest(quotient_path)}}


def orthogonal_polynomials(moments, degree):
    polys, norms = [[F(1)]], [moments[0]]
    for j in range(degree):
        alpha = inner([F(0)] + polys[j], polys[j], moments) / norms[j]
        following = [F(0)] + polys[j]
        for i, coefficient in enumerate(polys[j]):
            following[i] -= alpha * coefficient
        if j:
            beta = norms[j] / norms[j-1]
            for i, coefficient in enumerate(polys[j-1]):
                following[i] -= beta * coefficient
        norm = inner(following, following, moments)
        assert norm > 0
        polys.append(following)
        norms.append(norm)
    # Directly audit the whole Gram matrix, independent of the recurrence.
    for i in range(degree + 1):
        for j in range(degree + 1):
            assert inner(polys[i], polys[j], moments) == (norms[i] if i == j else 0)
    return polys, norms


def kernel_polynomial(polys, norms, center):
    degree = len(polys) - 1
    return [sum((evaluate(polys[j], center) / norms[j] * polys[j][i]
                 for j in range(i, degree + 1)), F(0)) for i in range(degree + 1)]


def endpoint_certificate(polys, norms, moments, orientation, endpoint):
    original = kernel_polynomial(polys, norms, orientation * endpoint)
    poly = [coefficient * orientation**j for j, coefficient in enumerate(original)]
    oriented_moments = [value * orientation**j for j, value in enumerate(moments)]
    norm = inner(poly, poly, oriented_moments)
    translated = shifted(poly, endpoint)
    assert all(value >= 0 for value in translated)
    assert translated[0] > 0 and translated[0]**2 > norm
    assert norm == translated[0]
    return {"orientation": orientation, "endpoint": str(endpoint),
            "polynomial_coefficients_low_to_high": rational_list(poly),
            "shifted_coefficients_low_to_high": rational_list(translated),
            "sum_of_squares_over_cosets": str(norm),
            "positive_exclusion_margin": str(translated[0]**2 - norm)}


def uniqueness_certificate(polys, norms, moments):
    center, left, right, threshold = F(43), F(42), F(4381, 100), F(3, 5)
    poly = kernel_polynomial(polys, norms, center)
    norm = inner(poly, poly, moments)
    adjusted = poly.copy()
    adjusted[0] -= threshold
    coefficients = bernstein(adjusted, left, right)
    assert all(value > 0 for value in coefficients)
    assert norm < 2 * threshold**2
    return {"center": str(center), "interval": [str(left), str(right)],
            "pointwise_polynomial_lower_bound": str(threshold),
            "polynomial_coefficients_low_to_high": rational_list(poly),
            "bernstein_coefficients_of_polynomial_minus_lower_bound": rational_list(coefficients),
            "sum_of_squares_over_cosets": str(norm),
            "two_coset_exclusion_margin": str(2 * threshold**2 - norm),
            "maximum_cosets_in_interval": 1}


def atan_partial(x, last):
    return sum(((-1)**j * x**(2*j+1) / (2*j+1) for j in range(last+1)), F(0))


def floor_at(value, denominator):
    return F(value.numerator * denominator // value.denominator, denominator)


def ceil_at(value, denominator):
    return -floor_at(-value, denominator)


def cos_partial(x, last):
    term, total = F(1), F(1)
    for j in range(1, last+1):
        term *= -x*x / ((2*j-1)*(2*j))
        total += term
    return total


def period_one_interval(p, subgroup):
    # Machin's identity is proved in the accompanying note; these are
    # alternating-series bounds and outward decimal rounding.
    fifth, two_thirty_ninth = F(1, 5), F(1, 239)
    pi_lo_raw = 16 * atan_partial(fifth, 25) - 4 * atan_partial(two_thirty_ninth, 6)
    pi_hi_raw = 16 * atan_partial(fifth, 24) - 4 * atan_partial(two_thirty_ninth, 7)
    assert pi_lo_raw < pi_hi_raw
    pi_lo = floor_at(pi_lo_raw, 10**30)
    pi_hi = ceil_at(pi_hi_raw, 10**30)
    assert 3 < pi_lo < pi_hi < F(22, 7)
    lower, upper = F(0), F(0)
    for h in subgroup:
        a = min(int(h), p-int(h))
        x_lo, x_hi = 2*pi_lo*a/p, 2*pi_hi*a/p
        # Both angles are below pi, so cosine is decreasing. The lower
        # Taylor truncation has odd index; the upper has even index.
        assert x_hi < pi_lo
        lower += cos_partial(x_hi, 21)
        upper += cos_partial(x_lo, 20)
    assert 43 < lower < upper < F(4381, 100)
    rounded_lower = floor_at(lower, 10**18)
    rounded_upper = ceil_at(upper, 10**18)
    assert rounded_upper - rounded_lower <= F(2, 10**18)
    return {"pi_lower": str(pi_lo), "pi_upper": str(pi_hi),
            "atan_truncations_at_one_fifth": {"lower_last_index": 25, "upper_last_index": 24},
            "atan_truncations_at_one_239th": {"lower_last_index": 7, "upper_last_index": 6},
            "cosine_truncations": {"lower_last_index": 21, "upper_last_index": 20},
            "raw_rational_lower": str(lower), "raw_rational_upper": str(upper),
            "outward_rounded_lower": str(rounded_lower), "outward_rounded_upper": str(rounded_upper),
            "lower_decimal_numerator_at_1e18": rounded_lower.numerator * (10**18 // rounded_lower.denominator),
            "upper_decimal_numerator_at_1e18": rounded_upper.numerator * (10**18 // rounded_upper.denominator)}


def main():
    source = Path(__file__).resolve()
    root = source.parents[1]
    subgroup, moments, count_audit = signed_moments(root)
    print("All 25 signed moment orders agree with independent exact carry counts.", flush=True)
    polys, norms = orthogonal_polynomials(moments, 12)
    alphas = [inner([F(0)] + polys[j], polys[j], moments) / norms[j] for j in range(12)]
    betas = [norms[j] / norms[j-1] for j in range(1, 13)]
    assert all(abs(value) <= 3 * (j+1) for j, value in enumerate(alphas))
    assert all(value**2 <= 64 * (j+1) for j, value in enumerate(alphas))
    assert all(value <= 64 * (j+1) for j, value in enumerate(betas))
    mean = moments[1] / moments[0]
    variance = moments[2] / moments[0] - mean**2
    fourth = (moments[4] - 4*mean*moments[3] + 6*mean**2*moments[2]
              - 4*mean**3*moments[1] + mean**4*moments[0]) / moments[0]
    assert mean == alphas[0] and variance == betas[0]
    assert fourth == betas[0] * (betas[0] + betas[1] + (alphas[1]-alphas[0])**2)
    assert betas[1] < 2*betas[0] and fourth > 3*variance**2
    upper = endpoint_certificate(polys, norms, moments, 1, F(4381, 100))
    lower = endpoint_certificate(polys, norms, moments, -1, F(26))
    unique = uniqueness_certificate(polys, norms, moments)
    interval = period_one_interval(count_audit["p"], subgroup)
    result = {"source_sha256": digest(source),
              "dependency_hashes": {"experiments/quotient_moments.py": digest(root/'experiments/quotient_moments.py'),
                                    "experiments/paley_exact.py": digest(root/'experiments/paley_exact.py')},
              "scope": "Finite exact maximum and unique maximizing coset, not an asymptotic Paley proof; not Lean.",
              "signed_moment_audit": count_audit,
              "orthogonal_polynomials": [{"degree": j, "coefficients": rational_list(poly),
                                          "squared_norm": str(norms[j])} for j, poly in enumerate(polys)],
              "recurrence_coefficients": {"alpha_j_for_j_0_through_11": rational_list(alphas),
                                          "beta_j_for_j_1_through_12": rational_list(betas),
                                          "finite_sqrt_alpha_bound_A": 1,
                                          "additional_finite_linear_alpha_bound": 3,
                                          "finite_beta_bound_B": 1,
                                          "uniform_coefficient_bound_proved": False},
              "central_fourth_moment_decomposition": {"variance": str(variance),
                                                       "fourth_moment": str(fourth),
                                                       "positive_gaussian_excess": str(fourth-3*variance**2),
                                                       "beta2_less_than_twice_beta1": True},
              "upper_endpoint_certificate": upper, "lower_endpoint_certificate": lower,
              "unique_large_coset_certificate": unique, "period_one_interval": interval,
              "conclusion": {"M_equals_period_one": True, "maximizing_frequency_set": "H=<2>",
                             "all_other_cosets_absolute_period_less_than": 42,
                             "global_period_interval": ["-26", "4381/100"]}}
    path = root / 'results/period_polynomial_bounds.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    for label in ['lower', 'upper']:
        digits = interval[f'{label}_decimal_numerator_at_1e18']
        print(f"M {label}: {digits // 10**18}.{digits % 10**18:018d}")
    print("Only the subgroup coset has a period at least 42; all other absolute periods are below 42.")
    print(path)


if __name__ == '__main__':
    main()
