#!/usr/bin/env python3
"""Exact arithmetic for the official profile's trace-zero obstruction.

This checks a field and its trace, not the prize's full reduction or a
Paley conjecture. All acceptance checks use integers only.
"""
from collections import Counter
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from cyclotomic_norm_audit import determinant


def trim(poly):
    poly = list(poly)
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def remainder(poly, modulus, p):
    poly = trim([x % p for x in poly])
    modulus = trim(modulus)
    while len(poly) >= len(modulus):
        shift = len(poly)-len(modulus)
        coefficient = poly[-1]*pow(modulus[-1], -1, p) % p
        for j, value in enumerate(modulus):
            poly[j+shift] = (poly[j+shift]-coefficient*value) % p
        poly = trim(poly)
    return poly


def multiply(left, right, modulus, p):
    answer = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i+j] += a*b
    return remainder(answer, modulus, p)


def power(value, exponent, modulus, p):
    answer = [1]
    while exponent:
        if exponent & 1:
            answer = multiply(answer, value, modulus, p)
        value = multiply(value, value, modulus, p)
        exponent >>= 1
    return answer


def subtract(left, right, p):
    return trim([((left[j] if j < len(left) else 0)
                  -(right[j] if j < len(right) else 0)) % p
                 for j in range(max(len(left), len(right)))])


def gcd(left, right, p):
    while right:
        left, right = right, remainder(left, right, p)
    inverse = pow(left[-1], -1, p)
    return [x*inverse % p for x in left]


def trace(value, degree, modulus, p):
    result = [0]*degree
    for _ in range(degree):
        for j, coefficient in enumerate(value):
            result[j] = (result[j]+coefficient) % p
        value = power(value, p, modulus, p)
    return trim(result)


def official_profile():
    p, degree, n = 2**31-2**24+1, 6, 2**18
    assert all(p % d for d in range(2, isqrt(p)+1))
    modulus, theta = [1, 0, 0, 1, 0, 0, 1], [0, 1]
    generator = 0x6C4A8A45
    assert pow(generator, n, p) == 1
    assert pow(generator, n//2, p) == p-1
    frobenius = [theta]
    for _ in range(degree):
        frobenius.append(power(frobenius[-1], p, modulus, p))
    # Rabin's irreducibility criterion for degree 6: proper prime divisors 2,3.
    assert frobenius[6] == theta
    gcds = [gcd(modulus, subtract(frobenius[j], theta, p), p) for j in [2, 3]]
    assert gcds == [[1], [1]]
    trace_basis = [trace([0]*j+[1], degree, modulus, p) for j in range(degree)]
    assert trace_basis == [[6], [], [], [p-3], [], []]
    # Independent definition of trace: diagonal sum of multiplication matrix.
    def matrix_trace(value):
        result = 0
        for j in range(degree):
            column = multiply(value, [0]*j+[1], modulus, p)
            result += column[j] if j < len(column) else 0
        return result % p
    assert [matrix_trace([0]*j+[1]) for j in range(degree)] == [6, 0, 0, p-3, 0, 0]
    gram = [[matrix_trace(power(theta, i+j, modulus, p))
             for j in range(degree)] for i in range(degree)]
    gram_determinant = determinant(gram) % p
    assert gram_determinant == (-3**9) % p != 0
    assert trace(theta, degree, modulus, p) == [] and theta != []
    assert p < n**4//4 and p**6 > n**4
    assert (p**5-1) % n == 0
    return {
        "p": p, "extension_degree": degree, "extension_size": p**degree,
        "n": n, "base_subgroup_generator": generator,
        "modulus_coefficients_low_to_high": modulus,
        "theta_frobenius_orbit_then_repeat": frobenius,
        "Rabin_gcds_at_degrees_2_and_3": gcds,
        "trace_of_power_basis": trace_basis,
        "trace_pairing_Gram_matrix": gram,
        "trace_pairing_determinant_mod_p": gram_determinant,
        "traces_independently_verified_by_multiplication_matrices": True,
        "nonzero_trace_zero_witness": theta,
        "nonzero_trace_zero_frequency_count": p**5-1,
        "trace_zero_nonzero_multiplicative_cosets_of_H": (p**5-1)//n,
        "prime_field_coset_count": (p-1)//n,
        "maximum_extension_period": n,
        "maximum_justification": "Trace(theta*h)=h*Trace(theta)=0 for every base-field h",
        "prime_below_quartic_window": True, "extension_size_above_quartic_window": True,
        "scope": "Field arithmetic and obstruction only; no prize theorem verified"}


def small_exhaustive():
    p, degree, n = 5, 2, 4
    modulus = [3, 0, 1]  # X^2-2 over F_5.
    assert all((x*x-2) % p for x in range(p))
    field = [[a, b] for a in range(p) for b in range(p)]
    subgroup = [1, 2, 3, 4]
    periods, trace_histogram = Counter(), Counter()
    for value in field:
        tr = trace(value, degree, modulus, p)
        scalar = tr[0] if tr else 0
        assert tr == trim([2*value[0] % p])
        trace_histogram[scalar] += 1
        exponents = Counter()
        for h in subgroup:
            image = trace(multiply(value, [h], modulus, p), degree, modulus, p)
            exponent = image[0] if image else 0
            assert exponent == h*scalar % p
            exponents[exponent] += 1
        if not scalar:
            assert exponents == {0: n}
            period = n
        else:
            assert exponents == {1: 1, 2: 1, 3: 1, 4: 1}
            period = -1  # Sum of all nontrivial fifth roots.
        if any(value):
            periods[period] += 1
    assert periods == {4: 4, -1: 20}
    assert trace_histogram == {j: 5 for j in range(5)}
    counts = [1]+[0]*(p-1)
    moments = []
    for r in range(1, 9):
        counts = [sum(counts[(x-h) % p] for h in subgroup) for x in range(p)]
        energy = sum(value**2 for value in counts)
        q_base = (p*energy-n**(2*r))//n
        q_extension = (p**2*energy-n**(2*r))//n
        assert q_base == 1
        assert q_extension == 5*q_base+4*n**(2*r-1)
        assert sum(count*value**(2*r) for value, count in periods.items()) == n*q_extension
        moments.append({"r": r, "energy": energy,
                        "prime_coset_moment": q_base, "extension_coset_moment": q_extension})
    return {"p": p, "degree": degree, "n": n,
            "all_field_elements_checked": len(field),
            "nonzero_frequency_period_histogram": dict(periods), "moments": moments}


if __name__ == "__main__":
    source = Path(__file__).resolve()
    root = source.parents[1]
    archive = root/"sources/official-prize-2026-09-04"
    manifest = json.loads((archive/"manifest.json").read_text())
    for group in [manifest]+manifest["dependencies"]:
        for entry in group["files"]:
            assert sha256((archive/entry["path"]).read_bytes()).hexdigest() == entry["sha256"]
    result = {"status": "passed", "source_sha256": sha256(source.read_bytes()).hexdigest(),
              "determinant_helper_sha256": sha256((source.parent/"cyclotomic_norm_audit.py").read_bytes()).hexdigest(),
              "source_archive_manifest_sha256": sha256((archive/"manifest.json").read_bytes()).hexdigest(),
              "official_profile": official_profile(), "small_exhaustive": small_exhaustive()}
    (root/"results/extension_trace_audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"],
                      "official_p": result["official_profile"]["p"],
                      "official_n": result["official_profile"]["n"],
                      "official_maximum_extension_period": result["official_profile"]["maximum_extension_period"],
                      "trace_basis": result["official_profile"]["trace_of_power_basis"],
                      "small_cases": result["small_exhaustive"]["all_field_elements_checked"]}, indent=2))
