#!/usr/bin/env python3
"""Integer certificates for research/riesz-tail-route.md.

Checks finite group-ring products and a strict counterexample to the
stronger assertion Z_+(t)<=1. Does not prove the Paley target.
"""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import isqrt
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]


def group(p, n, generator):
    assert all(p % divisor for divisor in range(2, isqrt(p) + 1))
    assert n % 2 == 0 and pow(generator, n, p) == 1
    # All orders used here are powers of two.
    assert n & (n - 1) == 0 and pow(generator, n // 2, p) == p - 1
    values = sorted({pow(generator, j, p) for j in range(n)})
    assert len(values) == n
    half = sorted({min(x, p - x) for x in values})
    assert len(half) * 2 == n
    return values, half


def riesz_coefficients(frequencies, p, sign=1):
    """Product of (8 + sign*X^h + sign*X^-h), modulo X^p-1.

Dividing by 8^d yields product(1 + sign*cos(2*pi*b*h/p)/4)
at a p-th root of unity. No numerical trigonometry is used.
"""
    coeffs = {0: 1}
    for h in frequencies:
        following = defaultdict(int)
        for index, count in coeffs.items():
            following[index] += 8 * count
            following[(index + h) % p] += sign * count
            following[(index - h) % p] += sign * count
        coeffs = {index: count for index, count in following.items() if count}
    assert sum(coeffs.values()) == (8 + 2 * sign) ** len(frequencies)
    return coeffs


def dissociated_subset(domain, p):
    selected, sums = [], {0}
    for x in domain:
        shifted = {(s + x) % p for s in sums}
        if sums.isdisjoint(shifted):
            selected.append(x)
            sums |= shifted
    assert len(sums) == 2 ** len(selected) <= p
    assert all(not sums.isdisjoint({(s + x) % p for s in sums})
               for x in domain if x not in selected)
    return selected


def small_example(p, n):
    generator = next(g for g in range(2, p)
                     if pow(g, n, p) == 1 and pow(g, n // 2, p) == p - 1)
    domain, half = group(p, n, generator)
    independent = dissociated_subset(domain, p)
    normalization = 8 ** len(independent)
    for sign in [1, -1]:
        coeffs = riesz_coefficients(independent, p, sign)
        assert coeffs.get(0, 0) == normalization
    # Exhaust all signed relations in the smaller independent sets.
    zero_relations = sum(sum(e * x for e, x in zip(eps, independent)) % p == 0
                         for eps in product([-1, 0, 1], repeat=len(independent)))
    assert zero_relations == 1
    products = []
    for sign in [1, -1]:
        coeffs = riesz_coefficients(half, p, sign)
        denominator = 8 ** len(half)
        principal = F((8 + 2 * sign) ** len(half), denominator)
        average = F(coeffs.get(0, 0), denominator)
        centered = (p * average - principal) / (p - 1)
        assert centered > 0
        # The entire product is invariant under multiplication by H.
        assert all(count == coeffs.get(index * generator % p, 0)
                   for index, count in coeffs.items())
        products.append({"sign": sign, "constant_coefficient": coeffs.get(0, 0),
                         "denominator": denominator, "all_frequency_mean": str(average),
                         "principal_value": str(principal), "nonprincipal_mean": str(centered),
                         "nonprincipal_mean_exceeds_one": centered > 1,
                         "coefficient_sha256": sha256(json.dumps(sorted(coeffs.items())).encode()).hexdigest()})
    return {"p": p, "n": n, "generator": generator, "half_domain": half,
            "maximal_dissociated_subset": independent,
            "dissociated_dimension": len(independent),
            "distinct_subset_sums": 2 ** len(independent),
            "signed_assignments_checked": 3 ** len(independent),
            "riesz_products_at_t_one_quarter": products}


def quartic_counterexample():
    p, n, generator = 67403009, 128, 64701253
    domain, half = group(p, n, generator)
    assert n ** 4 <= 4 * p <= 4 * n ** 4
    signed_pairs = defaultdict(list)
    for i, j in combinations(range(len(half)), 2):
        for a, b in product([-1, 1], repeat=2):
            signed_pairs[(a * half[i] + b * half[j]) % p].append(((i, a), (j, b)))
    four_relations = set()
    for value, pairs in signed_pairs.items():
        for left in pairs:
            for right in signed_pairs.get(-value % p, []):
                if len({i for i, _ in left + right}) == 4:
                    four_relations.add(tuple(sorted(left + right)))
    three_relations = set()
    for i, x in enumerate(half):
        for sign in [-1, 1]:
            for pair in signed_pairs.get(-sign * x % p, []):
                if i not in {j for j, _ in pair}:
                    three_relations.add(tuple(sorted(pair + ((i, sign),))))
    assert len(four_relations) == 128 and len(three_relations) == 0
    for relation in four_relations:
        assert sum(sign * half[i] for i, sign in relation) % p == 0
    # Independent parametrization by multiplying a known zero quadruple.
    representative = [1, 1074964, 1550267, 64777777]
    orbit = {tuple(sorted(x * h % p for x in representative)) for h in domain}
    recovered = {tuple(sorted(sign * half[i] % p for i, sign in relation))
                 for relation in four_relations}
    assert sum(representative) == p and recovered == orbit and len(orbit) == n
    t = F(1, 4)
    # The degree-four term of the zero coefficient is 128*(t/2)^4=8*t^4.
    # Every omitted coefficient is nonnegative at positive t.
    mean_lower = (p * (1 + 8 * t ** 4) - (1 + t) ** (n // 2)) / (p - 1)
    assert mean_lower > F(201, 200)
    assert (p + 32) * 4 ** 64 > 32 * 5 ** 64
    return {"p": p, "n": n, "generator": generator, "half_domain": half,
            "squarefree_signed_zero_counts_degrees_0_to_4": [1, 0, 0, 0, 128],
            "degree_four_riesz_coefficient": 8,
            "signed_quadruples": sorted(four_relations),
            "all_quadruples_in_one_H_orbit": True,
            "positive_t": str(t), "nonprincipal_mean_lower": str(mean_lower),
            "strict_lower_exceeds_201_over_200": True,
            "scope": "Only Z_+(1/4)<=1 is refuted; exp(C*n*t^2) and Paley remain open."}


def main():
    result = {"status": "passed; main goal unproved", "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "small_cases": [small_example(97, 32), small_example(2017, 8), small_example(17393, 16)],
              "quartic_counterexample": quartic_counterexample()}
    destination = ROOT / "results/riesz_tail_audit.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"],
                      "small_dissociated_dimensions": [r["dissociated_dimension"] for r in result["small_cases"]],
                      "quarter_scale_mean_lower": result["quartic_counterexample"]["nonprincipal_mean_lower"],
                      "output": str(destination)}, indent=2))


if __name__ == "__main__":
    main()
