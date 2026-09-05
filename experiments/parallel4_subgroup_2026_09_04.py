"""Bounded exact audits for the growing-class theorem; no prime scan.

The analytic prime-count asymptotic is supplied by Thorner--Zaman, not
tested numerically here. Only the fixed certificates below are evaluated.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import comb, factorial, isqrt, log, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "research/parallel4-subgroup-2026-09-04.md"
PRIOR = ROOT / "results/cyclotomic_norm_audit.json"
CERTIFICATES = [
    (5, 4), (17, 8), (41, 8), (17, 16), (97, 16), (113, 16),
    (193, 16), (257, 16), (337, 16), (2113, 32),
    (73, 4), (2017, 8), (65521, 16),
    (6700417, 64), (67403009, 128), (17189277697, 512),
]


def prime(p):
    """Exact trial division, only for the fixed inputs (largest < 2^35)."""
    return p >= 2 and (p == 2 or (p % 2 == 1 and all(
        p % d for d in range(3, isqrt(p) + 1, 2))))


def generator(p, n):
    # A bounded search for an element, not a search for primes.
    for a in range(2, min(p, 1002)):
        g = pow(a, (p - 1) // n, p)
        if pow(g, n // 2, p) == p - 1:
            return g
    raise AssertionError("Fixed seed bound exhausted")


def pair_counts(a, b, p):
    return Counter((x + y) % p for x in a for y in b)


def square_norm(counts):
    return sum(c * c for c in counts.values())


def determinant(matrix):
    """Fraction-free integer determinant with exact divisions checked."""
    a = [row[:] for row in matrix]
    d, sign, previous = len(a), 1, 1
    for j in range(d - 1):
        if not a[j][j]:
            pivot = next((i for i in range(j + 1, d) if a[i][j]), None)
            if pivot is None:
                return 0
            a[j], a[pivot] = a[pivot], a[j]
            sign = -sign
        pivot = a[j][j]
        for i in range(j + 1, d):
            for h in range(j + 1, d):
                numerator = pivot * a[i][h] - a[i][j] * a[j][h]
                assert numerator % previous == 0
                a[i][h] = numerator // previous
        for i in range(j + 1, d):
            a[i][j] = 0
        previous = pivot
    return sign * a[-1][-1]


def polynomial(n):
    """Integer coefficients of P_n, from root-filter power sums."""
    d = n // 2 - 1
    powers = [0]
    for r in range(1, d + 1):
        twice = n * sum(comb(n * r, n * j) for j in range(r + 1)) - 2**(n * r)
        assert twice % 2 == 0
        powers.append(twice // 2)
    elementary = [1]
    for r in range(1, d + 1):
        numerator = sum((-1)**(i - 1) * elementary[r - i] * powers[i]
                        for i in range(1, r + 1))
        assert numerator % r == 0
        elementary.append(numerator // r)
    return [(-1)**(d - i) * elementary[d - i] for i in range(d + 1)]


def evaluate(coefficients, x, modulus=None):
    value = 0
    for c in reversed(coefficients):
        value = value * x + c
        if modulus is not None:
            value %= modulus
    return value


POLYNOMIALS = {n: polynomial(n) for n in (4, 8, 16, 32, 64)}


def multiset_orbit_audit(h, p, expected_extra):
    """Independent enumeration, capped at subgroup order 16."""
    n = len(h)
    assert n <= 16
    intrinsic, extra = 0, 0
    orbit_members = {}
    orbit_types = {}
    for values in combinations_with_replacement(sorted(h), 4):
        if sum(values) % p:
            continue
        counts = Counter(values)
        orderings = factorial(4) // prod(factorial(c) for c in counts.values())
        if all(counts[x] == counts[-x % p] for x in counts):
            intrinsic += orderings
            continue
        extra += orderings
        images = {tuple(sorted(a * x % p for x in values)) for a in h}
        assert len(images) == n
        representative = min(images)
        orbit_members.setdefault(representative, set()).add(values)
        orbit_types[representative] = tuple(sorted(counts.values(), reverse=True))
    assert intrinsic == 3 * n * n - 3 * n
    assert extra == expected_extra
    assert all(len(v) == n for v in orbit_members.values())
    histogram = Counter(orbit_types.values())
    assert set(histogram) <= {(3, 1), (2, 1, 1), (1, 1, 1, 1)}
    assert extra == n * (4 * histogram[(3, 1)] + 12 * histogram[(2, 1, 1)]
                         + 24 * histogram[(1, 1, 1, 1)])
    return {"intrinsic": intrinsic, "extra": extra,
            "free_orbits_by_multiplicity": {str(k): v for k, v in histogram.items()}}


def tower(p, n):
    assert prime(p) and (p - 1) % n == 0
    g = generator(p, n)
    rows = []
    previous_D = 0
    for j in range(1, n.bit_length()):
        s = 2**j
        gs = pow(g, n // s, p)
        h = [pow(gs, a, p) for a in range(s)]
        assert len(set(h)) == s and pow(gs, s // 2, p) == p - 1
        r_h = pair_counts(h, h, p)
        energy = square_norm(r_h)
        D = energy - (3 * s * s - 3 * s)
        assert D >= 0
        row = {"parent_order": s, "energy": energy, "D": D}
        if s == 2:
            assert D == 0
        else:
            k = s // 2
            child, other = h[::2], h[1::2]
            r_k = pair_counts(child, child, p)
            r_mixed = pair_counts(child, other, p)
            B = square_norm(r_mixed)
            T = sum(c * r_mixed[x] for x, c in r_k.items())
            primitive_values = [pow(1 + pow(gs, a, p), s, p)
                                for a in range(1, s // 2, 2)]
            fibers = Counter(primitive_values)
            C = sum(comb(c, 2) for c in fibers.values())
            assert B == s * sum(c * c for c in fibers.values())
            assert B == k * k + 2 * s * C
            assert D == 2 * previous_D + 6 * (B - k * k) + 8 * T
            assert D >= 12 * s * C
            boundary = pow(2, s, p)
            F_mod = prod((boundary - pow(1 + pow(gs, a, p), s, p)) % p
                         for a in range(1, s // 2)) % p
            h_set = set(h)
            repeated_pairs = sum((2 - a) % p in h_set for a in h)
            assert (repeated_pairs > 1) == (F_mod == 0)
            if s in POLYNOMIALS:
                assert evaluate(POLYNOMIALS[s], boundary, p) == F_mod
            epsilon = int(3 % p in h_set)
            assert (repeated_pairs - 1 - 2 * epsilon) % 2 == 0
            u = (repeated_pairs - 1 - 2 * epsilon) // 2
            assert u >= 0
            residual = D - s * (4 * epsilon + 12 * u)
            assert residual >= 0 and residual % (24 * s) == 0
            if F_mod != 0:
                assert D % (24 * s) == 0
                assert D == 0 or D >= 24 * s
            row.update({"B": B, "T": T, "primitive_collision_mass": C,
                        "primitive_fiber_histogram": dict(Counter(fibers.values())),
                        "F_mod_p": F_mod, "repeated_pair_count": repeated_pairs,
                        "epsilon": epsilon, "u_orbits": u,
                        "v_orbits": residual // (24 * s)})
        if s <= 16:
            row["independent_multiset_audit"] = multiset_orbit_audit(h, p, D)
        rows.append(row)
        previous_D = D
    total_C = sum(row.get("primitive_collision_mass", 0) for row in rows)
    weighted_T = sum((n // row["parent_order"]) * row.get("T", 0) for row in rows)
    assert rows[-1]["D"] == 12 * n * total_C + 8 * weighted_T
    eta = sum((Fraction(row["D"], row["parent_order"]**2) for row in rows),
              Fraction(0))
    assert sum((Fraction(row.get("primitive_collision_mass", 0), row["parent_order"])
                for row in rows), Fraction(0)) <= eta / 12
    assert total_C <= eta * n / 12
    if rows[-1]["D"] == 0:
        assert all(row["D"] == 0 for row in rows)
        assert all(row.get("T", 0) == row.get("primitive_collision_mass", 0) == 0
                   for row in rows)
    return {"p": p, "parent_order": n, "generator": g,
            "in_quartic_window": n**4 // 4 <= p <= n**4,
            "total_primitive_collision_mass": total_C,
            "normalized_tower_defect": str(eta), "levels": rows}


def norm_audit(case, tower_lookup):
    """Rebuild integer norms; use the existing bounded factor certificate."""
    n, d = case["n"], case["n"] // 2
    polynomial_counts = Counter()
    for exponents in product(range(n), repeat=3):
        coefficients = [0] * d
        for e in (0,) + exponents:
            coefficients[e % d] += 1 if e < d else -1
        polynomial_counts[tuple(coefficients)] += 1
    intrinsic = n * polynomial_counts[(0,) * d]
    assert intrinsic == 3 * n * n - 3 * n
    histogram = Counter()
    for coefficients, multiplicity in polynomial_counts.items():
        if not any(coefficients):
            continue
        matrix = [[0] * d for _ in range(d)]
        for j in range(d):
            for i, c in enumerate(coefficients):
                matrix[(i + j) % d][j] += c * (-1 if i + j >= d else 1)
        norm = abs(determinant(matrix))
        assert norm > 0
        histogram[norm] += multiplicity
    assert histogram == Counter({int(k): v for k, v in case["nonzero_norm_histogram"].items()})
    R = (n**4 - intrinsic) // n
    assert sum(histogram.values()) == R
    norm_product = prod(norm**multiplicity for norm, multiplicity in histogram.items())
    factored_product = prod(int(p)**v for p, v in case["norm_product_prime_factorization"].items())
    assert norm_product == factored_product
    # Denominator-cleared AM-GM: this is an exact integer comparison.
    assert norm_product**2 * R**(d * R) <= (4 * n**3)**(d * R)
    required_product = 1
    for row in case["splitting_primes_with_extra_relations"]:
        p = row["p"]
        D = tower_lookup[(p, n)]["levels"][-1]["D"]
        assert D == row["extra_zero_quadruples"] and D % 2 == 0
        required_product *= p**(D // 2)
    assert norm_product % required_product == 0
    return {"n": n, "normalized_tuples": n**3, "nonintrinsic_normalized_tuples": R,
            "integer_AM_GM_passed": True, "required_prime_product_divides": True,
            "norm_product_bits": norm_product.bit_length()}


def main():
    towers = [tower(p, n) for p, n in CERTIFICATES]
    lookup = {(row["p"], row["parent_order"]): row for row in towers}
    repeated_witness = lookup[(6700417, 64)]["levels"][-1]
    assert repeated_witness["D"] == 768 == 12 * 64
    assert repeated_witness["F_mod_p"] == 0 and repeated_witness["u_orbits"] == 1
    assert repeated_witness["v_orbits"] == 0
    distinct_witness = lookup[(17189277697, 512)]["levels"][-1]
    assert distinct_witness["D"] == 12288 == 24 * 512
    assert distinct_witness["F_mod_p"] != 0 and distinct_witness["v_orbits"] == 1
    old_norm = json.loads(PRIOR.read_text())
    norm_cases = [norm_audit(case, lookup) for case in old_norm["cases"]]
    boundaries = []
    for n, coefficients in POLYNOMIALS.items():
        value = evaluate(coefficients, 2**n)
        height_exponent = (n + 1) * (n // 2 - 1)
        assert 0 < value < 2**height_exponent
        boundaries.append({"n": n, "F_hex": hex(value), "F_bits": value.bit_length(),
                           "strict_binary_height_exponent": height_exponent})
    epsilon_analytic = Fraction(1, 12)
    theta = Fraction(7, 12)
    assert theta + epsilon_analytic == Fraction(2, 3)
    assert 4 * (theta + epsilon_analytic) == Fraction(8, 3) < 3
    margins = []
    for m in range(2, 33):
        n = 2**m
        assert sum(4**j for j in range(1, m + 1)) == (4 * n * n - 4) // 3
        # (3N^3/2)^3 >= (N^(8/3))^3, using integers only.
        assert 27 * n**9 >= 8 * n**8
        assert n**4 // 4 > 2 and (n**4 // 4) & (n**4 // 4 - 1) == 0
        margins.append({"m": m, "all_level_square_sum": (4 * n * n - 4) // 3})
    # Leading constant ratios; logarithms only decorate these exact fractions.
    assert Fraction(1, 24) * Fraction(1, 4) / Fraction(3, 8) == Fraction(1, 36)
    assert Fraction(1, 12) * Fraction(1, 4) / Fraction(3, 8) == Fraction(1, 18)
    sources = [Path(__file__).resolve(), NOTE, PRIOR,
               ROOT / "experiments/cyclotomic_norm_audit.py",
               ROOT / "research/cyclotomic-prime-average.md",
               ROOT / "research/parallel2-subgroup-2026-09-04.md"]
    result = {
        "status": "passed",
        "scope": "Exact finite identities checked; asymptotic proof uses Thorner--Zaman Corollary 3.1/(3.2). No prime scan or factorization search.",
        "sources_sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                           for path in sources},
        "analytic_source": "https://arxiv.org/html/2108.10878#S3.SS1",
        "analytic_specialization": {"q": "N=2^m", "rad_q": 2, "a": 1, "x": "N^4",
                                    "h": "3N^4/4", "fixed_epsilon": "1/12",
                                    "prime_count_leading_constant": "3/8",
                                    "exact_good_proportion_lower_bound": "1-log(2)/36",
                                    "exact_good_proportion_decimal": 1 - log(2) / 36,
                                    "collision_mean_upper_bound": "log(2)/18",
                                    "finite_threshold_computed": False},
        "prime_certificates": towers, "small_norm_audits": norm_cases,
        "boundary_products": boundaries,
        "rational_exponent_and_dyadic_sum_checks": margins,
    }
    output = ROOT / "results/parallel4_subgroup_2026_09_04.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "fixed_prime_certificates": len(towers),
                      "dyadic_levels": sum(len(row["levels"]) for row in towers),
                      "small_norm_orders": [row["n"] for row in norm_cases],
                      "repeated_quartic_witness": repeated_witness,
                      "four_distinct_quartic_witness": distinct_witness,
                      "sources_sha256": result["sources_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
