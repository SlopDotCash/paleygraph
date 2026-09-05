"""Bounded exact audits for centered sixth moments and fixed amplification.

No prime scan, numerical asymptotic inference, or integer factorization.
The trilinear estimate is a stated external theorem, not tested here.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import factorial, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "research/parallel5-subgroup-2026-09-04.md"
FIXED_CASES = [
    (5, 4, 4), (13, 4, 4), (17, 4, 4), (17, 8, 4), (41, 8, 4),
    (73, 4, 4), (1153, 8, 4), (2017, 8, 4), (18433, 16, 4),
    (65521, 16, 4), (262657, 32, 4), (6700417, 64, 3),
]


def prime(p):
    return p >= 2 and (p == 2 or (p % 2 == 1 and all(
        p % d for d in range(3, isqrt(p) + 1, 2))))


def generator(p, n):
    for a in range(2, min(p, 1002)):
        g = pow(a, (p - 1) // n, p)
        if pow(g, n // 2, p) == p - 1:
            return g
    raise AssertionError("Bounded generator search exhausted")


def intrinsic(n, r):
    def multiply(a, b):
        return [sum((a[i] * b[j - i] for i in range(j + 1)), Fraction(0))
                for j in range(r + 1)]
    coefficients = [Fraction(1)] + [Fraction(0)] * r
    power = [Fraction(1, factorial(j)**2) for j in range(r + 1)]
    exponent = n // 2
    while exponent:
        if exponent & 1:
            coefficients = multiply(coefficients, power)
        power = multiply(power, power)
        exponent //= 2
    result = factorial(2 * r) * coefficients[r]
    assert result.denominator == 1
    return result.numerator


def convolve_step(counts, h, p):
    result = Counter()
    for x, c in counts.items():
        for y in h:
            result[(x + y) % p] += c
    return result


def fixed_moments(p, n, depth):
    assert prime(p) and (p - 1) % n == 0
    g = generator(p, n)
    h = [pow(g, j, p) for j in range(n)]
    assert len(set(h)) == n
    convolution = [Counter({0: 1})]
    moments = []
    for r in range(1, depth + 1):
        counts = convolve_step(convolution[-1], h, p)
        convolution.append(counts)
        assert sum(counts.values()) == n**r
        assert all(counts[x] == counts[-x % p] == counts[g * x % p]
                   for x in counts)
        energy = sum(c * c for c in counts.values())
        centered_numerator = p * energy - n**(2 * r)
        assert centered_numerator >= 0 and centered_numerator % n == 0
        # Exactly p^2 times the variance of the r-step count, including zeros.
        variance_numerator = sum((p * c - n**r)**2 for c in counts.values())
        variance_numerator += (p - len(counts)) * n**(2 * r)
        assert variance_numerator == p * centered_numerator
        # First coordinate normalized to 1; remaining 2r-1 split r-1 and r.
        normalized = sum(c * counts[(-1 - x) % p]
                         for x, c in convolution[r - 1].items())
        assert n * normalized == energy
        tr = intrinsic(n, r)
        assert energy >= tr
        moments.append({"r": r, "E_r": energy, "intrinsic": tr,
                        "extra_relations": energy - tr,
                        "principal_term": str(Fraction(n**(2 * r), p)),
                        "centered_energy": str(Fraction(centered_numerator, p)),
                        "centered_numerator": centered_numerator,
                        "coset_moment_Q_r": centered_numerator // n,
                        "sum_support_size": len(counts)})
    for i in range(1, len(moments) - 1):
        q0, q1, q2 = (moments[j]["coset_moment_Q_r"] for j in (i - 1, i, i + 1))
        assert q1 * q1 <= q0 * q2
    witness = None
    if (p, n) == (262657, 32):
        assert moments[1]["E_r"] == 2976 == intrinsic(n, 2)
        assert moments[2]["E_r"] == 458240
        assert moments[2]["intrinsic"] == 446720
        by_sum = defaultdict(list)
        exponent = {x: j for j, x in enumerate(h)}
        for triple in combinations_with_replacement(sorted(h), 3):
            value = sum(triple) % p
            for previous in by_sum[value]:
                word = triple + tuple(-x % p for x in previous)
                counts = Counter(word)
                if any(counts[x] != counts[-x % p] for x in counts):
                    assert sum(word) % p == 0
                    assert all((-x % p) not in counts for x in counts)
                    witness = {"values": word,
                               "exponents": [exponent[x] for x in word],
                               "sum_mod_p": 0, "has_opposite_pair": False,
                               "multiplicities": sorted(counts.values(), reverse=True)}
                    break
            if witness is not None:
                break
            by_sum[value].append(triple)
        assert witness is not None
    return {"p": p, "n": n, "generator": g, "depth": depth,
            "in_quartic_window": n**4 // 4 <= p <= n**4,
            "moments": moments, "nonintrinsic_six_word": witness}


def determinant(matrix):
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
            for k in range(j + 1, d):
                value = pivot * a[i][k] - a[i][j] * a[j][k]
                assert value % previous == 0
                a[i][k] = value // previous
        for i in range(j + 1, d):
            a[i][j] = 0
        previous = pivot
    return sign * a[-1][-1]


def norm_budget_certificate(n, order, cases):
    """All normalized tuples, capped at 32768; no factorization."""
    d, r = n // 2, order // 2
    assert n**(order - 1) <= 32768
    coefficient_counts = Counter()
    for exponents in product(range(n), repeat=order - 1):
        coefficients = [0] * d
        for a in (0,) + exponents:
            coefficients[a % d] += 1 if a < d else -1
        coefficient_counts[tuple(coefficients)] += 1
    T = intrinsic(n, r)
    assert n * coefficient_counts[(0,) * d] == T
    histogram = Counter()
    nonintrinsic_polynomials = []
    for coefficients, multiplicity in coefficient_counts.items():
        if not any(coefficients):
            continue
        matrix = [[0] * d for _ in range(d)]
        for j in range(d):
            for i, c in enumerate(coefficients):
                matrix[(i + j) % d][j] += c * (-1 if i + j >= d else 1)
        norm = abs(determinant(matrix))
        assert 0 < norm <= order**d
        histogram[norm] += multiplicity
        nonintrinsic_polynomials.append((coefficients, multiplicity))
    R = (n**order - T) // n
    assert sum(histogram.values()) == R
    norm_product = prod(norm**multiplicity for norm, multiplicity in histogram.items())
    # Squared AM--GM with no logarithmic or floating-point comparisons.
    assert norm_product**2 * R**(d * R) <= (order * n**(order - 1))**(d * R)
    required_product, checked_primes = 1, []
    for case in cases:
        if case["n"] != n or case["depth"] < r:
            continue
        p, g = case["p"], case["generator"]
        roots = [pow(g, j, p) for j in range(1, n, 2)]
        nullities = sum(multiplicity * sum(
            sum(c * pow(root, j, p) for j, c in enumerate(coefficients)) % p == 0
            for root in roots) for coefficients, multiplicity in nonintrinsic_polynomials)
        extra = case["moments"][r - 1]["extra_relations"]
        assert n * nullities == d * extra
        required_product *= p**nullities
        checked_primes.append(p)
    assert norm_product % required_product == 0
    return {"n": n, "word_order": order, "normalized_tuples": n**(order - 1),
            "intrinsic_count": T, "nonintrinsic_normalized_tuples": R,
            "norm_product_bits": norm_product.bit_length(),
            "AM_GM_integer_comparison_passed": True,
            "prime_nullity_and_divisibility_checks": checked_primes}


def select_fractional(values):
    """Audit the elementary dyadic selection lemma using exact fractions."""
    a = sum(values, Fraction(0)) / len(values)
    assert 0 < a <= 1 and all(0 <= f <= 1 for f in values)
    ell, power = 1, Fraction(1)
    while power < 2 / a:
        power *= 2
        ell += 1
    groups = defaultdict(list)
    for f in values:
        if f < a / 2:
            continue
        t = Fraction(1)
        while f <= t / 2:
            t /= 2
        groups[t].append(f)
    assert len(groups) <= ell
    t, selected = max(groups.items(), key=lambda kv: sum(kv[1]))
    assert t >= a / 2
    assert all(t / 2 < f <= t for f in selected)
    assert len(selected) >= a * len(values) / (2 * ell * t)
    return {"mean": str(a), "dyadic_level": str(t), "input_count": len(values),
            "selected_count": len(selected), "interval_count_bound": ell}


def ledger():
    saving = Fraction(1, 24)
    supports = {"X": 3 - 6 * saving, "Y": 3 - 18 * saving, "Z": 2 - 36 * saving}
    assert supports == {"X": Fraction(11, 4), "Y": Fraction(9, 4), "Z": Fraction(1, 2)}
    assert all(v > 0 for v in supports.values())
    # The three support estimates contribute L^-2, L^-2, and L^-1.
    assert 2 + 2 + 1 == 5
    # Delta^6 Delta_1^4 Delta_3^3; Delta_3>=c Delta_2^2;
    # Delta_2>=c Delta_1^3; Delta_1>=c Delta^3.
    delta_power = 6 + 3 * (4 + 3 * (2 * 3))
    assert delta_power == 72
    n_power = Fraction(1) + Fraction(4 - 7, 72)
    log_power = Fraction(5, 72) + Fraction(1, 36)
    assert n_power == Fraction(23, 24) < Fraction(2849, 2880)
    assert log_power == Fraction(7, 72)
    higher = []
    for r in range(3, 17):
        exponent = 3 + Fraction(23, 24) * (2 * r - 6)
        logs = 1 + Fraction(7, 72) * (2 * r - 6)
        assert exponent == Fraction(23 * r - 33, 12)
        assert logs == 1 + Fraction(7 * (r - 3), 36)
        higher.append({"r": r, "centered_exponent": str(exponent),
                       "log_exponent": str(logs), "raw_norm_average_exponent": 2 * r - 3})
    assert higher[1]["centered_exponent"] == "59/12"
    assert higher[1]["log_exponent"] == "43/36"
    for m in range(1, 33):
        n = 2**m
        assert sum(8**j for j in range(1, m + 1)) == (8 * n**3 - 8) // 7
        for r in (2, 3, 4):
            T = intrinsic(n, r)
            if r == 2:
                assert T == 3 * n * n - 3 * n
            elif r == 3:
                assert T == 15 * n**3 - 45 * n * n + 40 * n
            else:
                assert T == 105 * n**4 - 630 * n**3 + 1435 * n * n - 1155 * n
    return {"deletion_support_exponents": {k: str(v) for k, v in supports.items()},
            "delta_power": delta_power, "log_denominator_power": 5,
            "amplitude_exponent": str(n_power), "amplitude_log_exponent": str(log_power),
            "higher_centered_moment_bounds": higher,
            "all_level_cube_sum_checked_through_m": 32}


def main():
    cases = [fixed_moments(*case) for case in FIXED_CASES]
    norms = [norm_budget_certificate(n, order, cases)
             for n, order in [(4, 6), (8, 6), (4, 8)]]
    selection_examples = [
        [Fraction(0), Fraction(1)],
        [Fraction(j, 16) for j in range(17)],
        [Fraction(1, 2**j) for j in range(40)],
        [Fraction(0)] * 99 + [Fraction(1, 2**80)],
        [Fraction(1, 3)] * 7 + [Fraction(2, 5)] * 9,
    ]
    selections = [select_fractional(values) for values in selection_examples]
    proof_ledger = ledger()
    dependencies = [
        Path(__file__).resolve(), NOTE,
        ROOT / "research/parallel4-subgroup-2026-09-04.md",
        ROOT / "research/cyclotomic-prime-average.md",
        ROOT / "research/analytic-bounds-and-amplification.md",
        ROOT / "sources/analytic-bounds-2026-09-04/di-benedetto-et-al-2003.06165v1.html",
    ]
    result = {
        "status": "passed",
        "scope": "Exact finite moment/norm identities and selection algebra. Analytic inputs are cited theorems; no prime scan or factorization.",
        "sources_sha256": {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                           for path in dependencies},
        "analytic_sources": [
            "https://arxiv.org/html/2108.10878#S3.SS1",
            "https://arxiv.org/html/2003.06165#S4",
            "https://arxiv.org/html/1604.08469v4#S1.SS3"],
        "fixed_moment_certificates": cases, "norm_budget_certificates": norms,
        "rational_selection_certificates": selections, "exponent_ledger": proof_ledger,
    }
    output = ROOT / "results/parallel5_subgroup_2026_09_04.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    counterexample = next(case for case in cases if case["p"] == 262657)
    print(json.dumps({
        "status": "passed", "fixed_prime_certificates": len(cases),
        "moment_levels": sum(len(case["moments"]) for case in cases),
        "bounded_norm_word_orders": [(row["n"], row["word_order"]) for row in norms],
        "counterexample": counterexample,
        "maximum_period_exponent": proof_ledger["amplitude_exponent"],
        "maximum_period_log_exponent": proof_ledger["amplitude_log_exponent"],
        "sources_sha256": result["sources_sha256"],
    }, indent=2))


if __name__ == "__main__":
    main()
