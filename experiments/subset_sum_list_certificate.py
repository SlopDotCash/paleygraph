#!/usr/bin/env python3
"""Exact subset-sum and RS list certificates; standard library only.

Proofs and limitations: research/subset-sums-and-lists.md.
This is not a proof of either Paley conjecture or a Lean prize submission.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, isqrt
from pathlib import Path
import json

from mca_circuit_certificate import evaluate, mul, trim
from cyclotomic_norm_audit import determinant, multiplication_matrix


ROOT = Path(__file__).resolve().parents[1]


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def integer_record(value):
    magnitude = abs(value)
    result = {"sign": (value > 0) - (value < 0), "bits": magnitude.bit_length(),
              "magnitude_bytes_sha256": sha256(magnitude.to_bytes(
                  max(1, (magnitude.bit_length() + 7) // 8), "big")).hexdigest()}
    if magnitude.bit_length() < 200:
        result["value"] = value
    return result


def subgroup(p, n, generator):
    assert n > 1 and n & (n - 1) == 0 and (p - 1) % n == 0
    assert pow(generator, n, p) == 1 and pow(generator, n // 2, p) == p - 1
    domain, x = [], 1
    for _ in range(n):
        domain.append(x)
        x = x * generator % p
    assert x == 1 and len(set(domain)) == n
    return sorted(domain)


def count_subsets(domain, p):
    """All prefix counts; row j records sums of j DISTINCT chosen elements."""
    n = len(domain)
    initial = [[0] * p for _ in range(n + 1)]
    initial[0][0] = 1
    prefixes = [initial]
    for i, x in enumerate(domain):
        previous = prefixes[-1]
        current = [row[:] for row in previous]
        for j in range(1, i + 2):
            for c in range(p):
                current[j][c] += previous[j - 1][(c - x) % p]
        assert sum(current[j][c] for j in range(n + 1) for c in range(p)) == 2 ** (i + 1)
        prefixes.append(current)
    return prefixes


def recover_subset(domain, p, prefixes, size, total):
    assert prefixes[-1][size][total] > 0
    selected, j, c = [], size, total
    for i in range(len(domain), 0, -1):
        if prefixes[i - 1][j][c]:
            continue
        x = domain[i - 1]
        assert j > 0 and prefixes[i - 1][j - 1][(c - x) % p] > 0
        selected.append(x)
        j, c = j - 1, (c - x) % p
    assert j == c == 0 and len(selected) == size and sum(selected) % p == total
    return sorted(selected)


def vanishing(roots, p, a=1):
    poly = [1]
    for root in roots:
        poly = mul(poly, [-root % p] + [0] * (a - 1) + [1], p)
    return poly


def subset_bound(p, n, s, upper):
    assert 0 <= s <= n < p and upper >= 1
    t = min(s, n - s)
    total, error = comb(n, s), comb(upper + t - 1, t)
    numerator = total - (p - 1) * error
    lower = max(0, -(-numerator // p))
    return lower, {"n": n, "s": s, "t": t, "M_upper": upper,
                   "binomial_total": integer_record(total),
                   "binomial_error": integer_record(error),
                   "lower_numerator": integer_record(numerator),
                   "subset_count_lower": integer_record(lower),
                   "positive": numerator > 0,
                   "lower_log2_floor": lower.bit_length() - 1}


def small_case(p, n, generator, upper):
    assert all(p % d for d in range(2, isqrt(p) + 1))
    domain = subgroup(p, n, generator)
    # Full multiplicative group: all nonprincipal periods are exactly -1.
    # Otherwise use the classical Gauss-sum bound M <= sqrt(p).
    assert (n == p - 1 and upper == 1) or upper * upper >= p
    prefixes = count_subsets(domain, p)
    for size in range(n + 1):
        t = min(size, n - size)
        error = (p - 1) * comb(upper + t - 1, t)
        assert sum(prefixes[-1][size]) == comb(n, size)
        assert all(abs(p * count - comb(n, size)) <= error
                   for count in prefixes[-1][size])
        assert all(prefixes[-1][size][c] == prefixes[-1][n - size][(-c) % p]
                   for c in range(p))  # subgroup has sum zero
    k, s = n // 2, n // 2 + 1
    lower, bound = subset_bound(p, n, s, upper)
    assert lower > 0 and min(prefixes[-1][s]) >= lower
    witnesses = []
    for c in range(p):
        selected = recover_subset(domain, p, prefixes, s, c)
        qpoly = vanishing(selected, p)
        gamma = -c % p
        assert qpoly[-1] == 1 and qpoly[k] == gamma
        folded = [0] * (s + 1)
        folded[k], folded[s] = gamma, 1
        codeword = trim([(a - b) % p for a, b in zip(folded, qpoly)])
        assert len(codeword) <= k
        agreement = [x for x in domain
                     if evaluate(folded, x, p) == evaluate(codeword, x, p)]
        assert agreement == selected
        witnesses.append({"sum": c, "gamma": gamma, "roots": selected,
                          "codeword": codeword})
    # Independently enumerate all 11,440 subsets in the smaller example.
    exhaustive = None
    if p == 17:
        census = Counter(sum(xs) % p for xs in combinations(domain, s))
        assert [census[c] for c in range(p)] == prefixes[-1][s]
        zero_words = set()
        for xs in combinations(domain, s):
            if sum(xs) % p == 0:
                qpoly = vanishing(xs, p)
                codeword = tuple((-v) % p for v in qpoly[:k])
                assert qpoly[k] == 0
                zero_words.add(codeword)
        assert len(zero_words) == census[0]
        exhaustive = {"subsets": comb(n, s), "distinct_zero_sum_codewords": len(zero_words)}
    return {"p": p, "n": n, "k": k, "bound": bound,
            "zero_sum_count": prefixes[-1][s][0],
            "count_histogram": sorted(Counter(prefixes[-1][s]).items()),
            "all_sizes_and_sums_checked": (n + 1) * p,
            "mca_witnesses": witnesses, "exhaustive_check": exhaustive}


def small_root_lifts():
    p, n, k = 97, 32, 16
    domain = subgroup(p, n, 28)
    records = []
    for a in [1, 2, 4, 8]:
        image = sorted({pow(x, a, p) for x in domain})
        image_generator = pow(28, a, p)
        exponent = {pow(image_generator, j, p): j for j in range(len(image))}
        assert len(image) == n // a
        assert all(sum(pow(x, a, p) == z for x in domain) == a for z in image)
        size = k // a + 1
        prefixes = count_subsets(image, p)
        count = prefixes[-1][size][0]
        # One reconstructed exact witness suffices for a=1, whose exhaustive
        # subset space is large. Exhaust every image subset for a>1.
        candidates = ([recover_subset(image, p, prefixes, size, 0)] if count else [])
        if a > 1:
            candidates = [xs for xs in combinations(image, size) if sum(xs) % p == 0]
            assert len(candidates) == count
        words, norms = set(), []
        for xs in candidates:
            qpoly = vanishing(xs, p, a)
            assert len(qpoly) == k + a + 1 and qpoly[k] == 0
            vpoly = [0] * (k + a) + [1]
            upoly = trim([(x - y) % p for x, y in zip(vpoly, qpoly)])
            assert len(upoly) <= k and len(upoly) - 1 <= k - a
            roots = [x for x in domain if evaluate(qpoly, x, p) == 0]
            assert len(roots) == k + a
            assert set(roots) == {x for x in domain if pow(x, a, p) in xs}
            assert all(evaluate(vpoly, x, p) == evaluate(upoly, x, p) for x in roots)
            words.add(tuple(upoly))
            # Complementation gives an odd subset with t=N/2-1 elements.
            # Its nonzero cyclotomic norm must be divisible by p.
            if a > 1:
                complement = set(image) - set(xs)
                degree = len(image) // 2
                assert len(complement) % 2 == 1 and sum(complement) % p == 0
                coeffs = [0] * degree
                for x in complement:
                    j = exponent[x]
                    coeffs[j % degree] += 1 if j < degree else -1
                norm = abs(determinant(multiplication_matrix(coeffs)))
                assert 0 < norm <= len(complement) ** degree and norm % p == 0
                norms.append(norm)
        assert len(words) == len(candidates)
        records.append({"a": a, "image_size": len(image), "zero_sum_subsets": count,
                        "polynomials_checked": len(words), "all_candidates_checked": a > 1,
                        "complement_norm_histogram": sorted(Counter(norms).items()),
                        "codewords_sha256": digest(sorted(words))})
    assert any(r["a"] > 1 and r["zero_sum_subsets"] > 0 for r in records)
    return records


def energy_certificate(p, domain, independent=False):
    n = len(domain)
    assert p - 1 in domain
    kernel = Counter(pow((1 + x) % p, n, p) for x in domain if x != p - 1)
    assert sum(kernel.values()) == n - 1 and 0 not in kernel
    e2 = n * n + n * sum(v * v for v in kernel.values())
    q2, residue = divmod(p * e2 - n ** 4, n)
    assert residue == 0 and q2 > 0
    upper = isqrt(isqrt(q2))
    upper += upper ** 4 < q2
    assert (upper - 1) ** 4 < q2 <= upper ** 4
    # Independent normalized quadruple count: divide h1+h2=h3+h4 by h1.
    normalized = None
    if independent:
        members = set(domain)
        normalized = sum((1 + x - y) % p in members for x in domain for y in domain)
        assert n * normalized == e2
    return upper, {"energy_2": e2, "Q_2": q2, "M_upper_fourth_root": upper,
                   "kernel_histogram": sorted(Counter(kernel.values()).items()),
                   "kernel_sha256": digest(sorted(kernel.items())),
                   "independent_normalized_quadruples": normalized}


def official_profile():
    p, n, k = 2130706433, 262144, 131072
    assert all(p % d for d in range(2, isqrt(p) + 1))
    q = p ** 6
    full_generator = pow(3, (p - 1) // n, p)
    assert full_generator == 0x6C4A8A45
    records = []
    for a in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
        size = n // a
        generator = pow(full_generator, a, p)
        domain = subgroup(p, size, generator)
        upper, energy = energy_certificate(p, domain, independent=(a == 32))
        lower, bound = subset_bound(p, size, k // a + 1, upper)
        radius = Fraction(n - k - a, n)
        records.append({"a": a, "generator": generator, "energy": energy,
                        "subset_bound": bound, "radius": [radius.numerator, radius.denominator],
                        "lower_exceeds_q": lower > q})
        print(f"a={a}: M<={upper}, positive={lower > 0}, lower log2 floor={lower.bit_length()-1}", flush=True)
    selected = next(r for r in records if r["a"] == 32)
    assert selected["radius"] == [4095, 8192]
    assert selected["subset_bound"]["lower_log2_floor"] == 8154
    assert selected["lower_exceeds_q"]
    assert 2 ** 128 < p ** 5  # monomial pair bad-scalar probability p/q < 2^-128
    return {"p": p, "q": q, "n": n, "scalar_dimension": k,
            "records": records, "selected_a": 32,
            "monomial_pair_bad_scalars_at_radius_1_2_minus_1_n": p,
            "monomial_pair_mca_below_2_neg_128": True,
            "meaning": "List lower bound excludes certifiedGammaError<=2^-128 for delta>=4095/8192. No winning-set lower bound or worst-case MCA upper bound is claimed."}


def main():
    archive = ROOT / "sources/official-prize-2026-09-04"
    manifest = json.loads((archive / "manifest.json").read_text())
    archived = manifest["files"] + [f for d in manifest["dependencies"] for f in d["files"]]
    for item in archived:
        assert sha256((archive / item["path"]).read_bytes()).hexdigest() == item["sha256"]
    result = {"status": "Exact finite certificates; full goal unproved; no novelty or Lean claim",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "polynomial_helper_sha256": sha256((ROOT / "experiments/mca_circuit_certificate.py").read_bytes()).hexdigest(),
              "norm_helper_sha256": sha256((ROOT / "experiments/cyclotomic_norm_audit.py").read_bytes()).hexdigest(),
              "official_manifest_sha256": sha256((archive / "manifest.json").read_bytes()).hexdigest(),
              "official_archive_files_verified": len(archived),
              "small_cases": [small_case(17, 16, 3, 1), small_case(97, 32, 28, 10)],
              "small_root_lifts": small_root_lifts(),
              "official_profile": official_profile()}
    destination = ROOT / "results/subset_sum_list_certificate.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(f"PASS: saved {destination}")


if __name__ == "__main__":
    main()
