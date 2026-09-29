#!/usr/bin/env python3
"""Exact finite checks for the Paley research notes; not an asymptotic proof.

Uses Python integers only. No numerical optimization, floating point thresholds,
or probabilistic primality tests participate in the verification.
"""

import argparse
import itertools
import json
from fractions import Fraction
from math import comb, isqrt, prod
from pathlib import Path


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def character(p):
    if not prime(p) or p == 2:
        raise ValueError("An odd prime is required")
    squares = {x * x % p for x in range(1, p)}
    return [0 if x == 0 else 1 if x in squares else -1 for x in range(p)]


def row_sums(p, chi, subset):
    return [sum(chi[(x - b) % p] for b in subset) for x in range(p)]


def moment(rows, r):
    return sum(x ** (2 * r) for x in rows)


def quartic_trace(p, chi, quad):
    return sum(prod(chi[(x - b) % p] for b in quad) for x in range(p))


def verify():
    cases = 0
    quartic_cases = 0
    gram_entries = 0
    for p, max_size in [(3, 3), (5, 5), (7, 7), (11, 11), (13, 13), (17, 5)]:
        chi = character(p)
        for a in range(p):
            for b in range(p):
                observed = sum(chi[(x - a) % p] * chi[(x - b) % p]
                               for x in range(p))
                assert observed == (p - 1 if a == b else -1)
                gram_entries += 1
        traces = {q: quartic_trace(p, chi, q)
                  for q in itertools.combinations(range(p), 4)}
        for n in range(max_size + 1):
            for subset in itertools.combinations(range(p), n):
                rows = row_sums(p, chi, subset)
                assert sum(rows) == 0
                assert moment(rows, 1) == p * n - n * n
                local_energy = sum(rows[a] ** 2 for a in subset)
                trace_sum = sum(traces[q] for q in itertools.combinations(subset, 4))
                rhs = (p * (3 * n * n - 2 * n) - 6 * n ** 3
                       + 14 * n * n - 9 * n - 6 * local_energy + 24 * trace_sum)
                assert moment(rows, 2) == rhs, (p, subset)
                # Direct tuple expansion independently checks the partition formula.
                if p <= 7 and n <= 4:
                    expanded = sum(
                        sum(prod(chi[(x - b) % p] for b in t) for x in range(p))
                        for t in itertools.product(subset, repeat=4))
                    assert expanded == rhs
                    quartic_cases += 1
                # All cardinalities of A: sorted rows optimize its signed sum exactly.
                ordered = sorted(rows)
                for m in range(1, p + 1):
                    best = max(abs(sum(ordered[:m])), abs(sum(ordered[-m:])))
                    # Centered spectral inequality, squared and cleared of denominators.
                    assert p * best * best <= m * (p - m) * n * (p - n)
                cases += 1
    return {"status": "passed", "subset_cases": cases,
            "gram_entries": gram_entries,
            "independent_direct_fourth_moment_expansions": quartic_cases,
            "scope": "Finite identities and baseline inequalities only"}


def extremum(p, m, n):
    chi = character(p)
    best = -1
    witness = None
    count = 0
    for subset in itertools.combinations(range(p), n):
        rows = row_sums(p, chi, subset)
        ordered = sorted(range(p), key=lambda x: (rows[x], x))
        for chosen in (ordered[:m], ordered[-m:]):
            value = sum(rows[a] for a in chosen)
            if abs(value) > best:
                best = abs(value)
                witness = {"A": sorted(chosen), "B": list(subset), "signed_sum": value}
        count += 1
    assert count == comb(p, n)
    assert abs(sum(chi[(a - b) % p] for a in witness["A"]
                   for b in witness["B"])) == best
    # For the smallest instances, independently enumerate both A and B.
    if p <= 7:
        brute = max(abs(sum(chi[(a - b) % p] for a in a_set for b in b_set))
                    for a_set in itertools.combinations(range(p), m)
                    for b_set in itertools.combinations(range(p), n))
        assert brute == best
    return {"p": p, "m": m, "n": n, "maximum_absolute_sum": best,
            "maximum_bias": str(Fraction(best, m * n)), "witness": witness,
            "B_subsets_exhausted": count,
            "all_A_covered_by_exact_sorting": True}


def spike_examples():
    output = []
    for p, n, r in [(101, 10, 4), (1009, 30, 4), (10007, 40, 4)]:
        chi = character(p)
        subset = sorted((-q) % p for q in range(1, p) if chi[q] == 1)[:n]
        rows = row_sums(p, chi, subset)
        assert rows[0] == n
        pairing_constant = prod(range(1, 2 * r, 2))
        gaussian_bound = pairing_constant * p * n ** r
        actual = moment(rows, r)
        output.append({"p": p, "n": n, "r": r, "B": subset,
                       "F_B_at_zero": rows[0], "moment": actual,
                       "one_row_lower_bound": n ** (2 * r),
                       "gaussian_pairing_bound": gaussian_bound,
                       "gaussian_pairing_bound_fails": actual > gaussian_bound,
                       "n_below_sqrt_p": n * n < p})
    for p, k, r in [(1009, 3, 2), (10007, 5, 2), (65537, 8, 2)]:
        chi = character(p)
        prescribed = list(range(k))
        subset = [b for b in range(p)
                  if all(chi[(a - b) % p] == 1 for a in prescribed)]
        n = len(subset)
        assert n > 0
        lower = k * n ** (2 * r)
        candidate = p * n ** r + n ** (2 * r)
        output.append({"p": p, "k": k, "n": n, "r": r,
                       "prescribed_rows": prescribed, "B": subset,
                       "moment_lower_bound": lower,
                       "single_spike_corrected_bound_C1": candidate,
                       "lower_bound_already_refutes_C1": lower > candidate,
                       "required_constant_at_least": str(Fraction(lower, candidate))})
    return output


def lean_witness_crosscheck():
    p = 1009
    subset = sorted({(-j * j) % p for j in range(1, 31)})
    rows = row_sums(p, character(p), subset)
    assert len(subset) == 30 and rows[0] == 30
    lower = 30 ** 8
    bound = 105 * p * 30 ** 4
    actual = moment(rows, 4)
    assert actual >= lower > bound
    return {"p": p, "B": subset, "n": 30, "r": 4,
            "F_B_at_zero": rows[0], "moment": actual,
            "gaussian_pairing_bound": bound, "one_row_lower_bound": lower}


def common_neighbor_counting_checks():
    cases = 0
    for p in [5, 7, 11, 13, 17, 19]:
        chi = character(p)
        d = (p - 1) // 2
        for k in range(1, min(d, 4) + 1):
            total = 0
            largest = 0
            for rows in itertools.combinations(range(p), k):
                count = sum(all(chi[(a - b) % p] == 1 for a in rows)
                            for b in range(p))
                total += count
                largest = max(largest, count)
            assert total == p * comb(d, k)
            assert largest * comb(p, k) >= total
            ratio = Fraction(comb(d, k), comb(p, k))
            lower = Fraction(1, 2 ** k) * (1 - Fraction(k * k, p - k + 1))
            assert ratio >= lower
            cases += 1
    return {"status": "passed", "p_k_cases": cases,
            "identity": "sum_U |N(U)| = p binom((p-1)/2,k)"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["verify", "extrema", "spikes", "all"])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = {"asymptotic_conjecture_proved": False}
    if args.mode in ("verify", "all"):
        result["verification"] = verify()
        result["common_neighbor_counting"] = common_neighbor_counting_checks()
    if args.mode in ("extrema", "all"):
        result["extrema"] = [extremum(p, m, n) for p, m, n in
                             [(5, 2, 2), (7, 3, 3), (13, 3, 3), (17, 4, 4),
                              (29, 4, 4), (29, 6, 4), (37, 5, 4)]]
    if args.mode in ("spikes", "all"):
        result["spikes"] = spike_examples()
        result["lean_witness_crosscheck"] = lean_witness_crosscheck()
    serialized = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
        print(f"Wrote {args.output}")
        if "verification" in result:
            print(json.dumps(result["verification"]))
        for item in result.get("extrema", []):
            print(f"p={item['p']} m={item['m']} n={item['n']} max={item['maximum_absolute_sum']}")
        for item in result.get("spikes", []):
            print(json.dumps({k: v for k, v in item.items() if k != "B"}))
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
