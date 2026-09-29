#!/usr/bin/env python3
"""Exact checks for biased-row moment obstructions; not a Paley proof.

All combinatorics and character moments use Python integers. Decimal values
of the asymptotic variational coefficient are explicitly numerical evidence.
"""
from decimal import Decimal, localcontext
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json

from paley_exact import character, prime

ROOT = Path(__file__).resolve().parents[1]


def coefficient(r):
    with localcontext() as context:
        context.prec = 90
        one, two = Decimal(1), Decimal(2)

        def rate(t):
            return ((one + t) * (one + t).ln()
                    + (one - t) * (one - t).ln()) / two

        def derivative_sign(t):
            atanh = ((one + t).ln() - (one - t).ln()) / two
            return 2 * r * rate(t) - t * atanh

        lo, hi = Decimal('0.5'), one - Decimal('1e-70')
        assert derivative_sign(lo) > 0 > derivative_sign(hi)
        for _ in range(250):
            mid = (lo + hi) / two
            if derivative_sign(mid) > 0:
                lo = mid
            else:
                hi = mid
        t = (lo + hi) / two
        d = t ** (2 * r) / rate(t)
        return {
            'r': r,
            'numerical_optimizer_t': str(t),
            'numerical_D_r': str(d),
            'numerical_necessary_LM_C_r': str((one - one / r) * d),
        }


def main():
    records = []
    subset_checks = moment_checks = 0
    for p, max_k in [(5, 2), (7, 3), (13, 5), (17, 6), (29, 4)]:
        assert prime(p)
        chi = character(p)
        matrix = [[chi[(a - b) % p] for b in range(p)] for a in range(p)]
        d = (p - 1) // 2
        for k in range(1, min(max_k, d) + 1):
            js = range(k // 2 + 1, k + 1)
            totals = {j: 0 for j in js}
            best = {j: (-1, None, None) for j in js}
            for u in combinations(range(p), k):
                counts = []
                for b in range(p):
                    if b not in u:
                        counts.append((b, sum(matrix[a][b] == 1 for a in u)))
                for j in js:
                    good = [b for b, count in counts if count >= j]
                    totals[j] += len(good)
                    if len(good) > best[j][0]:
                        best[j] = (len(good), u, good)
                subset_checks += 1
            for j in js:
                expected_total = p * sum(comb(d, i) * comb(d, k - i)
                                         for i in range(j, k + 1))
                assert totals[j] == expected_total
                denominator = comb(p, k)
                n = (expected_total + denominator - 1) // denominator
                maximum, u, good = best[j]
                assert maximum * denominator >= expected_total
                bset = good[:n]
                rows = [sum(matrix[a][b] for b in bset) for a in range(p)]
                signed_rectangle = sum(rows[a] for a in u)
                assert signed_rectangle >= (2 * j - k) * n
                moments = {}
                for r in (2, 3, 4):
                    moment = sum(x ** (2 * r) for x in rows)
                    local_moment = sum(rows[a] ** (2 * r) for a in u)
                    assert k ** (2 * r - 1) * local_moment >= signed_rectangle ** (2 * r)
                    assert k ** (2 * r - 1) * moment >= (2 * j - k) ** (2 * r) * n ** (2 * r)
                    moments[str(2 * r)] = moment
                    moment_checks += 1
                records.append({
                    'p': p, 'k': k, 'j': j,
                    'total_qualifying_columns_over_all_row_sets': expected_total,
                    'row_set_count': denominator,
                    'guaranteed_n': n,
                    'maximum_qualifying_column_count': maximum,
                    'witness_U': list(u), 'witness_B': bset,
                    'signed_rectangle_sum': signed_rectangle,
                    'exact_moments': moments,
                })

    hypergeometric_checks = 0
    for p in [5, 7, 13, 17, 29, 41, 97, 257, 1009]:
        assert prime(p)
        d = (p - 1) // 2
        for k in range(1, min(d, 25) + 1):
            for j in range(k + 1):
                # P(exactly j plus signs and no zeros)
                # >= (1-k^2/p) * Binomial(k,1/2)[j].
                left = p * 2 ** k * comb(d, j) * comb(d, k - j)
                right = (p - k * k) * comb(k, j) * comb(p, k)
                assert left >= right
                hypergeometric_checks += 1

    report = {
        'claim': 'Biased-row double-counting lower bounds; no upper bound or full conjecture proof',
        'exact_subset_checks': subset_checks,
        'exact_counting_identities': len(records),
        'exact_moment_lower_bound_checks': moment_checks,
        'exact_hypergeometric_lower_bound_checks': hypergeometric_checks,
        'finite_records': records,
        'coefficient_values_are_numerical_not_certified': [coefficient(r) for r in range(2, 9)],
        'source_sha256': {
            path: sha256((ROOT / path).read_bytes()).hexdigest()
            for path in [
                'experiments/parallel_classical_2026_09_04.py',
                'research/parallel-classical-2026-09-04.md',
            ]
        },
    }
    target = ROOT / 'results/parallel_classical_2026_09_04.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ['finite_records', 'coefficient_values_are_numerical_not_certified']}, indent=2))


if __name__ == '__main__':
    main()
