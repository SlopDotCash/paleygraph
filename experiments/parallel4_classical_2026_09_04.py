#!/usr/bin/env python3
"""Exact transfers to structured polynomial-size sets; no Paley proof.

The Burgess and rank-two asymptotic bounds are imported theorems. Finite
checks verify the algebraic transfer and exponent arithmetic, not those
theorems' unspecified constants.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import isqrt
from pathlib import Path
import json

from paley_exact import character, row_sums

ROOT = Path(__file__).resolve().parents[1]


def check_subset(p, chi, subset, elliptic):
    n = len(subset)
    rows = row_sums(p, chi, subset)
    m2 = sum(x * x for x in rows)
    m4 = sum(x ** 4 for x in rows)
    maximum = max(map(abs, rows), default=0)
    assert m2 == p * n - n * n
    local = sum(rows[x] ** 2 for x in subset)
    falling = n * (n - 1) * (n - 2) * (n - 3)
    remainder = (p * (3 * n * n - 2 * n) - 6 * n ** 3 + 14 * n * n
                 - 9 * n - 6 * local + Fraction(3 * falling, p - 2))
    assert 0 <= remainder <= 3 * p * n * n
    assert m4 <= maximum * maximum * m2
    direction = m4 - remainder
    assert abs(direction) <= max(3 * p * n * n, maximum * maximum * m2)
    v = [0] * p
    for a, b, c, d in permutations(subset, 4):
        denominator = (c - a) * (d - b) % p
        lam = (c - b) * (d - a) * pow(denominator, -1, p) % p
        v[lam] += chi[denominator]
    center = Fraction(falling, (p - 2) * (p - 3))
    literal_pairing = sum((v[t] - center * elliptic[t]) * elliptic[t]
                          for t in range(2, p))
    assert literal_pairing == direction
    return {'n': n, 'M2': m2, 'M4': m4, 'maximum_row': maximum,
            'directional_pairing': str(direction), 'remainder': str(remainder)}


def main():
    exhaustive = 0
    exhaustive_fields = []
    for p in (5, 7, 11, 13, 17):
        chi = character(p)
        elliptic = [sum(chi[x * (x - 1) * (x - t) % p] for x in range(p))
                    for t in range(p)]
        count = 0
        for n in range(isqrt(p) + 1):
            for subset in combinations(range(p), n):
                check_subset(p, chi, subset, elliptic)
                count += 1
        exhaustive += count
        exhaustive_fields.append({'p': p, 'subsets': count})

    progression_records = []
    perturbation_checks = 0
    for p, h1, h2 in [(17, 2, 2), (41, 2, 3), (97, 3, 3), (257, 4, 4)]:
        chi = character(p)
        e = [sum(chi[x * (x - 1) * (x - t) % p] for x in range(p))
             for t in range(p)]
        n = h1 * h2
        assert n * n <= p
        for scale, shift in [(1, 0), (2, p - 1), (3, 7)]:
            progression = tuple(sorted({(shift + scale * (i + (h1 + 1) * j)) % p
                                        for i in range(h1) for j in range(h2)}))
            assert len(progression) == n  # Proper rank-two parametrization.
            for x in range(p):
                translated = {(x - b) % p for b in progression}
                assert len(translated) == n
            record = check_subset(p, chi, progression, e)
            record.update({'p': p, 'kind': 'proper_rank_two_GAP',
                           'shape': [h1, h2], 'scale': scale, 'shift': shift})
            progression_records.append(record)

        for scale, shift in [(1, 0), (2, p - 1)]:
            progression = {(shift + scale * i) % p for i in range(n)}
            original_rows = row_sums(p, chi, tuple(progression))
            original_max = max(map(abs, original_rows))
            edited = set(progression)
            edited.remove(min(edited))
            edited.add(next(x for x in range(p) if x not in progression))
            k = len(edited ^ progression)
            record = check_subset(p, chi, tuple(sorted(edited)), e)
            new_rows = row_sums(p, chi, tuple(edited))
            assert all(abs(a - b) <= k for a, b in zip(original_rows, new_rows))
            assert record['maximum_row'] <= original_max + k
            assert record['M4'] <= (original_max + k) ** 2 * (p * n - n * n)
            record.update({'p': p, 'kind': 'perturbed_arithmetic_progression',
                           'symmetric_difference': k,
                           'original_maximum_row': original_max,
                           'meets_finite_edit_budget': p * k ** 30 <= n ** 30})
            progression_records.append(record)
            perturbation_checks += 1

    thresholds = []
    for r in range(2, 101):
        threshold = Fraction(r * r + r + 1, 2 * r * (r + 2))
        gap = Fraction((2 * r - 5) * (r - 3), 30 * r * (r + 2))
        assert threshold - Fraction(13, 30) == gap
        assert gap >= 0 and (gap == 0) == (r == 3)
        if r <= 8:
            thresholds.append({'r': r, 'threshold': str(threshold)})
    exponent_records = []
    for beta in map(Fraction, ['1/120', '1/60', '1/30', '1/20']):
        alpha = Fraction(13, 30) + beta
        eta = beta / 3
        assert alpha < Fraction(1, 2)
        row_saving = alpha / 3 - Fraction(1, 9) - eta
        assert row_saving == Fraction(1, 30)
        directional_saving = alpha + Fraction(1, 15) - Fraction(1, 2)
        assert directional_saving == beta
        baseline_saving = 2 * alpha - Fraction(1, 2)
        assert baseline_saving >= beta
        exponent_records.append({'beta': str(beta), 'minimum_size_exponent': str(alpha),
                                 'row_saving': str(row_saving),
                                 'directional_saving': str(directional_saving)})

    report = {
        'claim': 'Exact transfer checks for polynomial-size structured-class corollaries; imported asymptotic constants are not inferred from finite tests',
        'exhaustive_remainder_and_directional_cases': exhaustive,
        'exhaustive_fields': exhaustive_fields,
        'progression_cases': len(progression_records),
        'perturbation_cases': perturbation_checks,
        'progression_records': progression_records,
        'integer_Burgess_parameters_checked': 99,
        'thresholds': thresholds,
        'exponent_records': exponent_records,
        'primary_source': {
            'url': 'https://arxiv.org/html/2509.07765v1',
            'inputs': ['Introduction: classical Burgess estimate', 'Theorem 1.1: proper rank-two GAPs'],
            'scope': 'Theorems imported; no numerical value assigned to delta(epsilon)',
        },
        'source_sha256': {path: sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in ['research/parallel4-classical-2026-09-04.md',
                                       'experiments/parallel4_classical_2026_09_04.py']},
    }
    (ROOT / 'results/parallel4_classical_2026_09_04.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'progression_records'}, indent=2))


if __name__ == '__main__':
    main()
