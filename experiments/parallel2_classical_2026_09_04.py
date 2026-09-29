#!/usr/bin/env python3
"""Exact elliptic fourth-moment reduction and a failed sufficient L2 route.

This verifies identities and finite examples, not a Paley proof. The
asymptotic refutation is proved in the accompanying mathematical note.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
from random import Random
import json

from paley_exact import character, prime, quartic_trace, row_sums

ROOT = Path(__file__).resolve().parents[1]


def falling4(n):
    return n * (n - 1) * (n - 2) * (n - 3)


def crossratio(p, chi, a, b, c, d):
    denominator = (c - a) * (d - b) % p
    lam = (c - b) * (d - a) * pow(denominator, -1, p) % p
    assert lam not in (0, 1)
    return lam, chi[denominator]


def elliptic(p, chi):
    return {lam: sum(chi[x * (x - 1) * (x - lam) % p] for x in range(p))
            for lam in range(2, p)}


def vector_from_quads(p, chi, subset):
    vector = Counter()
    for quad in permutations(subset, 4):
        lam, sign = crossratio(p, chi, *quad)
        vector[lam] += sign
    return vector


def rational(value):
    value = Fraction(value)
    return {'numerator': value.numerator, 'denominator': value.denominator}


def main():
    pointwise_checks = subset_checks = 0
    field_records = []
    rng = Random(9042026)
    for p in (5, 7, 13, 17, 29):
        chi = character(p)
        e = elliptic(p, chi)
        h = p * p - 2 * p - 3
        assert sum(value * value for value in e.values()) == h
        quad_vectors = {}
        full_vector = Counter()
        for quad in combinations(range(p), 4):
            literal = quartic_trace(p, chi, quad)
            vector = Counter()
            for ordered in permutations(quad):
                lam, sign = crossratio(p, chi, *ordered)
                assert literal == sign * e[lam] - 1
                vector[lam] += sign
                pointwise_checks += 1
            quad_vectors[quad] = vector
            full_vector.update(vector)
        assert all(full_vector[lam] == p * (p - 1) * e[lam] for lam in e)

        if p <= 13:
            subsets = [s for n in range(p + 1) for s in combinations(range(p), n)]
        else:
            subsets = [tuple(range(p)), (), (0,), (0, 1)]
            subsets += [tuple(sorted(rng.sample(range(p), rng.randrange(2, 10))))
                        for _ in range(40)]
        for subset in subsets:
            n = len(subset)
            vector = Counter()
            for quad in combinations(subset, 4):
                vector.update(quad_vectors[quad])
            rows = row_sums(p, chi, subset)
            m4 = sum(x ** 4 for x in rows)
            local = sum(rows[a] ** 2 for a in subset)
            pairing = sum(vector[lam] * value for lam, value in e.items())
            base = p * (3 * n * n - 2 * n) - n ** 4 + 3 * n * n - 3 * n - 6 * local
            assert m4 == base + pairing
            centering = Fraction(falling4(n), (p - 2) * (p - 3))
            centered = {lam: vector[lam] - centering * value for lam, value in e.items()}
            rterm = (p * (3 * n * n - 2 * n) - 6 * n ** 3 + 14 * n * n
                     - 9 * n - 6 * local + Fraction(3 * falling4(n), p - 2))
            assert m4 == rterm + sum(centered[lam] * value for lam, value in e.items())
            assert rterm <= 3 * p * n * n
            wenergy = sum(value * value for value in centered.values())
            centered_pairing = m4 - rterm
            assert centered_pairing ** 2 <= h * wenergy
            assert max(0, m4 - 3 * p * n * n) ** 2 <= p * p * wenergy
            subset_checks += 1
        field_records.append({'p': p, 'subsets_checked': len(subsets),
                              'elliptic_norm_squared': h})

    interval_records = []
    for n, p in [(8, 3361), (12, 18481), (16, 960961)]:
        assert prime(p)
        chi = character(p)
        assert p % 4 == 1
        assert all(chi[t % p] == 1 for t in range(-(n - 1), n) if t)
        subset = tuple(range(1, n + 1))
        vector = vector_from_quads(p, chi, subset)
        assert sum(vector.values()) == falling4(n)
        assert all(value > 0 for value in vector.values())
        energy = sum(value * value for value in vector.values())
        m = (n + 1) // 2
        collision_floor = m * m * (m - 1) * (m - 2) * (m - 3)
        assert energy >= collision_floor
        rows = row_sums(p, chi, subset)
        m4 = sum(x ** 4 for x in rows)
        local = sum(rows[a] ** 2 for a in subset)
        base = p * (3 * n * n - 2 * n) - n ** 4 + 3 * n * n - 3 * n - 6 * local
        pairing_from_identity = m4 - base
        c = Fraction(falling4(n), (p - 2) * (p - 3))
        h = p * p - 2 * p - 3
        centered_energy = energy - 2 * c * pairing_from_identity + c * c * h
        assert centered_energy >= 0
        interval_records.append({
            'n': n, 'p': p,
            'all_nonzero_interval_differences_are_squares': True,
            'uncentered_vector_norm_squared': energy,
            'translation_collision_lower_bound': collision_floor,
            'centered_vector_norm_squared_via_exact_moment_identity': rational(centered_energy),
            'centered_norm_squared_divided_by_n4': rational(centered_energy / n ** 4),
            'exact_fourth_moment': m4,
            'meets_asymptotic_p_ge_n8_condition': p >= n ** 8,
        })

    report = {
        'claim': 'Exact elliptic reduction and finite collision checks; CS-input is disproved by the accompanying asymptotic proof, not by these finite tests',
        'pointwise_fractional_linear_checks': pointwise_checks,
        'subset_fourth_moment_and_centering_checks': subset_checks,
        'fields': field_records,
        'interval_examples': interval_records,
        'source_sha256': {
            path: sha256((ROOT / path).read_bytes()).hexdigest()
            for path in [
                'research/parallel2-classical-2026-09-04.md',
                'experiments/parallel2_classical_2026_09_04.py',
            ]
        },
    }
    (ROOT / 'results/parallel2_classical_2026_09_04.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
