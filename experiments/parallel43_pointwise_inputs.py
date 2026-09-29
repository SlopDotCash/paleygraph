#!/usr/bin/env python3
"""Exact collision lifts, a quartic witness, and the energy/triangle bound.

Finite certificates only; no full Paley proof or uniform energy saving.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import json

from mixed_period_collisions import run, subgroup
from parallel41_independent_verify import triangle

ROOT = Path(__file__).resolve().parents[1]


def product_counts(modulus, n, generator):
    H = {pow(generator, j, modulus) for j in range(n)}
    assert len(H) == n and pow(generator, n, modulus) == 1
    R = sorted((h - 1) % modulus for h in H if h != 1)
    buckets = defaultdict(lambda: [0, 0])
    for a, b in combinations_with_replacement(R, 2):
        buckets[a * b % modulus][int(a != b)] += 1
    # Counts for square-square, square-distinct-pair, and pair-pair.
    parts = [sum(d * (d - 1) for d, o in buckets.values()),
             4 * sum(d * o for d, o in buckets.values()),
             4 * sum(o * (o - 1) for d, o in buckets.values())]
    B = sum((d + 2 * o) ** 2 for d, o in buckets.values())
    X = B - 2 * (n - 1) ** 2 + (n - 1)
    assert sum(parts) == X >= 0
    return dict(X=X, B=B, contributions_by_distinct_entries=parts,
                distinct_product_values=len(buckets))


def lift_sequence(p, n, generator):
    subgroup(p, n, generator)
    g, modulus, sequence, parts = generator, p, [], []
    while True:
        counts = product_counts(modulus, n, g)
        X = counts['X']
        if sequence:
            assert X <= sequence[-1]
        sequence.append(X)
        parts.append(counts['contributions_by_distinct_entries'])
        if X == 0:
            break
        # Explicit bound only guards this finite checker, not the theorem.
        assert len(sequence) < 12
        residual = (pow(g, n, modulus * p) - 1) // modulus
        correction = -residual * pow(n * pow(g, n - 1, p) % p, -1, p) % p
        g += correction * modulus
        modulus *= p
        assert pow(g, n, modulus) == 1
    return dict(p=p, n=n, generator=generator, X_by_precision=sequence,
                types_by_precision=parts, valuation_P=sum(sequence),
                first_zero_precision=len(sequence), all_passed=True)


def witness():
    p, n, g = 278177, 32, 160164
    H = subgroup(p, n, g)
    collision = run(p, n, g)
    tri = triangle(p, n)
    a = Counter(pow((h - 1) % p, n, p) for h in H if h != 1)
    assert Counter(a.values()) == {1: 1, 2: 15}
    assert n ** 4 <= 4 * p <= 4 * n ** 4
    assert collision['additive_energy'] == 3 * n ** 2 - 3 * n == 2976
    assert collision['row_zero_ordered_triples'] == 0
    assert tri['X'] == tri['X_dist'] == 72
    assert product_counts(p, n, g)['contributions_by_distinct_entries'] == [0, 0, 72]
    assert len(collision['all_cells_above_two']) == 12
    for row in collision['all_cells_above_two']:
        u, v = row['coset_power_labels']
        T = sum(z * z * a[u * c % p] ** 2 * a[v * c % p] ** 2 for c, z in a.items())
        assert len({1, u, v}) == 3 and row['size'] == 3 and T == 144
        row['T_b'] = T
    assert tri['W'] == 1849344 and tri['high'] == 165888
    # Independently verify the single-coset cardinality hypothesis and
    # normalization used with the published incidence lemma.
    x = collision['all_cells_above_two'][0]['members'][0]
    Q1, Q2 = {x * h % p for h in H}, {(x + 1) * h % p for h in H}
    ratio_image = {(a0 * pow(c, -1, p) % p, b0 * pow(c, -1, p) % p)
                   for a0 in Q1 for b0 in Q2 for c in H}
    assert len(ratio_image) == n ** 2
    incidence = sum((y - x0) % p in H for x0 in Q1 for y in Q2)
    assert incidence == n * 3
    return dict(**collision, positive_a_histogram=dict(Counter(a.values())),
                triangle=tri, source_hypothesis_check=dict(ratio_image=n ** 2,
                required_cardinality=n ** 2, incidence=incidence), all_passed=True)


def main():
    prior = json.loads((ROOT / 'results/parallel42_norm_budget_2026_09_06.json').read_text())
    lifts = []
    for case in prior['cases']:
        for exceptional in case['all_split_exceptional_primes']:
            row = lift_sequence(exceptional['p'], case['n'], exceptional['generator'])
            assert row['valuation_P'] == exceptional['valuation_P']
            assert row['X_by_precision'][0] == exceptional['X']
            lifts.append(row)
    quartic = witness()
    larger_lifts = [lift_sequence(p, n, g) for p, n, g in [
        (278177, 32, 160164), (6700417, 64, 2), (67403009, 128, 64701253)]]
    tests = []
    for p, n in [(97, 8), (353, 16), (1153, 8)]:
        g = next(h for h in range(2, p) if pow(h, n, p) == 1 and pow(h, n // 2, p) != 1)
        coll = run(p, n, g)
        tri = triangle(p, n)
        rho_max = max([2] + [row['size'] for row in coll['all_cells_above_two']])
        E_star = coll['additive_energy'] - n ** 2
        assert tri['W'] * n ** 2 <= rho_max * E_star ** 3
        tests.append(dict(p=p, n=n, W=tri['W'], E_star=E_star,
                          rho_max=rho_max, exact_bound_passed=True))
    E_star = quartic['additive_energy'] - quartic['n'] ** 2
    assert quartic['triangle']['W'] * quartic['n'] ** 2 <= 3 * E_star ** 3
    tests.append(dict(p=278177, n=32, W=quartic['triangle']['W'],
                      E_star=E_star, rho_max=3, exact_bound_passed=True))
    assert Q(2, 3) - 2 + 3 * Q(7, 3) == Q(17, 3)
    assert 4 - Q(7, 3) == Q(5, 3)
    assert Q(63, 31) - Q(5, 3) == Q(34, 93)
    result = dict(scope='Exact finite and prime-power identities; ordinary uniform proof remains separate.',
                  quartic_zero_diagonal_witness=quartic,
                  independent_norm_valuation_comparisons=lifts,
                  larger_lift_certificates=larger_lifts, energy_triangle_checks=tests,
                  exceptional_prime_exponent='5/3', previous_exponent='63/31',
                  exponent_improvement='34/93', all_passed=True)
    dest = ROOT / 'results/parallel43_pointwise_inputs_2026_09_06.json'
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_passed=True, norm_valuation_comparisons=len(lifts),
                          quartic_witness=dict(p=278177, n=32, X=72),
                          larger_lifts=larger_lifts, output=str(dest))))


if __name__ == '__main__':
    main()
