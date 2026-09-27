#!/usr/bin/env python3
"""Finite checks for supergroup cosets and an abstract interval obstruction."""
from collections import Counter
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

from mixed_period_collisions import prime, primitive_root, subgroup

ROOT = Path(__file__).resolve().parents[1]


def supergroup_case():
    p, n, g = 215535361, 128, 25525303
    H = set(subgroup(p, n, g))
    ambient = primitive_root(p)
    a = Counter(pow((h - 1) % p, n, p) for h in H if h != 1)
    groups = {}
    cap_checks = []
    for d in [1, 2, 3, 5, 6, 15]:
        assert (p - 1) % (n * d) == 0
        gen = pow(ambient, (p - 1) // (n * d), p)
        S = {pow(gen, j, p) for j in range(n * d)}
        assert len(S) == n * d and H <= S
        groups[d] = S
        for scale in [1, 2, 3, 5]:
            field_coset = {scale * x % p for x in S}
            labels = {pow(x, n, p) for x in field_coset}
            assert len(labels) == d
            mass = sum(a[c] for c in labels)
            intersection = sum((h - 1) % p in field_coset for h in H if h != 1)
            assert mass == intersection
            cap_checks.append(dict(quotient_size=d, scale=scale, mass=mass,
                                   mass_cubed_over_n2d=str(F(mass ** 3, n * n * d)),
                                   mitkin_size_condition=n * n * d >= 33 ** 3,
                                   mitkin_field_condition=(n * n * d) ** 2 < p ** 3))
    incidence_checks = []
    # The source hypotheses hold literally for these larger groups.
    for scales in [(1, 1, 2), (1, 2, 3), (2, 3, 5)]:
        ds = (3, 3, 15)
        sets = [{scale * x % p for x in groups[d]} for scale, d in zip(scales, ds)]
        Z, X, Y = sets
        assert groups[3] <= groups[15]
        z0 = scales[0]
        local_count = sum((x + z0) % p in Y for x in X)
        direct = sum((y - x) % p in Z for x in X for y in Y)
        assert direct == len(Z) * local_count and direct > 0
        size_product = len(X) * len(Y)
        assert size_product >= 33 ** 3 and size_product ** 2 < p ** 3
        incidence_checks.append(dict(quotient_sizes=list(ds), scales=list(scales),
            direct_incidence=direct, normalized_intersection=local_count,
            scaling_multiplicity=len(Z), source_hypotheses_passed=True))
    return dict(p=p, n=n, ambient_generator=ambient, mass_checks=cap_checks,
                nested_incidence_checks=incidence_checks, all_passed=True)


def interval_model(t):
    assert t >= 2 and t & (t - 1) == 0
    n, A, d = t ** 5, t ** 3, t
    fillers = (n - 2 - A * d) // 2
    B = 2 * fillers ** 2 + 1
    S = {B * j + j * j + 10 * d + 1 for j in range(1, fillers + 1)}
    start = max(n ** 3, 4 * max(S) + 1)
    modulus = start | 1
    while not prime(modulus):
        modulus += 2
    assert modulus > 2 * max(S)
    I = set(range(1, d + 1))
    assert not I & S and d + 1 not in I | S
    a = {x: A for x in I} | {x: 2 for x in S} | {d + 1: 1}
    assert 0 not in a
    assert sum(a.values()) == n - 1
    assert sum(v % 2 for v in a.values()) == 1
    A2 = sum(v * v for v in a.values())
    assert A2 ** 20 <= 3 ** 20 * n ** 29
    # Exact Sidon check for the filler, with no modular wrap in pair sums.
    sums = Counter(x + y for x in S for y in S)
    assert sum(v * v for v in sums.values()) == 2 * fillers ** 2 - fillers
    corr = Counter()
    for x, vx in a.items():
        for y, vy in a.items():
            corr[(x - y) % modulus] += vx * vy
    K = sum(v * v for v in corr.values())
    assert 3 * K >= 2 * n ** 3 and K <= 180 * n ** 3
    # All values are powers of two; these are their exact dyadic levels.
    levels = Counter(a.values())
    Q = max(v ** 3 * count for v, count in levels.items())
    assert Q == n ** 2
    # A prime-order group has only singleton and whole-group cosets.
    assert max(a.values()) ** 3 <= n ** 2
    assert sum(a.values()) ** 3 <= n ** 2 * modulus
    width = d // 4
    targets = [(u, v) for u in range(-width, width + 1)
               for v in range(-width, width + 1)]
    assert len(targets) * 4 >= d * d
    minimum = None
    for u, v in targets:
        T = sum(z ** 2 * a.get((c + u) % modulus, 0) ** 2
                * a.get((c + v) % modulus, 0) ** 2 for c, z in a.items())
        assert 2 * T >= A ** 6 * d
        minimum = T if minimum is None else min(minimum, T)
    # ||T||_(3/2)^3 >= min(T)^3 * number_of_targets^2.
    lower_cube = minimum ** 3 * len(targets) ** 2
    assert 128 * lower_cube >= t ** 61
    return dict(t=t, n=n, amplitude=A, interval_size=d, prime_group_order=modulus,
                filler_size=fillers, mass=sum(a.values()), A2=A2,
                K=K, Q=Q, tested_target_pairs=len(targets),
                norm_cube_lower_bound=lower_cube,
                subgroup_coset_caps_passed=True, is_actual_subgroup_difference=False,
                all_passed=True)


def matrix_identities(p, n):
    assert prime(p) and (p - 1) % n == 0
    m, g = (p - 1) // n, primitive_root(p)
    q = pow(g, n, p)
    labels = {pow(q, j, p): j for j in range(m)}
    C = [Counter() for _ in range(m)]
    for x in range(1, p - 1):
        C[labels[pow(x, n, p)]][labels[pow(x + 1, n, p)]] += 1
    a = C[0]
    square = [Counter() for _ in range(m)]
    for i, row in enumerate(C):
        for j, v in row.items():
            for k, w in C[j].items():
                square[i][k] += v * w
    for i in range(m):
        for j in range(m):
            lhs = n * (i == j) - n * (i == j == 0)
            lhs += sum(v * C[(i - k) % m][(j - k) % m] for k, v in a.items())
            assert square[i][j] == lhs
    for k in range(m):
        inner = sum(v * C[(i - k) % m][(j - k) % m]
                    for i, row in enumerate(C) for j, v in row.items())
        assert inner == n * (n - 1) + (p - 2 * n) * (k == 0)
    for weights in [dict(a), {i: v * v for i, v in a.items()}, {0: 1, 1: -1}]:
        matrix_norm = sum(sum(v * C[(i - k) % m][(j - k) % m]
                         for k, v in weights.items()) ** 2
                         for i in range(m) for j in range(m))
        assert matrix_norm == n * (n - 1) * sum(weights.values()) ** 2 + (p - 2 * n) * sum(v * v for v in weights.values())
    return dict(p=p, n=n, quotient_order=m, squared_matrix_entries=m * m,
                shifted_inner_products=m, weighted_frobenius_checks=3, all_passed=True)


def main():
    actual = supergroup_case()
    models = [interval_model(t) for t in [2, 4]]
    matrices = [matrix_identities(p, n) for p, n in [(97, 8), (353, 16), (1153, 8)]]
    assert F(5, 3) + 4 == F(17, 3)
    assert 1 - F(4, 3) == -F(1, 3)
    assert F(1, 3) - F(4, 3) == -1
    assert 2 * F(23, 60) + 2 * F(49, 20) == F(17, 3)
    assert F(3, 5) + F(52, 15) == F(61, 15)
    result = dict(scope='Exact structural identities and abstract functional obstruction; no full Paley proof.',
                  actual_supergroup_checks=actual, abstract_interval_models=models,
                  matrix_identity_checks=matrices,
                  all_passed=True)
    dest = ROOT / 'results/parallel44_structured_levels_2026_09_06.json'
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_passed=True, actual_mass_checks=len(actual['mass_checks']),
        nested_incidence_checks=len(actual['nested_incidence_checks']),
        abstract_models=[{k: m[k] for k in ['n', 'prime_group_order', 'A2', 'K', 'Q']} for m in models],
        output=str(dest))))


if __name__ == '__main__':
    main()
