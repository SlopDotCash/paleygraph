#!/usr/bin/env python3
"""Exact checks of triangle edge-coset decomposition, with explicit scope."""
from collections import Counter
from math import isqrt
from pathlib import Path
import json


def check(p, generator, n):
    assert p > 2 and all(p % d for d in range(2, isqrt(p) + 1))
    h = {pow(generator, j, p) for j in range(n)}
    assert len(h) == n and pow(generator, n, p) == 1 and p - 1 in h
    f = Counter((a - b) % p for a in h for b in h)
    f.pop(0)
    labels, representatives = {}, []
    for x in sorted(f):
        if x in labels:
            continue
        representatives.append(x)
        for z in h:
            y = x * z % p
            assert y in f and f[y] == f[x]
            labels[y] = x
    assert len(labels) == len(f) == n * len(representatives)
    moments = {k: sum(v ** k for v in f.values()) for k in range(2, 7)}
    kappa = f.get(1, 0)
    weights = {x: v * v for x, v in f.items()}
    by_distinct = Counter({1: 0, 2: 0, 3: 0})
    p_ab = p_ac = p_bc = 0
    for a, w_a in weights.items():
        for b, w_b in weights.items():
            c = (a - b) % p
            w_c = weights.get(c, 0)
            if not w_c:
                continue
            weight = w_a * w_b * w_c
            la, lb, lc = labels[a], labels[b], labels[c]
            by_distinct[len({la, lb, lc})] += weight
            p_ab += weight * (la == lb)
            p_ac += weight * (la == lc)
            p_bc += weight * (lb == lc)
    p_h = sum(f[a] ** 4 * f.get((1 - z) * a % p, 0) ** 2
              for z in h if z != 1 for a in f)
    p_q = n * sum(f[a] ** 4 * f[d] * f.get(a * d % p, 0) ** 2
                  for a in representatives for d in representatives)
    assert p_ab == p_ac == p_bc == p_h == p_q
    assert n * p_ab <= moments[3] * moments[4]
    assert by_distinct[1] == kappa * moments[6]
    repeated = by_distinct[1] + by_distinct[2]
    assert repeated == 3 * p_ab - 2 * kappa * moments[6]
    assert n * repeated <= 3 * moments[3] * moments[4]
    full = sum(by_distinct.values())
    assert full <= moments[3] ** 2
    return {
        'p': p, 'generator': generator, 'n': n,
        'difference_support_size': len(f),
        'nonzero_difference_moments': moments,
        'kappa': kappa,
        'triangle_weight_by_number_of_distinct_edge_cosets': dict(by_distinct),
        'one_equal_edge_pair_weight': p_ab,
        'repeated_edge_coset_weight': repeated,
        'full_triangle_weight': full,
        'full_triangle_weight_Young_upper': moments[3] ** 2,
        'pair_bound_numerator': moments[3] * moments[4],
        'pair_bound_denominator': n,
    }


def main():
    rows = [check(*case) for case in [
        (5, 4, 2), (13, 8, 4), (17, 2, 8), (97, pow(5, 12, 97), 8),
        (97, pow(5, 4, 97), 24), (1153, 75, 8), (6700417, 2, 64),
    ]]
    out = Path(__file__).resolve().parents[1] / 'results/parallel39_edge_cosets_2026_09_06.json'
    out.write_text(json.dumps({
        'scope': 'exact finite verification of identities and inequalities; no uniform off-region bound',
        'full_goal_proved': False, 'rows': rows,
    }, indent=2) + '\n')
    for row in rows:
        print('p', row['p'], 'n', row['n'],
              'weights', row['triangle_weight_by_number_of_distinct_edge_cosets'],
              'P', row['one_equal_edge_pair_weight'])
    print(out)


if __name__ == '__main__':
    main()
