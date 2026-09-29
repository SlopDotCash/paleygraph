#!/usr/bin/env python3
"""Exact checks for repeated-incidence row capacity and two abstract models."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
from pathlib import Path
import json

from mixed_period_collisions import prime, primitive_root

ROOT = Path(__file__).resolve().parents[1]


def ceil_sqrt(x):
    y = isqrt(x)
    return y + (y * y < x)


def actual(p, n):
    assert prime(p) and n % 2 == 0 and (p - 1) % n == 0
    m, g = (p - 1) // n, primitive_root(p)
    labels = {pow(g, n * i, p): i for i in range(m)}
    C = [Counter() for _ in range(m)]
    for x in range(1, p - 1):
        C[labels[pow(x, n, p)]][labels[pow(x + 1, n, p)]] += 1
    a = [C[0][i] for i in range(m)]
    h = [sum(c * (c - 1) for c in row.values()) for row in C]
    X = sum(c * (c - 1) * (c - 2) for row in C for c in row.values())
    for i in range(m):
        corr = sum(a[j] * a[(j + i) % m] for j in range(m))
        assert h[i] == corr - (n - 1) * (i == 0)
    assert sum(h) == (n - 1) * (n - 2)
    sets = {frozenset(range(m))}
    if m <= 12:
        sets.update(frozenset(i for i in range(m) if mask >> i & 1)
                    for mask in range(1, 1 << m))
    else:
        for size in [1, 2, 3]:
            sets.update(map(frozenset, combinations(range(m), size)))
        sets.update(frozenset(i for i in range(m) if h[i] > r) for r in set(h) | {0})
    checks, nonzero = 0, 0
    for D in sets:
        if not D:
            continue
        d = len(D)
        R = max((h[i] for i in range(m) if i not in D), default=0)
        B = max(0, sum(h[i] for i in D) - n * d * max(1, ceil_sqrt(R)))
        XD = sum(C[i][j] * (C[i][j] - 1) * (C[i][j] - 2) for i in D for j in D)
        assert B ** 3 <= 6 * d * d * XD * XD <= 6 * d * d * X * X
        checks += 1
        nonzero += B > 0
    return dict(p=p, n=n, m=m, in_sparse_source_range=n * n < p,
                X=X, subset_checks=checks, positive_lower_certificates=nonzero, passed=True)


def threshold_check(h, n, X):
    """Sufficient integer verification of every total-X subset inequality.

    At fixed outside threshold R, the strongest subset is {h>R}.
    A floor square root makes the tested positive numerator larger.
    """
    distribution = Counter(h.values())
    count = sum(c for value, c in distribution.items() if value > 0)
    mass = sum(value * c for value, c in distribution.items())
    assert max(h.values()) < n * n
    tests = 0
    largest_ratio = F(0)
    for R in sorted(set(distribution) | {0}):
        if R > 0:
            count -= distribution[R]
            mass -= R * distribution[R]
        if count:
            numerator = max(0, mass - n * count * max(1, isqrt(R)))
            ratio = F(numerator ** 3, 6 * count * count * X * X)
            assert ratio <= 1
            largest_ratio = max(largest_ratio, ratio)
            tests += 1
    return dict(thresholds=tests, largest_conservative_squared_ratio=str(largest_ratio), passed=True)


def abstract(t, kind):
    n, A, d = t ** 5, t ** 3, t
    l = (n - 2 - A * d) // 2
    I = set(range(1, d + 1))
    if kind == 'sidon':
        B = 2 * l * l + 1
        filler = {B * j + j * j + 10 * d + 1 for j in range(1, l + 1)}
    else:
        filler = set(range(10 * n + 1, 10 * n + 1 + l))
    m = n ** 3 + 1
    while not prime(m):
        m += 2
    assert m > 2 * max(filler)
    assert not I & filler and d + 1 not in I | filler
    a = {x: A for x in I} | {x: 2 for x in filler} | {d + 1: 1}
    assert 0 not in a and sum(a.values()) == n - 1
    assert sum(x % 2 for x in a.values()) == 1
    corr = Counter()
    for x, vx in a.items():
        for y, vy in a.items():
            corr[(x - y) % m] += vx * vy
    h = dict(corr)
    h[0] -= n - 1
    assert min(h.values()) >= 0 and sum(h.values()) == (n - 1) * (n - 2)
    K = sum(x * x for x in corr.values())
    A2 = sum(x * x for x in a.values())
    levels = Counter(a.values())
    Q = max(value ** 3 * count for value, count in levels.items())
    assert K <= 180 * n ** 3 and A2 <= 3 * t ** 7 and Q == n * n
    assert max(a.values()) ** 3 <= n * n and (n - 1) ** 3 <= n * n * m
    width = d // 4
    values = [sum(value ** 2 * a.get((x + u) % m, 0) ** 2 * a.get((x + v) % m, 0) ** 2
                  for x, value in a.items())
              for u in range(-width, width + 1) for v in range(-width, width + 1)]
    norm_cube_lower = min(values) ** 3 * len(values) ** 2
    assert 128 * norm_cube_lower >= t ** 61
    row = dict(kind=kind, t=t, n=n, m=m, A2=A2, K=K, Q=Q,
               norm_cube_lower_bound=norm_cube_lower, actual_subgroup_realization=False)
    if kind == 'sidon':
        D = {i % m for i in range(-(d - 1), d)}
        mass = sum(h.get(i, 0) for i in D)
        outside = max(v for i, v in h.items() if i not in D)
        assert mass == A * A * d * d + n - 2 - 2 * A
        assert outside <= 6 * A + 8
        row.update(core_size=len(D), core_mass=mass, max_outside_row_mass=outside)
    else:
        X = 192 * n * n
        B = 2 * (n - 1) ** 2 - (n - 1) + X
        assert X % 6 == 0 and K <= n * B
        assert X >= 3 * sum(v * (v - 1) * (v - 2) for v in a.values())
        row.update(model_X_budget=X, all_subset_total_X_capacity=threshold_check(h, n, X))
    return row


def symbolic_sidon(t):
    n, A, d = t ** 5, t ** 3, t
    count, mass = 2 * d - 1, A * A * d * d + n - 2 - 2 * A
    B = mass - n * count * ceil_sqrt(6 * A + 8)
    assert B * 2 >= t ** 8
    lower_squared = F(B ** 3, 6 * count * count)
    if t >= 4096:
        assert lower_squared > (192 * n * n) ** 2
    return dict(t=t, n=n, core_mass=mass, valid_integer_numerator=B,
                lower_X_squared_over_n4=str(lower_squared / n ** 4),
                excludes_192n2_budget=lower_squared > (192 * n * n) ** 2,
                scope='Closed-form necessary inequality; no huge group enumerated.')


def main():
    actuals = [actual(p, n) for p, n in [(7, 6), (13, 2), (97, 8), (257, 16), (353, 16)]]
    models = [abstract(t, kind) for t in [2, 4] for kind in ['sidon', 'interval']]
    symbolic = [symbolic_sidon(t) for t in [1024, 4096]]
    result = dict(scope='Necessary row-capacity inequality, exclusion of the Sidon-filled family, and a surviving abstract family; no uniform bound improved.',
                  actual_cases=actuals, abstract_cases=models, symbolic_checks=symbolic,
                  all_passed=True)
    path = ROOT / 'results/parallel46_row_capacity_2026_09_06.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_passed=True, actual_subset_checks=sum(r['subset_checks'] for r in actuals),
        positive_actual_certificates=sum(r['positive_lower_certificates'] for r in actuals),
        abstract_models=len(models), symbolic_checks=len(symbolic), output=str(path))))


if __name__ == '__main__':
    main()
