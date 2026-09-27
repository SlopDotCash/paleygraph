#!/usr/bin/env python3
"""Exact fixed-conductor norm-budget certificates for n=4,8,16.

No field-sized scan or general uniform estimate. The product is stored by
its exact prime factorization, using conjugate-pair product differences.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import isqrt
from pathlib import Path
import json

from parallel41_shell_real_recovery import multiply, norm_adjugate
from parallel40_independent_verify import determinant

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def factor(value):
    assert value >= 1
    remaining = value
    output = []
    d = 2
    while d * d <= remaining:
        exponent = 0
        while remaining % d == 0:
            remaining //= d
            exponent += 1
        if exponent:
            output.append((d, exponent))
        d = 3 if d == 2 else d + 2
    if remaining > 1:
        output.append((remaining, 1))
    restored = 1
    for prime, exponent in output:
        restored *= prime ** exponent
    assert restored == value
    return tuple(output)


def excess(p, n):
    assert (p - 1) % n == 0
    assert all(p % d for d in range(2, isqrt(p) + 1))
    g = next(g for x in range(2, p)
             if (g := pow(x, (p - 1) // n, p)) != 1
             and pow(g, n // 2, p) != 1)
    H = {pow(g, j, p) for j in range(n)}
    assert len(H) == n and p - 1 in H
    R = {(h - 1) % p for h in H if h != 1}
    ratios = Counter(x * pow(y, -1, p) % p for x in R for y in R)
    B = sum(v * v for v in ratios.values())
    X = B - 2 * (n - 1) ** 2 + (n - 1)
    assert X >= 0 and X % 6 == 0
    return X, g


def conductor(n):
    N, m = n // 2, n - 1
    shifted = {}
    for a in range(1, n):
        v = [0] * N
        v[0] = -1
        v[a % N] += -1 if a >= N else 1
        shifted[a] = v
    products = []
    for a, d in combinations_with_replacement(range(1, n), 2):
        products.append((multiply(shifted[a], shifted[d]), 1 if a == d else 2))
    norm_valuations = Counter()
    pair_count = tuple_count = determinant_checks = coefficient_square_sum = 0
    max_norm = 1
    distinct_norms = set()
    for (left, wl), (right, wr) in combinations(products, 2):
        v = [a - b for a, b in zip(left, right)]
        value = norm_adjugate(v)[0]
        assert 0 < value <= 6 ** N
        if pair_count % 137 == 0 or n <= 8:
            matrix = [[v[(i - j) % N] * (-1 if i < j else 1)
                       for j in range(N)] for i in range(N)]
            assert determinant(matrix) == value
            determinant_checks += 1
        weight = 2 * wl * wr
        coefficient_square_sum += weight * sum(x * x for x in v)
        for prime, exponent in factor(value):
            norm_valuations[prime] += weight * exponent
        pair_count += 1
        tuple_count += weight
        max_norm = max(max_norm, value)
        distinct_norms.add(value)
    M = m ** 4 - 2 * m ** 2 + m
    assert tuple_count == M
    second_moment = 8 * n * n * m * m - 2 * n ** 4
    assert coefficient_square_sum == second_moment
    # Product of all norms = |P_n|^N, since P_n is a rational integer.
    assert all(e % N == 0 for e in norm_valuations.values())
    valuations = {p: e // N for p, e in sorted(norm_valuations.items())}
    split = []
    for p, exponent in valuations.items():
        if (p - 1) % n:
            continue
        X, g = excess(p, n)
        assert 0 < X <= exponent
        split.append(dict(p=p, generator=g, X=X, valuation_P=exponent,
                          strict_valuation_loss=exponent-X))
    # Check additional split primes independently, including primes outside
    # the exceptional set and beyond every possible exceptional prime.
    tested = 0
    upper = min(max_norm + 2 * n, 5000)
    for p in range(n + 1, upper + 1, n):
        if not all(p % d for d in range(2, isqrt(p) + 1)):
            continue
        X, _ = excess(p, n)
        assert (X > 0) == (p in valuations)
        tested += 1
    # Exact integer form of the global archimedean bound.
    P_abs = 1
    for p, e in valuations.items():
        P_abs *= p ** e
    assert P_abs <= 6 ** M
    assert P_abs ** 2 * M ** M <= second_moment ** M
    return dict(n=n, degree=N, tuple_count=M, unordered_product_count=len(products),
                product_pair_count=pair_count, independent_determinants=determinant_checks,
                distinct_norms=len(distinct_norms), max_individual_norm=max_norm,
                tuple_second_moment=second_moment, sharp_budget_log_argument=[second_moment, M],
                P_absolute_prime_factorization={str(p):e for p,e in valuations.items()},
                all_split_exceptional_primes=split, extra_split_primes_checked=tested,
                all_checks_passed=True)


def main():
    cases = [conductor(n) for n in (4, 8, 16)]
    # A separate actual-subgroup counterexample to every constant-factor
    # domination of off-diagonal rich mass by diagonal rich mass.
    from parallel41_independent_verify import triangle
    witness = triangle(353, 16)
    assert witness['X'] == witness['X_dist'] == 72 and witness['high'] == 82944
    p, n, g = 353, 16, 304
    H = {pow(g, j, p) for j in range(n)}
    assert len(H) == n
    a = Counter(pow((h - 1) % p, n, p) for h in H if h != 1)
    histogram = Counter(a.values())
    assert histogram == {1: 1, 2: 7}
    rho = Counter((pow(x, n, p), pow(x + 1, n, p)) for x in range(1, p - 1))
    rich = []
    for (u, v), count in sorted(rho.items()):
        if count < 3:
            continue
        T = sum(value ** 2 * a[u * c % p] ** 2 * a[v * c % p] ** 2
                for c, value in a.items())
        assert count == 3 and T == 144 and len({1, u, v}) == 3
        rich.append(dict(u=u, v=v, rho=count, T_b=T))
    assert len(rich) == 12 and n * sum(cell['rho'] * cell['T_b'] for cell in rich) == witness['high']
    X, _ = excess(p, n)
    assert X == witness['X']
    witness.update(generator=g, positive_a_histogram=dict(histogram), rich_cells=rich,
                   quartic_window_n4_over4_to_n4=(n ** 4 <= 4 * p <= 4 * n ** 4))
    result = dict(scope='Complete norm-product factorizations for three fixed conductors; no all-conductor or all-prime proof.',
                  formula='sum_(p=1 mod n) X(H_p) log p <= [(n-1)^4-2(n-1)^2+(n-1)] log 6',
                  sharper_formula='sum X log p <= (M/2) log([8n²(n−1)²−2n⁴]/M)',
                  cases=cases, zero_diagonal_witness=witness, all_passed=True)
    dest = ROOT / 'results/parallel42_norm_budget_2026_09_06.json'
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_passed=True, cases=[dict(n=c['n'], exceptions=len(c['all_split_exceptional_primes']),
                         tuple_count=c['tuple_count'], determinants=c['independent_determinants']) for c in cases],
                         result=str(dest)), indent=2))


if __name__ == '__main__':
    main()
