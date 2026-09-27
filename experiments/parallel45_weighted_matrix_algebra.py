#!/usr/bin/env python3
"""Exact finite checks for the full weighted coset multiplication algebra."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

from mixed_period_collisions import prime, primitive_root

ROOT = Path(__file__).resolve().parents[1]


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def matvec(A, v):
    return [sum(value * v[j] for j, value in row.items()) for row in A]


def matmul(A, B):
    rows = []
    for row in A:
        out = Counter()
        for j, x in row.items():
            for k, y in B[j].items():
                out[k] += x * y
        rows.append(Counter({j: x for j, x in out.items() if x}))
    return rows


def trace(A):
    return sum(row[i] for i, row in enumerate(A))


def example(p, n):
    assert prime(p) and n % 2 == 0 and n * n < p and (p - 1) % n == 0
    m, g = (p - 1) // n, primitive_root(p)
    labels = {pow(g, n * i, p): i for i in range(m)}
    C = [Counter() for _ in range(m)]
    for x in range(1, p - 1):
        C[labels[pow(x, n, p)]][labels[pow(x + 1, n, p)]] += 1
    a = [C[0][j] for j in range(m)]

    def L(u):
        result = [Counter() for _ in range(m)]
        for k, weight in enumerate(u):
            if weight:
                for i, row in enumerate(C):
                    for j, value in row.items():
                        result[(i + k) % m][(j + k) % m] += weight * value
        return [Counter({j: x for j, x in row.items() if x}) for row in result]

    weights = [[int(i == 0) for i in range(m)],
               [int(i == 1) for i in range(m)], a,
               [x * x for x in a],
               [int(i == 0) - int(i == 1) for i in range(m)]]
    matrices = [L(u) for u in weights]
    ones = L([1] * m)
    for i, j in product(range(m), repeat=2):
        assert ones[i][j] == n - int(i == j)
    product_checks = 0
    for (u, U), (v, V) in product(zip(weights, matrices), repeat=2):
        uv, vu = matvec(U, v), matvec(V, u)
        assert uv == vu
        UV, VU, Luv = matmul(U, V), matmul(V, U), L(uv)
        for i, j in product(range(m), repeat=2):
            assert UV[i][j] == n * dot(u, v) * int(i == j) - n * u[i] * v[j] + Luv[i][j]
            assert UV[i][j] - VU[i][j] == n * (v[i] * u[j] - u[i] * v[j])
        product_checks += 1
    traces, frame_checks = [], []
    alpha, beta = n * (n - 1), p - 2 * n
    D = alpha * m + beta
    for u, U in zip(weights, matrices):
        S, Q = sum(u), dot(u, u)
        c = matvec(U, u)
        tau, R = dot(u, c), dot(c, c)
        # Independent field convolution, including the zero coefficient.
        field = {x: u[labels[pow(x, n, p)]] for x in range(1, p)}
        field = {x: value for x, value in field.items() if value}
        convolution = Counter()
        for x, vx in field.items():
            for y, vy in field.items():
                convolution[(x + y) % p] += vx * vy
        assert convolution[0] == n * Q
        for x in range(1, p):
            assert convolution[x] == c[labels[pow(x, n, p)]]
        assert sum(value * convolution[(-x) % p] for x, value in field.items()) == n * tau
        assert sum(value * value for value in convolution.values()) == n * n * Q * Q + n * R
        U2 = matmul(U, U)
        U3, U4 = matmul(U2, U), matmul(U2, U2)
        assert trace(U3) == n * n * (n - 1) * S ** 3 + (p - 3 * n) * tau
        assert trace(U4) == n ** 3 * (n - 1) * S ** 4 + n * (p - 2 * n) * Q ** 2 + (p - 4 * n) * R
        q = [sum(u[i] * u[j] * C[(i - k) % m][(j - k) % m]
                 for i, j in product(range(m), repeat=2)) for k in range(m)]
        Z = n * S * S - Q
        assert sum(q) == Z and dot(u, q) == tau
        coefficients = [F(x, beta) - F(alpha * Z, beta * D) for x in q]
        projection = L(coefficients)
        residual = [[F(u[i] * u[j]) - projection[i][j] for j in range(m)] for i in range(m)]
        residual_norm = sum(x * x for row in residual for x in row)
        centered_norm = sum((F(x) - F(Z, m)) ** 2 for x in q)
        assert Q ** 2 == F(Z * Z, m * D) + centered_norm / beta + residual_norm
        for k in range(m):
            assert sum(residual[i][j] * C[(i - k) % m][(j - k) % m]
                       for i, j in product(range(m), repeat=2)) == 0
        main = F(S * Z, m)
        bound_square = (Q - F(S * S, m)) * beta * (Q ** 2 - F(Z * Z, m * D))
        assert (tau - main) ** 2 <= bound_square
        traces.append(dict(sum_weights=S, norm_squared=Q, triangle=tau,
                           cubic=trace(U3), quartic=trace(U4)))
        frame_checks.append(dict(residual_norm_squared=str(residual_norm),
                                 centered_bound_square=str(bound_square), passed=True))
    return dict(p=p, n=n, m=m, algebra_weight_pairs=product_checks,
                algebra_entries_checked=product_checks * m * m,
                trace_checks=traces, exact_projection_checks=frame_checks,
                field_convolution_checks=len(weights),
                all_passed=True)


def main():
    cases = [example(p, n) for p, n in [(13, 2), (97, 8), (257, 16), (353, 16)]]
    # Power comparison in the quartic window: Q=A4 <= n^(8/3).
    assert 1 + 2 + F(3, 2) * F(8, 3) == 7
    assert 7 > F(86, 15) > F(17, 3)
    # The shifted-product lower bound at the candidate scale is insufficient.
    assert F(4, 5) * F(11, 9) == F(44, 45) < 1 < F(6, 5)
    # A hypothetical cubic product bound would exclude that scale.
    assert 3 * F(4, 5) - (1 + F(6, 5)) == F(1, 5)
    result = dict(scope='Exact weighted matrix identities and source-power comparisons; no new uniform bound.',
                  cases=cases, all_passed=True)
    target = ROOT / 'results/parallel45_weighted_matrix_algebra_2026_09_06.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_passed=True, cases=len(cases),
        algebra_entries_checked=sum(c['algebra_entries_checked'] for c in cases),
        trace_identities=2 * sum(len(c['trace_checks']) for c in cases),
        exact_projection_checks=sum(len(c['exact_projection_checks']) for c in cases),
        output=str(target))))


if __name__ == '__main__':
    main()
