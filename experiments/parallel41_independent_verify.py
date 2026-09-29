#!/usr/bin/env python3
"""Bounded, exact pass41 checks; finite evidence is not the uniform proof."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import isqrt
from pathlib import Path
import json

from parallel40_independent_verify import determinant

ROOT = Path(__file__).resolve().parents[1]


def energy(values, p):
    counts = Counter(x * pow(y, -1, p) % p for x in values for y in values)
    return sum(v * v for v in counts.values())


def triangle(p, n):
    g = next(x for x in range(2, p)
             if pow(x, n, p) == 1 and pow(x, n // 2, p) != 1)
    H = {pow(g, j, p) for j in range(n)}
    assert len(H) == n and p - 1 in H and n * n < p
    f = Counter((x - y) % p for x in H for y in H)
    a = Counter(pow((h - 1) % p, n, p) for h in H if h != 1)
    assert all(f[x] == a[pow(x, n, p)] for x in range(1, p))
    rho = Counter((pow(x, n, p), pow(x + 1, n, p)) for x in range(1, p - 1))
    R0 = {(h - 1) % p for h in H if h != 1}
    B = energy(R0, p)
    X = B - (2 * (n - 1) ** 2 - (n - 1))
    assert X == sum(v * (v - 1) * (v - 2) for v in rho.values())
    fall3 = lambda v: v * (v - 1) * (v - 2)
    X_dist = sum(fall3(v) for (u, w), v in rho.items() if len({1, u, w}) == 3)
    assert X_dist == X - 3 * sum(fall3(v) for v in a.values()) + 2 * fall3(a[1])
    assert sum(fall3(v) for v in a.values()) <= X
    assert 2 * sum(v ** 3 for v in a.values()) <= 8 * (n - 1) + 9 * X
    for (u, v), count in rho.items():
        assert rho[v, u] == count
        assert rho[pow(u, -1, p), v * pow(u, -1, p) % p] == count
    b = {c: v * v for c, v in a.items()}
    T = Counter()
    for c, d, e in product(b, repeat=3):
        ci = pow(c, -1, p)
        T[d * ci % p, e * ci % p] += b[c] * b[d] * b[e]
    W = sum(f[x] ** 2 * f[y] ** 2 * f[(x - y) % p] ** 2
            for x in f if x for y in f if y and x != y)
    assert W == n * sum(v * rho[k] for k, v in T.items())
    A2 = sum(b.values())
    assert sum(T.values()) == A2 ** 3
    assert sum(v for (u, w), v in T.items() if len({1, u, w}) == 3) == (
        A2 ** 3 - 3 * A2 * sum(v ** 4 for v in a.values())
        + 2 * sum(v ** 6 for v in a.values()))
    S = Counter()
    for x, y in product(a, repeat=2):
        S[x * pow(y, -1, p) % p] += a[x] * a[y]
    K = sum(v * v for v in S.values())
    assert K <= n * B
    low = n * sum(v * rho[k] for k, v in T.items() if rho[k] <= 2)
    high = W - low
    assert low <= 2 * n * A2 ** 3
    assert 2 * sum(v ** 3 for v in rho.values() if v >= 3) <= 9 * X
    # A rational lower enclosure of sum T^(3/2) gives a conservative RHS.
    scale = 10 ** 12
    s_lower = Q(sum(isqrt(v ** 3 * scale ** 2) for v in T.values()), scale)
    assert 2 * high ** 3 <= 9 * n ** 3 * X * s_lower ** 2
    levels = {}
    for c, v in a.items():
        levels.setdefault(1 << (v.bit_length() - 1), set()).add(c)
    correlations = {t: Counter(x * pow(y, -1, p) % p for x in d for y in d)
                    for t, d in levels.items()}
    for t, d in levels.items():
        Ri = {x for x in R0 if pow(x, n, p) in d}
        Ei = energy(Ri, p)
        assert t ** 4 * sum(v * v for v in correlations[t].values()) <= n * Ei
        assert Ei <= B
    blocks = 0
    for ti, tj, tk in product(levels, repeat=3):
        di, dj, dk = levels[ti], levels[tj], levels[tk]
        m = Counter((x * pow(z, -1, p) % p, y * pow(z, -1, p) % p)
                    for x in di for y in dj for z in dk)
        Pd = len(di) * len(dj) * len(dk)
        PT = ti * tj * tk
        L2sq = sum(v * v for v in m.values())
        assert sum(m.values()) == Pd
        assert L2sq == sum(v * correlations[tj][q] * correlations[tk][q]
                          for q, v in correlations[ti].items())
        assert L2sq ** 3 * PT ** 4 <= (n * B) ** 3 * Pd
        blocks += 1
    return dict(p=p, n=n, W=W, low=low, high=high, B=B, X=X,
                X_dist=X_dist, K=K, level_blocks=blocks, identities_passed=True)


def shell():
    data = json.loads((ROOT / 'results/parallel41_shell_real_recovery_2026_09_06.json').read_text())
    p, N = data['p'], data['N']
    def norm(v):
        return determinant([[v[(i - j) % N] * (-1 if i < j else 1)
                             for j in range(N)] for i in range(N)])
    def mul(v, w):
        return [sum(v[(i - j) % N] * (-1 if i < j else 1) * w[j]
                    for j in range(N)) for i in range(N)]
    def bar(v):
        return [v[0]] + [-v[N - j] for j in range(1, N)]
    h, c, F, f = (data[k] for k in ('h', 'c_p_over_h', 'F', 'recovered_f'))
    assert h == bar(h) and mul(h, c) == [p] + [0] * (N - 1)
    assert norm(h) == p * p and norm(f) == p * 1217
    assert mul(h, bar(F)) == [p * v for v in f]
    assert mul(c, bar(f)) == F
    for key, expected in [('F', 1217), ('reverse_raw_F', 641),
                          ('reverse_centered_F', 178771841)]:
        assert norm(data[key]) == p ** (N - 1) * expected
    roots = [pow(2, j, p) for j in range(1, 64, 2)]
    evaluate = lambda v, x: sum(y * pow(x, j, p) for j, y in enumerate(v)) % p
    assert [x for x in roots if not evaluate(h, x)] == [2, pow(2, -1, p)]
    assert not evaluate(f, pow(2, -1, p))
    for i in (1, 7, 17, 31):
        # Test the lattice bijection beyond the displayed actual coset vector.
        f0 = [0] * N
        f0[i] = 1
        f0[0] -= pow(pow(2, -1, p), i, p)
        F0 = mul(c, bar(f0))
        assert mul(h, bar(F0)) == [p * v for v in f0]
        assert all(not evaluate(F0, x) for x in roots if x != pow(2, -1, p))
    unit = [1, 1] + [0] * (N - 3) + [-1]
    assert norm(unit) == 1
    h_work, f_work = h[:], f[:]
    unit_energies = []
    for h_expected, f_expected in [(19, 13), (53, 43), (305, 293)]:
        assert norm(h_work) == p * p and norm(f_work) == 1217 * p
        assert sum(v * v for v in h_work) == h_expected
        assert sum(v * v for v in f_work) == f_expected
        unit_energies.append([h_expected, f_expected])
        h_work, f_work = mul(unit, h_work), mul(unit, f_work)
    assert 5 ** 16 - 2 ** 37 == 15148937153 > 0
    return dict(norm_method='32-by-32 integer Bareiss determinants',
                recovered_cofactor=1217, centering_defects=[641, 178771841],
                additional_lattice_vectors=4, unit_energies=unit_energies, all_passed=True)


def smoothing(p, V, t):
    chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)]
    F = [sum(chi[(x - v) % p] for v in V) for x in range(p)]
    M = {j: sum(v ** j for v in F) for j in range(1, 7)}
    assert M[1] == 0
    h6 = d6 = d2 = 0
    for shifts in product(range(p), repeat=t):
        for x in range(p):
            h = sum(F[(x - s) % p] for s in shifts)
            d = t * F[x] - h
            h6 += h ** 6
            d6 += d ** 6
            d2 += d * d
    actual = Q(h6, p ** t * t ** 6)
    predicted = (Q(M[6], t ** 5) + Q(15 * (t - 1) * M[4] * M[2], p * t ** 5)
                 + Q(10 * (t - 1) * M[3] ** 2, p * t ** 5)
                 + Q(15 * (t - 1) * (t - 2) * M[2] ** 3, p ** 2 * t ** 5))
    assert actual == predicted
    assert Q(d6, p ** t * t ** 6) >= M[6]
    assert Q(d2, p ** t * t ** 2) == Q(t + 1, t) * M[2]
    n = len(V)
    if n ** 4 <= p and t ** 5 >= n and t >= 2:
        assert actual <= Q(325, 32) * p * n ** 3
    return dict(p=p, V=V, t=t, third_moment=M[3], shift_tuples=p ** t,
                smoothed_sixth=str(actual), all_passed=True)


def quartic_excess(p, n, g, expected):
    assert all(p % d for d in range(2, isqrt(p) + 1))
    H = {pow(g, j, p) for j in range(n)}
    assert len(H) == n and p - 1 in H
    R0 = {(h - 1) % p for h in H if h != 1}
    a = Counter(pow(x, n, p) for x in R0)
    B = energy(R0, p)
    X = B - 2 * (n - 1) ** 2 + (n - 1)
    fall3 = lambda v: v * (v - 1) * (v - 2)
    diagonal = 3 * sum(fall3(v) for v in a.values()) - 2 * fall3(a[1])
    assert [X, diagonal, X - diagonal] == expected
    return dict(p=p, n=n, X=X, diagonal=diagonal, X_dist=X-diagonal,
                method='Direct shifted-ratio counts and the proved diagonal identity')


def main():
    assert Q(1) + Q(2, 3) + Q(3, 5) + Q(52, 15) == Q(86, 15)
    assert Q(1, 3) + Q(1, 5) + 3 == Q(53, 15)
    assert Q(1, 3) + Q(1, 5) == Q(8, 15)
    assert (Q(86, 15) + 9) / 6 == Q(221, 90)
    assert (Q(8, 15) + 1) / 6 == Q(23, 90)
    assert Q(221, 90) - Q(49, 20) == Q(1, 180)
    assert Q(8, 5) + Q(31, 15) * Q(61, 31) == Q(17, 3)
    result = dict(triangles=[triangle(p, n) for p, n in
                            [(17, 4), (41, 4), (97, 8), (193, 8), (257, 16), (1153, 32)]],
                  quartic_excess=[quartic_excess(*args) for args in
                      [(6700417, 64, 2, [114, 78, 36]),
                       (215535361, 128, 25525303, [1740, 1056, 684]),
                       (67403009, 128, 64701253, [720, 216, 504])]],
                  shell=shell(), smoothing=[smoothing(*x) for x in
                      [(13, [0, 1, 3], 3), (41, [0, 1], 2), (17, [0, 1], 2)]],
                  scope='Finite exact checks plus rational exponent audit; no uniform theorem certification.',
                  all_passed=True)
    dest = ROOT / 'results/parallel41_independent_verification_2026_09_06.json'
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
