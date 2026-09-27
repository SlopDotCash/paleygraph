#!/usr/bin/env python3
"""Exact finite checks of reflection smoothing and its paired recovery cost."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json


def check(p, V, A, r):
    assert p % 4 == 1
    n, m, t = len(V), len(A), 2 ** r
    chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)]
    F = [sum(chi[(x - v) % p] for v in V) for x in range(p)]
    indicator = [int(x in A) for x in range(p)]
    M = {j:sum(v ** j for v in F) for j in range(1, 7)}
    assert M[1] == 0
    S = sum(F[a] for a in A)
    b6 = pair = a2 = 0
    for centers in product(range(p), repeat=r):
        fn, an = F[:], indicator[:]
        for c in centers:
            fn = [fn[x] + fn[(2 * c - x) % p] for x in range(p)]
            an = [an[x] + an[(2 * c - x) % p] for x in range(p)]
        assert sum(fn) == 0 and sum(an) == t * m
        b6 += sum(v ** 6 for v in fn)
        pair += sum(x * y for x, y in zip(fn, an))
        a2 += sum(v * v for v in an)
    samples = p ** r
    B = Q(b6, samples * t ** 6)
    assert Q(pair, samples * t * t) == Q(S, t)
    assert Q(a2, samples * t * t) == Q(m * m, p) + Q(m, t) - Q(m * m, p * t)
    assert B >= Q(M[6], t ** 5)
    if r == 1:
        assert B == Q(M[6], 32) + Q(15 * M[4] * M[2], 32 * p) + Q(5 * M[3] ** 2, 16 * p)
    recovery = t ** 5 * (Q(1, m) + Q(t - 1, p)) * B / n ** 6
    direct = Q(M[6], m * n ** 6)
    assert recovery >= direct * (1 + Q(m * (t - 1), p))
    assert Q(S, m * n) ** 6 <= recovery
    if n ** 4 <= p:
        bound = Q(M[6], t ** 5) + Q(10 * p * n ** 3, t) * (1 - Q(1, t ** 4))
        assert B <= bound
        if t >= 2 and t ** 5 >= n:
            assert B <= Q(325, 32) * p * n ** 3
    return dict(p=p, n=n, V=V, A=A, reflections=r, averaging_size=t,
                center_tuples=samples, initial_pairing=S, initial_third_moment=M[3],
                expected_sixth=str(B), direct_sixth_certificate=str(direct),
                reflection_recovery_certificate=str(recovery), all_passed=True)


def main():
    specs = [(13,[0,1,3],[0,1,2,4,8],1), (13,[0,1,3],[0,1,2,4,8],2),
             (17,[0,1],[0,1,2,4,8],1), (17,[0,1],[0,1,2,4,8],2),
             (41,[0,1],[0,2,3,5,7,11],2), (97,[0,1,3],[0,1,2,4,8,16],1)]
    cases = [check(*s) for s in specs]
    result = dict(scope='Finite reflection-projection identities and a specific recovery-certificate comparison; no original Paley bound.',
                  cases=cases, center_tuples=sum(c['center_tuples'] for c in cases), all_passed=True)
    dest = Path(__file__).resolve().parents[1] / 'results/parallel42_classical_affine_2026_09_06.json'
    dest.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(all_passed=True, center_tuples=result['center_tuples'], result=str(dest))))


if __name__ == '__main__':
    main()
