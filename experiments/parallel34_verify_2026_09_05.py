#!/usr/bin/env python3
"""Exact finite checks of the opposite-pair transform, not a Paley proof."""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, isqrt, prod
from pathlib import Path
from time import perf_counter
import json

ROOT = Path(__file__).resolve().parents[1]
CHECKS = Counter()


def check(ok, kind, context=None):
    if not ok:
        raise AssertionError((kind, context))
    CHECKS[kind] += 1


def rat(q):
    q = Fraction(q)
    return dict(numerator=q.numerator, denominator=q.denominator)


def falling(n, r):
    return prod(range(n-r+1, n+1)) if r <= n else 0


def gaussian(n, s):
    return prod(range(1, 2*s, 2))*n**s


def cycle(M, j):
    if not j:
        return 1
    assert M >= 2*j >= 2
    return comb(M-j, j)+comb(M-j-1, j-1)


def forward(N, s, t):
    return factorial(2*s)//factorial(2*t)*comb(N-2*t, s-t)


def inverse(N, s, t):
    return (-1)**(s-t)*factorial(2*s)//factorial(2*t)*cycle(N-2*t, s-t)


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def subgroup(p, n):
    assert prime(p) and (p-1) % n == 0 and n % 2 == 0
    for g in range(2, p):
        h = pow(g, (p-1)//n, p)
        H = sorted({pow(h, j, p) for j in range(n)})
        if len(H) == n:
            return H
    raise AssertionError('missing generator')


def distinct_dp(A, p, limit):
    """One element at a time: coefficients of product(1+z X^a)."""
    rows = [{0: 1}]+[{} for _ in range(limit)]
    for count, a in enumerate(A, 1):
        for r in range(min(count, limit), 0, -1):
            target = rows[r]
            for b, v in rows[r-1].items():
                key = (a+b) % p
                target[key] = target.get(key, 0)+v
    return rows


def opposite_free_dp(A, p, limit):
    """One opposite class at a time, with choices none, +a, or -a."""
    representatives = [a for a in A if a < (-a) % p]
    assert len(representatives)*2 == len(A)
    rows = [{0: 1}]+[{} for _ in range(limit)]
    for count, a in enumerate(representatives, 1):
        for r in range(min(count, limit), 0, -1):
            target = rows[r]
            for b, v in rows[r-1].items():
                for key in [(b+a) % p, (b-a) % p]:
                    target[key] = target.get(key, 0)+v
    return rows


def all_word_counts(A, p, limit):
    row = {0: 1}
    zeros = [1]
    for r in range(1, limit+1):
        next_row = {}
        for b, v in row.items():
            for a in A:
                key = (b+a) % p
                next_row[key] = next_row.get(key, 0)+v
        row = next_row
        check(sum(row.values()) == len(A)**r, 'full word totals')
        zeros.append(row.get(0, 0))
    return zeros


def case(p, A, limit, label, is_subgroup):
    n = len(A)
    N = n//2
    assert 0 not in A and len(set(A)) == n and set(A) == {(-a) % p for a in A}
    assert limit <= N and limit < p and limit % 2 == 0
    I = distinct_dp(A, p, limit)
    O = opposite_free_dp(A, p, limit)
    E = all_word_counts(A, p, limit)
    for r in range(limit+1):
        check(sum(I[r].values())*factorial(r) == falling(n, r), 'distinct totals')
        check(sum(O[r].values())*factorial(r) == 2**r*falling(N, r), 'opposite-free totals')
    iz = [I[2*s].get(0, 0)*factorial(2*s) for s in range(limit//2+1)]
    oz = [O[2*s].get(0, 0)*factorial(2*s) for s in range(limit//2+1)]
    Q = [Fraction(iz[s])-Fraction(falling(n, 2*s), p) for s in range(limit//2+1)]
    B = [Fraction(oz[s])-Fraction(2**(2*s)*falling(N, 2*s), p) for s in range(limit//2+1)]
    rows = []
    for s in range(limit//2+1):
        r = 2*s
        check(iz[s] == sum(forward(N, s, t)*oz[t] for t in range(s+1)), 'forward raw')
        check(oz[s] == sum(inverse(N, s, t)*iz[t] for t in range(s+1)), 'inverse raw')
        check(Q[s] == sum(forward(N, s, t)*B[t] for t in range(s+1)), 'forward centered')
        check(B[s] == sum(inverse(N, s, t)*Q[t] for t in range(s+1)), 'inverse centered')
        row = dict(s=s, r=r, energy=E[r], distinct=iz[s], opposite_free=oz[s],
                   centered_distinct=rat(Q[s]), centered_opposite_free=rat(B[s]),
                   gaussian=gaussian(n, s))
        if s >= 2:
            mark = comb(r, 2)
            check(E[r-2]**s <= E[r]**(s-1), 'normalized Lyapunov')
            check(E[r]*p >= n**r, 'principal lower bound')
            check(E[r] >= n**s, 'Jensen lower bound')
            check(iz[s]-oz[s] <= mark*n*E[r-2], 'distinct opposite union bound')
            repeated = E[r]-iz[s]
            check(n*repeated**2 <= mark**2*E[r]**2, 'symmetric-set repeated fraction')
            # For combined exclusion, the repeated part is exact here and
            # the opposite part is a union bound, so no roots are evaluated.
            check(oz[s] >= E[r]-repeated-mark*n*E[r-2], 'combined exclusion')
            row.update(repeated=repeated, marked_opposite_upper=mark*n*E[r-2])
        rows.append(row)
    check(oz[1] == 0 and B[1] == -Fraction(n*(n-2), p), 'degree two base')
    return dict(p=p, n=n, label=label, subgroup=is_subgroup, elements=A, rows=rows)


def main():
    started = perf_counter()
    cases = []
    for p, n, limit in [(17, 16, 8), (97, 16, 8), (193, 16, 8),
                         (193, 32, 12), (257, 16, 8), (257, 32, 12),
                         (257, 64, 12), (33713, 16, 8)]:
        cases.append(case(p, subgroup(p, n), limit, 'multiplicative subgroup', True))
    A = sorted({x % 257 for a in [1, 2, 3, 5, 8, 13, 21, 34, 55, 89] for x in [a, -a]})
    check(any(a*b % 257 not in A for a in A for b in A), 'non-subgroup witness')
    cases.append(case(257, A, 10, 'symmetric non-subgroup', False))

    cycles = []
    for M in range(2, 17):
        counts = Counter()
        for mask in range(1 << M):
            shifted = ((mask << 1) | (mask >> (M-1))) & ((1 << M)-1)
            if not mask & shifted:
                counts[mask.bit_count()] += 1
        for j in range(M//2+1):
            check(counts[j] == cycle(M, j), 'cycle enumeration', (M, j))
        cycles.append(dict(M=M, independent_set_counts=dict(sorted(counts.items()))))

    coefficient_rows = []
    for N in list(range(1, 33))+[64, 128]:
        depth = min(N//2, 32)
        for s in range(depth+1):
            forward_norm = 0
            inverse_norm = 0
            for t in range(s+1):
                bound = comb(s, t)*gaussian(2*N, s)//gaussian(2*N, t)
                check(forward(N, s, t) <= bound, 'forward Gaussian coefficient')
                check(abs(inverse(N, s, t)) <= bound, 'inverse Gaussian coefficient')
                forward_norm += forward(N, s, t)*gaussian(2*N, t)
                inverse_norm += abs(inverse(N, s, t))*gaussian(2*N, t)
            for u in range(s+1):
                check(sum(inverse(N, s, j)*forward(N, j, u) for j in range(u, s+1))
                      == int(s == u), 'triangular inverse identity', (N, s, u))
            check(max(forward_norm, inverse_norm) <= 2**s*gaussian(2*N, s), 'Gaussian row norm')
        coefficient_rows.append(dict(N=N, maximum_s=depth))

    # An explicit sufficient threshold for every eligible n >= 2^35.
    repeated_upper = Fraction(1, 2**17)
    opposite_upper = Fraction(1, 128)
    excluded_upper = 45*(repeated_upper+opposite_upper)
    check(excluded_upper < Fraction(1, 2), 'ten-term threshold')
    check((2**35)**2 >= (2**17)**4, 'square-root threshold')
    check(2**35 == 128**5, 'fifth-root threshold')
    inputs = ['research/parallel33-centered-distinct-moments-2026-09-05.md',
              'experiments/parallel34_verify_2026_09_05.py',
              'results/parallel34_prior_state_2026_09_05.json']
    result = dict(status='all exact finite checks passed; full targets unproved',
                  created_at_utc=datetime.now(timezone.utc).isoformat(),
                  duration_seconds=perf_counter()-started,
                  input_sha256={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
                  check_counts=dict(CHECKS), cases=cases, cycle_enumeration=cycles,
                  coefficient_parameter_ranges=coefficient_rows,
                  ten_term_threshold=dict(n_at_least=2**35, p_at_most='n^4',
                       excluded_fraction_upper=rat(excluded_upper), conclusion='O_10 >= n^6/2',
                       scope='For any eligible symmetric subset of the nonzero prime field; no prime existence claim.'),
                  scope='Finite checks of independent counting implementations and exact bounds. The uniform proofs are ordinary mathematics; no Lean or separate-author verification.')
    (ROOT/'results/parallel34_verification_2026_09_05.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],checks=sum(CHECKS.values()),
          cases=len(cases),even_order_cases=sum(len(x['rows']) for x in cases),
          duration_seconds=result['duration_seconds'],check_counts=dict(CHECKS)),indent=2))


if __name__ == '__main__':
    main()
