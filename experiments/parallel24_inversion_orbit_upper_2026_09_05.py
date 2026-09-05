#!/usr/bin/env python3
"""Exact finite-field checks for the pass24 unsigned inversion theorem."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, isqrt, prod
from pathlib import Path
from random import Random
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()

def check(condition, label, multiplicity=1):
    assert condition, label
    COUNTS[label] += multiplicity

def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))

def characters(p):
    v = np.full(p, -1, dtype=np.int64)
    v[0] = 0
    for x in range(1, p):
        v[x*x % p] = 1
    return v

def bh(D, h, p):
    sums = [sum(t) % p for t in combinations_with_replacement(D, h)]
    return len(sums) == len(set(sums))

def qbound(k, h):
    return k+(2*h-2)*comb(comb(k+h-1, h), 2)

def ubound(p, k, r):
    return prod(range(1, 2*r, 2))*p*p*k**r+(2*r-1)**2*p*k**(2*r)

def case(p, C, expand=False):
    k = len(C)
    check(p*p*k**6 <= np.iinfo(np.int64).max, 'signed_int64_accumulation_safe')
    chi = characters(p)
    X = np.arange(p)
    A = chi[(X[:, None]-np.array(C)) % p]
    F = A.sum(axis=1)
    Q = A @ A.T
    moments = {r: int(np.sum(Q**(2*r))) for r in (1, 2, 3)}
    original = {r: int(np.sum(F**(2*r))) for r in (1, 2, 3)}
    check(original[1] == p*k-k*k, 'second_moment')
    for r in moments:
        check(moments[r] <= ubound(p, k, r), 'two_row_upper')
    pole_rows = []
    for z in range(p):
        if z in C:
            continue
        D = [pow((c-z) % p, -1, p) for c in C]
        w = chi[(np.array(C)-z) % p]
        AD = chi[(X[:, None]-np.array(D)) % p]
        FD, FW = AD.sum(axis=1), AD @ w
        nz = X[X != z]
        inv = np.array([pow(int(x-z) % p, -1, p) for x in nz])
        check(np.array_equal(FD[inv], chi[(nz-z) % p]*Q[nz, z]),
              'unsigned_inversion_coordinates', p-1)
        check(FD[0] == F[z] and Q[z, z] == k, 'unsigned_boundaries')
        check(np.array_equal(FW[inv], chi[-1]*chi[(nz-z) % p]*F[nz]),
              'signed_inversion_coordinates', p-1)
        check(FW[0] == chi[-1]*k, 'signed_boundary')
        row = {'z': z, 'bh': {h: bh(D, h, p) for h in (2, 3) if p > h}, 'M': {}}
        for r in moments:
            md, mw = int(np.sum(FD**(2*r))), int(np.sum(FW**(2*r)))
            rhs = int(np.sum(Q[:, z]**(2*r)))
            check(md == rhs-k**(2*r)+int(F[z])**(2*r), 'unsigned_moment_identity')
            check(md <= rhs, 'unsigned_column_upper')
            check(mw == original[r]-int(F[z])**(2*r)+k**(2*r), 'signed_moment_identity')
            check(mw >= original[r], 'signed_moment_dominates_original')
            row['M'][r] = md
        pole_rows.append(row)
    certificates = []
    for h in (2, 3):
        if p <= h:
            continue
        good = [row for row in pole_rows if row['bh'][h]]
        qb = qbound(k, h)
        check(p-len(good) <= qb, 'bad_pole_count')
        if p <= qb:
            continue
        check(bool(good), 'good_pole_exists')
        for r in moments:
            U = ubound(p, k, r)
            check(min(row['M'][r] for row in good)*(p-qb) <= U, 'bounded_good_representative')
            check(sum(row['M'][r] for row in good) <= U, 'good_pole_sum')
            for lam in (2, 3):
                count = sum(row['M'][r]*(p-qb) <= lam*U for row in good)
                check(lam*count >= (lam-1)*(p-qb), 'positive_fraction_good_poles')
        certificates.append({'h': h, 'Q_h': qb, 'good_poles': len(good),
                             'best_M6': min(row['M'][3] for row in good),
                             'M6_bound': str(Fraction(ubound(p, k, 3), p-qb))})
    if expand:
        for r in (1, 2, 3):
            total, even_count = 0, 0
            for ct in product(range(k), repeat=2*r):
                value = int(np.sum(np.prod(A[:, ct], axis=1)))
                even = all(v % 2 == 0 for v in Counter(ct).values())
                if even:
                    check(abs(value) <= p, 'square_tuple_bound')
                    even_count += 1
                else:
                    check(value*value <= (2*r-1)**2*p, 'individual_Weil_bound')
                total += value*value
            check(total == moments[r], 'ordered_tuple_expansion')
            check(even_count <= prod(range(1, 2*r, 2))*k**r, 'matching_tuple_count')
    return {'p': p, 'C': list(C), 'original_M6': original[3], 'certificates': certificates}

def main():
    rng = Random(2026090524)
    cases = []
    for p in (5, 7, 11):
        for k in range(1, 4):
            for C in combinations(range(p), k):
                cases.append(case(p, C, expand=(C == tuple(range(k)))))
    for p in (13, 17):
        for k in (1, 2):
            for C in combinations(range(p), k):
                cases.append(case(p, C))
    for p, k in ((31, 3), (41, 3), (101, 3), (257, 4), (1297, 6)):
        check(prime(p), 'sample_primality')
        for C in (list(range(k)), sorted(rng.sample(range(p), k))):
            cases.append(case(p, C))
    check(Fraction(prod(range(1, 6, 2)), 1)/(1-Fraction(1, 4)) == 20,
          'sixth_moment_asymptotic_constant')
    for h in range(2, 11):
        check(Fraction(h-1, prod(range(1, h+1))**2) <= Fraction(1, 4),
              'boundary_constants_sample')
    inputs = ['research/parallel24-inversion-orbit-upper-2026-09-05.md',
              'experiments/parallel24_inversion_orbit_upper_2026_09_05.py',
              'research/parallel20-inversion-moments-2026-09-05.md',
              'research/sigma-sumproduct-2026-09-05.md',
              'sources/mcdonald-sahay-wyman-2210.03789v2.html']
    report = {'status': 'passed', 'checks': dict(COUNTS), 'cases': len(cases),
              'coverage': 'All sets of sizes 1..3 for p=5,7,11; all sizes 1..2 for p=13,17; ten larger interval/random cases; exact integer arithmetic.',
              'limitations': 'Finite checks do not prove asymptotics. The theorem is unsigned existence per inversion family and does not bound the original moment or transported signs.',
              'input_sha256': {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs},
              'larger_cases': cases[-10:]}
    out = ROOT/'results/parallel24_inversion_orbit_upper_2026_09_05.json'
    out.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'checks', 'cases')}))

if __name__ == '__main__':
    main()
