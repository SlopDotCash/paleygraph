#!/usr/bin/env python3
"""Exact integer audit of the quartic-star reformulation; no asymptotic bound test."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import comb, isqrt
from pathlib import Path
import json
import random
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/parallel21_positive_upper_review_2026_09_05.json'
COUNTS = Counter()

def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))

def chi_table(p):
    a = np.full(p, -1, dtype=np.int64)
    a[0] = 0
    a[(np.arange(1, p, dtype=np.int64)**2) % p] = 1
    return a

def dot(a, b):
    # Explicit integer bounds below keep every dot product strictly in int64.
    assert len(a) * max((abs(int(x)) for x in a), default=0) * max((abs(int(x)) for x in b), default=0) < 2**63
    return int(a @ b)

def transform(a, chi):
    p = len(chi)
    assert p * max((abs(int(x)) for x in a), default=0) < 2**63
    x = np.arange(p, dtype=np.int64)
    out = np.empty(p, dtype=np.int64)
    for start in range(0, p, 96):
        t = np.arange(start, min(start+96, p), dtype=np.int64)
        out[start:start+len(t)] = chi[(x[None, :]-t[:, None]) % p] @ a
    return out

def feature(Q, chi):
    x = np.arange(len(chi), dtype=np.int64)
    v = np.ones(len(chi), dtype=np.int64)
    for q in Q:
        v *= chi[(x-q) % len(chi)]
    return v

def sidon(C, p):
    sums = [(a+b) % p for a,b in combinations_with_replacement(C, 2)]
    return len(set(sums)) == len(sums)

def greedy_sidon(p, n, seed):
    rng = random.Random(seed)
    candidates = list(range(p))
    rng.shuffle(candidates)
    C = []
    for a in candidates:
        if sidon(C+[a], p):
            C.append(a)
            if len(C) == n:
                return tuple(sorted(C))
    raise AssertionError('greedy Sidon fixture did not complete')

def case(C, p, full_star=False, gram=False):
    assert prime(p)
    n = len(C)
    assert n >= 3 and len(set(C)) == n
    chi = chi_table(p)
    x = np.arange(p, dtype=np.int64)
    columns = [chi[(x-c) % p] for c in C]
    es = [np.ones(p, dtype=np.int64)] + [np.zeros(p, dtype=np.int64) for _ in range(6)]
    for j,col in enumerate(columns):
        for d in range(min(6, j+1), 0, -1):
            es[d] += col * es[d-1]
    F = sum(columns)
    N = np.full(p, n, dtype=np.int64)
    N[list(C)] -= 1
    assert np.array_equal(6*es[3], F**3-(3*N-2)*F)
    COUNTS['cubic_pointwise_checks'] += p
    assert all(p * comb(n, d) < 2**63 for d in range(7) if d <= n)
    T = {d:int(es[d].sum()) for d in [2,3,4,6]}
    assert T[2] == -comb(n,2)
    # Independent distinct-subset expansion, including the genuinely 6-root terms.
    for d in [3,4,6]:
        total = 0
        for Q in combinations(C,d):
            K = int(feature(Q,chi).sum())
            assert K*K <= (d-1)**2*p
            total += K
            COUNTS[f'direct_degree_{d}_correlations'] += 1
        assert total == T[d]
    U = transform(es[3], chi)
    E = dot(es[3],es[3])
    L = dot(U,U)
    assert L == p*E-T[3]**2
    assert int(U.sum()) == 0
    Z = sum(int(es[2][c])**2 for c in C)
    R = 6*(n-4)*T[4]-(n-2)*(n-3)*comb(n,2)+p*comb(n,3)-Z
    assert E == 20*T[6]+R
    assert sum(int(U[c]) for c in C) == 4*T[4]-n*comb(n-1,2)-sum(int(es[2][c]) for c in C)
    COUNTS['complete_energy_identities'] += 1
    COUNTS['overlap_ledgers'] += 1
    COUNTS['on_set_sum_identities'] += 1
    qlist = list(combinations(C,3))
    rows = range(p) if full_star else list(C)+list(t for t in range(p) if t not in C)[:4]
    for t in rows:
        actual = sum(int(feature(Q+(t,),chi).sum()) for Q in qlist)
        assert actual == int(U[t])
        COUNTS['direct_star_rows'] += 1
        if t in C:
            outside = sum(int(feature(Q+(t,),chi).sum()) for Q in combinations([c for c in C if c != t],3))
            assert int(U[t]) == outside-comb(n-1,2)-int(es[2][t])
            COUNTS['repeated_root_row_identities'] += 1
    if gram:
        # Every pair of cubic features, including all overlap sizes.
        fs = [feature(Q,chi) for Q in qlist]
        ts = [transform(f,chi) for f in fs]
        for i,Q in enumerate(qlist):
            for j,RQ in enumerate(qlist):
                I = set(Q)&set(RQ)
                D = tuple(sorted(set(Q)^set(RQ)))
                fD = feature(D,chi)
                inner = int(fD.sum())-sum(int(fD[a]) for a in I)
                rhs = p*inner-int(fs[i].sum())*int(fs[j].sum())
                assert dot(ts[i],ts[j]) == rhs
                COUNTS[f'quartic_family_gram_overlap_{len(I)}'] += 1
    Lin = sum(int(U[c])**2 for c in C)
    Lout = L-Lin
    M6 = sum(int(y)**6 for y in F)
    if n >= 6 and n**4 <= p:
        assert abs(R) <= p*n**3
        assert Lin <= p*p*n**3
        assert 9*T[3]**2 <= p*p*n**3
        # Two quantitative implications with the observed ratio substituted.
        assert Lout <= 20*p*T[6]+p*p*n**3
        assert 20*p*T[6] <= Lout+3*p*p*n**3
        assert T[6] >= -Fraction(p*n**3,20)
        COUNTS['slice_constant_audits'] += 1
    COUNTS['cases'] += 1
    return dict(p=p,n=n,C=list(C),sidon=sidon(C,p),T3=T[3],T4=T[4],T6=T[6],R=R,
                cubic_feature_energy=E,quartic_star_energy=L,
                outside_star_energy=Lout,on_set_star_energy=Lin,M6=M6,
                M6_over_pn3=str(Fraction(M6,p*n**3)),
                outside_energy_over_p2n3=str(Fraction(Lout,p*p*n**3)))

def main():
    # Exhaustive inputs in these two small fields, plus varied fields and sizes.
    for p in [5,7]:
        chi=chi_table(p)
        x=np.arange(p,dtype=np.int64)
        P=chi[(x[None,:]-x[:,None])%p]
        assert np.array_equal(P.T@P,p*np.eye(p,dtype=np.int64)-np.ones((p,p),dtype=np.int64))
        COUNTS['full_translate_matrix_checks'] += p*p
        for n in range(3,p+1):
            for C in combinations(range(p),n):
                case(C,p,full_star=True)
    rng=random.Random(210905)
    for p in [11,13,17,19,23,29,31,41]:
        for n in range(3,min(8,p)+1):
            for rep in range(2):
                C=tuple(sorted(rng.sample(range(p),n)))
                case(C,p,full_star=True,gram=(n==6 and rep==0 and p<=19))
    fixtures=[]
    for p in [1297,2411,4099,6563,10007]:
        n=isqrt(isqrt(p))
        C=greedy_sidon(p,n,p+21)
        assert n**4<=p<(n+1)**4
        fixtures.append(case(C,p))
    inputs={str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in [
        Path(__file__),
        ROOT/'research/parallel21-squarefree-moments-2026-09-05.md',
        ROOT/'research/parallel20-inversion-moments-2026-09-05.md',
    ]}
    result={'status':'exact finite identities passed; missing uniform upper bound remains unproved',
            'arithmetic':'integer numpy operations with checked int64 bounds; rational comparisons; no floating acceptance',
            'counts':dict(sorted(COUNTS.items())),
            'slice_fixtures':fixtures,'input_sha256':inputs}
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'counts':result['counts'],'result':str(OUT)},indent=2))

if __name__=='__main__':
    main()
