#!/usr/bin/env python3
"""Exact correlation identities and an abstract sign-kernel obstruction.

The modified sign functions are NOT asserted to arise from elliptic
Kummer sheaves or to satisfy the ambient Paley projection identities.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import isqrt, prod
from pathlib import Path
import json
import numpy as np

from parallel14_elliptic_model_2026_09_05 import Curve, quotient
from parallel10_word_aggregate_2026_09_05 import positive_definite_integer
from parallel3_conic_2026_09_04 import reduce_poly, cyclic_product

ROOT = Path(__file__).resolve().parents[1]


def parity_vector(table, f, shifts):
    count = Counter(shifts)
    odd = sorted(t for t,c in count.items() if c % 2)
    even = sorted(t for t,c in count.items() if c % 2 == 0)
    L = len(f)
    neg = [next(k for k in range(L) if table[h][k] == 0) for h in range(L)]
    raw = [prod(f[table[h][t]] for t in shifts) for h in range(L)]
    clean = [prod(f[table[h][t]] for t in odd) for h in range(L)]
    restored = clean.copy()
    for t in even:
        restored[neg[t]] -= clean[neg[t]]
    assert raw == restored
    if not odd:
        assert raw == [int(h not in {neg[t] for t in count}) for h in range(L)]
    return odd,even,raw,clean


def small_case(p,r):
    E = Curve(p,r)
    sums = [[E.index[E.add(P,Q)] for Q in E.points] for P in E.points]
    fs = [E.f(P) for P in E.points]
    two = [i for i in range(len(E.points)) if sums[i][i] == 0]
    qr,data = quotient(E,sums,fs,two)
    table,f,coords,n,d,M = data
    L = len(f)
    neg = [next(k for k in range(L) if table[h][k] == 0) for h in range(L)]
    if L <= 8:
        patterns = list(combinations_with_replacement(range(L),4))
    else:
        # Systematic repeated and nonrepeated patterns, bounded in number.
        patterns = [(0,t,t,t) for t in range(L)]
        patterns += [(0,0,t,t) for t in range(L)]
        patterns += [(0,t,(3*t+1) % L,(5*t+2) % L) for t in range(L)]
        patterns += [tuple((j*(t+1)+j*j) % L for j in range(2+2*(t%5))) for t in range(L)]
    parity_checks = 0
    for shifts in patterns:
        odd,even,raw,clean = parity_vector(table,f,shifts)
        if odd:
            assert sum(clean)**2 <= len(odd)**2*p
            assert abs(sum(raw)-sum(clean)) <= len(even)
        parity_checks += 1
    # Certify every character twist via the singular values of convolution.
    supports = [tuple([0])]
    supports += [(0,t) for t in range(1,L)]
    supports += [(0,t,(2*t+1) % L) for t in range(1,min(L,10))]
    supports += [(0,t,(2*t+1) % L,(3*t+2) % L) for t in range(1,min(L,10))]
    supports = sorted(set(tuple(sorted(set(s))) for s in supports))
    minors = 0
    for shifts in supports:
        odd,even,raw,clean = parity_vector(table,f,shifts)
        assert len(odd) == len(shifts) and not even and raw == clean
        A = np.array([[clean[table[h][neg[k]]] for k in range(L)]
                      for h in range(L)],dtype=object)
        positive_definite_integer(len(shifts)**2*p*np.eye(L,dtype=object)-A@A.T)
        minors += L
    # All four-point row correlations, including the exceptional torsion rows.
    mat = np.array(M,dtype=object)
    square = mat@mat
    strata = Counter()
    for s in range(L):
        for t in range(L):
            odd,even,raw,clean = parity_vector(table,f,[s,neg[s],t,neg[t]])
            assert int(square[s,t]) == sum(raw)
            if not odd:
                assert int(square[s,t]) == L-len(set([s,neg[s],t,neg[t]]))
            else:
                assert sum(clean)**2 <= len(odd)**2*p
            strata[f'odd_{len(odd)}_even_{len(even)}'] += 1
    # Fourier convolution identity for products of translates, in Z[zeta_n].
    fourier_entries = 0
    if L <= 8:
        chars = list(product(range(n),range(d)))
        exps = {a:[(a[0]*coords[h][0]+(n//d)*a[1]*coords[h][1]) % n
                   for h in range(L)] for a in chars}
        def transform(v):
            out = {}
            for a in chars:
                poly = [0]*n
                for h in range(L):
                    poly[-exps[a][h] % n] += v[h]
                out[a] = poly
            return out
        F = transform(f)
        for shifts in [(0,), (0,1), (0,0,1), (0,1,2,3), (0,0,1,1)]:
            current = {(0,0):[1]+[0]*(n-1)}
            for t in shifts:
                updated = {a:[0]*n for a in chars}
                for a,A in current.items():
                    for b in chars:
                        phase = exps[b][t]
                        shifted = [F[b][(j-phase) % n] for j in range(n)]
                        term = cyclic_product(A,shifted,n)
                        c = ((a[0]+b[0]) % n,(a[1]+b[1]) % d)
                        updated[c] = [x+y for x,y in zip(updated[c],term)]
                current = updated
            raw = parity_vector(table,f,shifts)[2]
            direct = transform(raw)
            for a in chars:
                assert reduce_poly(current[a],n) == reduce_poly([L**(len(shifts)-1)*v for v in direct[a]],n)
                fourier_entries += 1
    return dict(p=p,r=r,quotient_order=L,parity_patterns=parity_checks,
                twisted_supports=len(supports),twisted_bound_leading_minors=minors,
                square_matrix_entries=L*L,square_strata=dict(strata),
                cyclotomic_product_entries=fourier_entries)


def mul(E,P,n):
    out = None
    while n:
        if n & 1:
            out = E.add(out,P)
        P = E.add(P,P)
        n //= 2
    return out


def prime_factors(n):
    out = []
    t = 2
    while t*t <= n:
        if n % t == 0:
            out.append(t)
            while n % t == 0:
                n //= t
        t += 1
    if n > 1:
        out.append(n)
    return out


def large_quotient(p,r):
    """Use the group 2G isomorphic to G/G[2], avoiding a quadratic group table."""
    E = Curve(p,r)
    values = {}
    multiplicities = Counter()
    for P in E.points:
        Q = E.add(P,P)
        if Q in values:
            assert values[Q] == E.f(P)
        values[Q] = E.f(P)
        multiplicities[Q] += 1
    assert set(multiplicities.values()) == {4}
    L = len(values)
    factors = prime_factors(L)
    def order(P):
        out = L
        for q in factors:
            while out % q == 0 and mul(E,P,out//q) is None:
                out //= q
        return out
    def cycle(P,k):
        out,h = [],None
        for _ in range(k):
            out.append(h)
            h = E.add(h,P)
        assert h is None and len(set(out)) == k
        return out
    points = list(values)
    for g in points[1:]:
        n = order(g)
        d = L//n
        if n < d or n % d:
            continue
        cg = cycle(g,n)
        sg = set(cg)
        for b in points:
            if order(b) != d:
                continue
            cb = cycle(b,d)
            if sg.intersection(cb) == {None}:
                break
        else:
            continue
        break
    else:
        raise AssertionError('No product coordinates found')
    coords = {E.add(P,Q):(i,j) for i,P in enumerate(cg) for j,Q in enumerate(cb)}
    assert len(coords) == L and n*d == L and n % 2 == d % 2 == 0
    for P,(i,j) in coords.items():
        assert coords[E.add(P,g)] == ((i+1) % n,j)
        assert coords[E.add(P,b)] == (i,(j+1) % d)
    f = {coords[P]:v for P,v in values.items()}
    assert f[0,0] == 0 and all(abs(v) == 1 for h,v in f.items() if h != (0,0))
    assert all(f[-i % n,-j % d] == v for (i,j),v in f.items())
    return E,n,d,f


def planted_case(p,r):
    E,n,d,f = large_quotient(p,r)
    L = len(f)
    N = isqrt(p)+1
    assert (N-1)**2 < p < N*N and L >= 16*N
    k = min(N,n//4)
    ell = (N+k-1)//k
    assert n >= 8 and ell <= d
    B = {(i,j) for i in range(1,k+1) for j in range(ell)}
    add = lambda a,b:((a[0]+b[0]) % n,(a[1]+b[1]) % d)
    neg = lambda a:(-a[0] % n,-a[1] % d)
    sub = lambda a,b:add(a,neg(b))
    negB = {neg(a) for a in B}
    assert not B.intersection(negB)
    special = {h for h in f if add(h,h) == (0,0)}
    assert len(special) == 4 and not (B|negB).intersection(special)
    D = {add(a,b) for a in B for b in B}|{sub(a,b) for a in B for b in B}
    D |= {neg(a) for a in D}
    m = len(B)
    assert N <= m < 2*N and len(D) <= 12*m < 24*N
    changed = {h for h in D if h != (0,0) and f[h] == -1}
    modified = {h:(1 if h in changed else v) for h,v in f.items()}
    assert all(modified[neg(h)] == v for h,v in modified.items())
    assert sum(abs(modified[h]-f[h]) for h in f) == 2*len(changed)
    kernel = lambda a,b:modified[add(a,b)]*modified[sub(a,b)]
    paired = sorted(B)+sorted(negB)
    quadratic = sum(kernel(a,b) for a in paired for b in paired)
    assert quadratic == 4*m*(m-1)
    # All exceptional couplings and the exact singleton exceptional block.
    for a in special:
        for b in f:
            assert kernel(a,b) == (0 if a == b else 1)
    # All four translation symmetries on a deterministic grid of pairs.
    sample = sorted(f)[::max(1,L//31)]
    symmetry_checks = 0
    for t in special:
        for a in sample:
            for b in sample:
                assert kernel(add(a,t),add(b,t)) == kernel(a,b)
                symmetry_checks += 1
    # A rational lower bound on the centered normalized Rayleigh quotient.
    floor = isqrt(p)
    lower = Fraction(4*((m-1)*floor-m),3*N*floor)
    assert lower > 1
    assert 2*len(changed) <= 48*N and 48*N <= 49*floor
    # L1 telescoping certifies every frequency twist at once.
    patterns = [(sample[i],sample[(3*i+1) % len(sample)],sample[(5*i+2) % len(sample)])
                for i in range(min(16,len(sample)))]
    patterns += [(h,h,sample[0],sample[0]) for h in sample[:8]]
    perturbations = []
    for shifts in patterns:
        count = Counter(shifts)
        odd = [t for t,c in count.items() if c % 2]
        if not odd:
            before = [prod(f[add(h,t)] for t in shifts) for h in f]
            after = [prod(modified[add(h,t)] for t in shifts) for h in f]
            assert before == after
            perturbations.append(dict(odd=0,l1_difference=0,bound=0))
            continue
        before = [prod(f[add(h,t)] for t in odd) for h in f]
        after = [prod(modified[add(h,t)] for t in odd) for h in f]
        difference = sum(abs(a-b) for a,b in zip(before,after))
        bound = 2*len(odd)*len(changed)
        assert difference <= bound
        perturbations.append(dict(odd=len(odd),l1_difference=difference,bound=bound))
    return dict(p=p,r=r,elliptic_group_order=len(E.points),quotient_order=L,
                invariant_factors=[d,n],box_dimensions=[k,ell],box_size=m,
                modification_region_size=len(D),changed_signs=len(changed),
                paired_test_vector_support=2*m,kernel_quadratic_form=quadratic,
                centered_normalized_rayleigh='4*(m-1-m/sqrt(p))/sqrt(7*p)',
                rational_rayleigh_lower_bound=str(lower),
                rational_margin_numerator=lower.numerator-lower.denominator,
                exceptional_entries_checked=4*L,symmetry_checks=symmetry_checks,
                perturbation_certificates=perturbations,
                limitation='Modified even sign kernel with preserved square-root orders and larger constants; not a Paley graph, no Kummer realization or ambient projection identity claimed.')


def main():
    small = [small_case(p,r) for p,r in [(13,4),(29,5),(41,2),(97,2),(257,62)]]
    print(json.dumps(dict(small_cases=len(small))),flush=True)
    planted = []
    for p,r in [(10009,2),(65537,2)]:
        out = planted_case(p,r)
        planted.append(out)
        print(json.dumps({k:out[k] for k in ['p','invariant_factors','box_dimensions','box_size','changed_signs','rational_rayleigh_lower_bound']}),flush=True)
    inputs = ['research/parallel15-elliptic-correlations-2026-09-05.md',
              'experiments/parallel15_elliptic_correlations_2026_09_05.py',
              'research/parallel14-elliptic-model-2026-09-05.md',
              'experiments/parallel14_elliptic_model_2026_09_05.py',
              'experiments/parallel3_conic_2026_09_04.py',
              'experiments/parallel10_word_aggregate_2026_09_05.py',
              'sources/katz-gauss-kloosterman-monodromy.pdf',
              'sources/bekker-zarhin-1702.02255v2.html']
    totals = dict(small_curves=len(small),parity_patterns=sum(c['parity_patterns'] for c in small),
                  twisted_supports=sum(c['twisted_supports'] for c in small),
                  twisted_bound_leading_minors=sum(c['twisted_bound_leading_minors'] for c in small),
                  square_matrix_entries=sum(c['square_matrix_entries'] for c in small),
                  cyclotomic_product_entries=sum(c['cyclotomic_product_entries'] for c in small),
                  planted_models=len(planted),exceptional_entries=sum(c['exceptional_entries_checked'] for c in planted),
                  symmetry_checks=sum(c['symmetry_checks'] for c in planted),
                  perturbation_certificates=sum(len(c['perturbation_certificates']) for c in planted))
    report = dict(status='Exact checks passed. Higher twisted correlations and an abstract obstruction; full Paley and spectral-edge goals unproved.',
                  small_cases=small,planted_models=planted,totals=totals,
                  arithmetic='Integer finite fields and group coordinates, exact Bareiss positivity, cyclotomic polynomial identities and rational Rayleigh lower bounds; no floating-point acceptance.',
                  input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    dest = ROOT/'results/parallel15_elliptic_correlations_2026_09_05.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),totals=totals)),flush=True)


if __name__ == '__main__':
    main()
