#!/usr/bin/env python3
"""Exact checks of the three-anchor elliptic model; no asymptotic acceptance.

Run from any directory. Integer/cyclotomic arithmetic only. The cohomological
Fourier bound and the general identities are proved in the companion note;
these finite checks are additional evidence, not proofs of their quantifiers.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path
import json
import numpy as np

from parallel3_conic_2026_09_04 import reduce_poly, cyclic_product
from parallel10_word_aggregate_2026_09_05 import positive_definite_integer

ROOT = Path(__file__).resolve().parents[1]


class Curve:
    def __init__(self, p, r):
        assert p % 4 == 1 and all(p % d for d in range(2, isqrt(p)+1))
        self.p, self.r = p, r
        self.chi = [0]+[1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1, p)]
        assert 1 < r < p and self.chi[r] == self.chi[r-1] == 1
        self.roots = {}
        for y in range(p):
            self.roots.setdefault(y*y % p, []).append(y)
        self.points = [None]+[(x, y) for x in range(p)
                              for y in self.roots.get(x*(x-1)*(x-r) % p, [])]
        self.index = {P:i for i, P in enumerate(self.points)}
        self.cell = [x for x in range(p)
                     if all(self.chi[(x-a) % p] == 1 for a in [0, 1, r])]

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p, r = self.p, self.r
        x, y = P
        u, v = Q
        if x == u and (y+v) % p == 0:
            return None
        if P == Q:
            slope = (3*x*x-2*(1+r)*x+r)*pow(2*y, -1, p) % p
        else:
            slope = (v-y)*pow(u-x, -1, p) % p
        z = (slope*slope+1+r-x-u) % p
        return z, (slope*(x-z)-y) % p

    def neg(self, P):
        return None if P is None else (P[0], -P[1] % self.p)

    def rho(self, P):
        Q = self.add(P, P)
        return self.p if Q is None else Q[0]

    def f(self, P):
        return 0 if P is None else self.chi[P[1]]


def quotient(E, sums, fs, two):
    """G/G[2], with explicit product-cyclic coordinates and all group laws checked."""
    ids, reps = {}, []
    for i in range(len(E.points)):
        if i not in ids:
            c = len(reps)
            reps.append(i)
            for t in two:
                ids[sums[i][t]] = c
    L = len(reps)
    table = [[ids[sums[i][j]] for j in reps] for i in reps]
    values = [fs[i] for i in reps]
    for i, c in ids.items():
        assert fs[i] == values[c]
    cycles = []
    for a in range(L):
        orbit, h = [0], a
        while h:
            orbit.append(h)
            h = table[h][a]
        cycles.append(orbit)
    g = max(range(L), key=lambda a: len(cycles[a]))
    n = len(cycles[g])
    d = L//n
    b = next(b for b in range(L) if len(cycles[b]) == d
             and set(cycles[g]) & set(cycles[b]) == {0})
    assert n % d == 0
    coords = {table[u][v]:(i,j) for i,u in enumerate(cycles[g])
              for j,v in enumerate(cycles[b])}
    assert len(coords) == L
    assert sum(table[h][h] == 0 for h in range(L)) == 4
    for h in range(L):
        for k in range(L):
            assert coords[table[h][k]] == ((coords[h][0]+coords[k][0]) % n,
                                           (coords[h][1]+coords[k][1]) % d)
    neg = [next(k for k in range(L) if table[h][k] == 0) for h in range(L)]
    assert values[0] == 0 and all(abs(v) == 1 for v in values[1:])
    assert all(values[h] == values[neg[h]] for h in range(L))
    A = np.array([[values[table[h][neg[k]]] for k in range(L)]
                  for h in range(L)], dtype=object)
    assert np.array_equal(A, A.T)
    bits = positive_definite_integer(E.p*np.eye(L, dtype=object)-A@A)
    M = [[values[table[h][k]]*values[table[h][neg[k]]]
          for k in range(L)] for h in range(L)]
    return dict(order=L, invariant_factors=[d,n], fourier_bound_minors=L,
                final_minor_bits=bits), (table, values, coords, n, d, M)


def fourier_check(data):
    """Unnormalized DFT numerators, reduced in Z[zeta_n], for every entry."""
    table, values, coords, n, d, M = data
    L = len(table)
    chars = list(product(range(n), range(d)))
    exps = {a:[(a[0]*coords[h][0]+(n//d)*a[1]*coords[h][1]) % n
               for h in range(L)] for a in chars}
    coeffs = {}
    for a in chars:
        poly = [0]*n
        for h in range(L):
            poly[-exps[a][h] % n] += values[h]
        coeffs[a] = poly
        assert reduce_poly(poly, n) == reduce_poly([poly[-i % n] for i in range(n)], n)
    squares = [0]*n
    for a in chars:
        q = cyclic_product(coeffs[a], coeffs[a], n)
        squares = [x+y for x,y in zip(squares,q)]
    assert reduce_poly(squares, n) == [L*(L-1)]
    roots = {}
    for g in chars:
        roots.setdefault((2*g[0] % n, 2*g[1] % d), []).append(g)
    assert all(len(v) == 4 for v in roots.values())
    offdiagonal = None
    forbidden, checked = 0, 0
    for a in chars:
        for b in chars:
            lhs, rhs = [0]*n, [0]*n
            for s in range(L):
                for t in range(L):
                    lhs[(-exps[a][s]+exps[b][t]) % n] += M[s][t]
            gs = roots.get(((a[0]-b[0]) % n, (a[1]-b[1]) % d), [])
            for g in gs:
                gb = ((g[0]+b[0]) % n, (g[1]+b[1]) % d)
                poly = cyclic_product(coeffs[g], coeffs[gb], n)
                rhs = [x+y for x,y in zip(rhs,poly)]
            left, right = reduce_poly(lhs,n), reduce_poly(rhs,n)
            assert left == right, (a,b,left,right)
            checked += 1
            forbidden += not gs
            if a != b and left != [0] and offdiagonal is None:
                offdiagonal = dict(alpha=a,beta=b,numerator_coefficients=left,
                                   denominator=L,root_of_unity_order=n)
    return dict(entry_checks=checked,structural_zero_checks=forbidden,
                offdiagonal_witness=offdiagonal)


def one_case(p, r, do_fourier):
    E = Curve(p,r)
    G, C = E.points, E.cell
    size, m = len(G), len(C)
    sums = [[E.index[E.add(P,Q)] for Q in G] for P in G]
    neg = [E.index[E.neg(P)] for P in G]
    fs = [E.f(P) for P in G]
    rho = [E.rho(P) for P in G]
    two = [i for i in range(size) if sums[i][i] == 0]
    four = [i for i in range(size) if sums[i][i] in two]
    assert len(two) == 4 and len(four) == 16
    assert all(fs[i] == fs[neg[i]] for i in range(size))
    assert all(fs[sums[i][t]] == fs[i] for i in range(size) for t in two)
    fibers = Counter(rho)
    assert set(fibers) == set(C+[0,1,r,p])
    assert all(fibers[x] == (8 if x in C else 4) for x in fibers)
    assert size == 8*m+16
    triples = 0
    for x in C:
        inverse_points = set()
        for a,b,c in product(E.roots[x], E.roots[(x-1) % p], E.roots[(x-r) % p]):
            D = ((1-r)*a+r*b-c) % p
            assert D
            assert D == (a-b)*(a-c)*(b-c) % p
            y = r*(r-1)*pow(D,-1,p) % p
            t = (r+y*(a-b)) % p
            assert t == (a+b)*(a+c) % p and y == (a+b)*(a+c)*(b+c) % p
            assert (t,y) in E.index and E.rho((t,y)) == x
            assert (t*t-r)*pow(2*y,-1,p) % p == a
            assert (t*t-2*t+r)*pow(2*y,-1,p) % p == b
            assert (t*t-2*r*t+r)*pow(2*y,-1,p) % p == c
            inverse_points.add((t,y))
            triples += 1
        assert len(inverse_points) == 8
    M = [[fs[sums[i][j]]*fs[sums[i][neg[j]]] for j in range(size)]
         for i in range(size)]
    for i,x in enumerate(rho):
        for j,y in enumerate(rho):
            expected = 0 if x == y else 1 if p in (x,y) else E.chi[(x-y) % p]
            assert M[i][j] == expected
    # Non-normalized fiber indicators; Q^T M Q = D K D and M Q = Q K D.
    labels = C+[0,1,r,p]
    Q = np.array([[int(x == y) for y in labels] for x in rho],dtype=object)
    D = np.diag([fibers[x] for x in labels]).astype(object)
    K = np.array([[0 if x == y else 1 if p in (x,y) else E.chi[(x-y) % p]
                   for y in labels] for x in labels],dtype=object)
    mat = np.array(M,dtype=object)
    assert np.array_equal(Q.T@Q,D)
    assert np.array_equal(mat@Q,Q@K@D)
    assert np.array_equal(mat,Q@K@Q.T)
    # The complete border is 4(J4-I4) with constant coupling. Its other 3
    # eigenvalues are -4; the normalized constant column couples by 8 sqrt(2).
    assert np.array_equal(K[m:,m:],np.ones((4,4),dtype=object)-np.eye(4,dtype=object))
    assert all(K[i,j] == 1 for i in range(m) for j in range(m,m+4))
    good = [i for i in range(size) if i not in four]
    F = Q[np.ix_(good,range(m))]
    Mg = mat[np.ix_(good,good)]
    S = K[:m,:m]
    assert np.array_equal(F.T@F,8*np.eye(m,dtype=object))
    assert np.array_equal(F.T@Mg@F,64*S)
    assert np.array_equal(F.T@np.ones((len(good),len(good)),dtype=object)@F,
                          64*np.ones((m,m),dtype=object))
    assert np.array_equal(Mg,F@S@F.T)
    # Exact V4 operators, including fixed orbits and all four projectors.
    taus = [lambda x:x,lambda x:r*pow(x,-1,p) % p,
            lambda x:(x-r)*pow(x-1,-1,p) % p,
            lambda x:r*(x-1)*pow(x-r,-1,p) % p]
    index = {x:i for i,x in enumerate(C)}
    Ps = []
    for tau in taus:
        assert {tau(x) for x in C} == set(C)
        P = np.zeros((m,m),dtype=object)
        for x in C:
            P[index[tau(x)],index[x]] = 1
        assert np.array_equal(P@S,S@P)
        Ps.append(P)
    assert np.array_equal(Ps[1]@Ps[2],Ps[3])
    fixed = [int(np.trace(P)) for P in Ps]
    assert all(k <= 2 for k in fixed[1:])
    chars = [(1,a,b,a*b) for a,b in product([1,-1],repeat=2)]
    numerators, dims = [], []
    for signs in chars:
        A = sum(s*P for s,P in zip(signs,Ps))
        dim4 = sum(s*t for s,t in zip(signs,fixed))
        assert dim4 % 4 == 0 and dim4 >= 0
        assert np.array_equal(A@A,4*A)
        assert int(np.trace(A)) == dim4
        if signs != (1,1,1,1):
            assert np.array_equal(A@np.ones((m,m),dtype=object),np.zeros((m,m),dtype=object))
        numerators.append(A)
        dims.append(dim4//4)
    assert np.array_equal(sum(numerators),4*np.eye(m,dtype=object))
    assert all(np.array_equal(A@B,np.zeros((m,m),dtype=object))
               for i,A in enumerate(numerators) for B in numerators[i+1:])
    assert sum(dims) == m
    qreport,qdata = quotient(E,sums,fs,two)
    report = dict(p=p,r=r,group_order=size,cell_size=m,kernel_entries=size*size,
                  inverse_triples=triples,periodicity_checks=4*size,
                  full_border_checked=True,centered_compression_checked=True,
                  klein_fixed_counts=fixed[1:],klein_block_dimensions=dims,
                  quotient=qreport)
    if do_fourier:
        report['fourier'] = fourier_check(qdata)
    return report, (E,S,taus)


def reduced_certificate(data):
    E,S,taus = data
    old = json.loads((ROOT/'results/parallel12_bootstrap_obstruction_2026_09_05.json').read_text())['finite_paley_certificate']
    assert E.p == 257 and E.r == 62 and E.cell == old['cell']
    v = dict(zip(E.cell,old['integer_vector']))
    signs = (1,-1,1,-1)
    assert all(v[tau(x)] == s*v[x] for s,tau in zip(signs,taus) for x in E.cell)
    unused, reps = set(E.cell), []
    while unused:
        x = min(unused)
        orbit = {tau(x) for tau in taus}
        assert len(orbit) == 4
        reps.append(x)
        unused -= orbit
    W = np.zeros((len(E.cell),len(reps)),dtype=object)
    for j,x in enumerate(reps):
        for s,tau in zip(signs,taus):
            W[E.cell.index(tau(x)),j] = s
    B = np.array([[sum(s*E.chi[(x-tau(y)) % E.p] for s,tau in zip(signs,taus))
                   for y in reps] for x in reps],dtype=object)
    w = np.array([v[x] for x in reps],dtype=object)
    assert np.array_equal(W.T@W,4*np.eye(len(reps),dtype=object))
    assert np.array_equal(S@W,W@B)
    assert np.array_equal(W@w,np.array(old['integer_vector'],dtype=object))
    assert int(w@w) == 38 and int(w@B@w) == -406
    assert 812**2-19**2*1799 == 9905 > 0
    assert signs != (1,1,1,1)
    return dict(p=257,r=62,character=signs,orbit_representatives=reps,
                block=B.tolist(),integer_vector=w.tolist(),squared_norm=38,
                quadratic_form=-406,centered_term_zero=True,
                normalized_rayleigh='-812/(19*sqrt(1799))',strict_square_margin=9905,
                meaning='Smaller exact certificate for the previous finite outlier; no asymptotic counterexample.')


def main():
    cases = [(p,r) for p in [13,17,29,41] for r in range(2,p)
             if pow(r,(p-1)//2,p) == pow(r-1,(p-1)//2,p) == 1]
    for p in [73,97]:
        eligible = [r for r in range(2,p) if pow(r,(p-1)//2,p) == pow(r-1,(p-1)//2,p) == 1]
        cases.extend((p,r) for r in [eligible[0],eligible[len(eligible)//2],eligible[-1]])
    cases.append((257,62))
    reports, certificate = [], None
    for p,r in cases:
        row,data = one_case(p,r,(p,r) in [(13,4),(29,5),(41,2),(257,62)])
        reports.append(row)
        if (p,r) == (257,62):
            certificate = reduced_certificate(data)
        print(json.dumps(dict(p=p,r=r,group_order=row['group_order'],cell_size=row['cell_size'])),flush=True)
    assert certificate and reports[-1]['fourier']['offdiagonal_witness']
    inputs = ['research/parallel14-elliptic-model-2026-09-05.md',
              'experiments/parallel14_elliptic_model_2026_09_05.py',
              'experiments/parallel3_conic_2026_09_04.py',
              'experiments/parallel10_word_aggregate_2026_09_05.py',
              'research/parallel12-bootstrap-obstruction-2026-09-05.md',
              'results/parallel12_bootstrap_obstruction_2026_09_05.json',
              'research/parallel13-principal-budget-2026-09-05.md',
              'sources/katz-gauss-kloosterman-monodromy.pdf',
              'sources/bekker-zarhin-1702.02255v2.html']
    totals = dict(curve_cases=len(reports),kernel_entries=sum(r['kernel_entries'] for r in reports),
                  inverse_triples=sum(r['inverse_triples'] for r in reports),
                  periodicity_checks=sum(r['periodicity_checks'] for r in reports),
                  klein_projectors=4*len(reports),
                  fourier_bound_leading_minors=sum(r['quotient']['fourier_bound_minors'] for r in reports),
                  cyclotomic_fourier_entries=sum(r.get('fourier',{}).get('entry_checks',0) for r in reports),
                  fourier_offdiagonal_witnesses=sum(bool(r.get('fourier',{}).get('offdiagonal_witness')) for r in reports))
    report = dict(status='Exact elliptic representation and finite checks passed; no new asymptotic operator or clique bound.',
                  cases=reports,reduced_certificate=certificate,totals=totals,
                  arithmetic='Integer finite fields, object-integer matrices, exact Bareiss minors and cyclotomic polynomial identities. No floating-point acceptance.',
                  source_pdf_pages_visually_reviewed=[21,25,34,35],
                  input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    dest = ROOT/'results/parallel14_elliptic_model_2026_09_05.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),totals=totals)),flush=True)


if __name__ == '__main__':
    main()
