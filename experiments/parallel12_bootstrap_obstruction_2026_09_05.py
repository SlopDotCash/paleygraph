#!/usr/bin/env python3
"""Exact projection countermodels and a finite Paley Rayleigh certificate.

These checks do not prove or disprove an asymptotic Paley conjecture.
All acceptance decisions use integers, Fractions, and Q(sqrt(p)); no floats.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import isqrt
from pathlib import Path
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def qadd(x, y): return x[0] + y[0], x[1] + y[1]
def qscale(x, c): return x[0] * c, x[1] * c
def qsub(x, y): return qadd(x, qscale(y, -1))
def qmul(x, y, p): return x[0]*y[0] + p*x[1]*y[1], x[0]*y[1] + x[1]*y[0]


def qinv(x, p):
    den = x[0]**2 - p*x[1]**2
    assert den != 0
    return F(x[0], den), F(-x[1], den)


def signq(x, p):
    a, b = x
    sa, sb = int(a > 0) - int(a < 0), int(b > 0) - int(b < 0)
    if not sa: return sb
    if not sb or sa == sb: return sa
    difference = a*a - p*b*b
    assert difference != 0  # p is nonsquare and a,b are nonzero rationals.
    return sa if difference > 0 else sb


def madd(x, y): return x[0]+y[0], x[1]+y[1]
def msub(x, y): return x[0]-y[0], x[1]-y[1]
def mscale(x, c, p): return x[0]*c[0]+p*x[1]*c[1], x[0]*c[1]+x[1]*c[0]
def mmul(x, y, p): return x[0]@y[0]+p*(x[1]@y[1]), x[0]@y[1]+x[1]@y[0]
def transpose(x): return x[0].T, x[1].T
def trace(x): return np.trace(x[0]), np.trace(x[1])
def norm2(x, p): return np.sum(x[0]*x[0])+p*np.sum(x[1]*x[1]), 2*np.sum(x[0]*x[1])
def entry(x, i, j): return x[0][i, j], x[1][i, j]
def serial(x): return [str(q) for q in x]


def equal(x, y):
    return np.array_equal(x[0], y[0]) and np.array_equal(x[1], y[1])


def absolute_bound(x, bound, p):
    assert signq(qsub(bound, x), p) >= 0
    assert signq(qadd(bound, x), p) >= 0


def field(p):
    assert p > 2 and p % 4 == 1 and all(p % d for d in range(2, isqrt(p)+1))
    chi = [0] + [1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1, p)]
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    eye, jmat, zero = np.eye(p, dtype=object), np.ones((p, p), dtype=object), np.zeros((p, p), dtype=object)
    assert np.array_equal(s@s, p*eye-jmat)
    assert np.array_equal(s@jmat, zero)
    return chi, s, eye, jmat, zero


def plant_case(p, anchors, endpoint, depth, word_depth=0):
    chi, s, eye, jmat, zero = field(p)
    assert endpoint in (0, 1)
    assert all(chi[(x-y) % p] == 1 for x, y in combinations(anchors, 2))
    cell = [x for x in range(p) if all(chi[(x-z) % p] == 1 for z in anchors)]
    assert len(cell) >= 2
    x, y = cell[:2]
    v0 = np.zeros((p, 1), dtype=object)
    v0[x, 0], v0[y, 0] = 1, -1
    v = v0, np.zeros_like(v0)
    proj = (eye-jmat*F(1, p))*F(1, 2), s*F(1, 2*p)
    assert equal(mmul(proj, proj, p), proj)
    w = mmul(proj, v, p)
    q = entry(mmul(transpose(v), w, p), 0, 0)
    s2 = 2
    assert q == (F(1), F(-chi[(x-y) % p], p))
    assert signq(q, p) > 0 and signq(qsub((s2, 0), q), p) > 0
    old_direction = mscale(mmul(w, transpose(w), p), qinv(q, p), p)
    if endpoint:
        new_direction = mscale(mmul(v, transpose(v), p), (F(1, s2), 0), p)
    else:
        u = msub(mscale(w, (s2, 0), p), mscale(v, q, p))
        unorm = qscale(qmul(q, qsub((s2, 0), q), p), s2)
        assert entry(mmul(transpose(u), u, p), 0, 0) == unorm
        new_direction = mscale(mmul(u, transpose(u), p), qinv(unorm, p), p)
    modified = madd(msub(proj, old_direction), new_direction)
    assert equal(modified, transpose(modified))
    assert equal(mmul(modified, modified, p), modified)
    assert trace(modified) == trace(proj) == (F(p-1, 2), 0)
    ones = np.ones((p, 1), dtype=object), np.zeros((p, 1), dtype=object)
    assert equal(mmul(modified, ones, p), (np.zeros((p, 1), dtype=object), np.zeros((p, 1), dtype=object)))
    assert equal(mmul(modified, v, p), mscale(v, (endpoint, 0), p))
    assert all(entry(modified, i, z) == entry(proj, i, z) for i in range(p) for z in anchors)
    r = msub(mscale(proj, (2, 0), p), (eye, zero))
    rnew = msub(mscale(modified, (2, 0), p), (eye, zero))
    snew = mscale(madd(rnew, (jmat*F(1, p), zero)), (0, 1), p)
    assert equal(snew, transpose(snew))
    assert equal(mmul(snew, snew, p), (p*eye-jmat, zero))
    assert equal(mmul(snew, ones, p), (np.zeros((p, 1), dtype=object), np.zeros((p, 1), dtype=object)))
    assert all(entry(snew, i, z) == (s[i, z], 0) for i in range(p) for z in anchors)
    assert equal(mmul(rnew, rnew, p), (eye, zero))
    alphabet_difference = qsub(entry(snew, x, x), entry(snew, x, y))
    assert alphabet_difference == (0, 2*endpoint-1)
    # An irrational difference cannot arise from a zero diagonal and +/-1 entry.
    assert not (entry(snew, x, x) == (0, 0) and entry(snew, x, y) in ((1, 0), (-1, 0)))
    b, n = 2**len(anchors), 2**len(anchors)-1
    d = np.array([b*int(i in cell)-1 for i in range(p)], dtype=object)
    masks = [np.array([int(np.prod([chi[(i-z) % p] for k,z in enumerate(anchors) if mask >> k & 1])) for i in range(p)], dtype=object)
             for mask in range(1, b)]
    soft = sum(masks, np.zeros(p, dtype=object))
    for z in anchors: soft[z] -= b//2
    assert np.array_equal(soft, d)
    # U=sqrt(N) T avoids introducing a second quadratic extension.
    uold, unew = (r[0]*d, r[1]*d), (rnew[0]*d, rnew[1]*d)
    assert equal(mmul(unew, v, p), mscale(v, ((2*endpoint-1)*n, 0), p))
    oldpowers, newpowers = [(eye, zero)], [(eye, zero)]
    for k in range(1, 2*depth+1):
        oldpowers.append(mmul(oldpowers[-1], uold, p))
        newpowers.append(mmul(newpowers[-1], unew, p))
        # |tr(T'^k-T^k)| <= 4 k N^(k/2).
        absolute_bound(qsub(trace(newpowers[-1]), trace(oldpowers[-1])), (4*k*n**k, 0), p)
    energies = [qscale(norm2(z, p), F(1, n**j)) for j,z in enumerate(newpowers[:depth+1])]
    checks = []
    for j in range(1, depth+1):
        diffnorm = qscale(norm2(msub(newpowers[j], oldpowers[j]), p), F(1, n**j))
        assert signq(qsub((8*j*j*n**j, 0), diffnorm), p) >= 0
        assert signq(qsub(energies[j], (n**j, 0)), p) >= 0
        assert equal(mmul(newpowers[j], v, p), mscale(v, (((2*endpoint-1)*n)**j, 0), p))
        rows = newpowers[j-1][0][cell, :], newpowers[j-1][1][cell, :]
        row_energy = qscale(norm2(rows, p), F(1, n**(j-1)))
        predicted = qadd(qscale(energies[j-1], F(1, n)), qscale(row_energy, F(n*n-1, n)))
        assert predicted == energies[j]
        checks.append(dict(j=j,energy=serial(energies[j]),power_difference_norm_squared=serial(diffnorm)))
    word_checks = 0
    word_pairs = [(mask[:, None]*s, zero) for mask in masks]
    new_word_pairs = [(mask[:, None]*snew[0], mask[:, None]*snew[1]) for mask in masks]
    for k in range(1, word_depth+1):
        for word in product(range(n), repeat=k):
            original, altered = (eye, zero), (eye, zero)
            for z in word:
                original, altered = mmul(original, word_pairs[z], p), mmul(altered, new_word_pairs[z], p)
            bound = (4*k*p**(k//2), 0) if k % 2 == 0 else (0, 4*k*p**(k//2))
            absolute_bound(qsub(trace(altered), trace(original)), bound, p)
            word_checks += 1
    return dict(p=p,anchors=anchors,cell=cell,endpoint=endpoint,planted_support=[x,y],
                q=serial(q),rank=(p-1)//2,projection_and_character_square_identities=True,
                anchor_column_checks=p*len(anchors),alphabet_defect=dict(diagonal=serial(entry(snew,x,x)),off_diagonal=serial(entry(snew,x,y)),difference=serial(alphabet_difference)),
                depth=depth,power_trace_checks=2*depth,power_norm_checks=depth,energy_recurrences=depth,
                planted_power_checks=depth,word_trace_checks=word_checks,energies=checks)


def finite_paley_outlier():
    p, anchors = 257, [0, 1, 62]
    chi, sfull, _, _, _ = field(p)
    assert all(chi[(x-y) % p] == 1 for x,y in combinations(anchors,2))
    cell = [x for x in range(p) if all(chi[(x-z) % p] == 1 for z in anchors)]
    expected = [2,16,18,26,30,31,32,36,58,60,61,73,114,121,122,123,124,129,134,135,141,185,190,196,197,198,199,208,227,235,240,249]
    v = [-1,1,-3,3,-3,1,3,-1,-1,1,-1,-3,-3,-3,3,2,-2,2,1,1,2,-3,-2,1,-1,3,-2,2,-1,3,-2,3]
    assert cell == expected and len(v) == len(cell)
    ss = sum(x*x for x in v)
    q = sum(v[i]*v[j]*chi[(x-y) % p] for i,x in enumerate(cell) for j,y in enumerate(cell))
    total = sum(v)
    assert (ss,q,total) == (152,-1624,0)
    # Independent matrix-vector calculation of the same integer quadratic form.
    vnp = np.array(v,dtype=object)
    assert int(vnp@sfull[np.ix_(cell,cell)]@vnp) == q
    # Normalized Rayleigh quotient is -812/(19 sqrt(1799)).
    assert F(-4*q,ss) == F(812,19)
    assert p*7 == 1799
    edge_margin = 812**2 - 19**2*1799
    assert edge_margin == 9905 > 0
    growth_margin = 812**2*144**2 - 19**2*1799*145**2
    assert growth_margin == 17702209 > 0
    assert (F(9,8)+F(8,9))/2 == F(145,144)
    return dict(p=p,anchors=anchors,cell=cell,integer_vector=v,squared_norm=ss,
                character_quadratic_form=q,coordinate_sum=total,
                normalized_rayleigh_quotient='-812/(19 sqrt(1799))',edge_square_margin=edge_margin,
                hyperbolic_threshold='145/144 = (9/8 + 8/9)/2',growth_square_margin=growth_margin,
                consequence='The actual finite T has spectral radius >9/8, hence H_j>(81/64)^j for every j>=1.',
                scope='One finite prime. This does not refute asymptotic spectral-edge convergence or the Paley conjecture.')


def main():
    cases=[]
    for p,anchors,depth,word_depth in [(13,[0,1],6,3),(17,[0,1],3,0),(41,[0,1,2],2,0)]:
        for endpoint in (0,1):
            cases.append(plant_case(p,anchors,endpoint,depth,word_depth))
            print(json.dumps(dict(p=p,endpoint=endpoint,exact_checks='passed')),flush=True)
    inputs=['research/parallel12-bootstrap-obstruction-2026-09-05.md',
            'experiments/parallel12_bootstrap_obstruction_2026_09_05.py',
            'research/parallel11-full-energy-2026-09-05.md',
            'research/parallel2-spectral-transfer-2026-09-04.md']
    keys=('anchor_column_checks','power_trace_checks','power_norm_checks','energy_recurrences','planted_power_checks','word_trace_checks')
    report=dict(status='Exact finite checks passed; an elementary obstruction to a projection-only bootstrap, not a new asymptotic Paley upper bound.',
                arithmetic='Integer, rational, and exact quadratic-field arithmetic. Numerical discovery is not an acceptance condition.',
                projection_countermodels=cases,finite_paley_certificate=finite_paley_outlier(),
                totals={k:sum(c[k] for c in cases) for k in keys},
                input_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs})
    report['totals']['projection_countermodels']=len(cases)
    dest=ROOT/'results/parallel12_bootstrap_obstruction_2026_09_05.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),totals=report['totals'],finite_paley_certificate='passed')),flush=True)


if __name__ == '__main__': main()
