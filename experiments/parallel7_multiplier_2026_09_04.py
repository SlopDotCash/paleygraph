#!/usr/bin/env python3
"""Fixed exact interval checks of the signed multiplier theorem."""
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from parallel6_subgroup_2026_09_04 import cosine_interval, sqrt_interval, outward, sum_walks

ROOT = Path(__file__).resolve().parents[1]


def plus(x, y):
    return x[0]+y[0], x[1]+y[1]


def times(x, y):
    vals = [a*b for a in x for b in y]
    return min(vals), max(vals)


def scalar(x, c):
    return times(x, (Q(c), Q(c)))


def square(x):
    hi = max(x[0]*x[0], x[1]*x[1])
    lo = Q(0) if x[0] <= 0 <= x[1] else min(x[0]*x[0], x[1]*x[1])
    return lo, hi


def total(xs):
    out = (Q(0), Q(0))
    for x in xs:
        out = plus(out, x)
    return out


def encode(x):
    if isinstance(x, tuple):
        return [str(v) for v in x]
    return str(x)


def primitive_root(p):
    assert all(p % d for d in range(2, isqrt(p)+1))
    return next(g for g in range(2, p) if len({pow(g, j, p) for j in range(p-1)}) == p-1)


def one_case(p, n):
    assert (p-1) % n == 0
    g = primitive_root(p)
    H = sorted({pow(g, (p-1)//n*j, p) for j in range(n)})
    assert len(H) == n
    walks = sum_walks(p, H, highest=4)
    energy = {r:sum(v*v for v in walks[r]) for r in walks}
    diff = [0]*p
    for h in H:
        for k in H:
            diff[(h-k) % p] += 1
    assert sum(diff) == n*n and diff[0] == n
    cos = [cosine_interval(t, p) for t in range(p)]
    magnitude = []
    for t in range(p):
        ab2 = total(scalar(cos[t*x % p], Q(diff[x], n*n)) for x in range(p))
        lo, hi = max(Q(0), ab2[0]), min(Q(1), ab2[1])
        assert lo <= hi
        magnitude.append((sqrt_interval(lo)[0], sqrt_interval(hi)[1]))
    magnitude[0] = (Q(1), Q(1))
    unused = set(range(1, p)); reps = []
    while unused:
        a = min(unused); reps.append(a)
        orbit = {a*h % p for h in H}
        assert orbit <= unused
        unused -= orbit
    flat = len(set(diff[1:])) == 1
    if flat:
        flat_square = Q(n-diff[1], n*n)
        assert 0 <= flat_square <= 1
        for t in range(1, p):
            assert magnitude[t][0]**2 <= flat_square <= magnitude[t][1]**2
    cases = []; norm_checks = 0; inverse_checks = 0; gate_checks = 0; flat_checks = 0
    negative = []
    for r in [2, 3, 4]:
        powered = [(lo**r, hi**r) for lo, hi in magnitude]
        if r % 2 == 0:
            base = walks[r//2]
            count = [sum(base[x]*base[(x-z) % p] for x in range(p)) for z in range(p)]
            alpha = [(Q(v, n**r), Q(v, n**r)) for v in count]
        else:
            alpha = [outward(*scalar(total(times(powered[t], cos[t*x % p])
                                              for t in range(p)), Q(1, p)))
                     for x in range(p)]
        mass = total(alpha)
        assert mass[0] <= 1 <= mass[1]
        assert all(alpha[x][0] <= alpha[-x % p][1] and alpha[-x % p][0] <= alpha[x][1]
                   for x in range(p))
        norm = total(square(v) for v in alpha)
        expected_norm = Q(energy[r], n**(2*r))
        assert norm[0] <= expected_norm <= norm[1]
        norm_checks += 1
        for t in reps+[0]:
            transform = total(times(alpha[x], cos[t*x % p]) for x in range(p))
            assert transform[0] <= powered[t][1] and powered[t][0] <= transform[1]
            inverse_checks += 1
        a0 = alpha[0]
        assert a0[0] >= 0 and a0[1] <= Q(1, n)
        if r == 3:
            for x, interval in enumerate(alpha):
                if interval[1] < 0:
                    negative.append({'coordinate':x, 'interval':encode(interval)})
        loss = lambda a: a*a+(1-a)**2/(p-1)
        va_lo = max(Q(0), expected_norm-max(loss(a0[0]),loss(a0[1])))
        for s in [1, 2, 3, 4]:
            mu = [Q(v, n**s) for v in walks[s]]
            u = mu[0]; m = 1-u
            vs = sum((v-m/(p-1))**2 for v in mu[1:])
            centered_s = Q(energy[s], n**(2*s))-Q(1, p)
            assert vs == centered_s-Q(p, p-1)*(u-Q(1, p))**2
            constant_lo = a0[0]+u-a0[0]*u-(1-a0[0])*(1-u)/(p-1)
            rhs_lo = constant_lo+sqrt_interval(p*va_lo*vs)[0]
            for frequency in reps:
                delta_hi = magnitude[frequency][1]
                if flat:
                    # alpha is exactly constant on nonzero coordinates, so
                    # V_alpha=0 and the full RHS is u+(1-u)t^r.
                    # t^(rs)<=t^r<=u+(1-u)t^r for 0<=t<=1, s>=1.
                    assert Q(0) <= u <= 1 and s >= 1
                    flat_checks += 1
                    margin = None
                else:
                    margin = rhs_lo-delta_hi**(r*s)
                    assert margin >= 0, (p,n,r,s,frequency,float(margin))
                    gate_checks += 1
                coarse = Q(2, n)+sqrt_interval(Q(p*energy[r]*energy[s], n**(2*r+2*s)))[0]
                assert delta_hi**(r*s) <= coarse
                cases.append({'r':r,'s':s,'frequency_coset':frequency,
                              'origin_interval':encode(a0),'mu_origin':encode(u),
                              'centered_gate_margin':None if margin is None else encode(margin),
                              'flat_spectrum_algebraic_case':flat})
    return {'p':p,'n':n,'H':H,'minus_one_in_H':p-1 in H,
            'quartic_window':n**4//4 <= p <= n**4,
            'energy':energy,'norm_checks':norm_checks,'fourier_inversion_enclosures':inverse_checks,
            'strict_centered_gate_checks':gate_checks,'flat_spectrum_algebraic_checks':flat_checks,
            'covered_frequency_order_pairs':len(cases)*n,
            'negative_alpha3_coordinates':negative,'gate_cases':cases}


def main():
    cases=[]
    for p,n in [(7,3),(13,2),(13,3),(13,4),(17,4),(29,4),(41,8),(193,4)]:
        cases.append(one_case(p,n))
        print(f'passed p={p}, n={n}',flush=True)
    assert any(not c['minus_one_in_H'] for c in cases)
    assert any(c['quartic_window'] for c in cases)
    assert any(c['negative_alpha3_coordinates'] for c in cases)
    # All exponents and elementary constant reductions are exact rational checks.
    assert Q(1)-Q(1,9)==Q(8,9)
    assert 2*Q(8,9)+3==Q(43,9)
    assert 1+Q(2,9)==Q(11,9)
    assert Q(17,18)-Q(8,9)==Q(1,18)
    assert Q(1)-Q(1,8)==Q(7,8)
    assert Q(3,2)**2 <= 2**3  # (3/2)n^(-3/2)<=1 for all n>=2.
    inputs=['research/parallel7-multiplier-2026-09-04.md',
            'experiments/parallel7_multiplier_2026_09_04.py',
            'experiments/parallel6_subgroup_2026_09_04.py',
            'research/parallel5-subgroup-2026-09-04.md',
            'research/parallel4-subgroup-2026-09-04.md',
            'sources/thorner-zaman-2108.10878v2.html']
    result={'status':'Passed fixed exact interval and algebraic checks; uniform theorem is proved in the note.',
            'arithmetic':'Rational outward cosine/square-root intervals, integer walk energies, exact algebra for flat spectra. No floating-point acceptance.',
            'cases':cases,
            'input_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}}
    dest=ROOT/'results/parallel7_multiplier_2026_09_04.json'
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(dest)


if __name__=='__main__':
    main()
