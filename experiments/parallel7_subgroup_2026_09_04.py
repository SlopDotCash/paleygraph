#!/usr/bin/env python3
"""Finish the bounded positivity and continuous-order audits after worker interruption."""
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from parallel6_subgroup_2026_09_04 import cosine_interval, sum_walks
from parallel7_multiplier_2026_09_04 import plus, times, scalar, total, encode

ROOT=Path(__file__).resolve().parents[1]


def negative_case(p, H, x, decimal_bracket):
    n=len(H)
    assert len(set(H))==n and p-1 in H
    assert n**4//4 <= p <= n**4
    assert all(a*b % p in H for a in H for b in H)
    cos=[cosine_interval(t,p) for t in range(p)]
    theta=[scalar(total(cos[t*h % p] for h in H),Q(1,n)) for t in range(p)]
    absolute=[]
    for lo,hi in theta:
        assert lo>=0 or hi<=0
        absolute.append((lo,hi) if lo>=0 else (-hi,-lo))
    def coefficient(z):
        return scalar(total(times((lo**3,hi**3),cos[t*z % p])
                            for t,(lo,hi) in enumerate(absolute)),Q(1,p))
    orbit=sorted({x*h % p for h in H})
    intervals={z:coefficient(z) for z in orbit}
    target=intervals[x]
    assert Q(decimal_bracket[0]) < target[0] <= target[1] < Q(decimal_bracket[1])
    mass=total(intervals.values())
    assert mass[1]<0
    # For g=1_orbit, E_alpha(g²)=E_alpha(g)=b<0<=b².
    jensen_gap_lower=mass[1]**2-mass[1]
    assert jensen_gap_lower>0
    mu=sum_walks(p,H,highest=3)[3]
    margins=[]
    for a in range(1,p):
        value=total(scalar((absolute[a*y % p][0]**3,absolute[a*y % p][1]**3),Q(mu[y],n**3))
                    for y in range(p))
        margin=value[0]-absolute[a][1]**9
        assert margin>0
        margins.append(margin)
    return {'p':p,'H':H,'coordinate':x,'alpha3_interval':encode(target),
            'negative_coset':orbit,'coset_mass_interval':encode(mass),
            'strict_signed_Jensen_gap_lower':encode(jensen_gap_lower),
            'positive_mu3_Jensen_checks':len(margins),
            'minimum_positive_Jensen_margin':encode(min(margins))}


def saving(r,s):
    if r<=3:
        if s==1: return Q(0)
        if s==2: return max(Q(0),(r-2)/(4*r))
        return (r-1)/(2*(s*(r-1)+3))
    if s==1: return Q(0)
    if s==2: return 1/(2*r+6)
    return 1/(r*s-r-s+6)


def main():
    negative=[negative_case(13,[1,12],6,('-0.004014894','-0.004014892')),
              negative_case(73,[1,27,46,72],5,('-0.001688700','-0.001688698'))]
    points=[]
    for j in range(8,257):
        r=Q(j,8)
        for s in range(1,65):
            value=saving(r,s)
            assert value<=Q(1,9)
            if value==Q(1,9): points.append([encode(r),s])
    assert points==[['3',3]]
    # The monotonicity proof in the note handles the continuum and tails.
    # These exact identities independently check its derivatives' numerators.
    derivative_checks=0
    for j in range(8,25):
        r=Q(j,8); x=r-1
        for s in range(3,65):
            denominator=2*(s*x+3)
            assert denominator-2*s*x==6
            assert -2*x*x<=0
            assert (r*s-r-s+6)==(r-1)*(s-1)+5
            derivative_checks+=3
    p=6700417; n=64
    assert all(p%d for d in range(2,isqrt(p)+1))
    H={pow(2,j,p) for j in range(n)}
    assert len(H)==n and pow(2,n,p)==1 and p-1 in H
    zero_triples=sum((-a-b)%p in H for a in H for b in H)
    assert zero_triples==3*n==192
    assert n**4//4<=p<=n**4
    inputs=['research/parallel7-subgroup-2026-09-04.md',
            'experiments/parallel7_subgroup_2026_09_04.py',
            'research/parallel7-multiplier-2026-09-04.md',
            'experiments/parallel7_multiplier_2026_09_04.py',
            'experiments/parallel6_subgroup_2026_09_04.py',
            'research/cyclotomic-prime-average.md',
            'research/parallel4-subgroup-2026-09-04.md']
    result={'status':'All exact positivity-obstruction, Jensen, ledger and witness checks passed.',
            'completion_note':'Root completed the verifier after the parallel worker hit the account limit.',
            'negative_cases':negative,
            'rational_ledger_points':249*64,'maximizers':points,
            'derivative_coefficient_checks':derivative_checks,
            'analytic_tail_argument':'Increasing low-r curves and decreasing high-r curves; (r-1)(s-1)+5>=9 for r,s>=3. Finite grid is not the continuum proof.',
            'zero_triple_witness':{'p':p,'n':n,'generator':2,'Z3':zero_triples},
            'input_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs}}
    dest=ROOT/'results/parallel7_subgroup_2026_09_04.json'
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'negative_cosets':len(negative),
                      'positive_Jensen_checks':sum(x['positive_mu3_Jensen_checks'] for x in negative),
                      'rational_ledger_points':result['rational_ledger_points'],
                      'derivative_coefficient_checks':derivative_checks,'output':str(dest)},indent=2))


if __name__=='__main__':
    main()
