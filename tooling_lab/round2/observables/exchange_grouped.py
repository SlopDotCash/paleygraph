#!/usr/bin/env python3
"""Grouped character-row backend for exact all-input Johnson L2 spectra.

No input-set or field-element enumeration. Fixed degree truncates every
coefficient computation. Scope: p prime, p == 1 mod 4, n between d and p/2.
These are finite all-input second moments, never worst-case estimates.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb,isqrt
from pathlib import Path
import json
import time

from review_exchange import C,solve,eberlein,kernel_overlap_sums,exact_spectrum,encode

HERE=Path(__file__).resolve().parent


def power_terms(a,b,count,limits):
    """Coefficients of (1+abX+aY+bZ)^count inside the truncation box."""
    rx,ry,rz=limits
    out={}
    for i in range(rx+1):
        for j in range(ry+1):
            for k in range(rz+1):
                if i+j+k>count:continue
                sign=(a*b)**i*a**j*b**k
                if sign:
                    out[(i,j,k)]=sign*C(count,i)*C(count-i,j)*C(count-i-j,k)
    return out


def grouped_coefficient(groups,r,d):
    limits=(r,d-r,d-r)
    values={(0,0,0):1}
    for a,b,count in groups:
        terms=power_terms(a,b,count,limits)
        out={}
        for (i,j,k),v in values.items():
            for (x,y,z),w in terms.items():
                key=(i+x,j+y,k+z)
                if all(s<=t for s,t in zip(key,limits)):
                    out[key]=out.get(key,0)+v*w
        values={k:v for k,v in out.items() if v}
    return values.get(limits,0)


def grouped_kernel_sums(p,d=6):
    assert p>=5 and p%4==1 and all(p%q for q in range(2,isqrt(p)+1))
    t=(p-1)//4
    diagonal=[(0,0,1),(1,1,2*t),(-1,-1,2*t)]
    distinct=[(0,1,1),(1,0,1),(1,1,t-1),(1,-1,t),(-1,1,t),(-1,-1,t)]
    K=[p*grouped_coefficient(diagonal,r,d)+p*(p-1)*grouped_coefficient(distinct,r,d)
       for r in range(d+1)]
    return K,{'diagonal':diagonal,'distinct':distinct}


def falling(n,k):
    if k<0 or k>n:return 0
    out=1
    for i in range(k):out*=n-i
    return out


def containment_falling(p,n,j,r,d):
    """Fixed-size marked-point allocation, no factorials or large subsets."""
    u=2*d-r
    numerator=0
    for a in range(d-r+1):
        for b in range(d-r+1):
            numerator+=C(d-r,a)*C(d-r,b)*falling(n-j,u-a-b)*falling(j,a)*falling(j,b)
    return Q(numerator,falling(p,u))


def grouped_spectrum(p,n,d=6):
    assert d<=n<=p//2
    K,groups=grouped_kernel_sums(p,d)
    corr=[sum(Q(K[r])*containment_falling(p,n,j,r,d) for r in range(d+1))
          for j in range(d+1)]
    theta=[[eberlein(p,n,z,j) for z in range(d+1)] for j in range(d+1)]
    energies=solve(theta,corr)
    mean=Q(p*(-1)**(d//2)*C((p-1)//2,d//2)*C(n,d),C(p,d)) if d%2==0 else Q(0)
    variance=corr[0]-mean*mean
    assert energies[0]==mean*mean and all(e>=0 for e in energies)
    assert sum(energies)==corr[0]
    # Even-degree targets are invariant under all affine maps, forcing
    # levels one and two to vanish. Odd-degree targets change sign under
    # nonsquare dilation and can have a nonzero degree-two component.
    if d>=2 and d%2==0:assert energies[1]==energies[2]==0
    return {'p':p,'n':n,'degree':d,'kernel_overlap_sums':K,'row_groups':groups,
            'distance_correlations':corr,'energies':energies,'mean':mean,'variance':variance,
            'top_energy_ratio':energies[d]/variance if variance else None,
            'lower_energy_ratio':sum(energies[1:d])/variance if variance else None}


def fourth_root_floor(p):
    n=isqrt(isqrt(p))
    assert n**4<=p<(n+1)**4
    return n


def main():
    start=time.perf_counter()
    matches=[]
    # Replay all previous cases via independent grouped coefficients and the
    # new falling-factorial containment formula; compare every exact fraction.
    for p,n in [(13,6),(17,6),(61,6),(1297,6),(2437,7),(4129,8)]:
        old=exact_spectrum(p,n)
        new=grouped_spectrum(p,n)
        for key in ['kernel_overlap_sums','distance_correlations','energies','mean','variance']:
            assert old[key]==new[key]
        matches.append({'p':p,'n':n,'all_coefficients_correlations_and_energies_equal':True})
    odd_matches=[]
    for p,n,d in [(13,6,3),(13,6,5),(17,6,3)]:
        old=exact_spectrum(p,n,d)
        new=grouped_spectrum(p,n,d)
        for key in ['kernel_overlap_sums','distance_correlations','energies','mean','variance']:
            assert old[key]==new[key]
        assert new['energies'][2]>0
        odd_matches.append({'p':p,'n':n,'degree':d,'level_two_energy':new['energies'][2],
                            'all_coefficients_correlations_and_energies_equal':True})
    cases=[]
    for p,n in [(65537,16),(1000033,31),(6700417,50),(6700417,64)]:
        t=time.perf_counter()
        spec=grouped_spectrum(p,n)
        spec['critical_n_floor_p_quarter']=fourth_root_floor(p)
        spec['on_critical_size']=n==spec['critical_n_floor_p_quarter']
        spec['elapsed_seconds']=round(time.perf_counter()-t,6)
        cases.append(spec)
        print(p,n,'critical',spec['on_critical_size'],'high',float(spec['top_energy_ratio']),
              'lower',float(spec['lower_energy_ratio']),'seconds',spec['elapsed_seconds'],flush=True)
    out={'status':'exact finite all-input L2 spectra, not a worst-case certificate',
         'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'review_backend_sha256':sha256((HERE/'review_exchange.py').read_bytes()).hexdigest(),
         'prior_backend_matches':matches,'odd_degree_backend_matches':odd_matches,'cases':cases,
         'complexity_note':'For fixed degree, polynomial coefficient work is independent of p and n apart from integer arithmetic; prime validation uses trial division.',
         'elapsed_seconds':round(time.perf_counter()-start,3)}
    (HERE/'exchange_grouped_results.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    print('Saved grouped spectra in',out['elapsed_seconds'],'seconds')


if __name__=='__main__':main()
