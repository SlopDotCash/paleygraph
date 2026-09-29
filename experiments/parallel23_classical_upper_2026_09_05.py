#!/usr/bin/env python3
"""Exact finite audit of an almost-all Sidon sixth-moment upper bound.

This does not check or claim a uniform bound on exceptional sets.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import comb, isqrt
from pathlib import Path
import json
import random
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()

def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

def chitable(p):
    return [0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]

def rademacher(m, order=5):
    U=[1]+[0]*order
    for _ in range(m):
        U=[sum(comb(2*s,2*a)*U[a] for a in range(s+1)) for s in range(order+1)]
    return U

def Fvec(C,ch):
    p=len(ch)
    return [sum(ch[(x-c)%p] for c in C) for x in range(p)]

def moments(F,order=5):
    return [len(F)]+[sum(x**(2*s) for x in F) for s in range(1,order+1)]

def transform(f,ch):
    p=len(ch)
    return [sum(ch[(x-b)%p]*f[x] for x in range(p)) for b in range(p)]

def var(vals):
    mu=Q(sum(vals),len(vals))
    return Q(sum(v*v for v in vals),len(vals))-mu*mu

def variance_certificate(p,n):
    U=rademacher(n-1)
    A=36*U[5]+240*U[4]+472*U[3]+240*U[2]+36*U[1]
    B=U[4]+2*U[3]+U[2]
    rad=p*A*B
    z=isqrt(rad)
    ceil=z+(z*z<rad)
    bound=n*(p*A+225*B+30*ceil)
    gap=p*(30*n*n-16*n)
    return {'A':A,'B':B,'radicand':rad,'ceil_sqrt':ceil,
            'variance_upper_integer':bound,'mean_gap':gap,
            'all_set_bad_fraction_upper':str(min(Q(1),Q(bound,gap*gap))),
            'sidon_certificate_applicable':n>=6 and n**4<=p,
            'sidon_bad_fraction_upper':(str(min(Q(1),Q(8*bound,7*gap*gap)))
                                         if n>=6 and n**4<=p else None)}

def is_sidon(C,p):
    sums=[(a+b)%p for a,b in combinations_with_replacement(C,2)]
    return len(set(sums))==len(sums)

def small_exhaustive():
    reports=[]
    for p in [5,7,11,13]:
        assert prime(p)
        ch=chitable(p)
        # Independent total numbers of two-sum relations of both support sizes.
        ap=sum((2*a-b-c)%p==0 for a in range(p)
               for b,c in combinations([x for x in range(p) if x!=a],2))
        four=0
        for a,b,c,d in combinations(range(p),4):
            four+=sum((u+v-w-z)%p==0 for u,v,w,z in [(a,b,c,d),(a,c,b,d),(a,d,b,c)])
        assert ap==p*(p-1)//2
        assert four==p*(p-1)*(p-3)//8
        COUNTS['relation_count_fields']+=1
        for n in range(1,min(5,p-1)+1):
            sets=list(combinations(range(p),n))
            table={C:moments(Fvec(C,ch)) for C in sets}
            for s in range(1,6):
                avg=Q(sum(table[C][s] for C in sets),len(sets))
                assert avg<=p*rademacher(n)[s]
                COUNTS['balanced_population_moment_domination']+=1
            values=[table[C][3] for C in sets]
            actual_var=var(values)
            cond=[]
            # The variance inequality is also checked on unrelated slice functions.
            alt={C:sum((i+1)*(c+1)**2 for i,c in enumerate(C)) for C in sets}
            alt_cond=[]
            for AA in combinations(range(p),n-1):
                Aset=set(AA)
                candidates=[b for b in range(p) if b not in Aset]
                F=Fvec(AA,ch)
                M=moments(F)
                H=[6*y**5+20*y**3+6*y for y in F]
                V=[y**4+y**2 for y in F]
                QH=transform(H,ch)
                g=[QH[b]-15*V[b] for b in range(p)]
                assert sum(y*y for y in QH)==p*sum(y*y for y in H)-sum(H)**2
                COUNTS['full_character_transform_identities']+=1
                const=M[3]+15*M[2]+15*M[1]+p-1
                mvals=[];avals=[]
                for b in candidates:
                    C=tuple(sorted(AA+(b,)))
                    assert table[C][3]==const+g[b]
                    mvals.append(table[C][3]);avals.append(alt[C])
                    COUNTS['exact_insertion_identities']+=1
                cond.append(var(mvals));alt_cond.append(var(avals))
            factor=Q(n*(p-n+1),p)
            assert actual_var<=factor*sum(cond)/len(cond)
            assert var(list(alt.values()))<=factor*sum(alt_cond)/len(alt_cond)
            COUNTS['slice_poincare_checks']+=2
            cert=variance_certificate(p,n)
            assert actual_var<=cert['variance_upper_integer']
            # Check the unsquared-radical variance certificate exactly as well.
            lhs=actual_var/Q(n)-p*cert['A']-225*cert['B']
            assert lhs<=0 or lhs*lhs<=900*cert['radicand']
            COUNTS['sixth_variance_bound_checks']+=1
            bad=sum(table[C][3]>15*p*n**3 for C in sets)
            bad_sidon=sum(table[C][3]>15*p*n**3 and is_sidon(C,p) for C in sets)
            sidon=sum(is_sidon(C,p) for C in sets)
            collision=Q(n*(n-1)*(n-2)*(n+1),8*(p-2)) if n>=3 else Q(0)
            assert Q(len(sets)-sidon,len(sets))<=collision
            assert Q(bad,len(sets))<=Q(cert['variance_upper_integer'],cert['mean_gap']**2)
            COUNTS['exhaustive_slice_sets']+=len(sets)
            reports.append({'p':p,'n':n,'sets':len(sets),'sidon_sets':sidon,
                            'bad_moment_sets':bad,'bad_sidon_sets':bad_sidon,
                            'actual_variance':str(actual_var),'variance_certificate':cert})
    return reports

def fast_quartic_energy(C,p,ch):
    n=len(C)
    x=np.arange(p,dtype=np.int64)
    ch=np.asarray(ch,dtype=np.int64)
    F=np.zeros(p,dtype=np.int64)
    for c in C:F+=ch[(x-c)%p]
    N=np.full(p,n,dtype=np.int64);N[list(C)]-=1
    num=F**3-(3*N-2)*F
    assert np.all(num%6==0)
    e3=num//6
    assert p*max(abs(int(v)) for v in e3)<2**63
    U=np.empty(p,dtype=np.int64)
    for start in range(0,p,96):
        b=np.arange(start,min(start+96,p),dtype=np.int64)
        U[start:start+len(b)]=ch[(x[None,:]-b[:,None])%p]@e3
    E3=sum(int(v)**2 for v in e3)
    L=sum(int(v)**2 for v in U)
    assert L==p*E3-int(e3.sum())**2
    Lin=sum(int(U[c])**2 for c in C)
    M6=sum(int(v)**6 for v in F)
    M2=sum(int(v)**2 for v in F)
    assert M2==p*n-n*n
    assert 36*E3<=M6+9*n*n*M2
    if M6<=15*p*n**3:
        assert 3*(L-Lin)<=2*p*p*n**3
    COUNTS['actual_quartic_energy_checks']+=1
    return {'C':list(C),'M6':M6,'M6_over_pn3':str(Q(M6,p*n**3)),
            'outside_quartic_energy':L-Lin,'outside_energy_over_p2n3':str(Q(L-Lin,p*p*n**3))}

def sampled_slices():
    rng=random.Random(230905)
    reports=[]
    for p in [1297,2411,4099,6563,10007]:
        assert prime(p)
        n=isqrt(isqrt(p))
        assert n>=6 and n**4<=p<(n+1)**4
        collision=Q(n*(n-1)*(n-2)*(n+1),8*(p-2))
        assert collision<=Q(1,8)
        ch=chitable(p)
        chosen=[];trials=0;moment_bad=0;sidon_total=0
        for _ in range(64):
            C=tuple(sorted(rng.sample(range(p),n)))
            trials+=1
            if not is_sidon(C,p):continue
            sidon_total+=1
            F=Fvec(C,ch)
            M6=sum(x**6 for x in F)
            if M6>15*p*n**3:moment_bad+=1
            if len(chosen)<1:chosen.append(C)
            COUNTS['sampled_actual_sidon_sets']+=1
        energies=[fast_quartic_energy(C,p,ch) for C in chosen]
        reports.append({'p':p,'n':n,'sampling_trials':trials,'sampled_sidon_sets':sidon_total,
                        'sampled_bad_sidon_sets':moment_bad,'collision_union_bound':str(collision),
                        'probability_certificate':variance_certificate(p,n),
                        'quartic_energy_fixtures':energies})
    return reports

def large_probability_certificates():
    out=[]
    for n in [32,64,128,256,512]:
        p=n**4+1
        while not prime(p):p+=2
        assert isqrt(isqrt(p))==n
        cert=variance_certificate(p,n)
        bound=Q(cert['sidon_bad_fraction_upper'])
        out.append({'p':p,'n':n,'sidon_bad_fraction_upper':str(bound),
                    'scaled_coefficient_upper':str(bound*p/(n*n)),
                    'scope':'analytic probability certificate only; no field character evaluation'})
        COUNTS['large_prime_probability_certificates']+=1
    return out

def main():
    small=small_exhaustive()
    print('Small exhaustive insertion, transform, and variance checks passed.',flush=True)
    samples=sampled_slices()
    large=large_probability_certificates()
    inputs=['research/parallel21-positive-upper-review-2026-09-05.md',
            'research/parallel22-quartic-star-independent-review-2026-09-05.md',
            'research/parallel20-inversion-moments-2026-09-05.md',
            'research/parallel6-classical-2026-09-04.md',
            'experiments/parallel23_classical_upper_2026_09_05.py']
    result={'status':'exact finite checks passed; an almost-all Sidon bound is proved, while uniform exceptional-set control remains open',
            'arithmetic':'Python integers and Fraction; exact radical comparisons; int64 only for bounded character transforms with checked bounds',
            'counts':dict(sorted(COUNTS.items())),
            'small_exhaustive_slices':small,'actual_prime_slice_samples':samples,
            'large_analytic_probability_certificates':large,
            'input_sha256':{f:sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs}}
    target=ROOT/'results/parallel23_classical_upper_2026_09_05.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'counts':result['counts'],'result':str(target)},indent=2))

if __name__=='__main__':main()
