#!/usr/bin/env python3
"""Exact domain-preserving exchange spectrum of the sixth character U-statistic.

Known machinery: Johnson scheme spectral projectors / slice harmonic analysis.
Specialized implementation: character-row histograms evaluate a full distance
average without enumerating replacement sets. This diagnoses target content;
it must NEVER be used as a purported target-free predictor.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random
import numpy as np
from critical_audit import chi_table

HERE=Path(__file__).resolve().parent


def choose(n,k):
    return comb(n,k) if 0<=k<=n else 0


@lru_cache(None)
def elementary(a,b,degree=6):
    return tuple(sum((-1)**j*choose(b,j)*choose(a,k-j) for j in range(k+1)) for k in range(degree+1))


def distance_averages(c,p,degree=6):
    """E[T_degree(D) | |D minus C|=j], exactly, j=0,...,degree.

    Given r roots of a degree-t product in C and t-r outside C, the inclusion
    probability is (n-j)_r/(n)_r * (j)_(t-r)/(p-n)_(t-r).
    All-field character elementary coefficients are [t^k](1-t²)^((p-1)/2).
    Polynomial division supplies outside coefficients from the inside row.
    """
    n=len(c)
    assert degree<=n<=p//2 and len(set(c))==n
    chi=chi_table(p)
    f=chi[(np.arange(p)[:,None]-np.array(c)[None,:])%p].sum(axis=1)
    zs=np.zeros(p,dtype=int);zs[list(c)]=1
    hist=Counter((int((n-z+s)//2),int((n-z-s)//2)) for s,z in zip(f,zs))
    full=[0 if k%2 else (-1)**(k//2)*comb((p-1)//2,k//2) for k in range(degree+1)]
    mixed=[0]*(degree+1)
    for (a,b),count in hist.items():
        inside=elementary(a,b,degree);outside=[]
        for k in range(degree+1):
            outside.append(full[k]-sum(inside[i]*outside[k-i] for i in range(1,k+1)))
        for r in range(degree+1):mixed[r]+=count*inside[r]*outside[degree-r]
    averages=[]
    for j in range(degree+1):
        averages.append(sum((F(choose(n-j,r),comb(n,r))*F(choose(j,degree-r),comb(p-n,degree-r))*mixed[r]
                             for r in range(degree+1)),F(0)))
    assert averages[0].denominator==1
    return averages


def polynomial_times_linear(poly,a,b):
    out=[F(0)]*(len(poly)+1)
    for i,v in enumerate(poly):out[i]+=a*v;out[i+1]+=b*v
    return out


@lru_cache(None)
def projector_weights(p,n,degree=6):
    """Weights on distance averages for the degree-d Johnson projector."""
    denominator=n*(p-n)
    eigen=[1-F(d*(p-d+1),denominator) for d in range(degree+1)]
    # Distribution of the one-swap walk's distance from its initial set.
    walks=[[F(int(j==0)) for j in range(degree+1)]]
    for t in range(degree):
        nxt=[F(0)]*(degree+1)
        for j,weight in enumerate(walks[-1]):
            if not weight:continue
            up=F((n-j)*(p-n-j),denominator);down=F(j*j,denominator)
            nxt[j]+=weight*(1-up-down)
            if j<degree:nxt[j+1]+=weight*up
            else:assert up*weight==0
            if j:nxt[j-1]+=weight*down
        assert sum(nxt)==1
        walks.append(nxt)
    weights=[]
    for d in range(degree+1):
        poly=[F(1)]
        for e in range(degree+1):
            if e!=d:poly=polynomial_times_linear(poly,-eigen[e]/(eigen[d]-eigen[e]),1/(eigen[d]-eigen[e]))
        weights.append(tuple(sum(poly[t]*walks[t][j] for t in range(degree+1)) for j in range(degree+1)))
    return tuple(weights)


def components(c,p,degree=6):
    av=distance_averages(c,p,degree)
    values=[sum(w*v for w,v in zip(row,av)) for row in projector_weights(p,len(c),degree)]
    assert sum(values)==av[0]
    mean=F(p*(-1)**(degree//2)*comb((p-1)//2,degree//2)*comb(len(c),degree),comb(p,degree)) if degree%2==0 else F(0)
    assert values[0]==mean
    return values,av


def frac(v):return [v.numerator,v.denominator]


def tiny_verify():
    # Complete actual p=13,n=6 universe, independently form all distance sums.
    p,n=13,6;chi=chi_table(p);sets=list(combinations(range(p),n));values={}
    for c in sets:values[c]=int(np.prod(chi[(np.arange(p)[:,None]-np.array(c)[None,:])%p],axis=1).sum())
    for c in [sets[0],sets[531],sets[-1]]:
        counts=[0]*(n+1);sums=[0]*(n+1);cs=set(c)
        for d,v in values.items():j=n-len(cs.intersection(d));counts[j]+=1;sums[j]+=v
        direct=[F(s,k) for s,k in zip(sums,counts)]
        comp,averages=components(c,p)
        assert direct==averages
        # Independent local eigen-equation on each projected component.
        neighbor_sum=[F(0)]*7
        for a in c:
            for b in set(range(p))-cs:
                vs,_=components(tuple(sorted(cs-{a}|{b})),p)
                neighbor_sum=[x+y for x,y in zip(neighbor_sum,vs)]
        for d in range(7):assert neighbor_sum[d]==(n*(p-n)-d*(p-d+1))*comp[d]
    return {'full_universe_targets':len(sets),'exact_distance_profiles':3,'projector_neighbor_equations':21}


def main():
    validation=tiny_verify();cases=[]
    for p,n,count in [(13,6,None),(17,6,None),(61,6,512),(1297,6,512),(2437,7,512),(4129,8,512)]:
        # For affine-invariant statistics, uniform sets containing0,1 give the
        # same expectation as all n-sets. Tiny cases enumerate that family.
        exhaustive=count is None
        if exhaustive:sets=[(0,1)+cs for cs in combinations(range(2,p),n-2)]
        else:
            rng=random.Random(940000+p)
            sets=[(0,1)+tuple(sorted(rng.sample(range(2,p),n-2))) for _ in range(count)]
        vectors=[];examples=[];mean=None
        for i,c in enumerate(sets):
            comp,av=components(c,p);vectors.append(comp);mean=comp[0]
            if i<2:examples.append({'C':c,'components':[frac(x) for x in comp],'distance_averages':[frac(x) for x in av]})
        gram=[[sum(v[i]*v[j] for v in vectors)/len(vectors) for j in range(7)] for i in range(7)]
        if exhaustive:
            assert all(gram[i][j]==0 for i in range(7) for j in range(7) if i!=j)
        variance=sum(sum(v[1:])**2 for v in vectors)/len(vectors)
        case={'p':p,'n':n,'count':len(sets),'exhaustive_normalized':exhaustive,
              'mean':frac(mean),'mean_squared_components':[float(gram[i][i]) for i in range(7)],
              'observed_centered_second_moment':float(variance),
              'highest_component_energy_ratio':float(gram[6][6]/variance),
              'lower_components_energy_ratio':float(sum(gram[i][i] for i in range(1,6))/variance),
              'gram_matrix':[[frac(v) for v in row] for row in gram], 'examples':examples}
        cases.append(case)
        print(p,case['highest_component_energy_ratio'],case['lower_components_energy_ratio'],flush=True)
    result={'status':'exact pointwise projectors; L2 energy exact only for the two exhaustive cohorts; no worst-case bound',
            'validation':validation,'cases':cases,'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'exchange_spectrum_results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
