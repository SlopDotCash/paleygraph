#!/usr/bin/env python3
"""Independent exact audit and all-input L2 spectrum of the character U-statistic.

Does not import the generator or its algebra. Eberlein eigenvalues replace its
polynomial-in-the-walk projectors. Two row-pair generating functions replace
its sampled or normalized input-set cohorts.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from math import comb,factorial,isqrt
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent


def C(n,k):return comb(n,k) if 0<=k<=n else 0


def chi(p):
    assert p>=3 and all(p%d for d in range(2,isqrt(p)+1))
    sq={i*i%p for i in range(1,p)}
    return [0 if i==0 else (1 if i in sq else -1) for i in range(p)]


def eberlein(p,n,d,j):
    return Q(sum((-1)**t*C(d,t)*C(n-d,j-t)*C(p-n-d,j-t) for t in range(j+1)),
             C(n,j)*C(p-n,j))


def solve(matrix,rhs):
    a=[[Q(x) for x in row]+[Q(y)] for row,y in zip(matrix,rhs)]
    size=len(a)
    for c in range(size):
        pivot=next(i for i in range(c,size) if a[i][c])
        a[c],a[pivot]=a[pivot],a[c]
        v=a[c][c];a[c]=[x/v for x in a[c]]
        for i in range(size):
            if i!=c:
                v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[c])]
    return [row[-1] for row in a]


def coefficient(rows,r,d):
    """[X^r Y^(d-r) Z^(d-r)] prod(1+abX+aY+bZ); exact DP."""
    side=d-r
    arr=[[[0]*(side+1) for _ in range(side+1)] for _ in range(r+1)]
    arr[0][0][0]=1
    for a,b in rows:
        for i in range(r,-1,-1):
            for j in range(side,-1,-1):
                for k in range(side,-1,-1):
                    arr[i][j][k]+=(a*b*arr[i-1][j][k] if i else 0)+(a*arr[i][j-1][k] if j else 0)+(b*arr[i][j][k-1] if k else 0)
    return arr[r][side][side]


def kernel_overlap_sums(p,d):
    h=chi(p)
    diag=[(h[-x%p],h[-x%p]) for x in range(p)]
    off=[(h[-x%p],h[(1-x)%p]) for x in range(p)]
    K=[p*coefficient(diag,r,d)+p*(p-1)*coefficient(off,r,d) for r in range(d+1)]
    return K,dict(Counter(off))


def multinomial(n,parts):
    if any(x<0 for x in parts) or sum(parts)!=n:return 0
    out=1
    for x in parts:
        out*=comb(n,x);n-=x
    return out


def containment_probability(p,n,j,r,d):
    """P(S subset C,T subset D) for fixed |S|=|T|=d,|S cap T|=r,
    uniform ordered n-set pairs (C,D) at Johnson distance j.
    """
    u=2*d-r
    count=0
    for a in range(d-r+1):
        for b in range(d-r+1):
            common=n-j-u+a+b
            count+=C(d-r,a)*C(d-r,b)*multinomial(p-u,[common,j-a,j-b,p-n-j])
    return Q(count,C(p,n)*C(n,j)*C(p-n,j))


def exact_spectrum(p,n,d=6):
    assert d<=n<=p//2
    K,hist=kernel_overlap_sums(p,d)
    corr=[sum(Q(K[r])*containment_probability(p,n,j,r,d) for r in range(d+1))
          for j in range(d+1)]
    theta=[[eberlein(p,n,z,j) for z in range(d+1)] for j in range(d+1)]
    energies=solve(theta,corr)
    mean=Q(p*((-1)**(d//2))*C((p-1)//2,d//2)*C(n,d),C(p,d)) if d%2==0 else Q(0)
    assert energies[0]==mean*mean and all(x>=0 for x in energies)
    assert sum(energies)==corr[0]
    variance=corr[0]-mean*mean
    return {'p':p,'n':n,'degree':d,'kernel_overlap_sums':K,
            'distinct_row_pair_histogram':{str(k):v for k,v in hist.items()},
            'distance_correlations':corr,'energies':energies,'mean':mean,'variance':variance,
            'top_energy_ratio':energies[d]/variance if variance else None,
            'lower_energy_ratio':sum(energies[1:d])/variance if variance else None}


def direct_target(c,p,d,h):
    out=0
    for y in range(p):
        coeff=[1]+[0]*d
        for x in c:
            a=h[(y-x)%p]
            for j in range(d,0,-1):coeff[j]+=a*coeff[j-1]
        out+=coeff[d]
    return out


def exhaustive_audit(p=13,n=6,d=6):
    h=chi(p)
    sets=list(combinations(range(p),n));masks=[sum(1<<x for x in c) for c in sets]
    targets=[direct_target(c,p,d,h) for c in sets]
    sums=[0]*(n+1)
    for a,va in zip(masks,targets):
        for b,vb in zip(masks,targets):sums[n-(a&b).bit_count()]+=va*vb
    correlations=[Q(sums[j],C(p,n)*C(n,j)*C(p-n,j)) for j in range(d+1)]
    spec=exact_spectrum(p,n,d)
    assert correlations==spec['distance_correlations']
    if n==d:
        assert list(reversed(sums))==spec['kernel_overlap_sums']
    normalized=[i for i,c in enumerate(sets) if 0 in c and 1 in c]
    assert Q(sum(targets),len(sets))==Q(sum(targets[i] for i in normalized),len(normalized))
    assert Q(sum(x*x for x in targets),len(sets))==Q(sum(targets[i]**2 for i in normalized),len(normalized))
    # Pointwise projectors: invert the Eberlein distance eigenmatrix, without
    # invoking the generator's walk polynomials.
    theta=[[eberlein(p,n,z,j) for z in range(d+1)] for j in range(d+1)]
    checked=[]
    for index in [0,len(sets)//3,len(sets)-1]:
        cs=masks[index];distance_sums=[0]*(n+1)
        for mask,value in zip(masks,targets):distance_sums[n-(cs&mask).bit_count()]+=value
        averages=[Q(distance_sums[j],C(n,j)*C(p-n,j)) for j in range(d+1)]
        components=solve(theta,averages)
        assert sum(components)==targets[index]
        assert components[0]==spec['mean']
        checked.append({'C':sets[index],'components':components,'averages':averages})
    return {'p':p,'n':n,'degree':d,'all_inputs':len(sets),'ordered_input_pairs':len(sets)**2,
            'normalized_family_count':len(normalized),'pointwise':checked}


def encode(value):
    if isinstance(value,Q):return [value.numerator,value.denominator]
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v) for v in value]
    return value


def main():
    start=time.perf_counter()
    audits=[exhaustive_audit(13,6,6),exhaustive_audit(11,4,2)]
    saved=json.loads((HERE/'exchange_spectrum_results.json').read_text())
    spectra=[]
    for row in saved['cases']:
        spec=exact_spectrum(row['p'],row['n'])
        if row['exhaustive_normalized']:
            assert spec['energies']==[Q(*row['gram_matrix'][i][i]) for i in range(7)]
        spec['generator_sampled_top_ratio']=row['highest_component_energy_ratio']
        spectra.append(spec)
        print(spec['p'],spec['n'],'exact high',float(spec['top_energy_ratio']),
              'exact lower',float(spec['lower_energy_ratio']),flush=True)
    # Check generator's saved pointwise examples using complete distance
    # enumeration at p13; independently reconstructed components are above.
    checked=0
    for example in saved['cases'][0]['examples']:
        p,n,d=13,6,6;h=chi(p);c=set(example['C'])
        sums=[0]*7
        for D in combinations(range(p),n):sums[n-len(c.intersection(D))]+=direct_target(D,p,d,h)
        avg=[Q(sums[j],C(n,j)*C(p-n,j)) for j in range(7)]
        theta=[[eberlein(p,n,z,j) for z in range(7)] for j in range(7)]
        assert avg==[Q(*x) for x in example['distance_averages']]
        assert solve(theta,avg)==[Q(*x) for x in example['components']]
        checked+=1
    result={'status':'independent exact audit passed; exact all-input second moments added',
            'generator_sha256':sha256((HERE/'exchange_spectrum.py').read_bytes()).hexdigest(),
            'review_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'audits':audits,'saved_pointwise_examples_checked':checked,'spectra':spectra,
            'limits':['Target-dependent diagnostic, not a target-free predictor.',
                      'Exact global L2 energies do not bound the worst-case set.',
                      'Low-degree projectors only act as claimed on the degree-at-most-six subspace.',
                      'Normalized-set expectation argument requires affine invariance of the measured quantity.'],
            'elapsed_seconds':round(time.perf_counter()-start,3)}
    (HERE/'review_exchange_results.json').write_text(json.dumps(encode(result),indent=2)+'\n')
    print('Saved independent review in',result['elapsed_seconds'],'seconds')


if __name__=='__main__':main()
