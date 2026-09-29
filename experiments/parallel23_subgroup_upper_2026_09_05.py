#!/usr/bin/env python3
"""Actual subgroup product-ratio fibers and opposite-free six-term relations."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import factorial, isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def check(ok,name):
    assert ok,name
    COUNTS[name]+=1


def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))


def subgroup(p,n):
    check(prime(p) and (p-1)%n==0,'prime_and_divisibility')
    for a in range(2,p):
        z=pow(a,(p-1)//n,p)
        if pow(z,n//2,p)!=1:
            H=sorted(pow(z,j,p) for j in range(n))
            check(len(set(H))==n and pow(z,n,p)==1,'exact_dyadic_order')
            return H
    raise AssertionError('generator')


def weight(t):
    out=factorial(len(t))
    for v in Counter(t).values():out//=factorial(v)
    return out


def opposite_free(word,p):
    s=set(word)
    return all((-a)%p not in s for a in s)


def balanced_some_partition(word,p):
    # Each unordered 3+3 split has a unique half containing position zero.
    for other in combinations(range(1,6),2):
        I={0,*other}
        left=right=1
        for j,a in enumerate(word):
            if j in I:left=left*a%p
            else:right=right*a%p
        if (left+right)%p==0:return True
    return False


def normalized_six(H,p):
    n=len(H)
    pair_buckets=defaultdict(list)
    for a,b in combinations_with_replacement(H,2):
        pair_buckets[(1+a+b)%p].append((a,b,1 if a==b else 2))
    fibers=Counter();r6=balanced_fixed=balanced_union=0;witness=None
    inverses={h:pow(h,-1,p) for h in H}
    for triple in combinations_with_replacement(H,3):
        candidates=pair_buckets.get(-sum(triple)%p,())
        if not candidates:continue
        wt=weight(triple)
        triple_product=triple[0]*triple[1]%p*triple[2]%p
        inv_product=inverses[triple_product]
        for a,b,wp in candidates:
            word=(1,a,b,*triple)
            mass=n*wp*wt
            rho=-a*b*inv_product%p
            fibers[rho]+=mass
            COUNTS['normalized_weighted_six_terms']+=1
            if opposite_free(word,p):
                r6+=mass
                if rho==1:balanced_fixed+=mass
                any_balanced=balanced_some_partition(word,p)
                if any_balanced:balanced_union+=mass
                elif witness is None:witness=list(word)
                COUNTS['opposite_free_weighted_six_terms']+=1
    return fibers,r6,balanced_fixed,balanced_union,witness


def shifted_energy(H,p):
    R=[(h-1)%p for h in H if h!=1]
    prod_R=Counter(a*b%p for a in R for b in R)
    X=sum(v*v for v in prod_R.values())-(2*(len(H)-1)**2-(len(H)-1))
    prod_all=Counter((a-1)*(b-1)%p for a in H for b in H)
    w=Counter({(z-1)%p:c for z,c in prod_all.items()})
    check(sum(prod_all.values())==len(H)**2,'product_distribution_mass')
    check(sum(c*c for c in prod_all.values())==6*len(H)**2-9*len(H)+4+X,'zero_product_and_permutation_ledger')
    return X,prod_all,w


def coset_matrix(H,p):
    idx={};reps=[]
    for a in range(1,p):
        if a in idx:continue
        j=len(reps);reps.append(a)
        for h in H:idx[a*h%p]=j
    C=[Counter() for _ in reps]
    for x in range(1,p):
        if x==p-1:continue
        C[idx[x]][idx[(x+1)%p]]+=1
    k=C[0]
    C2={j:sum(v*k.get(s,0) for s,v in row.items()) for j,row in enumerate(C)}
    r2=Counter((a+b)%p for a in H for b in H)
    r3=Counter()
    for x,c in r2.items():
        for a in H:r3[(x+a)%p]+=c
    n=len(H)
    for j,a in enumerate(reps):
        check(r3[a]==n*(j==0)+C2[j],'coset_triple_count_identity')
    E3=sum(c*c for c in r3.values())
    ledger=n**3+2*n*n*sum(c*c for c in k.values())+n*n*k.get(0,0)**2+n*sum(c*c for c in C2.values())
    check(E3==ledger,'rooted_fourth_walk_ledger')
    return max(c for row in C for c in row.values()),sum(c*(c-1)*(c-2) for row in C for c in row.values())


def inspect(p,n):
    H=subgroup(p,n)
    check(-1%p in H and all(a*b%p in H for a in H for b in H),'actual_multiplicative_closure')
    fibers,R6,Bfix,Bunion,witness=normalized_six(H,p)
    X,prod,w=shifted_energy(H,p)
    Eshift=sum(c*c for c in prod.values())
    diag=6*n**3-9*n*n+4*n
    check(fibers[1]==n*Eshift,'equal_sum_product_shifted_energy_bijection_count')
    check(fibers[1]-diag==n*X,'nontrivial_equal_product_count')
    for rho in H:
        inv=pow(rho,-1,p)
        corr=sum(c*w.get(z*inv%p,0) for z,c in w.items())
        check(fibers[rho]==n*corr,'every_product_ratio_correlation')
        check(corr<=Eshift,'product_ratio_cauchy_bound')
    check(Bfix<=n*X and Bunion<=10*n*X,'opposite_free_balanced_upper_bounds')
    r2=Counter((a+b)%p for a in H for b in H)
    E2=sum(c*c for c in r2.values());E3=sum(fibers.values())
    T4=3*n*n-3*n;T6=15*n**3-45*n*n+40*n
    check(E3==T6+(15*n-60)*(E2-T4)+60*n*(r2[2]-1)-30*n*(3 in H)+R6,'prior_sixth_decomposition')
    if X==0:
        check(Bfix==Bunion==0,'circularity_removes_balanced_six')
        check(E2==T4,'circularity_intrinsic_fourth_energy')
        check(r2[1]==0,'circular_dyadic_no_zero_triple')
        check(E3<=2*n**4+n**3-4*n*n,'circular_explicit_fourth_power_bound')
        check(R6<=2*n**4-14*n**3+41*n*n-40*n,'circular_explicit_opposite_free_bound')
    if witness:
        check(sum(witness)%p==0 and opposite_free(witness,p) and not balanced_some_partition(witness,p),'actual_unbalanced_six_witness')
    maxC=directX=None
    if p<=2017:
        maxC,directX=coset_matrix(H,p)
        check(directX==X and (maxC<=2)==(X==0),'direct_cyclotomic_collision_identity')
    # Direct sum/product buckets give another independent check in small groups.
    if n<=16:
        bucket=Counter(((a+b+c)%p,a*b%p*c%p) for a,b,c in product(H,repeat=3))
        check(sum(v*v for v in bucket.values())==fibers[1],'direct_triple_sum_product_buckets')
    return {'p':p,'n':n,'generator_set':H,'quartic_window':n**4//4<=p<=n**4,
            'X':X,'shifted_multiplicative_energy_including_zero':Eshift,'E2':E2,'E3':E3,
            'permutation_triple_energy':diag,'T6':T6,'R6':R6,'fixed_partition_balanced_R6':Bfix,
            'some_partition_balanced_R6':Bunion,'unbalanced_R6':R6-Bunion,
            'unbalanced_witness_normalized_to_1':witness,'direct_max_cyclotomic_cell':maxC,
            'product_ratio_fibers':{str(r):fibers[r] for r in H}}


def main():
    cases=[inspect(p,n) for p,n in [(17,8),(97,8),(1049,8),(2017,8),(17393,16),
                                    (6700417,64),(7204033,64),(67403009,128),(1073748737,256)]]
    inputs=['experiments/parallel23_subgroup_upper_2026_09_05.py',
            'research/parallel21-subgroup-next-input-2026-09-05.md',
            'research/parallel22-subgroup-independent-review-2026-09-05.md',
            'research/mixed-periods-and-shifted-energy.md',
            'sources/mixed-periods-2026-09-04/shkredov-1504.04522.html']
    result={'status':'PASS','scope':'Actual subgroup partial upper bound; total opposite-free R6 remains uncontrolled at Gaussian order.',
            'counts':dict(COUNTS),'cases':cases,'input_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}}
    path=ROOT/'results/parallel23_subgroup_upper_2026_09_05.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','counts':dict(COUNTS),'cases':[{k:x[k] for k in ['p','n','X','E2','E3','R6','fixed_partition_balanced_R6','some_partition_balanced_R6','unbalanced_witness_normalized_to_1']} for x in cases]},indent=2))

if __name__=='__main__':main()
