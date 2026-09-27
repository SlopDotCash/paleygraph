#!/usr/bin/env python3
"""Constructor-free field/polynomial/arc verifier using integer determinants."""
from itertools import combinations
from math import comb,isqrt

if not __debug__:raise RuntimeError('Verification requires Python assertions; do not use -O.')


def determinant(rows):
    a=[list(v) for v in rows];n=len(a);sign=1;last=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        v=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=v*a[i][j]-a[i][k]*a[k][j]
                assert numerator%last==0
                a[i][j]=numerator//last
            a[i][k]=0
        last=v
    return sign*a[-1][-1]


def evaluate(poly,x,p):return sum(a*pow(x,j,p) for j,a in enumerate(poly))%p


def validate(c,max_queries=250000):
    assert c['schema']=='polynomial_arc_query_cover_v1'
    assert type(max_queries) is int and 1<=max_queries<=250000
    data=c['input'];p,n,k,s,d=(data[x] for x in ('p','n','k','s','dimension'))
    assert type(p) is int and 2<=p<=2**31-1 and all(p%i for i in range(2,isqrt(p)+1))
    assert all(type(v) is int for v in (n,k,s,d)) and 2<=d<=min(8,k) and k<=n<=min(p,4096) and d<=s<=n
    domain=data['domain'];origin=data['origin'];basis=data['basis']
    assert len(domain)==len(set(domain))==n and all(type(x) is int and 0<=x<p for x in domain)
    assert len(basis)==d and all(len(b)<=k and all(type(v) is int and 0<=v<p for v in b) for b in [origin,*basis])
    columns=[tuple(evaluate(b,x,p) for b in basis) for x in domain]
    blocks=c['blocks'];assert blocks and all(B and len(B)<=64 for B in blocks)
    assert all(B==sorted(set(B)) and all(type(i) is int and 0<=i<n for i in B) for B in blocks)
    assert sorted(i for B in blocks for i in B)==list(range(n))
    queries=sum(comb(len(B),d) for B in blocks if len(B)>=d)
    capacity=sum(min(len(B),d-1) for B in blocks)
    assert type(c['query_count']) is int and queries==c['query_count'] and 1<=queries<=max_queries
    assert type(c['maximum_query_avoiding_set']) is int and capacity==c['maximum_query_avoiding_set']<s
    zeros=[i for i,v in enumerate(columns) if not any(v)];assert c['zero_coordinates']==zeros
    checked=0
    for B in blocks:
        for Q in combinations(B,d):
            assert determinant([columns[i] for i in Q])%p!=0
            checked+=1
    assert checked==queries
    return {'status':'verified','p':p,'n':n,'k':k,'s':s,'dimension':d,'blocks':len(blocks),
            'queries_rank_checked':checked,'maximum_query_avoiding_set':capacity,
            'minimum_guaranteed_agreement':capacity+1,'zero_coordinates':zeros,
            'scope':'Complete d-subset queries inside checked arc blocks. Deterministic coverage for every received word inside the supplied affine space only.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv)!=2:raise SystemExit('Usage: arc_verifier.py certificate.json')
    print(json.dumps(validate(json.loads(Path(sys.argv[1]).read_text())),indent=2))
