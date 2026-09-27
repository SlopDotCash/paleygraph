#!/usr/bin/env python3
"""Constructor-free coverage verifier: Gaussian rank and literal local subsets."""
from itertools import combinations
from math import isqrt

if not __debug__:raise RuntimeError('Certificate verification requires Python assertions; do not use -O.')


def rank(vectors,p):
    rows=[list(v) for v in vectors];r=0
    for j in range(3):
        pivot=next((i for i in range(r,len(rows)) if rows[i][j]%p),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        inv=pow(rows[r][j]%p,-1,p);rows[r]=[a*inv%p for a in rows[r]]
        for i in range(r+1,len(rows)):
            a=rows[i][j];rows[i]=[(x-a*y)%p for x,y in zip(rows[i],rows[r])]
        r+=1
        if r==3:break
    return r


def evaluate(poly,x,p):return sum(a*pow(x,j,p) for j,a in enumerate(poly))%p


def validate(c):
    assert c['schema']=='affine_rank_three_query_cover_v1';data=c['input']
    p,n,k,s=(data[x] for x in ('p','n','k','s'));domain=data['domain'];basis=data['basis'];origin=data['origin']
    assert type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
    assert type(n) is int and type(k) is int and 3<=k<=n<=min(p,4096) and type(s) is int and 3<=s<=n
    assert len(domain)==n and len(set(domain))==n and all(type(x) is int and 0<=x<p for x in domain)
    assert len(basis)==3 and all(len(b)<=k and all(type(a) is int and 0<=a<p for a in b) for b in [origin,*basis])
    columns=[tuple(evaluate(b,x,p) for b in basis) for x in domain]
    assert rank(columns,p)==3
    records=c['blocks'];assert records
    assert all(row['coordinates'] and len(row['coordinates'])<=12 for row in records)
    assert sorted(i for row in records for i in row['coordinates'])==list(range(n))
    queries_checked=subsets_checked=capacity=0;witness=[]
    for row in records:
        block=row['coordinates'];assert block==sorted(set(block)) and all(type(i) is int and 0<=i<n for i in block)
        queries=row['queries'];assert queries==sorted(queries) and len({tuple(t) for t in queries})==len(queries)
        for t in queries:
            assert len(t)==3 and t==sorted(set(t)) and all(type(i) is int and i in block for i in t)
            assert rank([columns[i] for i in t],p)==3
        edges=[set(t) for t in queries];histogram=[];alpha=0
        for size in range(len(block)+1):
            count=0
            for S in combinations(block,size):
                subsets_checked+=1;selected=set(S)
                if all(not edge<=selected for edge in edges):count+=1
            histogram.append(count)
            if count:alpha=size
        assert row['no_query_subset_histogram']==histogram and row['no_query_capacity']==alpha
        W=row['no_query_witness'];assert W==sorted(set(W)) and len(W)==alpha and all(type(i) is int and i in block for i in W)
        assert all(not edge<=set(W) for edge in edges)
        queries_checked+=len(queries);capacity+=alpha;witness.extend(W)
    assert capacity==c['no_query_capacity']<s and queries_checked==c['query_count']
    # Independence numbers add for a disjoint union. This witness attains the
    # sum, so the coverage threshold is exact for THIS query family.
    assert len(set(witness))==capacity
    return {'status':'verified','p':p,'n':n,'k':k,'s':s,'blocks_checked':len(records),
            'queries_rank_checked':queries_checked,'literal_local_subsets_checked':subsets_checked,
            'maximum_query_avoiding_set':capacity,'query_avoiding_witness':sorted(witness),
            'minimum_guaranteed_agreement':capacity+1,
            'scope':'Deterministic query coverage in the declared affine space for every received word. No claim of a smallest query family or guaranteed construction on other spaces.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv)!=2:raise SystemExit('Usage: cover_verifier.py cover.certificate.json')
    print(json.dumps(validate(json.loads(Path(sys.argv[1]).read_text())),indent=2))
