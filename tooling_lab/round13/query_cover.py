#!/usr/bin/env python3
"""Small verified query families for one supplied affine rank-three space.

Partition/Turan coverage and coordinate interpolation are known mechanisms.
Search is bounded and may fail. A successful certificate has deterministic
coverage, independent of the randomness used to find its partition.
"""
from itertools import combinations
from math import isqrt
import random


def evaluate(poly,x,p):
    v=0
    for a in poly[::-1]:v=(v*x+a)%p
    return v


def determinant(a,b,c,p):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))%p


def inputs(data):
    p,n,k,s=(data[x] for x in ('p','n','k','s'));domain=data['domain'];basis=data['basis'];origin=data['origin']
    if not (type(p) is int and 2<=p<=2**31-1 and all(p%i for i in range(2,isqrt(p)+1))
            and type(n) is int and type(k) is int and 3<=k<=n<=min(p,4096) and type(s) is int and 3<=s<=n):raise ValueError('invalid field, length, degree or threshold')
    if len(domain)!=n or len(set(domain))!=n or any(type(x) is not int or not 0<=x<p for x in domain):raise ValueError('distinct canonical evaluation domain required')
    if len(basis)!=3 or any(len(b)>k or any(type(a) is not int or not 0<=a<p for a in b) for b in [origin,*basis]):raise ValueError('three canonical degree-bounded basis polynomials required')
    columns=[tuple(evaluate(b,x,p) for b in basis) for x in domain]
    return columns


def local_record(block,columns,p):
    queries=[list(t) for t in combinations(block,3) if determinant(*(columns[i] for i in t),p)]
    positions={v:i for i,v in enumerate(block)}
    masks=[sum(1<<positions[i] for i in t) for t in queries]
    histogram=[0]*(len(block)+1);witness=[];capacity=0
    for subset in range(1<<len(block)):
        if any(subset&q==q for q in masks):continue
        size=subset.bit_count();histogram[size]+=1
        if size>capacity:capacity=size;witness=[block[j] for j in range(len(block)) if subset>>j&1]
    return {'coordinates':list(block),'queries':queries,'no_query_capacity':capacity,
            'no_query_witness':witness,'no_query_subset_histogram':histogram}


def compile_partition(data,blocks):
    columns=inputs(data);n=data['n']
    if not blocks or any(not block or len(block)>12 or block!=sorted(set(block)) or any(type(i) is not int or not 0<=i<n for i in block) for block in blocks):raise ValueError('canonical nonempty blocks of size at most12 required')
    if sorted(i for block in blocks for i in block)!=list(range(n)):raise ValueError('blocks must partition all coordinates')
    rows=[local_record(block,columns,data['p']) for block in blocks]
    cap=sum(row['no_query_capacity'] for row in rows)
    if cap>=data['s']:raise ValueError('partition gives insufficient coverage')
    return {'schema':'affine_rank_three_query_cover_v1','input':data,'blocks':rows,
            'query_count':sum(len(row['queries']) for row in rows),'no_query_capacity':cap,
            'guarantee':'Every coordinate set of size at least s contains a listed independent triple. The queries recover every nearby member of this supplied affine space for every received word.'}


def find_partition(data,sizes,seed=1306,max_attempts=1000):
    columns=inputs(data);n=data['n'];p=data['p']
    if not sizes or any(type(size) is not int or not 1<=size<=12 for size in sizes) or sum(sizes)!=n:raise ValueError('block sizes must sum to n and lie in1..12')
    if type(seed) is not int or type(max_attempts) is not int or not 1<=max_attempts<=10000:raise ValueError('integer seed and bounded search required')
    rng=random.Random(seed);ids=list(range(n));best=n
    for attempt in range(1,max_attempts+1):
        rng.shuffle(ids);offset=0;rows=[]
        for size in sizes:
            block=sorted(ids[offset:offset+size]);offset+=size
            rows.append(local_record(block,columns,p))
        cap=sum(row['no_query_capacity'] for row in rows);best=min(best,cap)
        if cap<data['s']:
            return {'status':'found','seed':seed,'attempts':attempt,'best_capacity':best,
                    'certificate':{'schema':'affine_rank_three_query_cover_v1','input':data,'blocks':rows,
                                   'query_count':sum(len(row['queries']) for row in rows),'no_query_capacity':cap,
                                   'guarantee':'Every coordinate set of size at least s contains a listed independent triple. Coverage is deterministic once the certificate is checked; this is only the supplied affine space.'}}
    return {'status':'not_found','seed':seed,'attempts':max_attempts,'best_capacity':best,
            'scope':'Search exhaustion proves no impossibility statement.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv) not in (3,4):raise SystemExit('Usage: query_cover.py polynomial_input.json block_sizes.json [seed]')
    data=json.loads(Path(sys.argv[1]).read_text());sizes=json.loads(Path(sys.argv[2]).read_text())
    print(json.dumps(find_partition(data,sizes,int(sys.argv[3]) if len(sys.argv)==4 else 1306),separators=(',',':')))
