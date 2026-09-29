#!/usr/bin/env python3
"""Public complete-arc cover compiler with reserved zero coordinates."""
from itertools import combinations
from math import comb,isqrt
import random
from arc_profile import profile
from dimension_preflight import rank,horner

SCHEMA='polynomial_arc_query_cover_v1'


def inputs(data):
    p,n,k,s,d=(data[x] for x in ('p','n','k','s','dimension'))
    if not (type(p) is int and 2<=p<=2**31-1 and all(p%i for i in range(2,isqrt(p)+1))
            and all(type(v) is int for v in (n,k,s,d)) and 2<=d<=min(8,k) and k<=n<=min(p,4096) and d<=s<=n):raise ValueError('invalid prime field, dimension, degree, length or agreement')
    domain=data['domain'];basis=data['basis'];origin=data['origin']
    if len(domain)!=n or len(set(domain))!=n or any(type(x) is not int or not 0<=x<p for x in domain):raise ValueError('canonical distinct evaluation points required')
    if len(basis)!=d or any(len(b)>k or any(type(a) is not int or not 0<=a<p for a in b) for b in [origin,*basis]):raise ValueError('canonical degree-bounded affine basis required')
    padded=[b+[0]*(k-len(b)) for b in basis]
    if rank(padded,p)!=d:raise ValueError('basis dimension is deficient')
    columns=[tuple(horner(b,x,p) for b in basis) for x in domain]
    return columns


def certify_blocks(data,blocks,max_queries=250000):
    columns=inputs(data);p,n,d,s=(data[x] for x in ('p','n','dimension','s'))
    if type(max_queries) is not int or not 1<=max_queries<=250000:raise ValueError('query budget must lie in1..250000')
    if not blocks or any(not B or len(B)>64 or B!=sorted(set(B)) or any(type(i) is not int or not 0<=i<n for i in B) for B in blocks):raise ValueError('canonical nonempty blocks of size at most64 required')
    if sorted(i for B in blocks for i in B)!=list(range(n)):raise ValueError('blocks must partition every coordinate once')
    count=sum(comb(len(B),d) for B in blocks if len(B)>=d);capacity=sum(min(len(B),d-1) for B in blocks)
    if count>max_queries:raise ValueError('query budget exceeded')
    if capacity>=s:raise ValueError('block profile does not guarantee requested coverage')
    for B in blocks:
        for Q in combinations(B,d):
            if rank([columns[i] for i in Q],p)!=d:raise ValueError('block contains a dependent query')
    return {'schema':SCHEMA,'input':data,'blocks':blocks,'query_count':count,
            'maximum_query_avoiding_set':capacity,'zero_coordinates':[i for i,v in enumerate(columns) if not any(v)],
            'scope':'Every s-coordinate set contains a queried basis in this declared affine space. Every d-subset of a block of size at least d is queried. No cluster discovery or globally optimal cover is asserted.'}


def find_cover(data,seed=1406,max_attempts=64,max_queries=250000):
    columns=inputs(data);p,n,s,d=(data[x] for x in ('p','n','s','dimension'))
    if type(seed) is not int or type(max_attempts) is not int or not 1<=max_attempts<=10000:raise ValueError('integer seed and1..10000 attempts required')
    if type(max_queries) is not int or not 1<=max_queries<=250000:raise ValueError('query budget must lie in1..250000')
    zeros=[i for i,v in enumerate(columns) if not any(v)];ids=[i for i in range(n) if i not in zeros]
    if s-len(zeros)<d:
        witness=(zeros[:s]+ids[:max(0,s-len(zeros))])
        assert len(witness)==s and rank([columns[i] for i in witness],p)<d
        return {'status':'all_s_sets_coverage_impossible','rank_deficient_s_set':sorted(witness),
                'scope':'This s-set contains no coordinate basis. This obstructs covering ALL s-sets, not decoding a given received word.'}
    shape=profile(len(ids),s-len(zeros),d)
    if shape['status']!='profile_found':return {'status':'no_profile_within_block_cap','profile':shape}
    if shape['queries']>max_queries:return {'status':'query_budget_exceeded','profile':shape,'max_queries':max_queries}
    rng=random.Random(seed);queries=0
    for attempt in range(1,max_attempts+1):
        rng.shuffle(ids);blocks=[];offset=0;bad=False
        for size in shape['sizes']:
            B=sorted(ids[offset:offset+size]);offset+=size;blocks.append(B)
            for Q in combinations(B,d):
                queries+=1
                if rank([columns[i] for i in Q],p)<d:bad=True;break
            if bad:break
        if bad:continue
        blocks.extend([[i] for i in zeros])
        certificate=certify_blocks(data,blocks,max_queries)
        return {'status':'found','seed':seed,'attempts':attempt,'search_rank_queries':queries,'profile':shape,'certificate':certificate}
    return {'status':'not_found','seed':seed,'attempts':max_attempts,'search_rank_queries':queries,'profile':shape,
            'scope':'Bounded search exhaustion is not a nonexistence certificate.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv) not in (2,3):raise SystemExit('Usage: arc_certificate.py polynomial_input.json [seed]')
    print(json.dumps(find_cover(json.loads(Path(sys.argv[1]).read_text()),int(sys.argv[2]) if len(sys.argv)==3 else 1406),separators=(',',':')))
