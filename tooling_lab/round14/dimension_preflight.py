#!/usr/bin/env python3
"""General-dimensional arc certificates on declared enlargements of actual data.

Producer uses modular elimination and Horner evaluation. A separate review uses
integer Bareiss determinants and direct power evaluation. This preflight is not
yet the standalone general-dimensional pruning interface.
"""
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random
from time import perf_counter
from arc_profile import profile

HERE=Path(__file__).resolve().parent


def rank(rows,p):
    if not rows:return 0
    a=[list(v) for v in rows];r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]%p),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        inv=pow(a[r][j]%p,-1,p);a[r]=[v*inv%p for v in a[r]]
        for i in range(r+1,len(a)):
            c=a[i][j];a[i]=[(v-c*w)%p for v,w in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r


def horner(poly,x,p):
    a=0
    for v in poly[::-1]:a=(a*x+v)%p
    return a


def integer_determinant(rows):
    # Fraction-free integer elimination. No modular inverses or field row
    # reduction; divisibility is checked at every update.
    a=[list(v) for v in rows];n=len(a);sign=1;previous=1
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        v=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=v*a[i][j]-a[i][k]*a[k][j]
                assert numerator%previous==0
                a[i][j]=numerator//previous
            a[i][k]=0
        previous=v
    return sign*a[-1][-1]


def extend(data,d):
    k=data['k'];p=data['p'];basis=[b+[0]*(k-len(b)) for b in data['basis']];added=[]
    assert rank(basis,p)==len(basis)==3
    for j in range(k):
        if len(basis)==d:break
        monomial=[int(i==j) for i in range(k)]
        if rank(basis+[monomial],p)>len(basis):basis.append(monomial);added.append(j)
    assert len(basis)==d and rank(basis,p)==d
    return {**data,'basis':basis,'dimension':d},added


def construct(data,shape,seed=1406,max_attempts=64):
    p,n,d=data['p'],data['n'],data['dimension'];ids=list(range(n));rng=random.Random(seed)
    columns=[tuple(horner(b,x,p) for b in data['basis']) for x in data['domain']]
    checked=0
    for attempt in range(1,max_attempts+1):
        rng.shuffle(ids);blocks=[];offset=0;bad=False
        for size in shape['sizes']:
            block=sorted(ids[offset:offset+size]);offset+=size;blocks.append(block)
            for Q in combinations(block,d):
                checked+=1
                if rank([columns[i] for i in Q],p)<d:bad=True;break
            if bad:break
        if not bad:return {'status':'found','seed':seed,'attempts':attempt,'rank_queries_during_search':checked,'blocks':blocks}
    return {'status':'not_found','seed':seed,'attempts':max_attempts,'rank_queries_during_search':checked}


def review(data,blocks):
    p,n,d,s=data['p'],data['n'],data['dimension'],data['s']
    assert len(data['basis'])==d and all(len(b)<=data['k'] for b in data['basis'])
    assert len(data['domain'])==len(set(data['domain']))==n
    assert sorted(i for B in blocks for i in B)==list(range(n))
    assert all(B==sorted(set(B)) and len(B)<=64 for B in blocks)
    columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(b))%p for b in data['basis']) for x in data['domain']]
    checked=0
    for block in blocks:
        for Q in combinations(block,d):
            assert integer_determinant([columns[i] for i in Q])%p
            checked+=1
    capacity=sum(min(len(B),d-1) for B in blocks);assert capacity<s
    return {'status':'verified_preflight','dimension':d,'queries_rank_checked':checked,
            'maximum_query_avoiding_set':capacity,'minimum_guaranteed_agreement':capacity+1,
            'scope':'Complete arc-block coverage for this declared polynomial space, independently rechecked by integer determinants. General public input validation and generic decoding remain pending.'}


def main():
    sourcepath=HERE.parent/'round12/full_length.certificate.json';source=json.loads(sourcepath.read_text())['input']
    profiles=[profile(1024,410,d) for d in range(3,9)];rows=[]
    for d in (3,4,5):
        data,added=extend(source,d);shape=profiles[d-3]
        start=perf_counter();search=construct(data,shape);construction_seconds=perf_counter()-start
        row={'dimension':d,'added_monomial_degrees':added,'profile':shape,'search':search,'construction_seconds':construction_seconds}
        if search['status']=='found':
            start=perf_counter();checked=review(data,search['blocks']);review_seconds=perf_counter()-start
            assert checked['queries_rank_checked']==shape['queries'] and checked['maximum_query_avoiding_set']==shape['no_query_capacity']
            artifact=HERE/f'dimension{d}.preflight.json'
            artifact.write_text(json.dumps({'input':data,'blocks':search['blocks'],'review':checked},separators=(',',':'))+'\n')
            row.update({'review':checked,'review_seconds':review_seconds,'artifact':artifact.name})
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k not in ('profile','search') }|{'search_status':search['status'],'attempts':search['attempts']}),flush=True)
    out={'status':'preflight_passed','scope':'Exact block-profile calculation through dimension8; actual large declared-space arc certificates tested only at dimensions3,4,5. Added monomials deliberately enlarge the saved space; they are not newly discovered clusters. Generic decoder and standalone certificate interface remain pending.',
         'profiles':profiles,'cases':rows,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['dimension_preflight.py','arc_profile.py']},
         'input_sha256':{'../round12/full_length.certificate.json':sha256(sourcepath.read_bytes()).hexdigest()},
         'artifact_sha256':{r['artifact']:sha256((HERE/r['artifact']).read_bytes()).hexdigest() for r in rows if 'artifact' in r}}
    (HERE/'dimension_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')


if __name__=='__main__':main()
