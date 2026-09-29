#!/usr/bin/env python3
"""Literal integer partitions and independent determinant controls."""
from functools import lru_cache
from hashlib import sha256
from itertools import permutations,product
import json
from math import comb
from pathlib import Path
import random
from arc_profile import profile
from dimension_preflight import integer_determinant,rank

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def partitions(n,minimum=1):
    if n==0:return ((),)
    return tuple((v,)+tail for v in range(minimum,n+1) for tail in partitions(n-v,v))


def leibniz(matrix):
    d=len(matrix);result=0
    for sigma in permutations(range(d)):
        term=1
        for i in range(d):term*=matrix[i][sigma[i]]
        result+=(-1 if sum(sigma[i]>sigma[j] for i in range(d) for j in range(i+1,d))%2 else 1)*term
    return result


def main():
    tested=partitions_checked=0
    for n in range(3,25):
        all_parts=partitions(n);partitions_checked+=len(all_parts)
        for d in range(2,min(n,8)+1):
            for cap in (max(d,8),64):
                best={}
                for parts in all_parts:
                    if max(parts)>cap:continue
                    alpha=sum(min(b,d-1) for b in parts);queries=sum(comb(b,d) for b in parts if b>=d)
                    for s in range(max(d,alpha+1),n+1):best[s]=min(best.get(s,10**20),queries)
                for s in range(d,n+1):
                    got=profile(n,s,d,max_block=cap);tested+=1
                    if s in best:assert got['status']=='profile_found' and got['queries']==best[s]
                    else:assert got['status']=='no_profile_within_block_cap'
    matrices=0
    for values in product((-1,0,1),repeat=9):
        a=[list(values[i:i+3]) for i in (0,3,6)]
        expected=(a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))
        assert integer_determinant(a)==expected
        for p in (2,3,17):assert (expected%p!=0)==(rank(a,p)==3)
        matrices+=1
    rng=random.Random(1406);higher=[]
    for d,count in ((4,32),(5,32),(6,16),(7,4),(8,1)):
        for _ in range(count):
            a=[[rng.randrange(-65537,65538) for j in range(d)] for i in range(d)]
            det=integer_determinant(a);assert det==leibniz(a)
            assert (det%65537!=0)==(rank(a,65537)==d)
        higher.append({'dimension':d,'matrices':count})
    out={'status':'passed','literal_integer_partitions':partitions_checked,'profile_queries_checked':tested,
         'exhaustive_ternary_three_by_three_matrices':matrices,'finite_field_rank_comparisons':3*matrices,
         'higher_dimensional_leibniz_controls':higher,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['profile_controls.py','arc_profile.py','dimension_preflight.py']}}
    (HERE/'profile_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':main()
