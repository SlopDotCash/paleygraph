#!/usr/bin/env python3
"""Literal set-partition controls for the parallel-uniform model."""
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
from projective_profile import optimize,pack

HERE=Path(__file__).resolve().parent


def partitions(n):
    blocks=[]
    def visit(i):
        if i==n:
            yield tuple(tuple(B) for B in blocks);return
        for j in range(len(blocks)):
            blocks[j].append(i);yield from visit(i+1);blocks[j].pop()
        blocks.append([i]);yield from visit(i+1);blocks.pop()
    yield from visit(0)


def main():
    cases=literal=thresholds=0
    for m in (3,4):
        for capacities in product((1,2),repeat=m):
            for z in (0,1):
                n=sum(capacities)+z
                if n>8:continue
                groups=[];labels=[]
                for j,c in enumerate(capacities):groups.append(list(range(len(labels),len(labels)+c)));labels.extend([j]*c)
                zeros=list(range(len(labels),n));labels.extend([None]*z)
                bests={d:{} for d in (2,3)}
                for partition in partitions(n):
                    literal+=1
                    for d in (2,3):
                        if any(len(B)>=d and (any(labels[i] is None for i in B) or len({labels[i] for i in B})!=len(B)) for B in partition):continue
                        alpha=sum(min(len(B),d-1) for B in partition);queries=sum(comb(len(B),d) for B in partition if len(B)>=d)
                        for s in range(max(d,alpha+1),n+1):bests[d][s]=min(bests[d].get(s,10**9),queries)
                for d in (2,3):
                    for s in range(d,n+1):
                        result=optimize(list(capacities),z,s,d);thresholds+=1
                        if s not in bests[d]:assert result['status']=='no_parallel_feasible_profile';continue
                        assert result['status']=='parallel_profile_found' and result['best']['queries']==bests[d][s]
                        blocks=pack(groups,zeros,result)
                        assert sorted(i for B in blocks for i in B)==list(range(n))
                        assert all(len(B)<d or (None not in [labels[i] for i in B] and len({labels[i] for i in B})==len(B)) for B in blocks)
                cases+=1
    out={'status':'passed','scope':'Exact optimum among complete arc-block partitions for small parallel extensions of uniform matroids. Abstract component controls; no polynomial realization of every pattern or higher-rank arithmetic coverage is implied.',
         'capacity_patterns':cases,'literal_coordinate_set_partitions':literal,'threshold_queries':thresholds,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['parallel_controls.py','projective_profile.py']}}
    (HERE/'parallel_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':main()
