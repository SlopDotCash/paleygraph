#!/usr/bin/env python3
"""Compile a class-aware balanced arc cover, retaining failed query witnesses.

The profile lower bound is exact for the parallel-class relaxation. If a cover
at its cost is found, it also attains the optimum among complete-arc block
partitions subject to the same block-size cap. Exhaustion is not nonexistence.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
import sys
from projective_profile import optimize,pack

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round14'))
from arc_certificate import inputs,certify_blocks
from dimension_preflight import rank
from arc_verifier import validate


def det3(a,b,c,p):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))%p


def projective_groups(columns,p):
    groups=defaultdict(list);zeros=[]
    for i,v in enumerate(columns):
        pivot=next((a for a in v if a),None)
        if pivot is None:zeros.append(i);continue
        inv=pow(pivot,-1,p);groups[tuple(a*inv%p for a in v)].append(i)
    return list(groups.values()),zeros


def audit_queries(columns,blocks,p,d):
    count=0;bad=[]
    for j,B in enumerate(blocks):
        for Q in combinations(B,d):
            count+=1
            dependent=(not det3(*(columns[i] for i in Q),p)) if d==3 else rank([columns[i] for i in Q],p)<d
            if dependent:bad.append({'block':j,'coordinates':list(Q)})
    return {'queries_checked':count,'dependent_queries':bad}


def find_capacity_cover(data,seed=150906,max_attempts=64,max_queries=250000,max_block=64):
    columns=inputs(data);p=data['p'];d=data['dimension']
    if type(seed) is not int or type(max_attempts) is not int or not 1<=max_attempts<=10000:raise ValueError('integer seed and1..10000 attempts required')
    if type(max_queries) is not int or not 1<=max_queries<=250000:raise ValueError('query budget must lie in1..250000')
    groups,zeros=projective_groups(columns,p)
    shape=optimize([len(g) for g in groups],len(zeros),data['s'],d,max_block)
    out={'profile':shape,'seed':seed,'attempts':[],
         'scope':'Optimal profile under projective capacities; complete rank checks certify an actual cover if found. Bounded search failure is not a nonexistence result.'}
    if shape['status']!='parallel_profile_found':return dict(out,status='no_parallel_feasible_profile')
    if shape['best']['queries']>max_queries:return dict(out,status='query_budget_exceeded',max_queries=max_queries)
    rng=random.Random(seed)
    for attempt in range(max_attempts):
        if attempt:
            rng.shuffle(groups)
            for group in groups:rng.shuffle(group)
        blocks=pack(groups,zeros,shape);audit=audit_queries(columns,blocks,p,d)
        out['attempts'].append({'attempt_index':attempt,'blocks':blocks,**audit})
        assert audit['queries_checked']==shape['best']['queries']
        if audit['dependent_queries']:continue
        cert=certify_blocks(data,blocks,max_queries);review=validate(cert)
        assert review['queries_rank_checked']==shape['best']['queries']
        return dict(out,status='found',certificate=cert,standalone_certificate_review=review,
                    complete_arc_partition_optimum_attained=True)
    return dict(out,status='not_found',complete_arc_partition_optimum_attained=False)


def main():
    source=HERE/'root_class_results.json';data=json.loads(source.read_text())['input']
    result=find_capacity_cover(data)
    result['source_sha256']={f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['capacity_search.py','projective_profile.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_verifier.py']}
    result['input_sha256']={'root_class_results.json':sha256(source.read_bytes()).hexdigest()}
    (HERE/'capacity_search.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    if result['status']=='found':
        (HERE/'root_class.certificate.json').write_text(json.dumps(result['certificate'],separators=(',',':'))+'\n')
    print(json.dumps({'status':result['status'],'attempts':len(result['attempts']),
                      'dependent_queries_by_attempt':[len(r['dependent_queries']) for r in result['attempts']],
                      'profile_queries':result['profile']['best']['queries'],
                      'complete_arc_partition_optimum_attained':result['complete_arc_partition_optimum_attained']}),flush=True)


if __name__=='__main__':main()
