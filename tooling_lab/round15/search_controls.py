#!/usr/bin/env python3
"""Realized polynomial controls for capacity search and local exchange repair."""
from copy import deepcopy
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
from capacity_search import find_capacity_cover,projective_groups,inputs,det3
from projective_profile import optimize,pack
from exchange_repair import repair
from capacity_review import review_exchange

HERE=Path(__file__).resolve().parent


def interpolate(values,p):
    out=[0]*len(values)
    for x,y in enumerate(values):
        poly=[1];den=1
        for z in range(len(values)):
            if z==x:continue
            new=[0]*(len(poly)+1)
            for j,a in enumerate(poly):new[j]=(new[j]-z*a)%p;new[j+1]=(new[j+1]+a)%p
            poly=new;den=den*(x-z)%p
        scale=y*pow(den,-1,p)%p
        for j,a in enumerate(poly):out[j]=(out[j]+scale*a)%p
    return out


def expect_rejected(fn):
    try:fn()
    except (AssertionError,ValueError,KeyError,IndexError,StopIteration):return
    raise AssertionError('corruption unexpectedly accepted')


def main():
    rng=random.Random(150908);rows=[];corruptions=0
    # Distinct x columns guarantee no parallel classes; random y data generate
    # dependencies of rank two, so these controls exercise more than capacities.
    for case in range(32):
        values=[rng.randrange(7) for _ in range(7)]
        while all((values[x]-values[0])%7==x*(values[1]-values[0])%7 for x in range(7)):values=[rng.randrange(7) for _ in range(7)]
        data={'p':7,'n':7,'k':7,'s':5,'dimension':3,'domain':list(range(7)),'origin':[],'basis':[[1],[0,1],interpolate(values,7)]}
        columns=inputs(data);assert [v[2] for v in columns]==values
        groups,zeros=projective_groups(columns,7);shape=optimize(list(map(len,groups)),len(zeros),5,3);blocks=pack(groups,zeros,shape)
        result=repair(data,blocks,max_trials=64);review=review_exchange(data,blocks,result)
        restart=find_capacity_cover(data,seed=case,max_attempts=8)
        # Literal search across ALL 4+3 partitions supplies an existence oracle
        # at the optimum five-query profile. Success must match this oracle;
        # a bounded failure may occur even if a certificate exists.
        feasible=[]
        for A in combinations(range(7),4):
            B=[i for i in range(7) if i not in A]
            if all(det3(*(columns[i] for i in Q),7) for S in (A,B) for Q in combinations(S,3)):feasible.append([list(A),B])
        assert shape['best']['queries']==5
        if result['status']=='found' or restart['status']=='found':assert feasible
        row={'case':case,'values':values,'literal_optimal_covers':len(feasible),'repair_status':result['status'],'restart_status':restart['status'],'review':review};rows.append(row)
        if result['trials'] and result['status']=='found' and corruptions==0:
            bad=deepcopy(result);bad['trials'][0]['accepted']=not bad['trials'][0]['accepted'];expect_rejected(lambda:review_exchange(data,blocks,bad));corruptions+=1
            bad=deepcopy(result);bad['trials'][0]['local_queries_checked']+=1;expect_rejected(lambda:review_exchange(data,blocks,bad));corruptions+=1
            bad=deepcopy(result);bad['trials'][0]['removed_dependencies']=[];expect_rejected(lambda:review_exchange(data,blocks,bad));corruptions+=1
            bad=deepcopy(result);bad['blocks'][0][0]=bad['blocks'][0][-1];expect_rejected(lambda:review_exchange(data,blocks,bad));corruptions+=1
    # Six points on one line expose a genuine higher-flat obstruction. At s5,
    # the parallel relaxation is feasible but no complete all-s-set basis cover
    # exists. At s7, moving a singleton into the sole queried block repairs it.
    values=[0]*6+[1];base={'p':7,'n':7,'k':7,'s':5,'dimension':3,'domain':list(range(7)),'origin':[],'basis':[[1],[0,1],interpolate(values,7)]}
    columns=inputs(base);groups,zeros=projective_groups(columns,7);shape=optimize(list(map(len,groups)),0,5,3);blocks=pack(groups,zeros,shape)
    failure=repair(base,blocks,max_trials=64);assert failure['status']=='not_found';review_exchange(base,blocks,failure)
    assert all(not det3(*(columns[i] for i in Q),7) for Q in combinations(range(5),3))
    base7={**base,'s':7};shape7=optimize(list(map(len,groups)),0,7,3);blocks7=pack(groups,[],shape7)
    success=repair(base7,blocks7,max_trials=64);assert success['status']=='found';review_exchange(base7,blocks7,success)
    assert any(len(blocks7[r['swap'][2]])==1 for r in success['trials'] if r['accepted'])
    assert find_capacity_cover(base,max_queries=1)['status']=='query_budget_exceeded'
    expect_rejected(lambda:find_capacity_cover(base,max_attempts=0))
    expect_rejected(lambda:find_capacity_cover(base,max_queries=0))
    expect_rejected(lambda:repair(base,blocks,max_trials=0))
    assert corruptions==4
    out={'status':'passed','scope':'32 interpolated rank-three polynomial configurations, literal optimal-profile covers, independent full-block exchange replay, a higher-flat impossibility control, singleton repair and corrupt-evidence rejection.',
         'random_cases':rows,'literal_partitions_checked':32*35,'corrupt_exchange_histories_rejected':corruptions,
         'higher_flat_obstruction':{'input':base,'rank_deficient_s_set':list(range(5)),'queried_blocks':shape['best']['queried_blocks'],'line_size':6,'unqueried':shape['best']['unqueried_coordinates']},
         'singleton_repair':{'input':base7,'initial_blocks':blocks7,'result':success},
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['search_controls.py','capacity_search.py','projective_profile.py','exchange_repair.py','capacity_review.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_verifier.py']}}
    (HERE/'search_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'random_cases':len(rows),'literal_partitions_checked':out['literal_partitions_checked'],
                      'repair_successes':sum(r['repair_status']=='found' for r in rows),'restart_successes':sum(r['restart_status']=='found' for r in rows),
                      'corrupt_histories_rejected':corruptions,'higher_flat_obstruction_checked':True,'singleton_repair_checked':True}),flush=True)


if __name__=='__main__':main()
