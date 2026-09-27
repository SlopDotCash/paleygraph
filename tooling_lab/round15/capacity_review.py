#!/usr/bin/env python3
"""Separate integer review of capacity optimality, defects and repair history.

Imports neither profile optimization nor search/repair implementations. Input
evaluations use direct powers and rank tests use integer Bareiss determinants.
Each accepted swap is reviewed by scanning both entire affected blocks, not by
the producer's changed-triple update. All failed restart witnesses are checked;
unrecorded nonzero queries in failed restarts are not independently replayed.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round14'))
from arc_verifier import determinant,evaluate,validate


def defects(columns,blocks,p):
    return {(j,Q) for j,B in enumerate(blocks) for Q in combinations(B,3) if determinant([columns[i] for i in Q])%p==0}


def read_bad(rows):return {(r['block'],tuple(r['coordinates'])) for r in rows}


def review_exchange(data,original,result):
    p=data['p'];columns=[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']]
    blocks=[B[:] for B in original];bad=defects(columns,blocks,p);assert bad==read_bad(result['initial_audit']['dependent_queries'])
    local_tests=0;full_local_tests=0;accepted=0
    for record in result['trials']:
        a,i,b,j=record['swap'];assert a!=b and i in blocks[a] and j in blocks[b]
        after=[B[:] for B in blocks];after[a]=sorted([x for x in blocks[a] if x!=i]+[j]);after[b]=sorted([x for x in blocks[b] if x!=j]+[i])
        for h in (a,b):
            keys=[]
            for x in after[h]:
                v=columns[x];first=next(vv for vv in v if vv);keys.append(tuple(vv*pow(first,-1,p)%p for vv in v))
            assert len(keys)==len(set(keys))
        before_local={(h,Q) for h,Q in bad if h in (a,b)}
        # Full affected-block oracle, independently of the incremental update.
        after_local={(h,Q) for h in (a,b) for Q in combinations(after[h],3) if determinant([columns[x] for x in Q])%p==0}
        full_local_tests+=sum(comb(len(after[h]),3) for h in (a,b))
        removed={(h,tuple(T)) for h,T in record['removed_dependencies']};added={(h,tuple(T)) for h,T in record['added_dependencies']}
        assert removed=={(h,Q) for h,Q in before_local if (h==a and i in Q) or (h==b and j in Q)}
        assert added=={(h,Q) for h,Q in after_local if (h==a and j in Q) or (h==b and i in Q)}
        assert (before_local-removed)|added==after_local
        assert record['local_queries_checked']==comb(len(blocks[a])-1,2)+comb(len(blocks[b])-1,2)
        local_tests+=record['local_queries_checked']
        assert type(record['accepted']) is bool and record['accepted']==(len(after_local)<len(before_local))
        if record['accepted']:
            blocks=after;bad.difference_update(before_local);bad.update(after_local);accepted+=1
    assert blocks==result['blocks'] and bad==read_bad(result['final_audit']['dependent_queries'])
    assert local_tests==result['local_queries_checked'] and result['search_queries_including_initial_scan']==sum(comb(len(B),3) for B in original)+local_tests
    if result['status']=='found':
        assert not bad and result['certificate']['blocks']==blocks and result['certificate']['input']==data
        assert validate(result['certificate'])==result['standalone_certificate_review']
    else:assert result['status']=='not_found' and bad
    return {'initial_queries_fully_checked':sum(comb(len(B),3) for B in original),
            'affected_block_queries_fully_checked':full_local_tests,'accepted_swaps':accepted,
            'producer_local_queries_checked':local_tests,'final_certificate_queries_checked':result.get('certificate',{}).get('query_count',0)}


def main():
    names=['root_class_results.json','capacity_search.json','exchange_repair.json','root_class.certificate.json','repaired.certificate.json']
    pre,search,repair,c1,c2=[json.loads((HERE/f).read_text()) for f in names];data=pre['input'];n,s,d,p=(data[k] for k in ('n','s','dimension','p'))
    columns=[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']]
    groups=defaultdict(list);zeros=[]
    for i,v in enumerate(columns):
        first=next((a for a in v if a),None)
        if first is None:zeros.append(i)
        else:groups[tuple(a*pow(first,-1,p)%p for a in v)].append(i)
    assert sorted(sorted(g) for g in groups.values())==sorted(sorted(g) for g in pre['projective_classes'])
    assert not zeros and max(map(len,groups.values()))==62
    # Enumerate EVERY possible t and u, rather than selecting u by the
    # producer's closed formula. Convex balancing supplies the size lower bound.
    candidates=[];profiles_examined=0
    for t in range(1,n//d+1):
        for u in range(n+1):
            if t*(d-1)+u>=s or n-u<t*d:continue
            if u<len(zeros)+sum(max(len(g)-t,0) for g in groups.values()):continue
            q,r=divmod(n-u,t)
            if q+(r>0)>64:continue
            profiles_examined+=1;cost=(t-r)*comb(q,d)+r*comb(q+1,d);candidates.append((cost,t,u))
    minimum=min(candidates);assert minimum==(136155,33,29)
    best=search['profile']['best'];assert (best['queries'],best['queried_blocks'],best['unqueried_coordinates'])==minimum
    assert len(search['attempts'])==23 and search['attempts'][0]['blocks']==pre['allocated_blocks']
    witnesses=0
    for attempt in search['attempts']:
        blocks=attempt['blocks'];assert sorted(i for B in blocks for i in B)==list(range(n))
        assert sorted(map(len,blocks))==sorted(best['queried_block_sizes']+[1]*best['unqueried_coordinates'])
        assert attempt['queries_checked']==sum(comb(len(B),3) for B in blocks)==minimum[0]
        for row in attempt['dependent_queries']:
            Q=row['coordinates'];assert len(Q)==3 and Q==sorted(set(Q)) and set(Q)<=set(blocks[row['block']])
            assert determinant([columns[i] for i in Q])%p==0;witnesses+=1
    assert search['certificate']==c1 and search['attempts'][-1]['blocks']==c1['blocks']
    assert not search['attempts'][-1]['dependent_queries'] and validate(c1)==search['standalone_certificate_review']
    assert repair['certificate']==c2
    reviewed=review_exchange(data,pre['allocated_blocks'],repair)
    out={'status':'passed','scope':'Independent exhaustive t,u lower-bound enumeration, every reported failed-restart dependency, full initial/affected-block repair replay, and both final arc certificates. Nonzero tests in failed restart attempts are not independently replayed.',
         'feasible_t_u_profiles_checked':profiles_examined,'minimum_complete_arc_partition_queries':minimum[0],
         'failed_restart_witnesses_checked':witnesses,'restart_search_queries_reported':sum(r['queries_checked'] for r in search['attempts']),
         'repair_review':reviewed,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['capacity_review.py','../round14/arc_verifier.py']},
         'input_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in names}}
    (HERE/'capacity_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('_sha256')}),flush=True)


if __name__=='__main__':main()
