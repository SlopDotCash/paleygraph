#!/usr/bin/env python3
"""Strictly decreasing dependency repair with explicit local swap evidence.

Only triples containing an exchanged coordinate can change. The initial full
audit and saved deltas are enough to track all defects; a final complete check
and two existing certificate implementations still verify the result.
"""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from capacity_search import inputs,projective_groups,det3,audit_queries,certify_blocks,validate

HERE=Path(__file__).resolve().parent


def repair(data,blocks,max_trials=10000):
    if data['dimension']!=3:raise ValueError('this exchange prototype is rank three')
    if type(max_trials) is not int or not 1<=max_trials<=100000:raise ValueError('trial budget must lie in1..100000')
    columns=inputs(data);p=data['p'];groups,zeros=projective_groups(columns,p)
    labels={i:j for j,g in enumerate(groups) for i in g}
    blocks=[list(B) for B in blocks];initial=audit_queries(columns,blocks,p,3)
    bad={(r['block'],tuple(r['coordinates'])) for r in initial['dependent_queries']};trials=[];skipped=0;checks=0
    while bad and len(trials)<max_trials:
        moved=False
        for a,Q in sorted(bad):
            if moved:break
            for i in Q:
                if moved:break
                for b,B in enumerate(blocks):
                    if b==a:continue
                    if moved:break
                    for j in B:
                        if len(trials)>=max_trials:break
                        left=[x for x in blocks[a] if x!=i];right=[x for x in B if x!=j]
                        if i in zeros or j in zeros or labels[j] in {labels[x] for x in left} or labels[i] in {labels[x] for x in right}:
                            skipped+=1;continue
                        removed=sorted((h,T) for h,T in bad if (h==a and i in T) or (h==b and j in T))
                        added=[];checked=0
                        for h,ids,new in [(a,left,j),(b,right,i)]:
                            for pair in combinations(ids,2):
                                T=tuple(sorted((*pair,new)));checked+=1
                                if not det3(*(columns[x] for x in T),p):added.append((h,T))
                        checks+=checked;accept=len(added)<len(removed)
                        trials.append({'swap':[a,i,b,j],'removed_dependencies':[[h,list(T)] for h,T in removed],
                                       'added_dependencies':[[h,list(T)] for h,T in sorted(added)],
                                       'local_queries_checked':checked,'accepted':accept})
                        if accept:
                            bad.difference_update(removed);bad.update(added)
                            blocks[a]=sorted(left+[j]);blocks[b]=sorted(right+[i]);moved=True;break
                    if len(trials)>=max_trials:break
        if not moved:break
    final=audit_queries(columns,blocks,p,3)
    assert bad=={(r['block'],tuple(r['coordinates'])) for r in final['dependent_queries']}
    out={'status':'found' if not bad else 'not_found','initial_audit':initial,'trials':trials,
         'class_conflicting_moves_skipped':skipped,'local_queries_checked':checks,'final_audit':final,'blocks':blocks,
         'search_queries_including_initial_scan':initial['queries_checked']+checks,
         'scope':'Bounded strict local descent; every swap preserves block sizes and projective capacities. No complete repair algorithm or novelty claim.'}
    if not bad:
        cert=certify_blocks(data,blocks);out['certificate']=cert;out['standalone_certificate_review']=validate(cert)
    return out


def main():
    source=HERE/'root_class_results.json';pre=json.loads(source.read_text())
    result=repair(pre['input'],pre['allocated_blocks'])
    result['source_sha256']={f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['exchange_repair.py','capacity_search.py','projective_profile.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_verifier.py']}
    result['input_sha256']={source.name:sha256(source.read_bytes()).hexdigest()}
    (HERE/'exchange_repair.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    if result['status']=='found':(HERE/'repaired.certificate.json').write_text(json.dumps(result['certificate'],separators=(',',':'))+'\n')
    print(json.dumps({'status':result['status'],'initial_defects':len(result['initial_audit']['dependent_queries']),
                      'rank_tested_trials':len(result['trials']),'accepted_swaps':sum(r['accepted'] for r in result['trials']),
                      'local_queries_checked':result['local_queries_checked'],'search_queries_including_initial_scan':result['search_queries_including_initial_scan']}),flush=True)


if __name__=='__main__':main()
