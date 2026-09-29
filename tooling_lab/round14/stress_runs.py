#!/usr/bin/env python3
"""Threshold-tight and late-rejection words below the ambient Johnson bound."""
from hashlib import sha256
import json
from pathlib import Path
import random
from arc_pruning import prepare,run
from arc_review import review_run

HERE=Path(__file__).resolve().parent


def main():
    path=HERE/'d3_s128.certificate.json';c=json.loads(path.read_text());prepared=prepare(c);data=c['input'];p,n,d,s=(data[x] for x in ('p','n','dimension','s'));blocks=c['blocks']
    assert s*s<(data['k']-1)*n and c['maximum_query_avoiding_set']==s-1
    avoiding={i for B in blocks for i in B[:min(len(B),d-1)]}
    positive=next(j for j,B in enumerate(blocks) if len(B)>=d+1)
    unqueried=next(B[0] for B in blocks if len(B)==1)
    order=sorted((j for j in range(len(blocks)) if j!=positive),key=lambda j:(len(blocks[j]),j))
    last=next(j for j in order[::-1] if len(blocks[j])>=d)
    boundary=avoiding|{blocks[positive][d-1]}
    near_miss=boundary-{blocks[last][d-2]}
    guard_mismatch=(boundary-{unqueried})|{blocks[positive][d]}
    assert len(boundary)==len(guard_mismatch)==s and len(near_miss)==s-1
    rng=random.Random(140614);records=[];artifacts={}
    for name,support,should_accept in [('boundary_true',boundary,True),('late_near_miss',near_miss,False),('true_guard_mismatch',guard_mismatch,True)]:
        received=[a if i in support else (a+rng.randrange(1,p))%p for i,a in enumerate(prepared['origin'])]
        result=run(prepared,received)
        zero_index=next(i for i,r in enumerate(result['transcript']['candidates']) if r['parameters']==[0]*d)
        zero_decision=result['filter']['decisions'][zero_index]
        assert zero_decision['accepted']==should_accept
        if should_accept:assert next(o for o in result['filter']['outputs'] if o['parameters']==[0]*d)['agreement_support']==sorted(support)
        checked=review_run(c,result)
        artifact=HERE/f'stress_{name}.run.json';artifact.write_text(json.dumps(result,separators=(',',':'))+'\n');artifacts[artifact.name]=sha256(artifact.read_bytes()).hexdigest()
        row={'name':name,'known_origin_agreement':sorted(support),'known_origin_accepted':should_accept,
             'known_origin_query_occurrences':sum(v['count'] for v in result['transcript']['candidates'][zero_index]['positive_blocks']),
             'known_origin_coordinate_tests':len(zero_decision['tested_coordinates']),
             'known_origin_initial_upper':zero_decision['initial_agreement_upper'],
             'known_origin_final_upper':zero_decision['final_agreement_upper'],
             'run_file':artifact.name,'review':checked,
             'total_explicit_coordinate_evaluations':result['filter']['total_explicit_coordinate_evaluations'],
             'full_scan_coordinate_evaluations':result['filter']['full_scan_coordinate_evaluations']}
        records.append(row);print(json.dumps({k:v for k,v in row.items() if k!='known_origin_agreement'}),flush=True)
    assert records[0]['known_origin_query_occurrences']==records[1]['known_origin_query_occurrences']==1
    assert records[1]['known_origin_coordinate_tests']>n//2
    assert records[2]['known_origin_query_occurrences']==4
    out={'status':'passed','scope':'Deliberately planted boundary and late-rejection received words at agreement128, below the ambient Johnson threshold. All candidates and query solutions independently checked in the supplied space. No complete oracle outside that space and no uniform filtering speedup are claimed.',
         'certificate':path.name,'cases':records,'artifact_sha256':artifacts,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['stress_runs.py','arc_pruning.py','arc_review.py','arc_verifier.py']},
         'input_sha256':{path.name:sha256(path.read_bytes()).hexdigest()}}
    (HERE/'stress_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
