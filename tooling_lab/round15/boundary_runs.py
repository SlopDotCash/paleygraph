#!/usr/bin/env python3
"""Threshold-tight supplied-space decoding on the root-class polynomial input."""
from hashlib import sha256
import json
from pathlib import Path
import random
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round14'))
from arc_pruning import prepare,collect
from arc_filter_fast import filter_candidates
from arc_review import review_run


def main():
    path=HERE/'root_class.certificate.json';c=json.loads(path.read_text());prepared=prepare(c)
    data=c['input'];p,n,d,s=(data[x] for x in ('p','n','dimension','s'));blocks=c['blocks']
    assert (n,s,d)==(1024,96,3) and s*s<(data['k']-1)*n
    avoiding={i for B in blocks for i in B[:min(len(B),d-1)]}
    positive=next(j for j,B in enumerate(blocks) if len(B)>=d+1)
    order=sorted((j for j in range(len(blocks)) if j!=positive),key=lambda j:(len(blocks[j]),j))
    last=next(j for j in order[::-1] if len(blocks[j])>=d)
    guard=next(B[0] for B in blocks if len(B)==1)
    boundary=avoiding|{blocks[positive][d-1]}
    near=boundary-{blocks[last][d-2]}
    mismatch=(boundary-{guard})|{blocks[positive][d]}
    assert len(avoiding)==s-1 and len(boundary)==len(mismatch)==s and len(near)==s-1
    rng=random.Random(150907);records=[];artifacts={}
    for name,support,accept in [('boundary',boundary,True),('late_near_miss',near,False),('true_guard_mismatch',mismatch,True)]:
        received=[a if i in support else (a+rng.randrange(1,p))%p for i,a in enumerate(prepared['origin'])]
        start=perf_counter();transcript=collect(prepared,received);query_seconds=perf_counter()-start
        start=perf_counter();filtered=filter_candidates(prepared,transcript);filter_seconds=perf_counter()-start
        result={'transcript':transcript,'filter':filtered,'query_seconds':query_seconds,'filter_seconds':filter_seconds}
        zero=next(i for i,r in enumerate(transcript['candidates']) if r['parameters']==[0]*d);decision=filtered['decisions'][zero]
        assert decision['accepted']==accept
        if accept:assert next(o for o in filtered['outputs'] if o['parameters']==[0]*d)['agreement_support']==sorted(support)
        artifact=HERE/f'{name}.run.json';artifact.write_text(json.dumps(result,separators=(',',':'))+'\n')
        artifacts[artifact.name]=sha256(artifact.read_bytes()).hexdigest()
        print(json.dumps({'phase':'producer_done','case':name,'queries':transcript['queries_executed'],'candidates':len(transcript['candidates'])}),flush=True)
        checked=review_run(c,result)
        row={'name':name,'known_origin_agreement':sorted(support),'known_origin_accepted':accept,
             'known_origin_query_occurrences':sum(v['count'] for v in transcript['candidates'][zero]['positive_blocks']),
             'known_origin_coordinate_tests':len(decision['tested_coordinates']),
             'known_origin_initial_upper':decision['initial_agreement_upper'],'known_origin_final_upper':decision['final_agreement_upper'],
             'run_file':artifact.name,'review':checked,
             'total_explicit_coordinate_evaluations':filtered['total_explicit_coordinate_evaluations'],
             'full_scan_coordinate_evaluations':filtered['full_scan_coordinate_evaluations']}
        records.append(row);print(json.dumps({k:v for k,v in row.items() if k!='known_origin_agreement'}),flush=True)
    assert records[0]['known_origin_query_occurrences']==records[1]['known_origin_query_occurrences']==1
    assert records[1]['known_origin_coordinate_tests']>n//2
    assert records[2]['known_origin_query_occurrences']==4
    out={'status':'passed','scope':'Complete query replay and all candidate words checked within the supplied root-class space. No full ambient list claim at agreement96.',
         'cases':records,'artifact_sha256':artifacts,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['boundary_runs.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_review.py','../round14/arc_verifier.py']},
         'input_sha256':{path.name:sha256(path.read_bytes()).hexdigest()}}
    (HERE/'boundary_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
