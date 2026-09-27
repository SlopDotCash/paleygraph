#!/usr/bin/env python3
"""Decode a checked conditional cover while retaining every branch transcript."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from conditional_verifier import validate_conditional
from arc_pruning import prepare,collect,evaluate
from arc_filter_fast import filter_candidates

HERE=Path(__file__).resolve().parent


def run(c,on_branch=None):
    checked=validate_conditional(c);planned=c['plan'];data=planned['input'];p,n,k=(data[x] for x in ('p','n','k'))
    outside=planned['outside_coordinates'];received=planned['received'];w=planned['projective_direction'];pivot=planned['pivot']
    branches=[];outputs=[];lifts=[];seen=set();whole_tests=0
    for index,part in enumerate(c['parts']):
        prepared=prepare(part['certificate'])
        start=perf_counter();transcript=collect(prepared,[received[i] for i in outside]);query_seconds=perf_counter()-start
        start=perf_counter();filtered=filter_candidates(prepared,transcript);filter_seconds=perf_counter()-start
        result={'transcript':transcript,'filter':filtered,'query_seconds':query_seconds,'filter_seconds':filter_seconds}
        branches.append(result)
        if on_branch:on_branch(index,result)
        for j,row in enumerate(filtered['outputs']):
            if part['kind']=='light':
                params=list(row['parameters']);value=sum(a*b for a,b in zip(params,w))%p
                if value in part['excluded_values']:
                    lifts.append({'branch':index,'branch_output':j,'parameters':params,'status':'excluded_heavy_value','equation_value':value});continue
            else:
                params=[0]*3
                for t,a in zip(part['free_parameter_indices'],row['parameters']):params[t]=a
                params[pivot]=(part['value']-sum(w[t]*params[t] for t in part['free_parameter_indices']))%p
                value=part['value'];assert sum(a*b for a,b in zip(params,w))%p==value
            assert tuple(params) not in seen;seen.add(tuple(params))
            coefficients=[((data['origin'][j] if j<len(data['origin']) else 0)+sum(params[t]*(data['basis'][t][j] if j<len(data['basis'][t]) else 0) for t in range(3)))%p for j in range(k)]
            word=[evaluate(coefficients,x,p) for x in data['domain']];whole_tests+=n
            support=[i for i,(a,b) in enumerate(zip(word,received)) if a==b]
            assert [word[i] for i in outside]==row['codeword']
            accept=len(support)>=data['s']
            lifts.append({'branch':index,'branch_output':j,'parameters':params,'status':'accepted' if accept else 'insufficient_total_agreement','equation_value':value,'agreement_count':len(support)})
            if accept:outputs.append({'parameters':params,'coefficients':coefficients,'codeword':word,'agreement_support':support})
    outputs.sort(key=lambda r:r['parameters'])
    metrics={'query_systems_solved':sum(r['transcript']['queries_executed'] for r in branches),
             'branch_candidates_generated':sum(len(r['transcript']['candidates']) for r in branches),
             'branch_candidate_coordinate_evaluations':sum(r['filter']['candidate_coordinate_evaluations'] for r in branches),
             'branch_output_word_coordinate_evaluations':sum(r['filter']['output_word_coordinate_evaluations'] for r in branches),
             'final_whole_word_coordinate_evaluations':whole_tests,
             'branch_full_scan_coordinate_evaluations':sum(r['filter']['full_scan_coordinate_evaluations'] for r in branches)}
    metrics['total_explicit_coordinate_evaluations']=metrics['branch_candidate_coordinate_evaluations']+metrics['branch_output_word_coordinate_evaluations']+whole_tests
    assert metrics['query_systems_solved']==checked['queries_rank_checked']
    return {'schema':'conditional_arc_pruning_v1','certificate_review':checked,'branches':branches,'lifts':lifts,'outputs':outputs,'metrics':metrics,
            'scope':'Complete decoding within the certified supplied space. Counts distinguish query solving, adaptive branch tests, branch word materialization and final whole-word checks; certificate construction/validation costs are separate.'}


def main():
    records=[];bindings={};artifacts={}
    for name in ['boundary','late_near_miss','true_guard_mismatch']:
        path=HERE/f'{name}.conditional.json';c=json.loads(path.read_text());bindings[path.name]=sha256(path.read_bytes()).hexdigest()
        def emit(index,result):
            print(json.dumps({'phase':'branch_produced','case':name,'branch':index,'kind':c['parts'][index]['kind'],'queries':result['transcript']['queries_executed'],'candidates':len(result['transcript']['candidates']),'coordinate_evaluations':result['filter']['total_explicit_coordinate_evaluations']}),flush=True)
        result=run(c,emit);artifact=HERE/f'{name}.decoded.json';artifact.write_text(json.dumps(result,separators=(',',':'))+'\n');artifacts[artifact.name]=sha256(artifact.read_bytes()).hexdigest()
        row={'case':name,'file':artifact.name,'outputs':len(result['outputs']),'metrics':result['metrics']};records.append(row);print(json.dumps(row),flush=True)
    out={'status':'produced','scope':'All three conditional decodings produced; separate Cramer/full-word review required before acceptance.',
         'cases':records,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['conditional_pruning.py','conditional_verifier.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_verifier.py']},'input_sha256':bindings,'artifact_sha256':artifacts}
    (HERE/'conditional_decoding.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
