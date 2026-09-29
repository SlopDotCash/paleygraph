#!/usr/bin/env python3
"""Exhaustive query probabilities and failure-assumption controls."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from pruning_sampler import prepare,run,solve,trials_for
from pruning_verifier import replay,solve as gaussian_solve

HERE=Path(__file__).resolve().parent


def main():
    cpath=HERE/'hard_rank3.certificate.json';c=json.loads(cpath.read_text());prepared=prepare(c,46)
    sourcepath=HERE/'pruning_results.json';source=json.loads(sourcepath.read_text());records=source['small_runs'];columns=prepared['columns'];origin=prepared['origin_word'];p=c['input']['p'];rows=[];queries=0
    for record in records:
        target={tuple(r['parameters']):0 for r in record['run']['outputs']};received=record['run']['received']
        for B in combinations(range(c['input']['n']),3):
            matrix=[columns[i] for i in B];rhs=[(received[i]-origin[i])%p for i in B]
            parameters=solve(matrix,rhs,p);assert parameters==gaussian_solve(matrix,rhs,p);queries+=1
            if parameters in target:target[parameters]+=1
        for parameters,count in target.items():assert F(count,560)>=F(*prepared['verified']['success_probability_lower'])
        rows.append({'scalar':record['scalar'],'all_unordered_triples':560,'exact_recovery_counts':[{'parameters':list(parameters),'count':count} for parameters,count in sorted(target.items())]})
    record=next(r for r in records if len(r['run']['outputs'])==2);saved=record['run'];replayed=run(prepared,saved['received'],saved['queries'])
    assert all(replayed[key]==saved[key] for key in ('outputs','singular_queries','duplicate_candidates','distinct_candidates_evaluated'))
    bad=[]
    for label in ('missing_output','wrong_parameters','wrong_polynomial','bad_agreement','bad_query','missing_query','false_probability','false_trial_budget'):
        changed=deepcopy(saved)
        if label=='missing_output':changed['outputs'].pop()
        elif label=='wrong_parameters':changed['outputs'][0]['parameters'][0]=(changed['outputs'][0]['parameters'][0]+1)%p
        elif label=='wrong_polynomial':changed['outputs'][0]['coefficients'][0]=(changed['outputs'][0]['coefficients'][0]+1)%p
        elif label=='bad_agreement':changed['outputs'][0]['agreement_support'].pop()
        elif label=='bad_query':changed['queries'][0][0]=16
        elif label=='missing_query':changed['queries'].pop()
        elif label=='false_probability':changed['sampling_budget']['per_candidate_success_lower'][0]+=1
        else:changed['sampling_budget']['trials']-=1
        try:replay(c,prepared['verified'],changed)
        except (AssertionError,ValueError,IndexError):bad.append(label)
        else:raise AssertionError(('corruption accepted',label))
    zero_rejected=False
    try:trials_for(17,[0,1],40)
    except ValueError:zero_rejected=True
    assert zero_rejected and trials_for(17,[1,1],40)['trials']==1
    # An adversarial repeated triple is a replay, not independent sampling. It
    # can recover at most one candidate and must not inherit completeness.
    repeated=run(prepared,saved['received'],[saved['queries'][0]]*prepared['budget']['trials'])
    assert repeated['sampling_mode']=='replay_of_provided_queries' and len(repeated['outputs'])<2
    out={'status':'passed','exact_query_probability_cases':rows,'all_query_solutions_compared':queries,
         'deterministic_trace_replay_matched':True,'corrupted_traces_rejected':bad,'zero_probability_rejected':zero_rejected,
         'probability_one_trial_control_passed':True,'repeated_query_control_outputs':len(repeated['outputs']),
         'repeated_query_control_true_list_size':2,'repeated_queries_do_not_inherit_uniform_sampling_guarantee':True,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('pruning_controls.py','pruning_sampler.py','pruning_verifier.py','certificate_verifier.py')},
         'input_sha256':{'hard_rank3.certificate.json':sha256(cpath.read_bytes()).hexdigest(),'pruning_results.json':sha256(sourcepath.read_bytes()).hexdigest()}}
    (HERE/'pruning_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','all_query_solutions_compared':queries,'corrupted_traces_rejected':len(bad),'repeated_query_control_outputs':len(repeated['outputs'])}),flush=True)


if __name__=='__main__':main()
