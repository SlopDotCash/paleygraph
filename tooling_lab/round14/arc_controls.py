#!/usr/bin/env python3
"""Complete word censuses, public input controls and transcript counterexamples."""
from collections import Counter
from copy import deepcopy
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import subprocess,sys
from arc_certificate import certify_blocks,find_cover
from arc_verifier import validate,determinant,evaluate
from arc_pruning import prepare,run
from arc_review import review_run

HERE=Path(__file__).resolve().parent


def rejected(f):
    try:f()
    except (AssertionError,ValueError,TypeError,KeyError,IndexError):return
    raise AssertionError('invalid input or corruption accepted')


def main():
    cert=json.loads((HERE/'small.certificate.json').read_text())
    mutations={
        'composite_field':lambda c:c['input'].__setitem__('p',15),
        'duplicate_domain':lambda c:c['input']['domain'].__setitem__(0,c['input']['domain'][1]),
        'noncanonical_coefficient':lambda c:c['input']['basis'][0].__setitem__(0,17),
        'dependent_basis':lambda c:c['input']['basis'].__setitem__(1,c['input']['basis'][0]),
        'wrong_dimension':lambda c:c['input'].__setitem__('dimension',4),
        'degree_violation':lambda c:c['input'].__setitem__('origin',[0]*c['input']['k']+[1]),
        'missing_coordinate':lambda c:c['blocks'][0].pop(),
        'overlapping_blocks':lambda c:c['blocks'][1].__setitem__(0,c['blocks'][0][0]),
        'unsorted_block':lambda c:c['blocks'][0].reverse(),
        'false_zero_coordinates':lambda c:c.__setitem__('zero_coordinates',[]),
        'false_capacity':lambda c:c.__setitem__('maximum_query_avoiding_set',0),
        'false_query_count':lambda c:c.__setitem__('query_count',1),
        'threshold_at_capacity':lambda c:c['input'].__setitem__('s',c['maximum_query_avoiding_set']),
        'wrong_schema':lambda c:c.__setitem__('schema','unverified'),
    }
    for mutate in mutations.values():
        c=deepcopy(cert);mutate(c);rejected(lambda:validate(c))
    common=json.loads((HERE.parent/'round12/common_roots.certificate.json').read_text())['input'];common={**common,'dimension':3}
    common_search=find_cover(common);assert common_search['status']=='found';common_review=validate(common_search['certificate'])
    obstruction=find_cover({**common,'s':7});assert obstruction['status']=='all_s_sets_coverage_impossible'
    W=obstruction['rank_deficient_s_set'];cols=[tuple(evaluate(b,common['domain'][i],17) for b in common['basis']) for i in W]
    assert len(W)==7 and all(determinant(rows)%17==0 for rows in combinations(cols,3))
    limit=find_cover(cert['input'],max_queries=1);assert limit['status']=='query_budget_exceeded'
    exhausted=find_cover(cert['input'],seed=1410,max_attempts=1);assert exhausted['status']=='not_found'
    rejected(lambda:find_cover(cert['input'],max_attempts=0))
    rejected(lambda:find_cover(cert['input'],max_queries=0))
    rejected(lambda:certify_blocks(cert['input'],[list(range(16))]))
    rejected(lambda:validate(cert,max_queries=1))
    rejected(lambda:run(prepare(cert),[0]*15))
    rejected(lambda:run(prepare(cert),[17]*16))
    optimized=subprocess.run([sys.executable,'-O',str(HERE/'arc_verifier.py'),str(HERE/'small.certificate.json')],capture_output=True,text=True)
    assert optimized.returncode!=0 and 'do not use -O' in optimized.stderr
    census=[];total_queries=total_words=0;mutation_sample=None
    for d in (2,3,4,5):
        data={'p':5,'n':5,'k':d,'s':d,'dimension':d,'domain':list(range(5)),'origin':[],
              'basis':[[0]*j+[1] for j in range(d)]}
        c=certify_blocks(data,[list(range(5))]);prepared=prepare(c)
        words=[(params,[evaluate(params,x,5) for x in range(5)]) for params in product(range(5),repeat=d)]
        assert len({tuple(w) for _,w in words})==5**d
        exact_word_lookup={tuple(word):params for params,word in words}
        hist=Counter();queries=0
        for received in product(range(5),repeat=5):
            result=run(prepared,received)
            if d==5:
                # Agreement5 means literal word equality. The independently
                # enumerated complete code gives an exact lookup oracle.
                truth=[list(exact_word_lookup[received])] if received in exact_word_lookup else []
            else:truth=[list(params) for params,word in words if sum(a==b for a,b in zip(word,received))>=d]
            assert truth==[o['parameters'] for o in result['filter']['outputs']]
            review=review_run(c,result);queries+=review['queries_independently_replayed'];hist[len(truth)]+=1
            if d==3 and received==(0,0,0,1,1):mutation_sample=(c,result)
        total_queries+=queries;total_words+=3125
        census.append({'dimension':d,'received_words':3125,'space_members':5**d,'queries_replayed':queries,'list_size_histogram':sorted(hist.items())})
        print(json.dumps(census[-1]),flush=True)
    # A candidate may have real agreement in a block where it never appears.
    # It can also be a true boundary solution despite appearing only once.
    data={'p':11,'n':8,'k':3,'s':5,'dimension':3,'domain':list(range(8)),'origin':[], 'basis':[[1],[0,1],[0,0,1]]}
    boundary_cert=certify_blocks(data,[list(range(4)),list(range(4,8))]);received=[0 if i in (0,1,2,4,5) else 1 for i in range(8)]
    boundary=run(prepare(boundary_cert),received);boundary_review=review_run(boundary_cert,boundary)
    zero=next(r for r in boundary['transcript']['candidates'] if r['parameters']==[0,0,0]);assert zero['positive_blocks']==[{'block':0,'count':1,'support_mask':7}]
    assert next(o for o in boundary['filter']['outputs'] if o['parameters']==[0,0,0])['agreement_support']==[0,1,2,4,5]
    # Enumerate all11^3 ambient members for this separate boundary control.
    truth=[]
    for params in product(range(11),repeat=3):
        if sum(evaluate(params,x,11)==received[x] for x in range(8))>=5:truth.append(list(params))
    assert truth==[o['parameters'] for o in boundary['filter']['outputs']]
    sample_cert,sample=mutation_sample
    alters={
        'missing_candidate':lambda r:r['transcript']['candidates'].pop(),
        'wrong_query_total':lambda r:r['transcript'].__setitem__('queries_executed',0),
        'false_occurrence_count':lambda r:r['transcript']['candidates'][0]['positive_blocks'][0].__setitem__('count',999),
        'false_union_support':lambda r:r['transcript']['candidates'][0]['positive_blocks'][0].__setitem__('support_mask',0),
        'false_initial_bound':lambda r:r['filter']['decisions'][0].__setitem__('initial_agreement_upper',0),
        'false_final_bound':lambda r:r['filter']['decisions'][0].__setitem__('final_agreement_upper',0),
        'missing_output':lambda r:r['filter']['outputs'].pop(),
        'false_evaluation_count':lambda r:r['filter'].__setitem__('candidate_coordinate_evaluations',999),
    }
    for mutate in alters.values():
        r=deepcopy(sample);mutate(r);rejected(lambda:review_run(sample_cert,r))
    boundary_index=next(i for i,r in enumerate(boundary['transcript']['candidates']) if r['parameters']==[0,0,0])
    for kind in ('false_test_outcome','missing_required_test','false_acceptance'):
        changed=deepcopy(boundary);decision=changed['filter']['decisions'][boundary_index]
        if kind=='false_test_outcome':decision['tested_coordinates'][0][1]^=1
        elif kind=='missing_required_test':decision['tested_coordinates'].pop()
        else:decision['accepted']=False
        rejected(lambda:review_run(boundary_cert,changed))
    # This received word has a true zero polynomial with exactly3 agreements,
    # so exactly one query can discover it. Omitting that query loses a true
    # output; completeness is an essential hypothesis, not a counting detail.
    zero=next(r for r in sample['transcript']['candidates'] if r['parameters']==[0,0,0])
    assert zero['positive_blocks']==[{'block':0,'count':1,'support_mask':7}]
    out={'status':'passed','corrupted_certificates_rejected':list(mutations),'invalid_api_controls':6,'optimized_python_rejected':True,
         'common_root_certificate':common_search['certificate'],'common_root_review':common_review,
         'rank_deficient_s_set_obstruction':obstruction,'query_budget_control':limit,'bounded_search_control':exhausted,
         'complete_received_word_censuses':census,'total_received_words':total_words,'total_queries_replayed':total_queries,
         'boundary_absence_and_single_occurrence_control':{'certificate':boundary_cert,'run':boundary,'review':boundary_review,'ambient_members_checked':11**3},
         'corrupted_transcripts_and_filters_rejected':list(alters),'omitting_single_discovery_query_can_lose_true_output':True,
         'corrupted_boundary_test_traces_rejected':['false_test_outcome','missing_required_test','false_acceptance'],
         'malformed_array_access_is_a_rejection':True,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['arc_controls.py','arc_certificate.py','arc_profile.py','dimension_preflight.py','arc_verifier.py','arc_pruning.py','arc_review.py']},
         'input_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['small.certificate.json','../round12/common_roots.certificate.json']}}
    (HERE/'arc_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','received_words':total_words,'queries_replayed':total_queries,'corrupted_certificates':len(mutations),'corrupted_transcripts':len(alters)}),flush=True)


if __name__=='__main__':main()
