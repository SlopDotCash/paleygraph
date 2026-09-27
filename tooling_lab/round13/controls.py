#!/usr/bin/env python3
"""Exhaustive coverage/output controls, transformations and corrupt certificates."""
from collections import Counter
from copy import deepcopy
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
from query_cover import compile_partition,find_partition
from cover_verifier import validate,evaluate
from deterministic_pruning import prepare,run
from pruning_review import replay

HERE=Path(__file__).resolve().parent


def expect_rejected(operation):
    try:operation()
    except (AssertionError,ValueError,KeyError,TypeError):return
    raise AssertionError('corruption or invalid input accepted')


def main():
    small=json.loads((HERE/'small.certificate.json').read_text());large=json.loads((HERE/'full_length.certificate.json').read_text())
    review=validate(small);edges=[set(q) for b in small['blocks'] for q in b['queries']];hist=[0]*17
    for size in range(17):
        for subset in combinations(range(16),size):
            S=set(subset)
            if all(not e<=S for e in edges):hist[size]+=1
    assert max(i for i,c in enumerate(hist) if c)==review['maximum_query_avoiding_set']==10
    assert hist[11:]==[0]*6
    # The exact avoiding witness proves the threshold cannot be strengthened
    # for this particular query family.
    W=set(review['query_avoiding_witness']);assert len(W)==10 and all(not e<=W for e in edges)
    corruptions={
        'missing_coordinate':lambda c:c['blocks'][0]['coordinates'].pop(),
        'overlapping_blocks':lambda c:c['blocks'][1]['coordinates'].__setitem__(0,c['blocks'][0]['coordinates'][0]),
        'missing_query':lambda c:c['blocks'][0]['queries'].pop(),
        'duplicate_query':lambda c:c['blocks'][0]['queries'].append(c['blocks'][0]['queries'][0]),
        'query_outside_block':lambda c:c['blocks'][0]['queries'][0].__setitem__(0,c['blocks'][1]['coordinates'][0]),
        'false_local_capacity':lambda c:c['blocks'][0].__setitem__('no_query_capacity',0),
        'false_histogram':lambda c:c['blocks'][0]['no_query_subset_histogram'].__setitem__(0,0),
        'false_witness':lambda c:c['blocks'][0].__setitem__('no_query_witness',[]),
        'false_total_capacity':lambda c:c.__setitem__('no_query_capacity',0),
        'false_query_count':lambda c:c.__setitem__('query_count',0),
        'threshold_at_avoiding_capacity':lambda c:c['input'].__setitem__('s',c['no_query_capacity']),
        'dependent_basis':lambda c:c['input']['basis'].__setitem__(2,c['input']['basis'][0]),
        'duplicate_domain':lambda c:c['input']['domain'].__setitem__(1,c['input']['domain'][0]),
        'nonprime_field':lambda c:c['input'].__setitem__('p',15),
        'degree_violation':lambda c:c['input'].__setitem__('origin',[0]*c['input']['k']+[1]),
    }
    for mutate in corruptions.values():
        c=deepcopy(small);mutate(c);expect_rejected(lambda:validate(c))
    # A changed polynomial basis and origin preserve the same query family.
    transformed=deepcopy(small);d=transformed['input'];p=d['p'];k=d['k'];old=d['basis'];matrix=[[1,1,0],[0,1,1],[1,0,1]]
    d['basis']=[[sum(matrix[i][j]*(old[j][r] if r<len(old[j]) else 0) for j in range(3))%p for r in range(k)] for i in range(3)]
    d['origin']=[(a+1)%p for a in d['origin']]
    transformed_review=validate(transformed)
    assert transformed_review['queries_rank_checked']==10 and transformed_review['maximum_query_avoiding_set']==10
    # Coordinate permutations transport the actual polynomial domain and blocks.
    permuted=deepcopy(small);n=permuted['input']['n'];permuted['input']['domain'].reverse()
    for b in permuted['blocks']:
        b['coordinates']=sorted(n-1-i for i in b['coordinates'])
        b['queries']=sorted(sorted(n-1-i for i in q) for q in b['queries'])
        b['no_query_witness']=sorted(n-1-i for i in b['no_query_witness'])
    permuted_review=validate(permuted)
    assert permuted_review['maximum_query_avoiding_set']==10
    common=json.loads((HERE.parent/'round12/common_roots.certificate.json').read_text())['input']
    cc=compile_partition(common,[list(range(5)),list(range(5,10)),list(range(10,16))]);cr=validate(cc)
    assert cr['maximum_query_avoiding_set']==9 and cr['queries_rank_checked']==30
    # Guaranteed failure is a search result, never an impossibility assertion.
    failed=find_partition(small['input'],[2]*8,max_attempts=1)
    assert failed['status']=='not_found' and failed['best_capacity']==16
    expect_rejected(lambda:compile_partition(small['input'],[list(range(16))]))
    expect_rejected(lambda:find_partition(small['input'],[3]*5))
    expect_rejected(lambda:run(prepare(small),[0]*15))
    expect_rejected(lambda:run(prepare(small),[17]*16))
    # Full received-word census on the nontrivial 3-dimensional F5, n5 code.
    toy={'p':5,'n':5,'k':3,'s':3,'domain':list(range(5)),'origin':[], 'basis':[[1],[0,1],[0,0,1]]}
    tc=compile_partition(toy,[list(range(5))]);tp=prepare(tc)
    words=[(params,[evaluate(params,x,5) for x in range(5)]) for params in product(range(5),repeat=3)]
    assert len({tuple(w) for _,w in words})==125
    counts=Counter();query_replays=0
    for received in product(range(5),repeat=5):
        result=run(tp,received)
        truth=[list(params) for params,word in words if sum(a==b for a,b in zip(word,received))>=3]
        assert truth==[o['parameters'] for o in result['outputs']]
        checked=replay(tc,result);query_replays+=checked['queries_replayed'];counts[len(truth)]+=1
    # Corrupt execution outputs are rejected even when coverage is valid.
    actual=json.loads((HERE/'actual_results.json').read_text())['cases'][0]['runs'][5]['run']
    alters={
        'missing_output':lambda r:r['outputs'].pop(),
        'wrong_parameters':lambda r:r['outputs'][0]['parameters'].__setitem__(0,16),
        'wrong_polynomial':lambda r:r['outputs'][0]['coefficients'].__setitem__(0,16),
        'bad_agreement':lambda r:r['outputs'][0]['agreement_support'].pop(),
        'false_query_count':lambda r:r.__setitem__('queries_executed',1),
        'false_candidate_count':lambda r:r.__setitem__('distinct_candidates_evaluated',0),
    }
    for mutate in alters.values():
        r=deepcopy(actual);mutate(r);expect_rejected(lambda:replay(small,r))
    out={'status':'passed','literal_global_subsets_checked':65536,'small_query_avoiding_subset_histogram':hist,
         'exact_query_family_threshold_checked':11,'corrupted_certificates_rejected':list(corruptions),
         'transformation_controls':[transformed_review,permuted_review],
         'common_root_control':{'certificate':cc,'review':cr},'search_exhaustion_control':failed,
         'invalid_partition_and_received_controls':4,
         'toy_received_words_exhausted':3125,'toy_space_members':125,'toy_list_size_histogram':sorted(counts.items()),
         'toy_queries_independently_replayed':query_replays,'corrupted_runs_rejected':list(alters),
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['controls.py','query_cover.py','cover_verifier.py','deterministic_pruning.py','pruning_review.py']},
         'input_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['small.certificate.json','full_length.certificate.json','actual_results.json','../round12/common_roots.certificate.json']}}
    (HERE/'controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','literal_global_subsets_checked':65536,'toy_received_words_exhausted':3125,
                      'toy_queries_replayed':query_replays,'corrupt_certificates_rejected':len(corruptions),'corrupt_runs_rejected':len(alters)}),flush=True)


if __name__=='__main__':main()
