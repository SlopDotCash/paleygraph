#!/usr/bin/env python3
"""Exhaustive conditional decoding and routing on every F5-valued length5 word."""
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from conditional_cover import plan,compile_plan
from conditional_verifier import validate_conditional,evaluate
from arc_pruning import prepare,collect
from arc_filter_fast import filter_candidates

HERE=Path(__file__).resolve().parent


def decode(c):
    planned=c['plan'];data=planned['input'];p=data['p'];outside=planned['outside_coordinates'];pivot=planned['pivot'];w=planned['projective_direction'];received=planned['received']
    found=set()
    for part in c['parts']:
        prepared=prepare(part['certificate']);transcript=collect(prepared,[received[i] for i in outside]);filtered=filter_candidates(prepared,transcript)
        for row in filtered['outputs']:
            if part['kind']=='light':
                params=tuple(row['parameters'])
                if sum(a*b for a,b in zip(params,w))%p in part['excluded_values']:continue
            else:
                params=[0]*3
                for j,value in zip(part['free_parameter_indices'],row['parameters']):params[j]=value
                params[pivot]=(part['value']-sum(w[j]*params[j] for j in part['free_parameter_indices']))%p;params=tuple(params)
            word=[(evaluate(data['origin'],x,p)+sum(a*evaluate(b,x,p) for a,b in zip(params,data['basis'])))%p for x in data['domain']]
            if sum(a==b for a,b in zip(word,received))>=data['s']:found.add(params)
    return sorted(found)


def rejected(fn):
    try:fn()
    except (AssertionError,ValueError,KeyError,IndexError,StopIteration):return
    raise AssertionError('corruption accepted')


def main():
    data={'p':5,'n':5,'k':4,'s':4,'dimension':3,'domain':list(range(5)),'origin':[],'basis':[[1],[0,4,1],[0,0,4,1]]}
    all_params=list(product(range(5),repeat=3))
    all_words=[[(sum(a*evaluate(b,x,5) for a,b in zip(theta,data['basis'])))%5 for x in data['domain']] for theta in all_params]
    words=queries=0;query_hist={};branch_hist={};corruptions=0
    for values in product(range(5),repeat=5):
        c=compile_plan(plan(data,list(values)));assert c['status']=='found';checked=validate_conditional(c)
        truth=[theta for theta,word in zip(all_params,all_words) if sum(a==b for a,b in zip(word,values))>=4]
        assert decode(c)==truth;words+=1;queries+=checked['queries_rank_checked'];query_hist[c['query_count']]=query_hist.get(c['query_count'],0)+1
        branch_hist[len(c['parts'])]=branch_hist.get(len(c['parts']),0)+1
        # Literal routing for every qualifying polynomial, including all
        # values of theta.w that do not occur in the received histogram.
        for theta in truth:
            value=sum(a*b for a,b in zip(theta,c['plan']['projective_direction']))%5
            matches=[part for part in c['parts'] if (part['kind']=='fibre' and part['value']==value) or (part['kind']=='light' and value not in part['excluded_values'])]
            assert len(matches)==1
        if not corruptions:
            bad=deepcopy(c);bad['plan']['normalized_received_labels'][0][1]=(bad['plan']['normalized_received_labels'][0][1]+1)%5;rejected(lambda:validate_conditional(bad));corruptions+=1
            bad=deepcopy(c);bad['parts']=bad['parts'][1:];rejected(lambda:validate_conditional(bad));corruptions+=1
            bad=deepcopy(c);bad['parts'][0]['certificate']['input']['origin']=[1];rejected(lambda:validate_conditional(bad));corruptions+=1
            bad=deepcopy(c);bad['query_count']+=1;rejected(lambda:validate_conditional(bad));corruptions+=1
    assert words==3125 and corruptions==4
    out={'status':'passed','received_words_exhausted':words,'complete_space_members_per_word':125,
         'constituent_queries_checked':queries,'query_count_histogram':sorted(query_hist.items()),'branch_count_histogram':sorted(branch_hist.items()),
         'corrupt_conditional_certificates_rejected':corruptions,
         'scope':'Exhaustive F5 length5 decoding with a repeated projective class, all125 supplied-space polynomials per received word, exact branch routing and corrupt evidence rejection. Punctured polynomial reduction is exercised; large candidate decoding is still pending.',
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['toy_controls.py','conditional_cover.py','conditional_verifier.py','../round15/capacity_search.py','../round15/projective_profile.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_verifier.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_profile.py']}}
    (HERE/'toy_controls.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'}),flush=True)


if __name__=='__main__':main()
