#!/usr/bin/env python3
"""Test supported nonminimal cutoffs to exercise multi-branch decoding unions."""
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from conditional_cover import plan,compile_plan
from conditional_verifier import validate_conditional,evaluate
from toy_controls import decode,rejected

HERE=Path(__file__).resolve().parent


def main():
    data={'p':5,'n':5,'k':4,'s':4,'dimension':3,'domain':list(range(5)),'origin':[],'basis':[[1],[0,4,1],[0,0,4,1]]}
    params=list(product(range(5),repeat=3));words=[[(sum(a*evaluate(b,x,5) for a,b in zip(theta,data['basis'])))%5 for x in data['domain']] for theta in params]
    variants=queries=mixed=multi_fibre=corruptions=0
    for a,b,c in product(range(5),repeat=3):
        received=[a,b,c,c,c];planned=plan(data,received)
        truth=[theta for theta,word in zip(params,words) if sum(x==y for x,y in zip(word,received))>=4]
        for choice in planned['candidate_plans']:
            selected=deepcopy(planned);selected['best']=choice;cert=compile_plan(selected);reviewed=validate_conditional(cert)
            assert decode(cert)==truth;variants+=1;queries+=reviewed['queries_rank_checked']
            kinds=[part['kind'] for part in cert['parts']]
            if 'light' in kinds and 'fibre' in kinds:
                mixed+=1
                if corruptions==0:
                    bad=deepcopy(cert);bad['parts']=[p for p in bad['parts'] if p['kind']=='light'];rejected(lambda:validate_conditional(bad));corruptions+=1
                    bad=deepcopy(cert);next(p for p in bad['parts'] if p['kind']=='light')['excluded_values']=[];rejected(lambda:validate_conditional(bad));corruptions+=1
                    bad=deepcopy(cert);next(p for p in bad['parts'] if p['kind']=='fibre')['free_parameter_indices']=[2,1];rejected(lambda:validate_conditional(bad));corruptions+=1
            if kinds.count('fibre')>1:multi_fibre+=1
    assert (variants,mixed,multi_fibre,corruptions)==(250,25,100,3)
    out={'status':'passed','received_words':125,'supported_cutoffs_tested':variants,'mixed_light_and_fibre_covers':mixed,
         'multiple_fibre_covers':multi_fibre,'constituent_queries_checked':queries,'corrupt_routing_certificates_rejected':corruptions,
         'complete_space_members_per_word':125,'scope':'All supported cutoffs on125 structured received words, including deliberately nonminimal choices that exercise mixed residual/fibre and multiple-fibre decoding unions.',
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['routing_controls.py','toy_controls.py','conditional_cover.py','conditional_verifier.py','../round15/capacity_search.py','../round15/projective_profile.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_verifier.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_profile.py']}}
    (HERE/'routing_controls.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'}),flush=True)


if __name__=='__main__':main()
