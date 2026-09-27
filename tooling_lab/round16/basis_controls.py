#!/usr/bin/env python3
"""Exhaust affine parameter lifts with a nonfirst pivot and mixed direction."""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from conditional_cover import plan,compile_plan
from conditional_verifier import validate_conditional,evaluate
from toy_controls import decode

HERE=Path(__file__).resolve().parent


def main():
    records=[]
    cases=[('nonfirst_pivot',[[0,4,1],[1],[0,0,4,1]],[0,1,0]),
           ('mixed_direction',[[1,4,1],[2,0,4,1],[3,3,2]],[1,2,3])]
    for name,basis,expected_direction in cases:
        data={'p':5,'n':5,'k':4,'s':4,'dimension':3,'domain':list(range(5)),'origin':[4,1],'basis':basis}
        params=list(product(range(5),repeat=3));words=[[(evaluate(data['origin'],x,5)+sum(a*evaluate(b,x,5) for a,b in zip(theta,basis)))%5 for x in data['domain']] for theta in params]
        assert len(set(map(tuple,words)))==125
        count=queries=0
        for received in product(range(5),repeat=5):
            c=compile_plan(plan(data,list(received)));assert c['status']=='found' and c['plan']['projective_direction']==expected_direction
            reviewed=validate_conditional(c)
            truth=[theta for theta,word in zip(params,words) if sum(a==b for a,b in zip(word,received))>=4]
            assert decode(c)==truth;count+=1;queries+=reviewed['queries_rank_checked']
        records.append({'case':name,'received_words_exhausted':count,'complete_space_members_per_word':125,'queries_checked':queries,'direction':expected_direction})
        print(json.dumps(records[-1]),flush=True)
    out={'status':'passed','scope':'Exhaustive received words with affine origin, a nonfirst pivot and a projective direction with three nonzero entries. Checks the conditional parameter lift against complete supplied-space oracles.',
         'cases':records,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['basis_controls.py','toy_controls.py','conditional_cover.py','conditional_verifier.py','../round15/capacity_search.py','../round15/projective_profile.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_verifier.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_profile.py']}}
    (HERE/'basis_controls.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
