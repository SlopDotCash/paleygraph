#!/usr/bin/env python3
"""Full transcript/lift controls for the production conditional decoder."""
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from conditional_cover import plan,compile_plan
from conditional_pruning import run
from decoding_review import review_decoding,evaluate
from guard_profiles import guarded_cover
from toy_controls import rejected

HERE=Path(__file__).resolve().parent


def main():
    cases=[('base',[[1],[0,4,1],[0,0,4,1]],[]),
           ('nonfirst_pivot',[[0,4,1],[1],[0,0,4,1]],[4,1]),
           ('mixed_direction',[[1,4,1],[2,0,4,1],[3,3,2]],[4,1])]
    rows=[];corruptions=0
    for name,basis,origin in cases:
        data={'p':5,'n':5,'k':4,'s':4,'dimension':3,'domain':list(range(5)),'origin':origin,'basis':basis}
        params=list(product(range(5),repeat=3));words=[[(evaluate(origin,x,5)+sum(a*evaluate(b,x,5) for a,b in zip(theta,basis)))%5 for x in range(5)] for theta in params]
        count=queries=mixed=0
        for a,b,c in product(range(5),repeat=3):
            received=[a,b,c,c,c];planned=plan(data,received)
            truth=[list(theta) for theta,word in zip(params,words) if sum(x==y for x,y in zip(word,received))>=4]
            for choice in planned['candidate_plans']:
                selected=deepcopy(planned);selected['best']=choice;cert=compile_plan(selected);result=run(cert);checked=review_decoding(cert,result)
                assert [o['parameters'] for o in result['outputs']]==truth
                count+=1;queries+=result['metrics']['query_systems_solved'];mixed+=len(cert['parts'])>1
                if result['outputs'] and corruptions==0:
                    bad=deepcopy(result);bad['outputs'][0]['codeword'][0]=(bad['outputs'][0]['codeword'][0]+1)%5;rejected(lambda:review_decoding(cert,bad));corruptions+=1
                    bad=deepcopy(result);bad['outputs']=[];rejected(lambda:review_decoding(cert,bad));corruptions+=1
                    bad=deepcopy(result);bad['lifts'][0]['status']='insufficient_total_agreement';rejected(lambda:review_decoding(cert,bad));corruptions+=1
                    bad=deepcopy(result);bad['metrics']['final_whole_word_coordinate_evaluations']+=1;rejected(lambda:review_decoding(cert,bad));corruptions+=1
        rows.append({'case':name,'supported_cutoff_decodings':count,'query_systems_independently_replayed':queries,'multibranch_decodings':mixed,'complete_space_members_per_word':125})
    tiny={'p':5,'n':3,'k':2,'s':2,'dimension':2,'domain':[0,1,2],'origin':[],'basis':[[1],[0,1]]}
    assert guarded_cover(tiny)['status']=='no_guarded_profile'
    assert corruptions==4
    out={'status':'passed','cases':rows,'corrupt_decoding_artifacts_rejected':corruptions,'no_guarded_profile_control_passed':True,
         'scope':'Production decoder and independent complete branch reviewer tested against full125-polynomial oracles, with multibranch routing, affine origins, parameter changes and corrupted outputs/metrics. A tiny geometry correctly rejects an impossible singleton-reserving profile.',
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['decoding_controls.py','conditional_cover.py','conditional_verifier.py','conditional_pruning.py','decoding_review.py','guard_profiles.py','toy_controls.py','../round15/capacity_search.py','../round15/projective_profile.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_review.py','../round14/arc_verifier.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_profile.py']}}
    (HERE/'decoding_controls.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'}),flush=True)


if __name__=='__main__':main()
