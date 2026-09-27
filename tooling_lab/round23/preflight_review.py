#!/usr/bin/env python3
from fractions import Fraction as F
from hashlib import sha256
from math import floor,prod
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify,replay,rational
from metric_review import check_transform


def main():
    path=HERE/'preflight.json';summary=json.loads(path.read_text());inputs={path.name:sha256(path.read_bytes()).hexdigest()};cases=[]
    for row in summary['cases']:
        path=HERE/row['artifact'];data=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
        c=data['certificate'];prepared=certify(c);check_transform(c)
        bounds=[floor(2*min(rational(c['coordinate_radii'][j]),rational(c['conditioning']['cuts'][j]['radius']))) for j in c['visible_directions']]
        assert bounds==row['individual_difference_bounds'] and prod(2*b+1 for b in bounds)==row['difference_box_size']
        count=0
        for r in data['records']:
            out=replay(c,dict(r['known']),data['candidate_budget'],prepared);assert out==r['output']
            if out['status']=='complete' and r['kind']=='actual':assert r['anchor_scalar'] in [v['scalar'] for v in out['completions']]
            count+=out.get('scalar_candidates_checked',0)
        result={'name':data['name'],'queries_replayed':len(data['records']),'scalar_candidates_checked':count,
                'cuts_verified':len(c['erased']),'optimality_certificates':sum(x['optimality_certified'] for x in c['conditioning']['cuts']),
                'difference_box_size':row['difference_box_size'],'uniform_candidate_cap':c['universal_candidate_box_cap']}
        cases.append(result);print(json.dumps(result),flush=True)
    out={'status':'passed','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['preflight_review.py','../round22/conditioned_review.py','../round22/metric_review.py']}}
    (HERE/'preflight_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
