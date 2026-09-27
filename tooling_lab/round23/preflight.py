#!/usr/bin/env python3
"""Measure larger exact erasure boxes before choosing an enumeration method."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import floor,prod
from pathlib import Path
import json
import random
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from metric_erasure import change_basis
from conditioned_erasure import improve,decode,encode
from lattice_erasure import prepare


def main():
    source=HERE.parent/'round22'/'metric_consecutive32.json'; original=json.loads(source.read_text())['certificate']
    rng=random.Random(20260908); anchors=[rng.randrange(original['p']) for _ in range(8)]
    cases=[]
    for d in [33,36,40]:
        c=improve(change_basis(prepare(original['p'],original['g'],original['relation'],list(range(d)))))
        bounds=[floor(2*min(F(*c['coordinate_radii'][j]),F(*c['conditioning']['cuts'][j]['radius']))) for j in c['visible_directions']]
        records=[]; K=c['conditioning']['known_coordinates']
        for a in anchors:
            word=encode(c['p'],c['g'],c['relation'],a)[0]; known={j:word[j] for j in K}
            for edit in [False,True]:
                values=known.copy()
                if edit:
                    j=K[a%len(K)]; values[j]+=1 if values[j]<c['digit_bound'] else -1
                out=decode(c,values,50000)
                if not edit and out['status']=='complete': assert a in [x['scalar'] for x in out['completions']]
                records.append({'kind':'bounded_edit' if edit else 'actual','anchor_scalar':a,'known':list(map(list,values.items())),'output':out})
        name=f'consecutive{d}'; path=HERE/(name+'.json')
        path.write_text(json.dumps({'name':name,'certificate':c,'records':records,'candidate_budget':50000},separators=(',', ':'))+'\n')
        result={'name':name,'artifact':path.name,'erasures':d,'visible_directions':len(c['visible_directions']),
                'individual_difference_bounds':bounds,'difference_box_size':prod(2*b+1 for b in bounds),
                'difference_sign_representatives':(prod(2*b+1 for b in bounds)-1)//2,
                'uniform_candidate_cap':c['universal_candidate_box_cap'],
                'statuses':dict(Counter(r['output']['status'] for r in records)),
                'complete_counts':dict(Counter(str(r['output']['count']) for r in records if r['output']['status']=='complete')),
                'maximum_observed_box':max(r['output']['candidate_box_size'] for r in records),
                'scalar_candidates_checked':sum(r['output'].get('scalar_candidates_checked',0) for r in records)}
        cases.append(result);print(json.dumps(result),flush=True)
    out={'status':'produced','seed':20260908,'anchors':anchors,'cases':cases,
         'input_sha256':{'../round22/'+source.name:sha256(source.read_bytes()).hexdigest()},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['preflight.py','../round22/metric_erasure.py','../round22/conditioned_erasure.py','../round21/lattice_erasure.py']}}
    (HERE/'preflight.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__': main()
