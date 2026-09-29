#!/usr/bin/env python3
"""Public interface checks through dimension8 with unique root-bound oracles."""
from hashlib import sha256
import json
from pathlib import Path
from arc_certificate import certify_blocks,find_cover
from arc_verifier import validate,evaluate
from arc_pruning import prepare,run
from arc_review import review_run

HERE=Path(__file__).resolve().parent


def main():
    cases=[]
    for d in (6,7,8):
        data={'p':11,'n':d+2,'k':d,'s':d+1,'dimension':d,'domain':list(range(d+2)),
              'origin':[],'basis':[[0]*j+[1] for j in range(d)]}
        c=certify_blocks(data,[list(range(d+2))]);params=list(range(d));coeffs=params[:]
        received=[evaluate(coeffs,x,11) for x in data['domain']];received[-1]=(received[-1]+1)%11
        result=run(prepare(c),received);review=review_run(c,result)
        assert [o['parameters'] for o in result['filter']['outputs']]==[params]
        assert len(result['filter']['outputs'][0]['agreement_support'])==data['s']
        assert 2*data['s']-data['n']>data['k']-1
        cases.append({'name':f'dimension{d}','certificate':c,'certificate_review':validate(c),'run':result,'review':review,
                      'complete_oracle':'The planted word has s agreements. Two degree<k solutions would share at least2s-n>k-1 agreement roots, so the entire ambient list is unique.'})
    # Common factor X(X-1) gives two actual zero evaluation columns.
    data={'p':11,'n':11,'k':8,'s':10,'dimension':6,'domain':list(range(11)),'origin':[1],
          'basis':[[0]*j+[0,10,1] for j in range(6)]}
    search=find_cover(data);assert search['status']=='found';c=search['certificate'];params=list(range(6))
    coeffs=[((1 if i==0 else 0)+sum(params[j]*(data['basis'][j][i] if i<len(data['basis'][j]) else 0) for j in range(6)))%11 for i in range(8)]
    received=[evaluate(coeffs,x,11) for x in range(11)];received[-1]=(received[-1]+1)%11
    result=run(prepare(c),received);review=review_run(c,result)
    assert c['zero_coordinates']==[0,1] and [o['parameters'] for o in result['filter']['outputs']]==[params]
    assert len(result['filter']['outputs'][0]['agreement_support'])==10 and 2*data['s']-data['n']>data['k']-1
    cases.append({'name':'dimension6_common_roots','certificate':c,'certificate_review':validate(c),'run':result,'review':review,
                  'complete_oracle':'Full ambient unique list by2s-n>k-1.'})
    out={'status':'passed','cases':cases,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['higher_dimension_controls.py','arc_certificate.py','arc_profile.py','dimension_preflight.py','arc_verifier.py','arc_pruning.py','arc_review.py']}}
    (HERE/'higher_dimension_controls.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({'status':'passed','cases':len(cases),'dimensions':[r['certificate']['input']['dimension'] for r in cases]}),flush=True)


if __name__=='__main__':main()
