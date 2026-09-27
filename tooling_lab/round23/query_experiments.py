#!/usr/bin/env python3
"""Exercise a cyclic33 query with the verified refined cover."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
from cover_review import verify

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_erasure import decode,encode


def main():
    c=json.loads((HERE/'consecutive33.json').read_text())['certificate']
    proof=json.loads((HERE/'cover33_refined.json').read_text())['cover'];checked=verify(c,proof)
    assert checked['universal_unique_completion']
    N,p,g=c['N'],c['p'],c['g'];d=len(c['erased']);records=[];a=1234567
    word=encode(p,g,c['relation'],a)[0]
    for start in [0,1,7,31,32,47,57,63]:
        erased={(start+j)%N for j in range(d)};digits=[None if j in erased else v for j,v in enumerate(word)]
        rotated={j:digits[(j+start)%N]*(1 if j+start<N else -1) for j in range(d,N)}
        out=decode(c,rotated,50000);assert out['status']=='complete'
        actual=[]
        for item in out['completions']:
            scalar=item['scalar']*pow(g,-start,p)%p;original=encode(p,g,c['relation'],scalar)[0]
            actual.append({'scalar':scalar,'digits':original})
        out['completions']=sorted(actual,key=lambda x:x['scalar']);assert [v['scalar'] for v in actual]==[a]
        records.append({'start':start,'input':digits,'rotated_known':list(map(list,rotated.items())),'output':out})
    sample=records[6];(HERE/'erasure_input.json').write_text(json.dumps(sample['input'],indent=2)+'\n')
    out={'status':'produced','cyclic_queries':records,'cover_result':checked,
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['consecutive33.json','cover33_refined.json']},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['query_experiments.py','cover_review.py','../round22/conditioned_erasure.py']}}
    (HERE/'query_experiments.json').write_text(json.dumps(out,separators=(',', ':'))+'\n')
    print(json.dumps({'cyclic_queries':len(records),'scalar_candidates_checked':sum(r['output']['scalar_candidates_checked'] for r in records)}),flush=True)


if __name__=='__main__':main()
