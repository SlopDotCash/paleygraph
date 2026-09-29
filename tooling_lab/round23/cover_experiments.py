#!/usr/bin/env python3
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from box_cover import cover

HERE=Path(__file__).resolve().parent


def main():
    cases=[];inputs={}
    # First independent coverage control, then a larger input once preflight exists.
    for rel in ['../round22/metric_consecutive32.json','consecutive33.json','consecutive36.json','consecutive40.json']:
        path=HERE/rel;data=json.loads(path.read_text());c=data['certificate']
        inputs[rel]=sha256(path.read_bytes()).hexdigest();result=cover(c,cut_budget=64,node_budget=4000)
        name=f'cover{len(c["erased"])}';output=HERE/(name+'.json')
        output.write_text(json.dumps({'name':name,'source':rel,'cover':result},separators=(',', ':'))+'\n')
        row={'name':name,'artifact':output.name,'nodes':len(result['nodes']),'cuts':len(result['cuts']),
             'node_kinds':dict(Counter(n['kind'] for n in result['nodes'])),
             'continuous_witnesses':len(result['continuous_witness_cuts']),
             'unresolved_leaves':len(result['unresolved_leaves']),
             'universal_unique_completion':result['universal_unique_completion']}
        cases.append(row);print(json.dumps(row),flush=True)
    out={'status':'produced','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['cover_experiments.py','box_cover.py','../round22/conditioned_erasure.py']}}
    (HERE/'cover_summary.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
