#!/usr/bin/env python3
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
from box_cover import cover

HERE=Path(__file__).resolve().parent


def main():
    cases=[];inputs={}
    for name in ['p17_basic_erase2','p17_basic_erase4','p41_general_erase2','p97_general_erase2','squared_relation_erase3','scalar_p_squared_erase1','scalar_p_squared_erase3']:
        rel='../round22/small_'+name+'.json';path=HERE/rel;inputs[rel]=sha256(path.read_bytes()).hexdigest()
        c=json.loads(path.read_text())['certificate'];proof=cover(c,cut_budget=16,node_budget=500)
        output=HERE/('cover_'+name+'.json');output.write_text(json.dumps({'name':name,'source':rel,'cover':proof},separators=(',', ':'))+'\n')
        row={'name':name,'artifact':output.name,'nodes':len(proof['nodes']),'cuts':len(proof['cuts']),
             'node_kinds':dict(Counter(n['kind'] for n in proof['nodes'])),'unresolved_leaves':len(proof['unresolved_leaves']),
             'universal_unique_completion':proof['universal_unique_completion']}
        cases.append(row);print(json.dumps(row),flush=True)
    out={'status':'produced','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['small_cover.py','box_cover.py']}}
    (HERE/'small_cover_summary.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
