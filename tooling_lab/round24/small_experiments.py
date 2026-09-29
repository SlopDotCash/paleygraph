#!/usr/bin/env python3
from hashlib import sha256
from pathlib import Path
import json
from physical_constraints import compile_constraints
from physical_cover import cover

HERE=Path(__file__).resolve().parent


def main():
    cases=[];inputs={}
    for name in ['p17_basic_erase2','p17_basic_erase4','p41_general_erase2','p97_general_erase2','squared_relation_erase3','scalar_p_squared_erase1','scalar_p_squared_erase3']:
        rel='../round22/small_'+name+'.json';path=HERE/rel;c=json.loads(path.read_text())['certificate'];inputs[rel]=sha256(path.read_bytes()).hexdigest()
        constraints=compile_constraints(c);proof=cover(c,constraints,1000);artifact='small_'+name+'.json'
        (HERE/artifact).write_text(json.dumps({'name':name,'source':rel,'constraints':constraints,'cover':proof},separators=(',', ':'))+'\n')
        row={'name':name,'artifact':artifact,'cuts':len(constraints['cuts']),'nodes':len(proof['nodes']),
             'unresolved_leaves':len(proof['unresolved_leaves']),'unresolved_integer_points_up_to_sign':proof['unresolved_integer_points_up_to_sign'],
             'universal_unique_completion':proof['universal_unique_completion']}
        cases.append(row);print(json.dumps(row),flush=True)
    out={'status':'produced','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['small_experiments.py','physical_constraints.py','physical_cover.py']}}
    (HERE/'small_summary.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
