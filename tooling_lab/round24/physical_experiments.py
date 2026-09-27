#!/usr/bin/env python3
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
from physical_constraints import compile_constraints
from physical_cover import cover

HERE=Path(__file__).resolve().parent


def main():
    inputs={};cases=[]
    for name,rel in [('consecutive32','../round22/metric_consecutive32.json'),
                     ('consecutive33','../round23/consecutive33.json'),('consecutive36','../round23/consecutive36.json'),
                     ('consecutive40','../round23/consecutive40.json'),('alternating32','../round22/metric_alternating32.json'),
                     ('random32','../round22/metric_random32.json')]:
        path=HERE/rel;c=json.loads(path.read_text())['certificate'];inputs[rel]=sha256(path.read_bytes()).hexdigest()
        constraints=compile_constraints(c);proof=cover(c,constraints,4000)
        artifact=name+'.json';(HERE/artifact).write_text(json.dumps({'name':name,'source':rel,'constraints':constraints,'cover':proof},separators=(',', ':'))+'\n')
        row={'name':name,'artifact':artifact,'cuts':len(constraints['cuts']),
             'cut_kinds':dict(Counter(c['kind'] for c in constraints['cuts'])),
             'complete_continuous_projection':constraints['continuous_projection_complete'],
             'nodes':len(proof['nodes']),'narrowing_steps':sum(len(n['trace']) for n in proof['nodes']),
             'unresolved_leaves':len(proof['unresolved_leaves']),
             'unresolved_integer_points_up_to_sign':proof['unresolved_integer_points_up_to_sign'],
             'universal_unique_completion':proof['universal_unique_completion']}
        cases.append(row);print(json.dumps(row),flush=True)
    out={'status':'produced','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['physical_experiments.py','physical_constraints.py','physical_cover.py']}}
    (HERE/'physical_summary.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
