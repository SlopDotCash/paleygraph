#!/usr/bin/env python3
"""Locate information lost when scalar-invisible integer coordinates are relaxed."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def main():
    sources=['consecutive33.json','consecutive36.json','consecutive40.json',
             '../round22/metric_alternating32.json','../round22/metric_random32.json',
             '../round22/small_scalar_p_squared_erase1.json','../round22/small_scalar_p_squared_erase3.json']
    cases=[];inputs={}
    for rel in sources:
        path=HERE/rel;c=json.loads(path.read_text())['certificate'];inputs[rel]=sha256(path.read_bytes()).hexdigest()
        N,p,k,U,adj=c['N'],c['p'],c['k'],c['erased'],c['adj'];rows=[]
        for j,row in enumerate(c['basis']):
            v=[F(sum((adj[i-u] if i>=u else -adj[N+i-u])*x for u,x in zip(U,row)),k) for i in range(N)]
            integral=all(x.denominator==1 for x in v)
            congruent=integral and all((v[i]-v[0]*pow(c['g'],i,p))%p==0 for i in range(N))
            carry=(c['scalar_steps'][j]==0 and integral and all(x%p==0 for x in v))
            rows.append({'basis_index':j,'inverse_row':[[x.numerator,x.denominator] for x in v],
                         'integral':integral,'scalar_congruent':congruent,'integer_carry_direction':carry})
        cases.append({'source':rel,'rows':rows,'all_kernel_rows_integral_and_congruent':all(r['integral'] and r['scalar_congruent'] for r in rows),
                      'invisible_directions':len(c['invisible_directions']),
                      'integer_carry_directions':sum(r['integer_carry_direction'] for r in rows)})
    out={'status':'produced','scope':'Exact basis checks only. In the large tested inputs, invisible integer coordinates lift to p-multiple carry vectors. The continuous visible projection omits their integer values.',
         'cases':cases,'input_sha256':inputs,'source_sha256':{'lattice_information.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'lattice_information.json').write_text(json.dumps(out,separators=(',', ':'))+'\n')
    print(json.dumps([{k:v for k,v in c.items() if k!='rows'} for c in cases]),flush=True)


if __name__=='__main__':main()
