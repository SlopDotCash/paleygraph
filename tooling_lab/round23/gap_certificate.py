#!/usr/bin/env python3
"""An integral difference satisfying the geometric model but no actual pair."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def main():
    c=json.loads((HERE/'consecutive40.json').read_text())['certificate']
    orbit=json.loads((HERE/'orbit_large.json').read_text());p,N,f=c['p'],c['N'],c['relation'];m=p//2
    fixed=dict(orbit['fixed']);delta=orbit['delta']
    H=[fixed[i] if i in fixed else (delta*pow(c['g'],i,p)+m)%p-m for i in range(N)]
    assert all(abs(v)<=p-1 and (v-delta*pow(c['g'],i,p))%p==0 for i,v in enumerate(H))
    raw=[sum((f[i-j] if i>=j else -f[N+i-j])*H[j] for j in range(N)) for i in range(N)]
    assert all(v%p==0 for v in raw);D=[v//p for v in raw]
    assert all(D[i]==0 for i in c['conditioning']['known_coordinates'])
    z=[sum(F(*c['inverse'][i][j])*D[u] for i,u in enumerate(c['erased'])) for j in range(len(c['erased']))]
    assert all(v.denominator==1 for v in z)
    assert [int(z[j]) for j in c['visible_directions']]==orbit['target']
    assert max(map(abs,D))<=2*c['digit_bound'] and orbit['output']['count']==0
    out={'status':'produced','p':p,'N':N,'erased':c['erased'],'scalar_difference':delta,
         'centered_difference_lift':H,'digit_difference':D,'integer_lattice_coordinates':list(map(int,z)),
         'visible_target':orbit['target'],'max_digit_difference':max(map(abs,D)),
         'actual_pairs_with_this_visible_target':0,
         'scope':'This integral bounded lift passes the full erasure-lattice geometry and scalar congruences but cannot be an actual difference of two codewords agreeing outside the erasures. Exact orbit-pair count supplies the additional obstruction.',
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['consecutive40.json','orbit_large.json','orbit_review.json']},
         'source_sha256':{'gap_certificate.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'gap_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'delta':delta,'max_digit_difference':out['max_digit_difference'],'nonzero_digit_coordinates':sum(bool(x) for x in D)}),flush=True)


if __name__=='__main__':main()
