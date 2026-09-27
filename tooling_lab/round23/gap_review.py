#!/usr/bin/env python3
from hashlib import sha256
from pathlib import Path
from copy import deepcopy
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify


def verify(c,certificate):
    prepared=certify(c);N,p,f=c['N'],c['p'],c['relation']
    H,D,z=certificate['centered_difference_lift'],certificate['digit_difference'],certificate['integer_lattice_coordinates']
    assert len(H)==len(D)==N and len(z)==len(c['erased'])
    assert all(type(v)is int for v in H+D+z)
    delta=certificate['scalar_difference'];assert delta!=0
    assert all(abs(v)<=p-1 and (v-delta*pow(c['g'],i,p))%p==0 for i,v in enumerate(H))
    M=sp.Matrix([[f[(i-j)%N]*(1 if i>=j else -1) for j in range(N)] for i in range(N)])
    assert M*sp.Matrix(H)==p*sp.Matrix(D)
    assert all(D[i]==0 for i in c['conditioning']['known_coordinates'])
    assert sp.Matrix(1,len(z),z)*sp.Matrix(c['basis'])==sp.Matrix(1,len(z),[D[u] for u in c['erased']])
    assert [z[j] for j in c['visible_directions']]==certificate['visible_target']
    assert max(map(abs,D))==certificate['max_digit_difference']<=2*c['digit_bound']


def main():
    c=json.loads((HERE/'consecutive40.json').read_text())['certificate']
    data=json.loads((HERE/'gap_certificate.json').read_text());verify(c,data)
    orbit=json.loads((HERE/'orbit_large.json').read_text());review=json.loads((HERE/'orbit_review.json').read_text())
    assert data['visible_target']==orbit['target'] and data['scalar_difference']==orbit['delta']==review['large_scalar_difference']
    assert all(data['centered_difference_lift'][i]==v for i,v in orbit['fixed'])
    assert review['status']=='passed' and review['large_count']==orbit['output']['count']==data['actual_pairs_with_this_visible_target']==0
    assert review['input_sha256']['orbit_large.json']==sha256((HERE/'orbit_large.json').read_bytes()).hexdigest()
    rejected=[]
    for name in ['changed_lift','changed_lattice_coordinate']:
        changed=deepcopy(data)
        if name=='changed_lift':changed['centered_difference_lift'][0]+=1
        else:changed['integer_lattice_coordinates'][0]+=1
        try:verify(c,changed)
        except AssertionError:rejected.append(name)
        else:raise AssertionError('corrupt gap accepted')
    out={'status':'passed','integral_lift_checked':True,'full_integer_lattice_checked':True,
         'scalar_congruences_checked':c['N'],'erased_coordinates':40,'actual_pair_count':0,
         'scalar_difference':data['scalar_difference'],'max_digit_difference':data['max_digit_difference'],
         'corrupt_controls_rejected':rejected,
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['consecutive40.json','gap_certificate.json','orbit_large.json','orbit_review.json']},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['gap_review.py','../round22/conditioned_review.py']}}
    (HERE/'gap_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
