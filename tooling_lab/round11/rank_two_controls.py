#!/usr/bin/env python3
"""Basis, coordinate, common-root and interface controls for rank-two pinning."""
from copy import deepcopy
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random
from polynomial_cluster import certify, evaluate
from rank_two_pinning import compile_pinning
from independent_verifier import validate_certificate, determinant_partition, integer_optimum

HERE=Path(__file__).resolve().parent


def multiply(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
    return out


def main():
    rng=random.Random(110917);records=[];abstract=[];extremal=[]
    for name in ('hard_f17','large_f65537'):
        data=json.loads((HERE/f'{name}.input.json').read_text());base=certify(data)
        p=data['p'];n=len(data['domain']);target=base['pinning']['minimum_injective_pairs']
        padded=[b+[0]*(data['k']-len(b)) for b in data['basis']]
        for label,M in [('swap',((0,1),(1,0))),('shear',((1,3),(0,1))),('mix',((2,3),(5,7)))]:
            assert (M[0][0]*M[1][1]-M[0][1]*M[1][0])%p
            changed=deepcopy(data)
            changed['basis']=[[(M[j][0]*padded[0][i]+M[j][1]*padded[1][i])%p for i in range(data['k'])] for j in range(2)]
            cert=certify(changed);review=validate_certificate(cert);assert review['minimum']==target
            records.append({'case':name,'control':label,**review})
        changed=deepcopy(data);rng.shuffle(changed['domain']);cert=certify(changed)
        review=validate_certificate(cert);assert review['minimum']==target
        records.append({'case':name,'control':'coordinate_permutation',**review})
        # Independent nonzero coordinate rescaling preserves an evaluation matroid,
        # but need not preserve this Reed--Solomon degree bound. Test the abstract
        # column interface only and make no polynomial-subspace assertion here.
        columns=[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']]
        scales=[rng.randrange(1,p) for _ in columns]
        scaled=[tuple(a*v%p for v in col) for a,col in zip(scales,columns)]
        out=compile_pinning(p,scaled,data['s']);groups,zeros,count=determinant_partition(scaled,p)
        minimum,transitions=integer_optimum(list(map(len,groups)),len(zeros),data['s'])
        assert minimum==target==out['minimum_injective_pairs']
        abstract.append({'case':name,'minimum':minimum,'pairs_checked':count,'dynamic_program_transitions':transitions})
    # P times span(1,X): all zero coordinates are selected before an MDS residual.
    for p,n,k in [(5,5,2),(7,7,4),(17,16,8),(37,32,12)]:
        domain=list(range(n));P=[1]
        for x in domain[:k-2]:P=multiply(P,[(-x)%p,1],p)
        for s in range(k,n+1):
            data={'p':p,'k':k,'s':s,'domain':domain,'basis':[P,[0]+P]}
            cert=certify(data);review=validate_certificate(cert)
            assert review['minimum']==comb(s-k+2,2)
            extremal.append({'p':p,'n':n,'k':k,'s':s,'minimum':review['minimum']})
    bads=[];original=json.loads((HERE/'hard_f17.input.json').read_text())
    for label in ('composite','duplicate_domain','noncanonical_domain','large_degree','dependent_basis','threshold','bool_threshold','noncanonical_coefficient'):
        bad=deepcopy(original)
        if label=='composite':bad['p']=21
        elif label=='duplicate_domain':bad['domain'][0]=bad['domain'][1]
        elif label=='noncanonical_domain':bad['domain'][0]=-1
        elif label=='large_degree':bad['basis'][0]+=[0]*(bad['k']+1)
        elif label=='dependent_basis':bad['basis'][1]=bad['basis'][0][:]
        elif label=='threshold':bad['s']=len(bad['domain'])+1
        elif label=='bool_threshold':bad['s']=True
        else:bad['basis'][0][0]=bad['p']
        try:certify(bad)
        except ValueError:bads.append(label)
        else:raise AssertionError(('invalid input accepted',label))
    out={'status':'passed','scope':'Declared polynomial subspaces; rescaling controls are abstract linear evaluation spaces only.',
         'polynomial_controls':records,'abstract_scaling_controls':abstract,'common_root_extremal_cases':extremal,'invalid_inputs_rejected':bads,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('rank_two_controls.py','polynomial_cluster.py','rank_two_pinning.py','independent_verifier.py')},
         'input_sha256':{f'{name}.input.json':sha256((HERE/f'{name}.input.json').read_bytes()).hexdigest() for name in ('hard_f17','large_f65537')}}
    (HERE/'rank_two_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','polynomial_controls':len(records),'abstract_scaling_controls':len(abstract),'extremal_cases':len(extremal),'invalid_inputs_rejected':len(bads)}),flush=True)


if __name__=='__main__':main()
