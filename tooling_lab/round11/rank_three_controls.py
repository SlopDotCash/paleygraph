#!/usr/bin/env python3
"""Boundary and transformed polynomial controls with complete small censuses."""
from copy import deepcopy
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random
from rank_three_pinning import certify
from rank_three_verifier import validate,rank

HERE=Path(__file__).resolve().parent


def main():
    rng=random.Random(110617);rows=[];sets=0
    for p,n,k in [(5,3,3),(5,5,4),(7,6,5),(11,8,6),(17,10,8)]:
        domain=list(range(n))
        for label,f in [('conic',[0,0,1]),('cubic',[0,0,0,1])]:
            if len(f)>k:continue
            basis=[[1],[0,1],f];columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(b))%p for b in basis) for x in domain]
            if rank(columns,p)<3:continue
            triples=[sum(1<<i for i in t) for t in combinations(range(n),3) if rank([columns[i] for i in t],p)<3]
            minima=[comb(s,3) for s in range(n+1)]
            for A in range(1<<n):
                s=A.bit_count();value=comb(s,3)-sum(A&t==t for t in triples);minima[s]=min(minima[s],value);sets+=1
            for s in range(n+1):
                data={'p':p,'k':k,'s':s,'domain':domain,'basis':basis};review=validate(certify(data,100_000))
                assert review['interval']==[minima[s]]*2
                if label=='conic':assert minima[s]==comb(s,3)
                rows.append({'p':p,'n':n,'k':k,'s':s,'case':label,'minimum':minima[s]})
    data={'p':17,'k':8,'s':11,'domain':list(range(1,17)),'basis':[[1],[0,1],[0,0,1,0,0,0,1]]};target=147;transformed=[]
    for index in range(6):
        changed=deepcopy(data);rng.shuffle(changed['domain']);M=((1,index+1,2),(0,1,index+3),(0,0,1))
        padded=[b+[0]*(8-len(b)) for b in data['basis']]
        changed['basis']=[[sum(M[j][v]*padded[v][i] for v in range(3))%17 for i in range(8)] for j in range(3)]
        r=validate(certify(changed));assert r['interval']==[target]*2;transformed.append(r)
    rejected=[]
    for label in ('composite','duplicate_domain','noncanonical_domain','degree','dependent_basis','nonsimple','bad_threshold','bool_threshold'):
        bad=deepcopy(data)
        if label=='composite':bad['p']=25
        elif label=='duplicate_domain':bad['domain'][0]=bad['domain'][1]
        elif label=='noncanonical_domain':bad['domain'][0]=17
        elif label=='degree':bad['basis'][2]+=[0]*8
        elif label=='dependent_basis':bad['basis'][2]=[0,1]
        elif label=='nonsimple':bad['basis']=[[0,1],[0,0,1],[0,0,0,1]];bad['domain'][0]=0
        elif label=='bad_threshold':bad['s']=17
        else:bad['s']=True
        try:certify(bad)
        except ValueError:rejected.append(label)
        else:raise AssertionError(('invalid input accepted',label))
    out={'status':'passed','boundary_cases':rows,'literal_subsets':sets,'transformed_cases':transformed,'invalid_inputs_rejected':rejected,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('rank_three_controls.py','rank_three_pinning.py','rank_three_verifier.py')}}
    (HERE/'rank_three_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','boundary_cases':len(rows),'literal_subsets':sets,'transformed_cases':len(transformed),'invalid_inputs':len(rejected)}),flush=True)


if __name__=='__main__':main()
