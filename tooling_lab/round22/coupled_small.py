#!/usr/bin/env python3
"""Small exact-difference controls for the joint separation compiler."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import floor
from pathlib import Path
import json
from conditioned_erasure import discover, pair, encode

HERE=Path(__file__).resolve().parent


def main():
    cases=[]; inputs={}
    names=['p17_basic_erase2','p17_basic_erase4','p41_general_erase2','p97_general_erase2',
           'p97_general_erase4','squared_relation_erase3','scalar_p_squared_erase1','scalar_p_squared_erase3']
    for name in names:
        path=HERE/('small_'+name+'.json'); c=json.loads(path.read_text())['certificate']; inputs[path.name]=sha256(path.read_bytes()).hexdigest()
        N,p,f,U,V=c['N'],c['p'],c['relation'],c['erased'],c['visible_directions']; K=c['conditioning']['known_coordinates']
        rows=[[F(f[i-a] if i>=a else -f[N+i-a]) for a in range(N)] for i in K]
        Q=[[F(*v) for v in c['centered_coordinate_forms'][j]] for j in V]
        words=[encode(p,c['g'],f,a)[0] for a in range(p)]; actual={}; pair_count=0
        for a,b in combinations(range(p),2):
            if any(words[a][i]!=words[b][i] for i in K): continue
            t=tuple(sum(F(*c['inverse'][i][j])*(words[b][u]-words[a][u]) for i,u in enumerate(U)) for j in V)
            assert all(v.denominator==1 for v in t) and any(t)
            t=tuple(map(int,t)); pair_count+=1
            if next(v for v in t if v)<0: t=tuple(-v for v in t); a,b=b,a
            actual.setdefault(t,[a,b])
        bounds=[floor(2*min(F(*c['coordinate_radii'][j]),F(*c['conditioning']['cuts'][j]['radius']))) for j in V]
        other=[t for t in product(*(range(-b,b+1) for b in bounds)) if any(t) and next(v for v in t if v)>0 and t not in actual]
        chosen=[(t,actual[t]) for t in sorted(actual)[:6]]+[(t,None) for t in other[:6]]
        records=[]
        for t,anchors in chosen:
            pivot=next(i for i,v in enumerate(t) if v); others=[j for j in range(len(V)) if j!=pivot]
            q0=[v/F(t[pivot]) for v in Q[pivot]]
            extra=[[Q[j][a]-F(t[j],t[pivot])*Q[pivot][a] for a in range(N)] for j in others]
            cut=discover(q0,rows+extra); lam=[F(*v) for v in cut['lambda']]
            mu=[F(0)]*len(V); mu[pivot]=F(1,t[pivot])
            for j,value in zip(others,lam[len(K):]): mu[j]=-value; mu[pivot]+=F(t[j],t[pivot])*value
            norm=F(*cut['l1']); separated=F(p-1,p)*norm<1
            r={'target':list(t),'pivot':pivot,'certificate':cut,'direction':[pair(v) for v in mu],
               'actual_scalar_pair':anchors,'separated':separated}
            if cut['optimality_certified'] and norm and not separated:
                r['normalized_centered_difference']=[pair(F(*v)/norm) for v in cut['dual']]
            if anchors: assert not separated
            records.append(r)
        cases.append({'name':name,'source':path.name,'actual_shared_known_pairs':pair_count,
                      'distinct_signed_actual_differences':len(actual),'records':records})
        print(json.dumps({'name':name,'actual_pairs':pair_count,'records':len(records),'separated':sum(r['separated'] for r in records)}),flush=True)
    out={'status':'produced','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['coupled_small.py','conditioned_erasure.py']}}
    (HERE/'coupled_small.json').write_text(json.dumps(out,separators=(',', ':'))+'\n')


if __name__=='__main__': main()
