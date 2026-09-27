#!/usr/bin/env python3
"""Verify orbit intersections via actual pair differences, not interval filtering."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify,rational
from coupled_review import target_certificate,setup


def main():
    summarypath=HERE/'orbit_summary.json';summary=json.loads(summarypath.read_text());inputs={summarypath.name:sha256(summarypath.read_bytes()).hexdigest()}
    small=[]
    for case in summary['small_cases']:
        path=HERE/case['artifact'];data=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
        p,g,N=data['p'],data['g'],data['N'];m=p//2
        vectors=[[(a*pow(g,i,p)+m)%p-m for i in range(N)] for a in range(p)]
        checks=0
        for r in data['records']:
            delta=r['delta'];fixed=dict(r['fixed'])
            truth=[a for a in range(p) if all(vectors[a][i]-vectors[(a-delta)%p][i]==H for i,H in fixed.items())]
            assert r['output']['status']=='complete' and r['output']['count']==len(truth)
            assert sorted(r['output']['scalars_a'])==truth and not r['output']['truncated'];checks+=p
        small.append({'p':p,'queries_checked':len(data['records']),'direct_scalar_pair_checks':checks})
    path=HERE/summary['large_artifact'];large=json.loads(path.read_text());inputs[path.name]=sha256(path.read_bytes()).hexdigest()
    source=HERE/'consecutive40.json';c=json.loads(source.read_text())['certificate'];inputs[source.name]=sha256(source.read_bytes()).hexdigest();certify(c)
    cp=HERE/'cover40.json';cover=json.loads(cp.read_text())['cover'];inputs[cp.name]=sha256(cp.read_bytes()).hexdigest()
    record=cover['cuts'][large['source_cover_cut']];_,V,rows,Q=setup(c)
    target_certificate(c,record,rows,Q);assert large['target']==record['target']
    N,p,g,f,U=c['N'],c['p'],c['g'],c['relation'],c['erased'];m=p//2
    # Independently identify free directions by solving the inverse image
    # with the already saved exact basis image and verifying M*v=p*L_j.
    information=json.loads((HERE/'lattice_information.json').read_text())
    images=next(case['rows'] for case in information['cases'] if case['source']=='consecutive40.json')
    inputs['lattice_information.json']=sha256((HERE/'lattice_information.json').read_bytes()).hexdigest()
    M=sp.Matrix([[f[(i-j)%N]*(1 if i>=j else -1) for j in range(N)] for i in range(N)])
    free=[];physical=[]
    for j,item in enumerate(images):
        v=[rational(x) for x in item['inverse_row']];H=[0]*N
        for u,x in zip(U,c['basis'][j]):H[u]=p*x
        assert M*sp.Matrix(v)==sp.Matrix(H) and all(x.denominator==1 for x in v);physical.append(list(map(int,v)))
        if j in c['invisible_directions']:
            indices=[i for i,x in enumerate(v) if x];assert len(indices)==1 and abs(v[indices[0]])==p;free+=indices
    assert sorted(free)==large['free_coordinates'] and len(free)==len(set(free))
    full=[sum(t*physical[j][i] for j,t in zip(V,large['target'])) for i in range(N)]
    fixed={i:full[i] for i in range(N) if i not in free};assert list(map(list,fixed.items()))==large['fixed']
    delta=full[0]%p;assert delta==large['delta'] and delta!=0
    assert all((full[i]-delta*pow(g,i,p))%p==0 for i in range(N))
    seed=min(fixed,key=lambda i:(p-abs(fixed[i]),i));H=fixed[seed]
    # Enumerate centered F_b at the seed; form a=b+delta, and compare actual
    # centered differences at every other coordinate. No q-box inequalities.
    low=max(-m,-m-H);high=min(m,m-H);size=high-low+1
    assert size==large['output']['seed_candidates']==27904523
    assert (p-1)**2+m<np.iinfo(np.int64).max
    order=sorted(fixed,key=lambda i:(p-abs(fixed[i]),i));inverse=pow(g,-seed,p);count=0;retained=[]
    checks=0
    for start in range(low,high+1,131072):
        y=np.arange(start,min(high+1,start+131072),dtype=np.int64)
        b=(y*inverse)%p;a=(b+delta)%p;checks+=len(a)
        for i in order:
            power=pow(g,i,p);difference=(a*power+m)%p-(b*power+m)%p
            keep=difference==fixed[i];a=a[keep];b=b[keep]
            if not len(a):break
        count+=len(a)
        if len(retained)<100:retained.extend(map(int,a[:100-len(retained)]))
    assert count==large['output']['count'] and sorted(retained)==sorted(large['output']['scalars_a'])
    out={'status':'passed','small_cases':small,'large_seed_candidates_checked':checks,'large_count':count,
         'large_scalar_difference':delta,'fixed_coordinates':len(fixed),'free_integer_carry_coordinates':len(free),
         'scope':'Complete actual pair count for this one visible difference target, with independent scalar-pair arithmetic. This does not cover all remaining40-erasure differences.',
         'input_sha256':inputs,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['orbit_review.py','../round22/conditioned_review.py','../round22/coupled_review.py']}}
    (HERE/'orbit_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
