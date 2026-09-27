#!/usr/bin/env python3
"""Exact four-cell coarsening and cubic increment budgets.

No uniform saving for the increment or Paley proof is claimed.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]


def f3(c):return c*(c-1)*(c-2)


def shifted_excess(H,p):
    R=[(h-1)%p for h in H if h!=1]
    products=Counter(a*b%p for a in R for b in R)
    return sum(v*v for v in products.values())-(2*len(R)**2-len(R))


def coarsening(p,H,s):
    """Recover every cell of size>=2 from its ordered pair of points."""
    k=s//2
    nonidentity=[h for h in H if h!=1]
    points=set()
    for a in nonidentity:
        for b in nonidentity:
            if a!=b:
                x=(b-1)*pow((a-b)%p,-1,p)%p
                assert x not in (0,p-1)
                points.add(x)
    coarse=defaultdict(Counter)
    for x in points:
        fine=(pow(x,k,p),pow(x+1,k,p))
        label=(fine[0]*fine[0]%p,fine[1]*fine[1]%p)
        coarse[label][fine]+=1
    assert all(2<=sum(row.values()) for row in coarse.values())
    assert all(len(row)<=4 for row in coarse.values())
    assert sum(sum(row.values())*(sum(row.values())-1) for row in coarse.values())==(s-1)*(s-2)
    X_parent=sum(f3(sum(row.values())) for row in coarse.values())
    X_child=sum(f3(c) for row in coarse.values() for c in row.values())
    delta=0
    for row in coarse.values():
        c=list(row.values())
        total=sum(c)
        mixed2=3*sum(x*(x-1)*(total-x) for x in c)
        mixed3=6*sum(c[i]*c[j]*c[l] for i in range(len(c))
                     for j in range(i+1,len(c)) for l in range(j+1,len(c)))
        assert f3(total)-sum(map(f3,c))==mixed2+mixed3>=0
        delta+=mixed2+mixed3
    assert delta==X_parent-X_child
    if delta==0:
        assert all(sum(row.values())<=2 or len(row)==1 for row in coarse.values())
    return dict(X_parent=X_parent,X_child=X_child,delta=delta,
                recovered_cells=len(coarse),recovered_points=len(points),
                ordered_pair_parameters=(s-1)*(s-2))


def step(p,s,g,expected):
    H=[pow(g,j,p) for j in range(s)]
    K=H[::2]
    assert len(set(H))==s and pow(g,s//2,p)==p-1
    data=coarsening(p,H,s)
    assert data['X_parent']==shifted_excess(H,p)
    assert data['X_child']==shifted_excess(K,p)
    # Exact four children of every parent row-zero cell, including
    # the singleton cells absent from the pair reconstruction.
    blocks=defaultdict(Counter)
    for z in H:
        if z==p-1:continue
        first,second=pow(z,s//2,p),pow(z+1,s//2,p)
        blocks[second*second%p][(first,second)]+=1
    c_expected=Counter(pow(1+pow(g,j,p),s,p) for j in range(1,s//2,2))
    A_mass=Y=P=L=AB_mass=row_delta=0
    records=[]
    for label,block in blocks.items():
        u=min(second for _,second in block)
        a=block[(1,u)];b=block[(1,-u%p)]
        c=block[(p-1,u)];d=block[(p-1,-u%p)]
        assert c==d==c_expected[label]
        S=a+b
        inc=f3(S+2*c)-f3(a)-f3(b)-2*f3(c)
        expanded=3*a*b*(S-2)+6*c*c*(c-1)+6*c*S*(S+2*c-2)
        assert inc==expanded>=6*c*c*(c-1)+6*c*S
        A_mass+=c*(c-1)
        Y+=c*(c-2) if c>=3 else 0
        P+=c*c*(c-1)
        L+=c*S
        AB_mass+=a*b
        row_delta+=inc
        if inc:
            records.append(dict(label=label,children=[a,b,c,d],increment=inc))
    assert A_mass==expected['pair_collision_mass'] and Y==expected['Y']
    assert AB_mass==A_mass
    D_delta=6*A_mass+4*L
    assert Fraction(L,2)==Fraction(expected['T'],s)
    assert row_delta<=data['delta']
    assert data['delta']>=6*P>=36*Y
    assert 24*Y*Y<=s*data['delta']
    assert 2*data['delta']>=3*D_delta+6*A_mass
    if data['delta']==0:
        assert A_mass==Y==expected['T']==D_delta==0
    data.update(s=s,Y=Y,P=P,A=A_mass,row_zero_increment=row_delta,
                normalized_energy_increment=D_delta,nonzero_row_blocks=records,
                shifted_product_pair_checks=(s-1)**2+(s//2-1)**2,all_passed=True)
    return data


def main():
    source=ROOT/'results/parallel51_collision_transport_2026_09_06.json'
    old=json.loads(source.read_text())
    cache={};towers=[]
    for tower in old['towers']+[old['witness']]:
        p,N,g=tower['p'],tower['N'],tower['generator']
        steps=[]
        for level in tower['levels']:
            s=level['s'];key=(p,s)
            if key not in cache:
                cache[key]=step(p,s,pow(g,N//s,p),level)
            steps.append(cache[key])
        assert steps[0]['X_child']==0
        assert all(a['X_parent']==b['X_child'] for a,b in zip(steps,steps[1:]))
        assert sum(x['delta'] for x in steps)==steps[-1]['X_parent']
        for j in range(len(steps)):
            tail=steps[j:]
            dx=tail[-1]['X_parent']-tail[0]['X_child']
            assert dx==sum(x['delta'] for x in tail)
            assert dx>=6*sum(x['P'] for x in tail)
            assert dx>=36*sum(x['Y'] for x in tail)
            assert 2*dx>=3*sum(x['normalized_energy_increment'] for x in tail)+6*sum(x['A'] for x in tail)
        towers.append(dict(p=p,N=N,X_endpoint=steps[-1]['X_parent'],
                           increments=[dict(s=x['s'],delta=x['delta'],Y=x['Y']) for x in steps],
                           cutoffs_checked=len(steps),all_passed=True))
    witness=next(t for t in towers if t['p']==2144280833 and t['N']==256)
    assert witness['X_endpoint']==1260
    assert witness['increments'][-1]['delta']==0
    assert cache[(2144280833,128)]['X_parent']==1260
    frozen=sum(x['delta']==0 and x['X_parent']>0 for x in cache.values())
    assert Fraction(3,2)+Fraction(5,3)/2==Fraction(7,3)
    assert Fraction(3,2)+2/Fraction(2)==Fraction(5,2)
    result=dict(scope='Exact cubic increment budgets under incidence-cell coarsening; no uniform increment saving or full proof.',
                steps=[dict(p=p,**data) for (p,s),data in sorted(cache.items())],
                towers=towers,tower_cases=len(towers),distinct_step_cases=len(cache),
                tower_cutoff_checks=sum(t['cutoffs_checked'] for t in towers),
                reconstructed_parent_cells=sum(x['recovered_cells'] for x in cache.values()),
                ordered_pair_parameters=sum(x['ordered_pair_parameters'] for x in cache.values()),
                shifted_product_pair_checks=sum(x['shifted_product_pair_checks'] for x in cache.values()),
                zero_increment_steps_with_positive_excess=frozen,
                witness=witness,prior_certificate_sha256=sha256(source.read_bytes()).hexdigest(),all_passed=True)
    out=ROOT/'results/parallel52_cubic_increments_2026_09_06.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['tower_cases','distinct_step_cases','tower_cutoff_checks','reconstructed_parent_cells','ordered_pair_parameters','shifted_product_pair_checks','all_passed']}))


if __name__=='__main__':
    main()
