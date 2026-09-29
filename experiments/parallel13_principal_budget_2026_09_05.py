#!/usr/bin/env python3
"""Exact principal-rank budgets and full-energy checks; no proof of the open edge."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import importlib.util
import json

ROOT=Path(__file__).resolve().parents[1]


def module(name, path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


P9=module('p9','experiments/parallel9_all_degrees_2026_09_05.py')
P10=module('p10','experiments/parallel10_word_aggregate_2026_09_05.py')
P11=module('p11','experiments/parallel11_full_energy_2026_09_05.py')


def principal_levels(a, depth):
    """Sum local types over all words, using linearity instead of enumeration."""
    h,n=2**(a-1),2**a-1
    plus,minus=Counter({1:1}),Counter()
    rank,words=1,1
    output=[]
    for m in range(depth+1):
        assert rank==sum(k*v for k,v in (plus+minus).items())
        output.append(dict(depth=m,words=words,principal_rank_sum=rank,
                           plus=dict(plus),minus=dict(minus)))
        lengths=set(plus)|set(minus)
        gp=Counter({k:(h-1)*plus[k]+h*minus[k] for k in lengths})
        gm=Counter({k:h*plus[k]+(h-1)*minus[k] for k in lengths})
        next_rank=a*n*rank+n*words-(a+1)*sum(gp.values())
        next_plus=Counter({k+1:v for k,v in gm.items() if v})
        next_minus=Counter({k-1:v for k,v in gp.items() if k>1 and v})
        next_plus[1]=next_rank-sum(k*v for k,v in (next_plus+next_minus).items())
        assert next_plus[1]>=0
        rank,words,plus,minus=next_rank,n*words,next_plus,next_minus
    return output


class Quadratic:
    """Exact Q(lambda), with lambda^2=C lambda-N and lambda the larger root."""
    def __init__(self,a):
        self.a=a;self.h=2**(a-1);self.n=2**a-1
        self.c=F(self.h*(a*a+1),a)
    def add(self,x,y):return x[0]+y[0],x[1]+y[1]
    def scale(self,x,s):return x[0]*s,x[1]*s
    def sub(self,x,y):return self.add(x,self.scale(y,-1))
    def mul(self,x,y):
        return x[0]*y[0]-self.n*x[1]*y[1],x[0]*y[1]+x[1]*y[0]+self.c*x[1]*y[1]
    def inv(self,x):
        norm=x[0]*x[0]+self.c*x[0]*x[1]+self.n*x[1]*x[1]
        assert norm
        return F(x[0]+self.c*x[1],norm),F(-x[1],norm)
    def sign(self,x):
        # a+b lambda=(a+b*C/2)+(b/2)*sqrt(C^2-4N).
        a,b=F(x[0])+F(x[1])*self.c/2,F(x[1])/2
        if not a:return int(b>0)-int(b<0)
        if not b:return int(a>0)-int(a<0)
        if (a>0)==(b>0):return int(a>0)-int(a<0)
        d=a*a-(self.c*self.c-4*self.n)*b*b
        return (int(a>0)-int(a<0)) if d>0 else (int(b>0)-int(b<0)) if d<0 else 0


def potential_audit(a, depth=80):
    q=Quadratic(a);h,n=q.h,q.n
    lam=(F(0),F(1))
    u=q.scale(q.inv((-n,F(1))),(a-1)*n)
    ca=q.scale(q.sub(u,(F(1),F(0))),F(a,a-1))
    cb=q.scale(q.add(u,(F(a),F(0))),F(1,a-1))
    beta=q.sub(u,q.scale(cb,F(a-1,a)))
    assert q.sign(q.sub(lam,(n,0)))>0
    assert q.sign(q.sub((h*(a+1)-1,0),lam))>0
    assert q.sign(q.sub(u,(1,0)))>0 and q.sign(beta)>0
    def weight(k,sign):
        return q.sub(q.scale(u,k),q.scale(ca if sign==0 else cb,1-F(1,a**k)))
    assert weight(0,0)==weight(0,1)==(0,0) and weight(1,0)==(1,0)
    identities=0
    for k in range(1,depth+2):
        wa,wb=weight(k,0),weight(k,1)
        fa=(a-1)*n*k-((a+1)*h-a)
        fb=(a-1)*n*k-((a+1)*h-1)
        predicted_a=q.add((fa,0),q.add(q.scale(weight(k+1,0),h),q.scale(weight(k-1,1),h-1)))
        predicted_b=q.add((fb,0),q.add(q.scale(weight(k+1,0),h-1),q.scale(weight(k-1,1),h)))
        assert q.mul(lam,wa)==predicted_a and q.mul(lam,wb)==predicted_b
        assert q.sign(q.sub(wa,(k,0)))>=0 and q.sign(q.sub(wb,q.scale(beta,k)))>=0
        assert q.sign(q.sub(q.scale(u,k),wa))>=0 and q.sign(q.sub(q.scale(u,k),wb))>=0
        identities+=2
    levels=principal_levels(a,depth)
    potential=(F(1),F(0));lam_power=(F(1),F(0))
    values=[]
    for level in levels:
        m=level['depth'];d=level['principal_rank_sum']
        actual=(F(0),F(0))
        for k,count in level['plus'].items():actual=q.add(actual,q.scale(weight(k,0),count))
        for k,count in level['minus'].items():actual=q.add(actual,q.scale(weight(k,1),count))
        assert actual==potential
        assert q.sign(q.sub(q.scale(u,d),actual))>=0
        if m<10:values.append(d)
        potential=q.add(q.mul(lam,potential),(n**(m+1),0))
        lam_power=q.mul(lam,lam_power)
    return dict(a=a,depth=depth,principal_ranks_first_ten=values,
                recurrence_checks=depth+1,coefficient_identities=identities,
                quadratic_relation=f'lambda^2-({q.c})lambda+{n}=0',
                positive_weight_slope=[str(x) for x in u],minimum_minus_increment=[str(x) for x in beta])


def enumeration_audit(a,depth):
    expected=principal_levels(a,depth)
    measured=[dict(rank=0,plus=Counter(),minus=Counter(),infinity_plus=Counter(),infinity_minus=Counter(),words=0) for _ in range(depth+1)]
    checks=error_checks=0
    def visit(word,local,inventory):
        nonlocal checks,error_checks
        m=len(word);row=measured[m];d=P9.dim(local['y'])
        row['rank']+=d;row['words']+=1
        for (sign,k),number in local[0].items():row['plus' if sign==0 else 'minus'][k]+=number
        for (sign,k),number in local['infinity'].items():row['infinity_plus' if sign==0 else 'infinity_minus'][k]+=number
        # The raw error's conductor at y is zero, checked on its entire local inventory.
        if inventory is not None:
            error=inventory.copy()
            error[P10.freeze(local,a)]-=1
            assert min(error.values(),default=0)>=0
            _,conductors=P10.inventory_invariants(+error,a)
            assert conductors['y']==0
            for key,mult in (+error).items():
                if key[0]=='point':assert key[1]!='y'
                else:
                    component=P10.thaw(key,a)
                    assert P9.dim(component['y'])==P9.blocks(component['y'],0)
            error_checks+=1
        checks+=1
        if m==depth:return
        for mask in range(1,2**a):
            twisted,_=P9.mask_twist(local,mask,a)
            next_local,_=P9.mc(twisted,a)
            nxt=P10.inventory_step(inventory,mask,a) if inventory is not None and m<3 else None
            visit(word+(mask,),next_local,nxt)
    initial=P9.initial(a)
    visit((),initial,Counter({P10.freeze(initial,a):1}))
    for actual,expect in zip(measured,expected):
        assert actual['rank']==expect['principal_rank_sum'] and actual['words']==expect['words']
        assert actual['plus']==Counter(expect['plus']) and actual['minus']==Counter(expect['minus'])
        assert actual['infinity_minus']==actual['plus'] and actual['infinity_plus']==actual['minus']
    return dict(a=a,depth=depth,enumerated_words=checks,error_unramified_inventory_checks=error_checks,aggregate_level_checks=depth+1)


def improved_budget(p,a,m,j):
    n=2**a-1
    d=principal_levels(a,j)[-1]['principal_rank_sum']
    raw=(a*(2**(a-1)*(a+1)-1)**j-n**j)//(a-1) if a>=2 else j+1
    e=raw-d
    assert e>=0
    root=P11.sqrt_upper(F(1,p))
    alpha=(F(a*d*d,n**j)+1)*root
    v=1+alpha+F(2*a*d*e+e*e,p*n**j)
    u=p*v+F(4*(a+1)*raw*raw,n**j)+(a+1)*n**j
    old=P11.budget(p,a,m,j)
    full=(P11.sqrt_upper(u)+j*old['d0']*P11.sqrt_upper(n**(j-1)))**2
    return dict(d=d,e=e,raw=raw,v=v,u=u,bound=full,previous_bound=old['bound'])


def finite_energy_audit():
    original=json.loads((ROOT/'results/parallel11_full_energy_2026_09_05.json').read_text())
    rows=[]
    for case in original['cases']:
        p,a,m=case['p'],len(case['anchors']),len(case['cell'])
        for item in case['full_bounds']:
            j=item['j'];energy=tuple(F(x) for x in item['energy'])
            bd=improved_budget(p,a,m,j)
            assert P11.upperq(energy,p)<=bd['bound']
            rows.append(dict(p=p,anchors=case['anchors'],j=j,
                             principal_rank_sum=bd['d'],error_rank_sum=bd['e'],
                             upper_bound=str(bd['bound']),smaller_than_previous=bd['bound']<bd['previous_bound']))
    p=10009
    assert all(p%d for d in range(2,isqrt(p)+1))
    m=(p-5)//4
    bd=improved_budget(p,2,m,1)
    assert F(m*9+p-m,3)<=bd['bound']<bd['previous_bound']<2*p
    return dict(cases=rows,full_upper_bound_checks=len(rows),
                strict_improvements=sum(r['smaller_than_previous'] for r in rows),
                nonvacuous_example=dict(p=p,a=2,j=1,upper_bound_over_p=str(bd['bound']/p),
                                        previous_bound_over_p=str(bd['previous_bound']/p)))


def main():
    enumerations=[enumeration_audit(a,d) for a,d in [(1,12),(2,6),(3,4),(4,3),(5,2)]]
    potentials=[potential_audit(a,80) for a in range(2,11)]
    energies=finite_energy_audit()
    inputs=['research/parallel13-principal-budget-2026-09-05.md',
            'experiments/parallel13_principal_budget_2026_09_05.py',
            'research/parallel9-all-degrees-2026-09-05.md','research/parallel10-word-aggregate-2026-09-05.md',
            'research/parallel11-full-energy-2026-09-05.md','results/parallel11_full_energy_2026_09_05.json',
            'experiments/parallel9_all_degrees_2026_09_05.py','experiments/parallel10_word_aggregate_2026_09_05.py',
            'experiments/parallel11_full_energy_2026_09_05.py','experiments/parallel2_spectral_transfer_2026_09_04.py',
            'sources/katz-rigid-local-systems.pdf']
    totals=dict(enumerated_words=sum(r['enumerated_words'] for r in enumerations),
                error_unramified_inventory_checks=sum(r['error_unramified_inventory_checks'] for r in enumerations),
                aggregate_level_checks=sum(r['aggregate_level_checks'] for r in enumerations),
                exact_potential_recurrences=sum(r['recurrence_checks'] for r in potentials),
                exact_coefficient_identities=sum(r['coefficient_identities'] for r in potentials),
                full_upper_bound_checks=energies['full_upper_bound_checks'],strict_finite_improvements=energies['strict_improvements'])
    report=dict(status='Exact finite checks passed. The improved logarithmic range uses the written source-dependent proof; the long-depth target is unproved.',
                enumerations=enumerations,potentials=potentials,energies=energies,totals=totals,
                source_pages_visually_rechecked=[49,50,63,64,67,107,109,139],
                arithmetic='Integers, rational upper square roots, and exact quadratic fields; no floating-point acceptance.',
                input_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs})
    dest=ROOT/'results/parallel13_principal_budget_2026_09_05.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),totals=totals)),flush=True)


if __name__=='__main__':main()
