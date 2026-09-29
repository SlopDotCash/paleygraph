#!/usr/bin/env python3
"""Separate SymPy-resultant and polynomial-remainder review of carry evidence."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent
X=sp.Symbol('X')


def product(a,b,N):
    poly=sp.Poly.from_list(a[::-1],X,domain=sp.ZZ)*sp.Poly.from_list(b[::-1],X,domain=sp.ZZ)
    rem=poly.rem(sp.Poly(X**N+1,X,domain=sp.ZZ))
    return [int(rem.nth(j)) for j in range(N)]


def result_norm(a,N):return int(sp.resultant(sp.Poly(X**N+1,X,domain=sp.ZZ),sp.Poly.from_list(a[::-1],X,domain=sp.ZZ)))


def main():
    source=HERE/'norm_carry.json';produced=json.loads(source.read_text());reviews=[]
    for case in produced['cases']:
        p,n,N,g=(case[k] for k in ('p','n','N','g'));assert sp.isprime(p) and n==2*N and N>=4 and N&(N-1)==0
        assert pow(g,N,p)==p-1 and pow(g,n,p)==1 and case['p_ge_n_to_fourth']==(p>=n**4)
        unit=case['unit_coefficients'];assert unit==[1,1]+[0]*(N-3)+[-1] and result_norm(unit,N)==1
        u=pow(g,-1,p);m=sum(v*pow(u,j,p) for j,v in enumerate(unit))%p;assert m==case['multiplier']!=0
        cache={};ratios=[];carry_hist=Counter()
        for index,r in enumerate(case['records']):
            a=r['a'];b=r['next_a'];assert b==a*m%p and 1<=a<p
            expected=[]
            for j in range(N):
                value=a*pow(g,j,p)%p;expected.append(value if value<=p//2 else value-p)
            assert r['F']==expected
            raw=product(unit,r['F'],N);assert raw==r['raw_unit_image']
            centered=[v%p if v%p<=p//2 else v%p-p for v in raw];assert centered==r['next_F']
            assert all(centered[j]%p==b*pow(g,j,p)%p for j in range(N))
            carry=[(v-y)//p for v,y in zip(raw,centered)];assert carry==r['carry'] and all(x in (-1,0,1) for x in carry)
            for label,poly in [(a,r['F']),(b,centered)]:
                if label not in cache:
                    value=result_norm(poly,N);assert value>0 and value%p**(N-1)==0;cache[label]=value//p**(N-1)
            assert cache[a]==r['norm_defect'] and cache[b]==r['next_norm_defect']
            E=sum(x*x for x in r['F']);next_E=sum(x*x for x in centered);raw_E=sum(x*x for x in raw);dot=sum(x*y for x,y in zip(raw,carry));carry_E=sum(x*x for x in carry)
            assert [E,next_E,raw_E,dot,carry_E]==[r[k] for k in ('coefficient_energy','next_energy','raw_energy','raw_carry_pairing','carry_energy')]
            assert next_E==raw_E-2*p*dot+p*p*carry_E
            second=product(unit,centered,N);second_center=[v%p if v%p<=p//2 else v%p-p for v in second]
            second_carry=[(v-y)//p for v,y in zip(second,second_center)];two=[x+y for x,y in zip(product(unit,carry,N),second_carry)]
            assert second_center==r['second_F'] and second_carry==r['second_carry'] and two==r['two_step_carry'] and r['second_a']==b*m%p
            assert product(unit,raw,N)==[v+p*c for v,c in zip(second_center,two)]
            flag=next_E<E and cache[b]>cache[a];assert flag==r['radius_down_norm_up']
            if flag:ratios.append(r)
            carry_hist[carry_E]+=1
            if not case['complete_nonzero_scalar_census']:
                assert a==(case['start'] if index==0 else case['records'][index-1]['next_a'])
        if case['complete_nonzero_scalar_census']:assert [r['a'] for r in case['records']]==list(range(1,p))
        else:assert len(case['records'])==case['steps']
        assert len(cache)==case['distinct_norms_evaluated'] and len(ratios)==case['radius_down_norm_up_steps'] and sorted(carry_hist.items())==[tuple(x) for x in case['carry_support_histogram']]
        if ratios:
            winner=max(ratios,key=lambda r:Fraction(r['next_norm_defect'],r['norm_defect']));assert case['largest_observed_norm_ratio_during_radius_decrease']['a']==winner['a']
        else:assert case['largest_observed_norm_ratio_during_radius_decrease'] is None
        row={'case':case['name'],'p':p,'n':n,'steps_checked':len(case['records']),'independent_resultants':len(cache)+1,'radius_down_norm_up_steps':len(ratios),'quartic_window':case['p_ge_n_to_fourth']};reviews.append(row);print(json.dumps(row),flush=True)
    out={'status':'passed','scope':'Independent prime/order checks, polynomial products/remainders, integer resultants, all carries, two-step cocycles and energy/norm comparisons. Complete only for the stated small scalar census; other walks are bounded and no global arithmetic estimate is inferred.',
         'cases':reviews,'source_sha256':{'carry_review.py':sha256((HERE/'carry_review.py').read_bytes()).hexdigest()},'input_sha256':{source.name:sha256(source.read_bytes()).hexdigest()}}
    (HERE/'carry_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
