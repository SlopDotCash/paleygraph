#!/usr/bin/env python3
"""Independent rational reconstruction of the inventory-free enclosures."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb,isqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]


def digest(path):return sha256(path.read_bytes()).hexdigest()


def enc(x):return [x.numerator,x.denominator]


def at(poly,q):return sum(F(*c)*q**j for j,c in enumerate(poly))


def inclusion(polys,q,n):
    return sum(at(poly,q)*F(comb(n,k),comb(q,k)) for k,poly in enumerate(polys) if k<=n)


def rebuild(p,q,n):
    assert q%4==1 and q>=29 and all(q%d for d in range(2,isqrt(q)+1))
    H=isqrt(4*q);L=H*H;N=q-2;M2=(q-3)*(q+1)
    A0=F(q-5,4);B0=F(3*(q-1),4)
    mean=-q*comb((q-1)//2,3)*F(comb(n,6),comb(q,6))
    second=inclusion(p['second_moment_by_union_k'],q,n)
    variance=second-mean**2;assert variance>0
    baseline=inclusion(p['baseline_by_union_k'],q,n)
    polys=p['statistic_coefficients_by_union_k']
    coefficients=[inclusion([row[i] for row in polys],q,n) for i in range(7)]
    weight=F(comb(n,6),comb(q,6))
    a=15*q-85
    # f(x)=x^3-a*x^2 is concave on [0,L]. Its chord is a
    # lower bound; Jensen gives the upper bound using the exact M2.
    assert 6*L<=2*a and 0<=F(M2,N)<=L
    lower6=weight*F(q*(q-1),720)*M2*(L**2-a*L)
    upper6=weight*F(q*(q-1),720)*N*(F(M2,N)**3-a*F(M2,N)**2)
    assert lower6<=upper6
    expected6=[F(0)]*7
    expected6[4]=expected6[5]=-F(q*(q-1)*a,720)
    expected6[6]=F(q*(q-1),720)
    assert [at(poly,q) for poly in polys[6]]==expected6
    remainder=[c-weight*t for c,t in zip(coefficients,expected6)]
    ranges=[(-A0*H,A0*H),(F(0),A0*H**2),(-A0*H**3,A0*H**3),(-B0*H**3,B0*H**3),
            (F(0),A0*H**4),(F(0),B0*H**4),(F(0),F(N*H**6))]
    lower=baseline+lower6+sum(min(c*x,c*y) for c,(x,y) in zip(remainder,ranges))
    upper=baseline+upper6+sum(max(c*x,c*y) for c,(x,y) in zip(remainder,ranges))
    offset=-3*mean*second+2*mean**3
    return {'exact_mean':enc(mean),'exact_second_moment':enc(second),'exact_variance':enc(variance),
            'raw_third_interval':[enc(lower),enc(upper)],'central_third_interval':[enc(lower+offset),enc(upper+offset)],
            'exact_raw_width':enc(upper-lower),'standardized_width_squared':enc((upper-lower)**2/variance**3),
            'union6_arithmetic_interval':[enc(lower6),enc(upper6)],
            'tail_statistic_coefficients':[enc(x) for x in remainder],
            'statistic_ranges':[[enc(x),enc(y)] for x,y in ranges]}


def main():
    folder=LAB/'round6/critical_third'
    p=json.loads((folder/'q_polynomials.json').read_text())
    result=json.loads((folder/'inventory_free_results.json').read_text())
    assert digest(folder/'q_polynomials.json')==result['q_polynomials_sha256']
    assert digest(folder/'inventory_free.py')==result['source_sha256']
    assert digest(folder/'derive_scaling.py')==result['helper_sha256']
    cases=[]
    for row in result['cases']:
        expected=rebuild(p,row['q'],row['n'])
        for key,value in expected.items():assert row[key]==value,(row['q'],key)
        assert row['character_array_entries_computed']==row['trace_inventory_entries_required']==0
        cases.append({'q':row['q'],'n':row['n'],'all_exact_enclosure_fields_match':True,
                      'prime_verified_independently':True,'standardized_width_squared':expected['standardized_width_squared']})
    old=json.loads((LAB/'round5/third_moment/results.json').read_text())['scale_cases']
    controls=[]
    for row in old:
        if row['q']<29:continue
        expected=rebuild(p,row['q'],row['n'])
        low,high=map(lambda x:F(*x),expected['raw_third_interval'])
        assert low<=F(*row['third_moment'])<=high
        assert expected['exact_variance']==row['variance']
        controls.append({'q':row['q'],'n':row['n'],'contains_exact_saved_moment':True})
    output={'status':'passed','scope':'Independent rational reconstruction and concavity/Jensen audit; no production helper import or character array.',
            'large_cases':cases,'known_moment_controls':controls,
            'source_sha256':{str(path.relative_to(LAB)):digest(path) for path in
                (Path(__file__),folder/'inventory_free.py',folder/'derive_scaling.py',folder/'q_polynomials.json',
                 folder/'inventory_free_results.json',LAB/'round5/third_moment/results.json')}}
    (HERE/'enclosure_results.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Independent inventory-free audit passed: 2 large enclosures and 6 exact-moment controls.')


if __name__=='__main__':main()
