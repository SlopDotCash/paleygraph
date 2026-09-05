#!/usr/bin/env python3
"""Prime Paley global-moment enclosure from q,n alone; no character array."""
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt,sqrt
from pathlib import Path
import sys
import time
sys.dont_write_bytecode=True
from derive_scaling import decode,evaluate,inclusion,falling,enc

HERE=Path(__file__).resolve().parent


def enclose(q,n,polynomial_data):
    start=time.perf_counter()
    assert isinstance(q,int) and q>=29 and q%4==1 and 6<=n<=q
    checked=0
    for divisor in range(2,isqrt(q)+1):
        assert q%divisor,('not prime',q,divisor)
        checked+=1
    B=[decode(p) for p in polynomial_data['baseline_by_union_k']]
    C=[[decode(p) for p in row] for row in polynomial_data['statistic_coefficients_by_union_k']]
    V=[decode(p) for p in polynomial_data['second_moment_by_union_k']]
    baseline=inclusion(B,q,n)
    coefficients=[inclusion([row[i] for row in C],q,n) for i in range(7)]
    mean=-Fraction(falling(n,6),48*(q-2)*(q-4))
    second=inclusion(V,q,n);variance=second-mean*mean;assert variance>0
    H=isqrt(4*q);L=H*H;N=q-2;m2=(q-3)*(q+1)
    A0=Fraction(q-5,4);B0=Fraction(3*(q-1),4)
    ranges=[(-A0*H,A0*H),(0,A0*H**2),(-A0*H**3,A0*H**3),(-B0*H**3,B0*H**3),
            (0,A0*H**4),(0,B0*H**4),(0,N*H**6)]
    probability=Fraction(falling(n,6),falling(q,6));a=15*q-85
    assert 6*L-2*a<=0
    f=lambda x:x**3-a*x*x
    mu=Fraction(m2,N);assert 0<=mu<=L
    sixth_factor=probability*Fraction(q*(q-1),720)
    low6=sixth_factor*m2*(L*L-a*L)
    high6=sixth_factor*N*f(mu)
    # The free union6 part is exactly M6-(15q-85)M4 over720,
    # with the global ordered-row weight; all other union terms are retained.
    assert evaluate(C[6][6],q)==Fraction(q*(q-1),720)
    assert evaluate(C[6][4],q)==evaluate(C[6][5],q)==-Fraction(q*(q-1)*a,720)
    assert all(evaluate(C[6][i],q)==0 for i in range(4))
    tail=[c-probability*evaluate(C[6][i],q) for i,c in enumerate(coefficients)]
    lower=baseline+low6;upper=baseline+high6
    for coefficient,(lo,hi) in zip(tail,ranges):
        lower+=min(coefficient*lo,coefficient*hi)
        upper+=max(coefficient*lo,coefficient*hi)
    assert lower<=upper
    center=-3*mean*second+2*mean**3
    standardizer=float(variance)**1.5
    proxy=Fraction(q*falling(n,9),216)
    return {'q':q,'n':n,'degree':6,'prime_verified_by_trial_division':True,'trial_divisors_checked':checked,
            'character_array_entries_computed':0,'trace_inventory_entries_required':0,
            'exact_mean':enc(mean),'exact_second_moment':enc(second),'exact_variance':enc(variance),
            'raw_third_interval':[enc(lower),enc(upper)],'central_third_interval':[enc(lower+center),enc(upper+center)],
            'exact_raw_width':enc(upper-lower),'standardized_width_squared':enc((upper-lower)**2/variance**3),
            'standardized_width_float':sqrt(float((upper-lower)**2/variance**3)),
            'standardized_interval_float':[float(lower+center)/standardizer,float(upper+center)/standardizer],
            'finite_n_raw_proxy':enc(proxy),'finite_n_standardized_proxy':float(proxy)/float(Fraction(q*falling(n,6),720))**1.5,
            'limiting_standardized_proxy':40*sqrt(5/q),'hasse_integer_support_bound':H,
            'union6_arithmetic_interval':[enc(low6),enc(high6)],'tail_statistic_coefficients':[enc(x) for x in tail],
            'statistic_ranges':[[enc(Fraction(lo)),enc(Fraction(hi))] for lo,hi in ranges],
            'scope':'Unconditional global-moment enclosure for the stated prime Paley matrix using Hasse and conference identities. It uses no measured ordinary moments and gives no subset extremum or prize bound.',
            'elapsed_seconds':time.perf_counter()-start}


def main():
    path=HERE/'q_polynomials.json';data=json.loads(path.read_text())
    cases=[enclose(q,n,data) for q,n in((6700417,50),(2013265921,211))]
    # Controls reuse known exact scalars, not their inventories or arrays.
    known=json.loads((HERE.parents[1]/'round5/third_moment/results.json').read_text())
    controls=[]
    for row in known['scale_cases']:
        if row['q']<29:continue
        result=enclose(row['q'],row['n'],data)
        lo,hi=map(lambda x:Fraction(*x),result['raw_third_interval'])
        actual=Fraction(*row['third_moment']);assert lo<=actual<=hi
        assert result['exact_variance']==row['variance']
        controls.append({'q':row['q'],'n':row['n'],'contains_known_actual_third_moment':True,'exact_variance_matches':True})
    out={'scope':'Inventory-free prime Paley global-average application; no fixed ordinary moments are assumed.',
         'cases':cases,'known_scalar_controls':controls,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'q_polynomials_sha256':sha256(path.read_bytes()).hexdigest(),
         'helper_sha256':sha256((HERE/'derive_scaling.py').read_bytes()).hexdigest()}
    (HERE/'inventory_free_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{'q':row['q'],'n':row['n'],'standardized_interval':row['standardized_interval_float'],
                       'standardized_width':row['standardized_width_float'],'seconds':row['elapsed_seconds']} for row in cases],indent=2))


if __name__=='__main__':main()
