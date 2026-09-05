#!/usr/bin/env python3
"""Exact degree/exponent audit, finite replays, and Hasse sensitivity bounds."""
from fractions import Fraction
from hashlib import sha256
import json
from math import prod,isqrt,sqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]


def enc(x):return [x.numerator,x.denominator]
def decode(p):return [Fraction(*x) for x in p]
def evaluate(p,q):return sum((v*q**i for i,v in enumerate(p)),Fraction())


def subtract(a,b):
    p=[(a[j] if j<len(a) else Fraction())-(b[j] if j<len(b) else Fraction()) for j in range(max(len(a),len(b)))]
    while len(p)>1 and not p[-1]:p.pop()
    return p


def falling(n,k):
    if k>n:return 0
    return prod(range(n-k+1,n+1))


def inclusion(rows,q,n):
    probability=Fraction(1);total=Fraction()
    for k,p in enumerate(rows):
        if k>n:break
        if k:probability*=Fraction(n-k+1,q-k+1)
        total+=evaluate(p,q)*probability
    return total


def terms(polys,statistic_power=None):
    records=[]
    extra=Fraction() if statistic_power is None else 1+Fraction(statistic_power,2)
    for k,p in enumerate(polys):
        if not any(p):continue
        degree=len(p)-1
        records.append({'union_k':k,'q_degree':degree,'leading_coefficient':enc(p[-1]),
                        'critical_exponent':enc(Fraction(degree)-Fraction(3*k,4)+extra)})
    return records


def main():
    path=HERE/'q_polynomials.json';data=json.loads(path.read_text())
    B=[decode(p) for p in data['baseline_by_union_k']]
    C=[[decode(p) for p in row] for row in data['statistic_coefficients_by_union_k']]
    V=[decode(p) for p in data['second_moment_by_union_k']]
    atlas={'baseline':terms(B),'second_moment':terms(V),'statistics':{}}
    for i,name in enumerate(data['statistic_names']):atlas['statistics'][name]=terms([row[i] for row in C],int(name[-1]))
    d3=[subtract(row[2],row[3]) for row in C];d4=[subtract(row[4],row[5]) for row in C]
    atlas['fixed_ordinary_statistics']={'A2':terms([row[1] for row in C],2),'A3':terms(d3,3),'A4':terms(d4,4)}
    assert max(Fraction(*r['critical_exponent']) for r in atlas['baseline'])==Fraction(13,4)
    assert [r for r in atlas['baseline'] if Fraction(*r['critical_exponent'])==Fraction(13,4)]==[
        {'union_k':9,'q_degree':10,'leading_coefficient':[1,216],'critical_exponent':[13,4]}]
    assert max(Fraction(*r['critical_exponent']) for r in atlas['baseline'] if r['union_k']!=9)==Fraction(5,2)
    assert max(Fraction(*r['critical_exponent']) for r in atlas['second_moment'])==Fraction(5,2)
    assert V[6][-1]==Fraction(1,720) and len(V[6])==8
    assert max(Fraction(*r['critical_exponent']) for r in atlas['second_moment'] if r['union_k']!=6)==Fraction(7,4)
    assert max(Fraction(*r['critical_exponent']) for rows in atlas['statistics'].values() for r in rows)==Fraction(3,2)
    assert C[7][1]==[Fraction(),Fraction(134,3),Fraction(-140,3),Fraction(2)]
    assert d4[7]==[Fraction(),Fraction(1,3),Fraction(-1,3)]
    assert max(Fraction(*r['critical_exponent']) for rows in atlas['fixed_ordinary_statistics'].values() for r in rows)==Fraction(-1,4)
    finite=[]
    previous=json.loads((LAB/'round5/third_moment/results.json').read_text())
    summaries={}
    for target in previous['scale_cases']+previous['twins']:
        q,n=target['q'],target['n'];name=target.get('graph',f'Paley{q}')
        summary_path=LAB/'round5/trace_invariants'/f'seven_statistics_{name}.json'
        summary=json.loads(summary_path.read_text());summaries[str(summary_path.relative_to(LAB))]=sha256(summary_path.read_bytes()).hexdigest()
        s=list(summary['statistics'].values())
        assert list(summary['statistics'])==['mono_theta1','mono_theta2','mono_theta3','mixed_theta3','mono_theta4','mixed_theta4','total_theta6']
        coefficients=[inclusion([row[i] for row in C],q,n) for i in range(7)]
        exact=inclusion(B,q,n)+sum((a*b for a,b in zip(coefficients,s)),Fraction())
        second=inclusion(V,q,n)
        mean=-Fraction(falling(n,6),48*(q-2)*(q-4));variance=second-mean**2
        central=exact-3*mean*second+2*mean**3
        assert exact==Fraction(*target['third_moment'])
        assert second==Fraction(*target['second_moment']) and mean==Fraction(*target['mean'])
        assert variance==Fraction(*target['variance']) and central==Fraction(*target['central_third_moment'])
        proxy=Fraction(q*falling(n,9),216)
        row={'graph':name,'q':q,'n':n,'exact_first_three_moments_match':True,'raw_third_moment':enc(exact),
             'mean':enc(mean),'variance':enc(variance),'central_third_moment':enc(central),
             'finite_n_raw_proxy':enc(proxy),'variance_proxy':enc(Fraction(q*falling(n,6),720)),
             'exact_standardized_squared':enc(central**2/variance**3),
             'finite_n_standardized_proxy':float(proxy)/float(Fraction(q*falling(n,6),720))**1.5,
             'asymptotic_standardized_proxy':40*sqrt(5/q),
             'actual_standardized_float':target['standardized_third_float']}
        if proxy:row['raw_proxy_relative_error']=enc((exact-proxy)/proxy)
        if q>=101 and name.startswith('Paley'):
            H=isqrt(4*q);A0=Fraction(q-5,4);B0=Fraction(3*(q-1),4);N=q-2;m2=(q-3)*(q+1)
            # Triangle inequality with exact nonnegative even-moment ranges.
            ranges=[2*A0*H,A0*H**2,2*A0*H**3,2*B0*H**3,A0*H**4,B0*H**4,N*H**6]
            general=sum((abs(c)*bound for c,bound in zip(coefficients,ranges)),Fraction())
            # Sharpen union6 jointly using concavity of x^3-(15q-85)x^2.
            a=15*q-85;L=H*H;mu=Fraction(m2,N)
            assert 6*L-2*a<=0 and 0<=mu<=L
            f=lambda x:x**3-a*x*x
            range6=N*f(mu)-m2*(L*L-a*L)
            probability6=Fraction(falling(n,6),falling(q,6))
            bound6=probability6*Fraction(q*(q-1),720)*range6
            tailcoeff=[coefficients[i]-probability6*evaluate(C[6][i],q) for i in range(7)]
            refined=bound6+sum((abs(c)*bound for c,bound in zip(tailcoeff,ranges)),Fraction())
            # With all ordinary powers fixed, only one mono polynomial remains.
            a2=coefficients[1];a3=coefficients[2]-coefficients[3];a4=coefficients[4]-coefficients[5]
            mono_support=[t for t in range(-H,H+1) if t%8==(15-q)%8]
            vals=[a2*t*t+a3*t**3+a4*t**4 for t in mono_support]
            ordinary=A0*(max(vals)-min(vals))
            row['finite_hasse_sensitivity']={'all_arithmetic_triangle_bound':enc(general),
                 'all_arithmetic_concavity_bound':enc(refined),'union6_concavity_component':enc(bound6),
                 'fixed_ordinary_moments_range_bound':enc(ordinary),
                 'fixed_ordinary_standardized_bound_squared':enc(ordinary**2/variance**3),
                 'monochromatic_theta_support_points_evaluated':len(mono_support),
                 'range_minimum':enc(min(vals)),'range_maximum':enc(max(vals)),
                 'scope':'Differences between any two nonnegative Hasse-supported inventories with the stated shared constraints; no realizability claim.'}
        finite.append(row)
    lp=json.loads((LAB/'round6/trace_envelopes/results.json').read_text())
    envelope_checks=[]
    for target in lp['cases']:
        row=next(r for r in finite if(r['q'],r['n'])==(target['q'],target['n']))
        key='fixed_ordinary_moments_range_bound' if target['model']=='ordinary_powers_through_six' else 'all_arithmetic_concavity_bound'
        bound=Fraction(*row['finite_hasse_sensitivity'][key])
        witnesses=[Fraction(*v) for v in target['exact_relaxed_witnesses']]
        assert abs(witnesses[1]-witnesses[0])<=bound
        envelope_checks.append({'q':target['q'],'n':target['n'],'model':target['model'],'witness_difference':enc(witnesses[1]-witnesses[0]),
                                'independent_range_bound':enc(bound),'bound_to_witness_ratio':None if witnesses[0]==witnesses[1] else float(bound/(witnesses[1]-witnesses[0]))})
    out={'scope':'Exact polynomial exclusion and Hasse-supported global-moment sensitivity, not a pointwise bound.',
         'critical_scaling':'n=alpha*q^(1/4)+O(1), fixed alpha>0, prime q=1 mod4',
         'leading_laws':{'mean':'-alpha^6*q^(-1/2)/48','variance':'alpha^6*q^(5/2)/720',
                         'raw_and_central_third':'alpha^9*q^(13/4)/216','standardized_third':'40*sqrt(5)*q^(-1/2)'},
         'finite_n_proxy_errors':{'raw_and_central_third':'q*(n)_9/216 + O_alpha(q^(5/2))',
                                  'variance':'q*(n)_6/720 + O_alpha(q^(7/4))',
                                  'skewness':'40*sqrt(5)/sqrt(q) * (n)_9/(n)_6^(3/2) * (1+O_alpha(q^(-3/4)))'},
         'uniform_sensitivity_limsup':{'arbitrary_arithmetic':'limsup q^(-3/2)*|DeltaM3| <= alpha^6/24',
                                      'fixed_ordinary_moments':'limsup q^(1/4)*|DeltaM3| <= 3*alpha^7/4',
                                      'arbitrary_arithmetic_standardized':'limsup q^(9/4)*|DeltaSkew| <= 360*sqrt(5)*alpha^(-3)',
                                      'fixed_ordinary_standardized':'limsup q^4*|DeltaSkew| <= 6480*sqrt(5)*alpha^(-2)'},
         'degree_and_exponent_atlas':atlas,'finite_checks':finite,'envelope_witness_checks':envelope_checks,
         'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'q_polynomials_sha256':sha256(path.read_bytes()).hexdigest(),
         'input_sha256':{**summaries,str((LAB/'round5/third_moment/results.json').relative_to(LAB)):sha256((LAB/'round5/third_moment/results.json').read_bytes()).hexdigest(),
                         'round6/trace_envelopes/results.json':sha256((LAB/'round6/trace_envelopes/results.json').read_bytes()).hexdigest()}}
    (HERE/'scaling_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'finite_moment_cases':len(finite),'envelope_witness_checks':len(envelope_checks),'leading_laws':out['leading_laws']},indent=2))


if __name__=='__main__':main()
