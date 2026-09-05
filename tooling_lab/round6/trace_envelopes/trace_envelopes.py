#!/usr/bin/env python3
"""Exact certificates for a finite nonnegative trace-histogram relaxation.

HiGHS proposes primal/dual vectors. Only exact rational identities and checked
dual inequalities establish the exported bounds. This is not a graph
realizability test, nor a bound on a maximum over column sets.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import linprog
import sympy as sp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]


def encode(x):
    x = Fraction(x)
    return [x.numerator, x.denominator]


def evaluate(poly, x):
    value = Fraction()
    for c in reversed(poly):
        value = value*x+c
    return value


def hash_file(path):
    return sha256(path.read_bytes()).hexdigest()


def prime(q):
    return q >= 2 and all(q % d for d in range(2, isqrt(q)+1))


def problem(moment, folded, model):
    """Use scaled counts: w_m=N_m/q, w_u=N_mixed/(3q)."""
    q = moment['q']
    assert q == folded['q'] and q % 4 == 1 and prime(q)
    assert model in ('conference_only', 'ordinary_powers_through_six')
    assert moment['degree'] == 6
    bound, root_scale = isqrt(4*q), isqrt(q)+1
    points = []
    for mono in (True, False):
        residue = ((15-q) if mono else (q-7)) % 8
        points.extend((mono, theta) for theta in range(-bound, bound+1) if theta % 8 == residue)
    # At small q, Hasse permits formal values with negative sign-type counts.
    # Remove them using the actual bulk Walsh formulas, without importing the
    # moment compiler. The three boundary columns are accounted for in q-3.
    def nonnegative(mono, theta):
        e, tau = ((1,1,1),theta) if mono else ((1,1,-1),-theta)
        x,y,z = e
        for a in (-1,1):
            for b in (-1,1):
                for c in (-1,1):
                    v=q-3-a*(x+y)-b*(x+z)-c*(y+z)-a*b*(1+y*z)-a*c*(1+x*z)-b*c*(1+x*y)+a*b*c*tau
                    if v < 0 or v % 8:
                        return False
        return True
    points = [point for point in points if nonnegative(*point)]
    polys = {p['monochromatic']:[Fraction(*c) for c in p['polynomial_in_theta']]
             for p in moment['family_plans']}
    assert all(p['polynomial_scope']=='all theta by degree bound' for p in moment['family_plans'])
    objective = [q*q*(q-1)*(1 if mono else 3)*evaluate(polys[mono],theta) for mono,theta in points]
    constant = sum(Fraction(*moment[k]) for k in ('all_equal_rows_contribution','exactly_two_equal_rows_contribution'))
    rows = [[Fraction(int(mono)) for mono,theta in points],
            [Fraction(int(not mono)) for mono,theta in points],
            [Fraction(theta*(1 if mono else -1),root_scale) for mono,theta in points],
            [Fraction(theta**2*(1 if mono else 3),root_scale**2) for mono,theta in points]]
    rhs = [Fraction(q-5,4*q),Fraction(q-1,4*q),Fraction(2,q*root_scale),
           Fraction((q-3)*(q+1),q*root_scale**2)]
    labels = ['monochromatic count/q','one mixed class count/q',
              '(mono theta sum - mixed theta sum/3)/(q R)',
              'total theta squared sum/(q R^2)']
    if model == 'ordinary_powers_through_six':
        # Total zeroth and second moments are already present. Powers1,3,4,
        # 5,6 add the ordinary scalar inventory without its residue labels.
        for j in (1,3,4,5,6):
            value=sum(row['count']*row['theta']**j for row in folded['records'])
            rows.append([Fraction(theta**j*(1 if mono else 3),root_scale**j) for mono,theta in points])
            rhs.append(Fraction(value,q*root_scale**j))
            labels.append(f'total theta power {j}/(q R^{j})')
    actual = {(row['monochromatic'],row['theta']):Fraction(row['count'],q*(1 if row['monochromatic'] else 3))
              for row in folded['records']}
    assert set(actual) <= set(points)
    actual_vector = [actual.get(point,Fraction()) for point in points]
    assert all(sum(a*w for a,w in zip(row,actual_vector))==b for row,b in zip(rows,rhs))
    assert constant+sum(c*w for c,w in zip(objective,actual_vector))==Fraction(*moment['third_moment'])
    return {'q':q,'n':moment['n'],'model':model,'points':points,'rows':rows,'rhs':rhs,
            'objective':objective,'constant':constant,'labels':labels,'root_scale':root_scale,
            'actual_vector':actual_vector}


def quotient_objective(data):
    """Remove an exactly known constraint component before numerical search.

    Otherwise the residual of interest can be 10^18 times smaller than the
    objective entries. No precision setting repairs that conceptual loss.
    The removed component is added back to the exact dual certificate.
    """
    A=sp.Matrix(data['rows'])
    _,columns=A.rref()
    row_indices=list(range(len(data['rows'])))
    if len(columns)<len(row_indices):
        _,row_indices=A.T.rref()
    B=A.extract(row_indices,columns)
    target=sp.Matrix([data['objective'][i] for i in columns])
    values=B.T.inv()*target
    multiplier=[Fraction() for _ in data['rows']]
    for index,value in zip(row_indices,values):
        multiplier[index]=Fraction(value)
    residual=[c-sum(y*row[i] for y,row in zip(multiplier,data['rows']))
              for i,c in enumerate(data['objective'])]
    return multiplier,residual


def exact_primal(rows,rhs,solution):
    indices=[i for i,x in enumerate(solution) if x > 1e-11]
    if not indices:
        return None
    A=sp.Matrix([[sp.Rational(row[i].numerator,row[i].denominator) for i in indices] for row in rows])
    b=sp.Matrix([sp.Rational(x.numerator,x.denominator) for x in rhs])
    try:
        vector,parameters=A.gauss_jordan_solve(b)
    except ValueError:
        return None
    if parameters.rows:
        return None
    values=[Fraction(x) for x in vector]
    if any(x < 0 for x in values):
        return None
    if not all(sum(row[i]*v for i,v in zip(indices,values))==target for row,target in zip(rows,rhs)):
        return None
    return list(zip(indices,values))


def certify_upper(data, sign):
    rows,rhs,points=data['rows'],data['rhs'],data['points']
    c=[sign*x for x in data['objective']]
    base,residual=data['quotient']
    magnitude=max(abs(x) for x in residual)
    if magnitude:
        result=linprog(-np.array([float(sign*x/magnitude) for x in residual]),
                       A_eq=np.array([[float(x) for x in row] for row in rows]),
                       b_eq=np.array([float(x) for x in rhs]),bounds=(0,None),method='highs')
        assert result.success, result.message
        dual=[sign*b+Fraction(float(-x)).limit_denominator(10**10)*magnitude
              for b,x in zip(base,result.eqlin.marginals)]
        primal=exact_primal(rows,rhs,result.x)
        message=result.message
    else:
        dual=[sign*b for b in base]
        primal=[(i,x) for i,x in enumerate(data['actual_vector']) if x]
        message='Objective lies exactly in the constraint span; no numerical optimization needed.'
    # Repair the proposed dual by the two constant rows. Each point belongs
    # to exactly one class, so adding its largest deficit is always sound.
    repairs=[]
    for row_index,mono in enumerate((True,False)):
        deficit=max([Fraction()]+[ci-sum(y*row[i] for y,row in zip(dual,rows))
                                  for i,((m,theta),ci) in enumerate(zip(points,c)) if m==mono])
        dual[row_index]+=deficit
        repairs.append(deficit)
    slacks=[sum(y*row[i] for y,row in zip(dual,rows))-ci for i,ci in enumerate(c)]
    assert all(x>=0 for x in slacks)
    upper=sum(y*b for y,b in zip(dual,rhs))
    value=None
    if primal is not None:
        value=sum(c[i]*w for i,w in primal)
        assert value<=upper
    return {'side':'upper' if sign==1 else 'lower','sign':sign,
            'dual_multipliers':[encode(x) for x in dual],
            'class_constant_repairs':[encode(x) for x in repairs],
            'minimum_dual_slack':encode(min(slacks)),
            'signed_variable_upper_bound':encode(upper),
            'primal':None if primal is None else [{'point_index':i,'weight':encode(w)} for i,w in primal],
            'signed_primal_objective':None if value is None else encode(value),
            'exact_primal_dual_gap':None if value is None else encode(upper-value),
            'constraint_component_multipliers':[encode(sign*x) for x in base],
            'maximum_absolute_quotient_objective':encode(magnitude),
            'pointwise_dual_inequalities_checked':len(points),
            'highs_message':message}


def verify(data,certificate):
    sign=certificate['sign']
    dual=[Fraction(*x) for x in certificate['dual_multipliers']]
    assert len(dual)==len(data['rows'])
    for i,c in enumerate(data['objective']):
        assert sum(y*row[i] for y,row in zip(dual,data['rows']))>=sign*c
    upper=sum(y*b for y,b in zip(dual,data['rhs']))
    assert upper==Fraction(*certificate['signed_variable_upper_bound'])
    if certificate['primal'] is not None:
        weights={row['point_index']:Fraction(*row['weight']) for row in certificate['primal']}
        assert len(weights)==len(certificate['primal']) and all(0<=i<len(data['points']) and w>=0 for i,w in weights.items())
        assert all(sum(row[i]*w for i,w in weights.items())==b for row,b in zip(data['rows'],data['rhs']))
        value=sum(sign*data['objective'][i]*w for i,w in weights.items())
        assert value==Fraction(*certificate['signed_primal_objective'])<=upper
    return upper


def main():
    start=time.perf_counter()
    moment_file=LAB/'round5/trace_invariants/results.json'
    original_file=LAB/'round5/third_moment/results.json'
    moments=json.loads(moment_file.read_text())['comparisons']
    originals=json.loads(original_file.read_text())['scale_cases']
    output=[]
    for q,n in ((101,8),(1297,8),(65537,16),(1000033,31)):
        moment=next(x for x in moments if (x['q'],x['n'])==(q,n))
        original=next(x for x in originals if (x['q'],x['n'])==(q,n))
        inventory_file=LAB/'round5/trace_invariants'/f'inventory_Paley{q}.json'
        folded=json.loads(inventory_file.read_text())
        mean=Fraction(*original['mean']);second=Fraction(*original['second_moment'])
        variance=Fraction(*original['variance'])
        center=-3*mean*second+2*mean**3
        standardizer=float(variance)**1.5
        for model in ('conference_only','ordinary_powers_through_six'):
            data=problem(moment,folded,model)
            data['quotient']=quotient_objective(data)
            low,high=[certify_upper(data,sign) for sign in (-1,1)]
            lo=data['constant']-verify(data,low)
            hi=data['constant']+verify(data,high)
            exact=Fraction(*moment['third_moment'])
            assert lo<=exact<=hi
            witnesses=[]
            for cert in (low,high):
                if cert['signed_primal_objective'] is not None:
                    witnesses.append(data['constant']+cert['sign']*Fraction(*cert['signed_primal_objective']))
            row={'q':q,'n':n,'model':model,'support_points':[list(point) for point in data['points']],
                 'constraint_labels':data['labels'],'constraint_rhs':[encode(x) for x in data['rhs']],
                 'root_scale':data['root_scale'],'lower':encode(lo),'upper':encode(hi),
                 'actual':encode(exact),'lower_standardized_float':float(lo+center)/standardizer,
                 'upper_standardized_float':float(hi+center)/standardizer,
                 'actual_standardized_float':original['standardized_third_float'],
                 'width':encode(hi-lo),'width_relative_to_absolute_raw_third_moment':encode((hi-lo)/abs(exact)),
                 'standardized_width_squared':encode((hi-lo)**2/variance**3),
                 'central_lower':encode(lo+center),'central_upper':encode(hi+center),
                 'variance':encode(variance),'exact_relaxed_witnesses':[encode(x) for x in witnesses],
                 'two_distinct_relaxed_values_certified':len(witnesses)==2 and witnesses[0]!=witnesses[1],
                 'lower_certificate':low,'upper_certificate':high,
                 'inventory':str(inventory_file.relative_to(LAB)),
                 'inventory_sha256':hash_file(inventory_file)}
            output.append(row)
            print(f'{q=} {n=} {model}: standardized [{row["lower_standardized_float"]:.6g}, {row["upper_standardized_float"]:.6g}], actual {row["actual_standardized_float"]:.6g}; exact primal witnesses {len(witnesses)}',flush=True)
    record={'scope':'Exact rational bounds on a nonnegative Hasse-supported histogram relaxation for prime Paley global third moments; no graph-realizability, worst-case or prize claim.',
            'cases':output,'source_sha256':hash_file(Path(__file__)),
            'input_sha256':{str(p.relative_to(LAB)):hash_file(p) for p in (moment_file,original_file)},
            'elapsed_seconds':time.perf_counter()-start}
    (HERE/'results.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    main()
