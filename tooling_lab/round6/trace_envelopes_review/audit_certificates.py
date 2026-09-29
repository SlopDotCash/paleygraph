#!/usr/bin/env python3
"""Independent stdlib-only rational audit; imports no production model or LP.

This reconstructs support/constraints/objective from frozen trace polynomial
data. It calls no optimizer and does not rerun the old convolution backend.
"""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import isqrt,sqrt,isclose
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
EXPECTED_SOURCE='0315f4367914d1cad6ceef2e08647cc706357762df0c00c5b68175aab41ee6ed'
EXPECTED_RESULTS='10f8f77d25261422c9c0077a9cff4fe23843a8d08fe2b7688677ad3b3534689d'


def F(value):return Fraction(*value)
def enc(value):return [value.numerator,value.denominator]
def digest(path):return sha256(path.read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Fraction())


def reconstruct(case,moment,inventory):
    q=case['q'];assert q==moment['q']==inventory['q'] and q%4==1 and q>=13
    assert all(q%d for d in range(2,isqrt(q)+1))
    assert moment['n']==case['n'] and moment['degree']==6
    assert case['model'] in ('conference_only','ordinary_powers_through_six')
    # Independently invert the eight parity moments of the bulk signs.
    support=[]
    for family in (True,False):
        edges=(1,1,1) if family else (1,1,-1)
        boundary=((0,edges[0],edges[1]),(edges[0],0,edges[2]),(edges[1],edges[2],0))
        for theta in range(-isqrt(4*q),isqrt(4*q)+1):
            if theta%8!=((15-q) if family else(q-7))%8:continue
            tau=theta if family else -theta
            moments={0:q-3,7:tau}
            for mask in (1,2,4):
                coordinate=mask.bit_length()-1
                moments[mask]=-sum(row[coordinate] for row in boundary)
            for mask in (3,5,6):
                coordinates=[i for i in range(3) if mask&(1<<i)]
                moments[mask]=-1-sum(row[coordinates[0]]*row[coordinates[1]] for row in boundary)
            nums=[]
            for signs in product((-1,1),repeat=3):
                nums.append(sum(value*__import__('math').prod(signs[i] for i in range(3) if mask&(1<<i))
                                for mask,value in moments.items()))
            if all(v>=0 and v%8==0 for v in nums):support.append((family,theta))
    assert case['support_points']==[list(x) for x in support]
    R=isqrt(q)+1;assert case['root_scale']==R
    # Build unscaled *count* identities first, then substitute N_m=q*w_m,
    # N_mixed=3q*w_u and divide each equality by its stated row scale.
    functions=[lambda m,t:int(m),lambda m,t:int(not m),
               lambda m,t:Fraction(t if m else -t,R),
               lambda m,t:Fraction(t*t*(1 if m else 3),R*R)]
    rhs=[Fraction(q-5,4*q),Fraction(q-1,4*q),Fraction(2,q*R),Fraction((q-3)*(q+1),q*R*R)]
    if case['model']=='ordinary_powers_through_six':
        for degree in (1,3,4,5,6):
            functions.append(lambda m,t,d=degree:Fraction((1 if m else 3)*t**d,R**d))
            actual_power=sum(row['count']*row['theta']**degree for row in inventory['records'])
            rhs.append(Fraction(actual_power,q*R**degree))
    rows=[[fn(m,t) for m,t in support] for fn in functions]
    assert [enc(x) for x in rhs]==case['constraint_rhs']
    assert len(case['constraint_labels'])==len(rows)
    assert len(moment['family_plans'])==2
    polys={row['monochromatic']:[F(v) for v in row['polynomial_in_theta']] for row in moment['family_plans']}
    assert set(polys)=={True,False}
    assert all(row['polynomial_scope']=='all theta by degree bound' for row in moment['family_plans'])
    objective=[q*q*(q-1)*(1 if m else 3)*sum((co*t**degree for degree,co in enumerate(polys[m])),Fraction()) for m,t in support]
    constant=F(moment['all_equal_rows_contribution'])+F(moment['exactly_two_equal_rows_contribution'])
    actual_map={}
    for record in inventory['records']:
        point=(record['monochromatic'],record['theta']);assert point not in actual_map
        assert type(point[0]) is bool and isinstance(record['count'],int) and record['count']>0
        actual_map[point]=Fraction(record['count'],q*(1 if point[0] else 3))
    assert set(actual_map)<=set(support)
    actual=[actual_map.get(point,Fraction()) for point in support]
    assert [dot(row,actual) for row in rows]==rhs
    assert constant+dot(objective,actual)==F(moment['third_moment'])==F(case['actual'])
    return {'q':q,'n':case['n'],'points':support,'rows':rows,'rhs':rhs,'objective':objective,'constant':constant,'actual':actual}


def inspect_certificate(data,certificate):
    sign=certificate['sign'];assert sign in (-1,1)
    assert certificate['side']==('upper' if sign==1 else 'lower')
    A,b,c=data['rows'],data['rhs'],data['objective']
    y=[F(v) for v in certificate['dual_multipliers']];assert len(y)==len(A)
    slacks=[sum((yi*row[j] for yi,row in zip(y,A)),Fraction())-sign*cj for j,cj in enumerate(c)]
    assert min(slacks)>=0
    assert F(certificate['minimum_dual_slack'])==min(slacks)
    assert certificate['pointwise_dual_inequalities_checked']==len(c)
    bound=dot(y,b);assert bound==F(certificate['signed_variable_upper_bound'])
    repairs=[F(v) for v in certificate['class_constant_repairs']]
    assert len(repairs)==2 and min(repairs)>=0
    before=y[:];before[0]-=repairs[0];before[1]-=repairs[1]
    for family,index in ((True,0),(False,1)):
        deficits=[sign*c[j]-sum((yi*row[j] for yi,row in zip(before,A)),Fraction())
                  for j,(m,t) in enumerate(data['points']) if m==family]
        assert repairs[index]==max([Fraction(),*deficits])
    component=[F(v) for v in certificate['constraint_component_multipliers']]
    assert len(component)==len(A)
    quotient=[sign*cj-sum((yi*row[j] for yi,row in zip(component,A)),Fraction()) for j,cj in enumerate(c)]
    assert max(abs(x) for x in quotient)==F(certificate['maximum_absolute_quotient_objective'])
    assert dot(component,data['rhs'])+dot(quotient,data['actual'])==sign*dot(c,data['actual'])
    primal=certificate['primal'];weights=None;value=None
    if primal is not None:
        weights=[Fraction() for _ in c];seen=set()
        for row in primal:
            index=row['point_index'];assert isinstance(index,int) and 0<=index<len(c) and index not in seen
            seen.add(index);weights[index]=F(row['weight']);assert weights[index]>=0
        assert [dot(row,weights) for row in A]==b
        value=sign*dot(c,weights)
        assert value==F(certificate['signed_primal_objective']) and value<=bound
        assert bound-value==F(certificate['exact_primal_dual_gap'])
    else:
        assert certificate['signed_primal_objective'] is None and certificate['exact_primal_dual_gap'] is None
    return {'bound':bound,'primal_weights':weights,'signed_primal_value':value,
            'quotient_to_raw_maximum_ratio':max(abs(x) for x in quotient)/max(abs(x) for x in c)}


def histogram_powers(data,weights):
    return [sum((data['q']*(1 if m else 3)*w*t**j for (m,t),w in zip(data['points'],weights)),Fraction()) for j in range(7)]


def nullspace(matrix):
    A=[[Fraction(x) for x in row] for row in matrix];rank=0;pivots=[];width=len(A[0])
    for col in range(width):
        pivot=next((i for i in range(rank,len(A)) if A[i][col]),None)
        if pivot is None:continue
        A[rank],A[pivot]=A[pivot],A[rank]
        value=A[rank][col];A[rank]=[x/value for x in A[rank]]
        for i in range(len(A)):
            if i!=rank and A[i][col]:
                value=A[i][col];A[i]=[a-value*b for a,b in zip(A[i],A[rank])]
        pivots.append(col);rank+=1
    vectors=[]
    for free in set(range(width))-set(pivots):
        vector=[Fraction() for _ in range(width)];vector[free]=1
        for row,pivot in enumerate(pivots):vector[pivot]=-A[row][free]
        assert all(dot(row,vector)==0 for row in matrix);vectors.append(vector)
    return rank,vectors


def singleton_check(data):
    rank,basis=nullspace(data['rows']);assert len(basis)==1
    direction=basis[0];lower=[];upper=[]
    for w,v in zip(data['actual'],direction):
        if v>0:lower.append(-w/v)
        elif v<0:upper.append(-w/v)
    lo,hi=max(lower),min(upper);assert lo==hi==0
    return {'q':101,'equation_rank':rank,'support_size':len(data['points']),
            'equation_nullity':1,'null_direction':[enc(x) for x in direction],
            'support_points':[list(x) for x in data['points']],
            'actual_scaled_weights':[enc(x) for x in data['actual']],
            'zero_weight_constraints':[{'point':list(point),'direction':enc(v)}
                                       for point,w,v in zip(data['points'],data['actual'],direction) if w==0 and v],
            'feasible_parameter_interval':[enc(lo),enc(hi)],
            'conclusion':'At this finite support and these exact ordinary moments, nonnegativity leaves a unique feasible histogram.'}


def main():
    started=time.perf_counter()
    source=LAB/'round6/trace_envelopes/trace_envelopes.py';path=source.with_name('results.json')
    assert digest(source)==EXPECTED_SOURCE and digest(path)==EXPECTED_RESULTS
    saved=json.loads(path.read_text());assert saved['source_sha256']==EXPECTED_SOURCE
    for name,value in saved['input_sha256'].items():assert digest(LAB/name)==value
    moment_path=LAB/'round5/trace_invariants/results.json'
    original_path=LAB/'round5/third_moment/results.json'
    moments=json.loads(moment_path.read_text())['comparisons'];originals=json.loads(original_path.read_text())['scale_cases']
    reports=[];mutations=0;inequalities=0;equalities=0;witness=None;unique=None
    for case in saved['cases']:
        moment=next(x for x in moments if(x['q'],x['n'])==(case['q'],case['n']))
        original=next(x for x in originals if(x['q'],x['n'])==(case['q'],case['n']))
        inventory_path=LAB/case['inventory'];assert digest(inventory_path)==case['inventory_sha256']
        data=reconstruct(case,moment,json.loads(inventory_path.read_text()))
        low=inspect_certificate(data,case['lower_certificate']);high=inspect_certificate(data,case['upper_certificate'])
        inequalities+=2*len(data['points']);equalities+=2*len(data['rows'])
        lo,hi=data['constant']-low['bound'],data['constant']+high['bound']
        assert lo==F(case['lower']) and hi==F(case['upper'])
        actual=F(case['actual']);assert lo<=actual<=hi
        width=hi-lo;assert width==F(case['width'])
        mu,m2=F(original['mean']),F(original['second_moment']);variance=m2-mu*mu
        assert variance==F(original['variance'])==F(case['variance']) and variance>0
        center=-3*mu*m2+2*mu**3
        assert actual+center==F(original['central_third_moment'])
        assert (actual+center)**2/variance**3==F(original['standardized_third_squared'])
        assert lo+center==F(case['central_lower']) and hi+center==F(case['central_upper'])
        assert width**2/variance**3==F(case['standardized_width_squared'])
        assert width/abs(actual)==F(case['width_relative_to_absolute_raw_third_moment'])
        standardized=float(variance)**1.5
        assert isclose(case['lower_standardized_float'],float(lo+center)/standardized,rel_tol=1e-14,abs_tol=1e-14)
        assert isclose(case['upper_standardized_float'],float(hi+center)/standardized,rel_tol=1e-14,abs_tol=1e-14)
        assert isclose(case['actual_standardized_float'],float(actual+center)/standardized,rel_tol=1e-14,abs_tol=1e-14)
        values=[data['constant']-low['signed_primal_value'],data['constant']+high['signed_primal_value']]
        assert [enc(x) for x in values]==case['exact_relaxed_witnesses']
        distinct=values[0]!=values[1];assert distinct==case['two_distinct_relaxed_values_certified']
        assert lo<=values[0]<=values[1]<=hi
        powers=[histogram_powers(data,x['primal_weights']) for x in(low,high)]
        if case['model']=='ordinary_powers_through_six':
            assert powers[0]==powers[1]==histogram_powers(data,data['actual'])
            assert distinct==(case['q']!=101)
            if case['q']==101:unique=singleton_check(data)
            if case['q']==1297:
                witness={'q':1297,'n':8,'model':case['model'],
                         'count_convention':'Entry is N_monochromatic or N_mixed/3, NOT optimizer-scaled w. All are nonnegative rational counts.',
                         'histograms':[{'side':side,'records':[{'monochromatic':point[0],'theta':point[1],'count_per_edge_class':enc(data['q']*w)}
                                      for point,w in zip(data['points'],item['primal_weights']) if w]}
                                       for side,item in(('lower',low),('upper',high))],
                         'shared_ordinary_theta_powers_0_through_6':[enc(x) for x in powers[0]],
                         'third_moments':[enc(x) for x in values],'exact_difference':enc(values[1]-values[0]),
                         'scope':'Two feasible fractional histograms with identical ordinary moments; no graph-realizability claim.'}
        # These corruptions target inequalities, equalities, indexing, and metadata.
        for kind in ('negative_dual','negative_weight','positive_weight_perturbation','duplicate_index','wrong_bound','wrong_slack','wrong_quotient'):
            damaged=deepcopy(case['upper_certificate'])
            if kind=='negative_dual':damaged['dual_multipliers'][0]=[-10**150,1]
            elif kind=='negative_weight':damaged['primal'][0]['weight']=[-1,1]
            elif kind=='positive_weight_perturbation':damaged['primal'][0]['weight']=enc(F(damaged['primal'][0]['weight'])+1)
            elif kind=='duplicate_index':damaged['primal'].append(damaged['primal'][0])
            elif kind=='wrong_bound':damaged['signed_variable_upper_bound']=enc(F(damaged['signed_variable_upper_bound'])+1)
            elif kind=='wrong_slack':damaged['minimum_dual_slack']=enc(F(damaged['minimum_dual_slack'])+1)
            else:damaged['maximum_absolute_quotient_objective']=enc(F(damaged['maximum_absolute_quotient_objective'])+1)
            try:inspect_certificate(data,damaged)
            except AssertionError:mutations+=1
            else:raise AssertionError('corrupt certificate accepted: '+kind)
        reports.append({'q':case['q'],'n':case['n'],'model':case['model'],'support_points':len(data['points']),
                        'constraints':len(data['rows']),'actual_inventory_feasible':True,
                        'lower':enc(lo),'upper':enc(hi),'exact_width':enc(width),
                        'standardized_width_squared':enc(width**2/variance**3),
                        'standardized_width_float_from_exact_width':sqrt(float(width**2/variance**3)),
                        'exact_primal_witness_width':enc(values[1]-values[0]),
                        'primal_witness_width_float':float(values[1]-values[0]),
                        'two_distinct_relaxed_values_certified':distinct,
                        'quotient_to_raw_maximum_ratio':enc(high['quotient_to_raw_maximum_ratio'])})
        print(case['q'],case['model'],'passed',flush=True)
    assert len(reports)==8 and witness is not None and unique is not None
    (HERE/'rational_ambiguity_witness_1297.json').write_text(json.dumps(witness,indent=2)+'\n')
    (HERE/'q101_uniqueness.json').write_text(json.dumps(unique,indent=2)+'\n')
    out={'status':'passed','cases':reports,'pointwise_dual_inequalities':inequalities,'primal_equalities':equalities,
         'corrupted_certificates_rejected':mutations,'source_sha256':digest(Path(__file__)),
         'reviewed_source_sha256':EXPECTED_SOURCE,'reviewed_results_sha256':EXPECTED_RESULTS,
         'input_sha256':saved['input_sha256'],'runtime':{'executable':sys.executable,'version':sys.version},
         'elapsed_seconds':time.perf_counter()-started,
         'scope':'Independent model construction and exact saved-certificate audit using Python standard library; no optimizer, graph-realizability proof, old-backend replay, or prize bound.'}
    (HERE/'review_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
