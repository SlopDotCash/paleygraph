#!/usr/bin/env python3
"""Tiny exact row-triple bridge; never enumerates six-column subsets."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def mul(a,b):
    return ((a%7)*(b%7)-(a//7)*(b//7))%7+7*(((a%7)*(b//7)+(a//7)*(b%7))%7)


def sub(a,b):return (a%7-b%7)%7+7*((a//7-b//7)%7)


def power(a,k):
    out=1
    for _ in range(k):out=mul(out,a)
    return out


def e6_counts(positive,negative):
    return sum((-1)**j*(comb(negative,j) if j<=negative else 0)*
               (comb(positive,6-j) if 0<=6-j<=positive else 0) for j in range(7))


def e6_tau(nonzero,tau):
    assert abs(tau)<=nonzero and (nonzero+tau)%2==0
    value=e6_counts((nonzero+tau)//2,(nonzero-tau)//2)
    numerator=(tau**6-(15*nonzero-40)*tau**4+
               (45*nonzero**2-210*nonzero+184)*tau**2-
               15*nonzero*(nonzero-2)*(nonzero-4))
    assert numerator==720*value
    return value


def triple_types_from_invariants(q,edges,tau):
    # edges e01,e02,e12. The three marked columns have one zero each.
    e01,e02,e12=edges
    zero_patterns=[(0,e01,e02),(e01,0,e12),(e02,e12,0)]
    linear=[-e01-e02,-e01-e12,-e02-e12]
    quadratic=[-1-e02*e12,-1-e01*e12,-1-e01*e02]
    out=Counter(zero_patterns)
    for a,b,c in product((-1,1),repeat=3):
        numerator=q-3+a*linear[0]+b*linear[1]+c*linear[2]+a*b*quadratic[0]+a*c*quadratic[1]+b*c*quadratic[2]+a*b*c*tau
        assert numerator>=0 and numerator%8==0
        if numerator:out[a,b,c]+=numerator//8
    return out


def evaluate_case(name,sign,multiply,subtract,powfun,semilinear_exponent):
    q=len(sign)
    S=[[sign[subtract(x,y)] for y in range(q)] for x in range(q)]
    assert all(S[x][y]==S[y][x] for x in range(q) for y in range(q))
    assert all(sum(row)==0 for row in S)
    assert all(sum(S[x][z]*S[y][z] for z in range(q))==(q-1 if x==y else -1)
               for x in range(q) for y in range(q))
    normalization_checks=0
    for d in range(1,q):
        exponent=semilinear_exponent(d)
        de=powfun(d,exponent)
        inverse=next(x for x in range(1,q) if multiply(x,de)==1)
        image=[multiply(inverse,powfun(x,exponent)) for x in range(q)]
        assert sorted(image)==list(range(q)) and image[0]==0 and image[d]==1
        assert all(sign[image[x]]==sign[d]*sign[x] for x in range(q))
        normalization_checks+=1
    histogram=Counter();joint=Counter();normalized_sum=0;types_checked=0
    for t in range(2,q):
        values=[S[0][z]*S[1][z]*S[t][z] for z in range(q)]
        tau=sum(values);assert values.count(0)==3
        histogram[tau]+=1
        edges=(S[0][1],S[0][t],S[1][t])
        joint[(*edges,tau)]+=1
        normalized_sum+=e6_tau(q-3,tau)
        actual=Counter((S[0][z],S[1][z],S[t][z]) for z in range(q))
        assert actual==triple_types_from_invariants(q,edges,tau)
        types_checked+=1
    diagonal=-q*comb((q-1)//2,3)
    double=-3*q*(q-1)*comb((q-3)//2,3)
    distinct=q*(q-1)*normalized_sum
    numerator=diagonal+double+distinct
    # Separate unnormalized row-triple oracle, still no column six-sets.
    direct_distinct=0
    for a,b,c in combinations(range(q),3):
        values=[S[a][z]*S[b][z]*S[c][z] for z in range(q)]
        direct_distinct+=6*e6_counts(values.count(1),values.count(-1))
        actual=Counter((S[a][z],S[b][z],S[c][z]) for z in range(q))
        assert actual==triple_types_from_invariants(q,(S[a][b],S[a][c],S[b][c]),sum(values))
        types_checked+=1
    assert direct_distinct==distinct
    direct_double=0
    for a in range(q):
        for b in range(q):
            if a==b:continue
            values=[S[a][z]**2*S[b][z] for z in range(q)]
            direct_double+=3*e6_counts(values.count(1),values.count(-1))
    assert direct_double==double
    moment=Fraction(numerator,comb(q,6))
    return {'graph':name,'q':q,'normalized_parameters':q-2,
            'normalization_maps_verified':normalization_checks,
            'all_distinct_row_triple_type_checks':types_checked,
            'tau_histogram':sorted(histogram.items()),
            'tau_even_power_sums':{str(k):sum(count*tau**k for tau,count in histogram.items()) for k in (0,2,4,6)},
            'joint_edges_tau_histogram':[{'edges':list(key[:3]),'tau':key[3],'count':count} for key,count in sorted(joint.items())],
            'row_partition_numerators':{'all_equal':diagonal,'exactly_two_equal':double,'all_distinct':distinct},
            'sum_T6_cubed':numerator,'third_moment':[moment.numerator,moment.denominator],
            'six_column_subsets_enumerated':0}


def main():
    cases=[]
    for p in (13,17):
        sign=[0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1) for x in range(p)]
        cases.append(evaluate_case(f'Paley{p}',sign,lambda a,b:a*b%p,lambda a,b:(a-b)%p,
                                   lambda a,k:pow(a,k,p),lambda d:1))
    powers=[power(9,k) for k in range(48)]
    assert len(set(powers))==48 and power(9,48)==1
    logarithm={x:k for k,x in enumerate(powers)}
    for name,classes in (('Paley49',(0,2)),('Peisert49',(0,1))):
        sign=[0]+[1 if logarithm[x]%4 in classes else -1 for x in range(1,49)]
        exponent=(lambda d:1) if name=='Paley49' else (lambda d:1 if logarithm[d]%2==0 else 7)
        cases.append(evaluate_case(name,sign,mul,sub,power,exponent))
    assert cases[-2]['third_moment']==[44555,35673]
    assert cases[-1]['third_moment']==[6145,3243]
    known=Path(__file__).resolve().parents[2]/'round3'/'pair_type_twins'/'results.json'
    out={'status':'both twin third moments recovered without six-set enumeration','date':'2026-09-05',
         'cases':cases,'source_sha256':{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()},
         'acceptance_reference_sha256':sha256(known.read_bytes()).hexdigest(),
         'scope':'Tiny row-triple correctness preflight, not a general n>6 compiler or a historical novelty claim.'}
    (HERE/'third_moment_preflight.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{'graph':r['graph'],'third_moment':r['third_moment'],'tau_even_power_sums':r['tau_even_power_sums']} for r in cases],indent=2))


if __name__=='__main__':main()
