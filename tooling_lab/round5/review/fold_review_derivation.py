#!/usr/bin/env python3
"""Certified finite-degree reconstruction of the first folded class differences.

Raw E_k has total degree at most k in q,tau by Newton's recurrence with
linear sign power sums. At target(6,6,6), its tau degree is at most6 by the
separate logarithm argument. The grid sizes below follow those proven bounds.
"""
import sys
sys.dont_write_bytecode=True
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import time

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def interpolation(nodes,values):
    # Newton divided differences converted by nested multiplication. This
    # differs from the production adapter's Lagrange-basis construction.
    dd=list(map(Fraction,values))
    for j in range(1,len(nodes)):
        for i in range(len(nodes)-1,j-1,-1):dd[i]=(dd[i]-dd[i-1])/(nodes[i]-nodes[i-j])
    out=[dd[-1]]
    for i in range(len(nodes)-2,-1,-1):
        new=[Fraction() for _ in range(len(out)+1)]
        for k,v in enumerate(out):new[k]-=nodes[i]*v;new[k+1]+=v
        new[0]+=dd[i];out=new
    return out


def value(poly,x):return sum((c*x**i for i,c in enumerate(poly)),Fraction())


def main():
    started=time.perf_counter()
    source=HERE.parent/'third_moment'/'third_moment.py'
    compiler=load('fold_review_canonical_compiler',source)
    oracle=load('fold_review_literal_histogram',HERE/'direct_third_oracle.py')
    qs=list(range(49,89,4));histograms=[];plans=[]
    for q in qs:
        for mono in (False,True):
            edges=(1,1,1) if mono else (1,1,-1);sign=1 if mono else -1
            admissible=[]
            for tau in range(-(q-3),q-2):
                try:hist=oracle.triple_histogram(q,edges,tau)
                except AssertionError:continue
                admissible.append((sign*tau,hist))
            assert len(admissible)>=8
            picked=admissible[:8]
            plans.append((q,mono,len(histograms),[t for t,_ in picked]))
            histograms.extend(h for _,h in picked)
    rows,_=compiler.union_coefficients_batch(histograms,degree=6,max_union=8)
    theta_polys={};extra_checks=0
    for q,mono,index,nodes in plans:
        for k in (6,7,8):
            values=[row['coefficients'][k] for row in rows[index:index+8]]
            poly=interpolation(nodes[:7],values[:7])
            assert value(poly,nodes[7])==values[7];extra_checks+=1
            theta_polys[q,mono,k]=poly
    differences={(q,k):[a-b for a,b in zip(theta_polys[q,True,k],theta_polys[q,False,k])]
                 for q in qs for k in (6,7,8)}
    formulas={};q_extra=0
    for k in (6,7,8):
        bivariate=[]
        for theta_degree in range(7):
            # Total degree<=k bounds this coefficient's q degree<=k-j.
            q_degree=k-theta_degree
            samples=[differences[q,k][theta_degree] for q in qs]
            poly=interpolation(qs[:q_degree+1],samples[:q_degree+1])
            assert all(value(poly,q)==v for q,v in zip(qs,samples));q_extra+=len(qs)-q_degree-1
            bivariate.append([[v.numerator,v.denominator] for v in poly])
        formulas[str(k)]=bivariate
    assert all(v==[0,1] for p in formulas['6'] for v in p)
    assert all(v==[0,1] for k in ('7','8') for p in formulas[k][5:] for v in p)
    import sympy as sp
    Q,T=sp.symbols('q theta');expressions={}
    for k,rows in formulas.items():
        expression=sum(sp.Rational(*v)*Q**i*T**j for j,row in enumerate(rows) for i,v in enumerate(row))
        expressions[k]={'expanded':str(sp.expand(expression)),'factored':str(sp.factor(expression))}
    out={'status':'exact class-difference identities reconstructed under proved degree bounds',
         'definition':'Delta_k(q,theta)=coefficient[t1^6 t2^6 t3^6 z^k] mono minus mixed; mixed canonical tau=-theta',
         'synthetic_q_values':qs,'raw_coefficient_evaluations':len(histograms),
         'theta_extra_node_checks':extra_checks,'q_extra_node_checks':q_extra,
         'coefficient_arrays_theta_then_q':formulas,'expressions':expressions,
         'degree_justification':'Newton E_k uses products of at most k sign power sums linear in q and tau, so total degree<=k; the logarithm gives tau degree<=6.',
         'interpretation':'No graph realization is claimed for grid inventories. At n6 only k<=6 occurs and difference is zero. At n7 the difference is Delta7/C(q,7); at n8 it is8*Delta7/C(q,7)+Delta8/C(q,8).',
         'source_sha256':{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest(),
                         '../third_moment/third_moment.py':sha256(source.read_bytes()).hexdigest(),
                         '../third_moment/union_coefficients.cpp':sha256((source.parent/'union_coefficients.cpp').read_bytes()).hexdigest(),
                         'direct_third_oracle.py':sha256((HERE/'direct_third_oracle.py').read_bytes()).hexdigest()},
         'elapsed_seconds':time.perf_counter()-started}
    (HERE/'fold_review_derivation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='coefficient_arrays_theta_then_q'},indent=2))


if __name__=='__main__':main()
