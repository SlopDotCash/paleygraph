#!/usr/bin/env python3
"""Exact q-polynomial compiler with a formal degree bound, not fitted data.

At target degree(6,6,6), each q-dependent logarithm factor consumes >=2
target degrees. Thus the raw q degree is at most9. theta^5 vanishes and
theta^6 is universal, so six admissible nodes check each quartic remainder.
"""
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
import time
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'round5/trace_invariants'))
from folded_moment import compiler,polynomial_from_values,encode
from residue_moment import universal_sixth_union

spec=importlib.util.spec_from_file_location('critical_third_second',LAB/'round4/marked_moments/conditional_moments.py')
second=importlib.util.module_from_spec(spec);spec.loader.exec_module(second)


def trim(poly):
    poly=list(poly)
    while len(poly)>1 and not poly[-1]:poly.pop()
    return poly


def value(poly,q):return sum((c*q**j for j,c in enumerate(poly)),Fraction())


def add(*polys):
    out=[Fraction() for _ in range(max(map(len,polys)))]
    for poly in polys:
        for i,c in enumerate(poly):out[i]+=c
    return trim(out)


def mul(a,b):
    out=[Fraction() for _ in range(len(a)+len(b)-1)]
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]+=c*d
    return trim(out)


def scale(a,c):return trim([c*x for x in a])


def compile_raw(qs):
    histograms=[];plans=[];universal=universal_sixth_union()
    for q in qs:
        start=len(histograms);histograms.extend(compiler.repeated_histograms(q));families=[]
        for mono in (False,True):
            edges=(1,1,1) if mono else(1,1,-1)
            taus=compiler.admissible_tau_nodes(q,edges)[:6]
            assert len(taus)==6
            families.append((mono,len(histograms),[t if mono else-t for t in taus]))
            histograms.extend(compiler.distinct_histogram(q,edges,t) for t in taus)
        plans.append((q,start,families))
    raw,elapsed=compiler.union_coefficients_batch(histograms,6,max_union=18)
    output=[]
    for q,start,families in plans:
        assert raw[start+1]['coefficients']==raw[start+2]['coefficients']
        row={'q':q,'repeated_raw':[r['coefficients'] for r in raw[start:start+3]],'families':{}}
        for mono,index,nodes in families:
            polys=[]
            for k in range(19):
                vals=[Fraction(r['coefficients'][k])-universal[k]*t**6 for r,t in zip(raw[index:index+6],nodes)]
                poly=polynomial_from_values(nodes[:5],vals[:5])
                assert value(poly,nodes[-1])==vals[-1]
                polys.append(poly)
            row['families'][mono]=polys
        output.append(row)
    return output,{'coefficient_evaluations':len(histograms),'backend_seconds':elapsed,
                  'maximum_intermediate_bits':max(x['maximum_intermediate_bits'] for x in raw)}


def main():
    started=time.perf_counter();qs=list(range(81,162,8))
    raw,cost=compile_raw(qs);universal=universal_sixth_union()
    families={};repeated=[]
    for mono in (False,True):
        by_k=[]
        for k in range(19):
            by_j=[]
            for j in range(5):
                bound=min(k,(18-3*j)//2)
                nodes=qs[:bound+1]
                poly=trim(polynomial_from_values(nodes,[row['families'][mono][k][j] for row in raw[:bound+1]]))
                assert all(value(poly,q)==row['families'][mono][k][j] for q,row in zip(qs,raw))
                by_j.append(poly)
            by_k.append(by_j)
        families[mono]=by_k
    for index in range(3):
        by_k=[]
        for k in range(19):
            bound=min(k,9)
            poly=trim(polynomial_from_values(qs[:bound+1],[row['repeated_raw'][index][k] for row in raw[:bound+1]]))
            assert all(value(poly,q)==row['repeated_raw'][index][k] for q,row in zip(qs,raw))
            by_k.append(poly)
        repeated.append(by_k)
    B=[];stats=[];qq=[0,-1,1]
    for k in range(19):
        a,b=families[True][k],families[False][k]
        baseline=add(mul(a[0],[Fraction(-5,4),Fraction(1,4)]),mul(b[0],[Fraction(-3,4),Fraction(3,4)]),
                     scale(b[1],-6),mul(b[2],[-3,-2,1]))
        B.append(add(mul([0,1],repeated[0][k]),scale(mul(qq,repeated[1][k]),3),mul(qq,baseline)))
        raw_stats=[add(a[1],scale(b[1],3)),add(a[2],scale(b[2],-1)),a[3],b[3],a[4],b[4],[universal[k]]]
        stats.append([mul(qq,poly) for poly in raw_stats])
    # Independent-q checks use the other mod8 branch and remote large q.
    check_qs=[49,101,1297,65537,1000033]
    extra,extra_cost=compile_raw(check_qs)
    for row in extra:
        q=row['q']
        for k in range(19):
            for i in range(3):assert value(repeated[i][k],q)==row['repeated_raw'][i][k]
            for mono in (False,True):
                for j in range(5):assert value(families[mono][k][j],q)==row['families'][mono][k][j]
    # Degree12 target for the pair compiler permits q degree<=6; global
    # row weights increase this to at most8.
    second_values=[]
    for q in qs:
        diagonal=second.pair_coefficients(q,6,12,(0,()))
        pos=second.pair_coefficients(q,6,12,(1,()))
        neg=second.pair_coefficients(q,6,12,(-1,()))
        assert pos==neg
        second_values.append([q*a+q*(q-1)*b for a,b in zip(diagonal,pos)])
    second_polys=[]
    for k in range(13):
        poly=trim(polynomial_from_values(qs[:9],[row[k] for row in second_values[:9]]))
        assert all(value(poly,q)==row[k] for row,q in zip(second_values,qs))
        second_polys.append(poly)
    out={'scope':'Exact polynomial identities from proved finite degree bounds; synthetic q nodes are not asserted prime or graph-realizable.',
         'q_nodes':qs,'independent_q_checks':check_qs,'sample_cost':cost,'extra_check_cost':extra_cost,
         'raw_q_degree_bound':9,'theta_coefficient_q_degree_bound':'min(k,floor((18-3j)/2))',
         'statistic_names':['A1','A2','A3','B3','A4','B4','M6'],
         'baseline_by_union_k':[[encode(c) for c in p] for p in B],
         'statistic_coefficients_by_union_k':[[[encode(c) for c in p] for p in row] for row in stats],
         'second_moment_by_union_k':[[encode(c) for c in p] for p in second_polys],
         'family_quartics_by_union_k':{str(m):[[[encode(c) for c in p] for p in row] for row in families[m]] for m in(False,True)},
         'repeated_by_union_k':[[[encode(c) for c in p] for p in rows] for rows in repeated],
         'elapsed_seconds':time.perf_counter()-started,
         'source_sha256':{str(p.relative_to(LAB)):sha256(p.read_bytes()).hexdigest() for p in
             (Path(__file__),LAB/'round5/trace_invariants/folded_moment.py',LAB/'round5/trace_invariants/residue_moment.py',
              LAB/'round5/third_moment/third_moment.py',LAB/'round5/third_moment/union_coefficients.cpp',
              LAB/'round5/third_moment/union_coefficients',LAB/'round4/marked_moments/conditional_moments.py')}}
    (HERE/'q_polynomials.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'sample_cost':cost,'extra_check_cost':extra_cost,'elapsed_seconds':out['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
