#!/usr/bin/env python3
"""Independent literal-column checks and critical-weight audit of q-polynomials.

No production compiler or interpolation helper is imported. The earlier
independent literal-column oracle is reused for the three-row product.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb,factorial
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'round5/review'))
from review_union_backend import literal_columns


def encode(x):
    x=F(x)
    return [x.numerator,x.denominator]


def val(poly,q):return sum(F(*c)*q**j for j,c in enumerate(poly))


def paley(q):
    chi=[0]+[1 if pow(x,(q-1)//2,q)==1 else -1 for x in range(1,q)]
    return [[chi[(y-x)%q] for y in range(q)] for x in range(q)]


def pair_columns(columns):
    state={(0,0,0):1}
    for a,b in columns:
        nxt=state.copy()
        for (k,i,j),old in state.items():
            for x,y,c in ((1,0,a),(0,1,b),(1,1,a*b)):
                if c and k<12 and i+x<=6 and j+y<=6:
                    key=k+1,i+x,j+y
                    nxt[key]=nxt.get(key,0)+old*c
        state={key:v for key,v in nxt.items() if v}
    return [state.get((k,6,6),0) for k in range(13)]


def terms_by_weight(polys,extra=F(0)):
    records=[]
    for k,poly in enumerate(polys):
        for degree,c in enumerate(poly):
            if F(*c):records.append({'union':k,'q_degree':degree,'coefficient':c,'critical_weight':encode(F(degree)-F(3*k,4)+extra)})
    return sorted(records,key=lambda r:(F(*r['critical_weight']),r['union'],r['q_degree']),reverse=True)


def main():
    start=time.perf_counter()
    source=LAB/'round6/critical_third/q_polynomials.json'
    p=json.loads(source.read_text())
    for name,h in p['source_sha256'].items():assert sha256((LAB/name).read_bytes()).hexdigest()==h
    # Independently multiply K(z)^6; no imported universal-sixth helper.
    sixth=[F(1)]
    for _ in range(6):
        out=[F() for _ in range(len(sixth)+3)]
        for i,c in enumerate(sixth):
            for j,v in enumerate((0,1,-3,2)):out[i+j]+=c*v
        sixth=out
    sixth=[c/factorial(6) for c in sixth]
    rows=[];raw_equalities=0;second_equalities=0
    for q in (13,17,29):
        S=paley(q)
        assert all(sum(row)==0 for row in S)
        assert all(sum(S[a][x]*S[b][x] for x in range(q))==(q-1 if a==b else -1) for a in range(q) for b in range(q))
        pairs={s:next(t for t in range(1,q) if S[0][t]==s) for s in (-1,1)}
        triples=[(0,0,0),(0,0,pairs[-1]),(0,0,pairs[1])]
        distinct={}
        for a,b,c in combinations(range(q),3):
            edges=(S[a][b],S[a][c],S[b][c]);mono=len(set(edges))==1
            theta=edges[0]*edges[1]*edges[2]*sum(S[a][x]*S[b][x]*S[c][x] for x in range(q))
            distinct.setdefault((mono,theta),(a,b,c))
        # Every actually occurring family/theta type, including small fields
        # outside the q-interpolation sample interval.
        triples.extend(distinct.values())
        for a,b,c in triples:
            literal=literal_columns(list(zip(S[a],S[b],S[c])),6,min(q,18))
            literal += [0]*(19-len(literal))
            if a==b==c:
                expected=[val(poly,q) for poly in p['repeated_by_union_k'][0]]
                tag='all_equal'
            elif a==b:
                index=1 if S[a][c]==-1 else 2
                expected=[val(poly,q) for poly in p['repeated_by_union_k'][index]]
                tag='two_equal'
            else:
                edges=(S[a][b],S[a][c],S[b][c]);mono=len(set(edges))==1
                theta=edges[0]*edges[1]*edges[2]*sum(S[a][x]*S[b][x]*S[c][x] for x in range(q))
                expected=[sum(val(poly,q)*theta**j for j,poly in enumerate(by_j))+sixth[k]*theta**6
                          for k,by_j in enumerate(p['family_quartics_by_union_k'][str(mono)])]
                tag=f'{mono=},{theta=}'
            assert list(map(F,literal))==expected,(q,(a,b,c))
            raw_equalities+=19
            rows.append({'q':q,'triple':[a,b,c],'type':tag,'exact_union_coefficients':19})
        same=pair_columns(list(zip(S[0],S[0])))
        different=[pair_columns(list(zip(S[0],S[pairs[s]]))) for s in (-1,1)]
        assert different[0]==different[1]
        second=[q*a+q*(q-1)*b for a,b in zip(same,different[0])]
        assert list(map(F,second))==[val(poly,q) for poly in p['second_moment_by_union_k']]
        second_equalities+=13
        print(f'q{q}: {len(triples)} literal triple types and all second-moment union coefficients passed',flush=True)

    baseline=terms_by_weight(p['baseline_by_union_k'])
    second=terms_by_weight(p['second_moment_by_union_k'])
    assert baseline[0]=={'union':9,'q_degree':10,'coefficient':[1,216],'critical_weight':[13,4]}
    assert F(*baseline[1]['critical_weight'])==F(5,2)
    assert second[0]=={'union':6,'q_degree':7,'coefficient':[1,720],'critical_weight':[5,2]}
    assert F(*second[1]['critical_weight'])==F(7,4)
    statistics=[]
    for i,(name,j) in enumerate(zip(p['statistic_names'],(1,2,3,3,4,4,6))):
        polys=[row[i] for row in p['statistic_coefficients_by_union_k']]
        terms=terms_by_weight(polys,F(1)+F(j,2))
        assert F(*terms[0]['critical_weight'])<=F(3,2)
        statistics.append({'statistic':name,'Hasse_mass_weight':encode(F(1)+F(j,2)),'maximum_term':terms[0]})
    # Fixed ordinary moments leave only A2,A3,A4 differences after A1 is
    # fixed by M1. This is independent of the production sensitivity code.
    fixed=[]
    for j in (2,3,4):
        a_index,b_index={2:(1,None),3:(2,3),4:(4,5)}[j]
        polys=[]
        for row in p['statistic_coefficients_by_union_k']:
            a=[F(*x) for x in row[a_index]]
            b=[F(*x) for x in row[b_index]] if b_index is not None else []
            poly=[(a[t] if t<len(a) else F())-(b[t] if t<len(b) else F()) for t in range(max(len(a),len(b)))]
            polys.append([encode(x) for x in poly])
        terms=terms_by_weight(polys,F(1)+F(j,2))
        assert F(*terms[0]['critical_weight'])<=-F(1,4)
        fixed.append({'statistic':f'A{j}','maximum_term':terms[0]})
    assert fixed[0]['maximum_term']['coefficient']==[2,1]
    assert fixed[2]['maximum_term']['coefficient']==[-1,3]
    # Combinatorial constants: each of the nine used columns appears in
    # precisely two of three six-sets; each pair-specific class has size3.
    assert factorial(3)**3==216 and factorial(6)==720
    assert F(720)**3/F(216)**2==8000  # (40 sqrt(5))^2
    out={'status':'passed','scope':'Independent literal polynomial checks and exact critical-weight exclusion; no production compiler import.',
         'literal_triple_cases':len(rows),'raw_union_equalities':raw_equalities,
         'second_moment_union_equalities':second_equalities,'records':rows,
         'baseline_top_terms':baseline[:8],'second_moment_top_terms':second[:8],
         'arithmetic_statistic_top_weights':statistics,'fixed_ordinary_top_weights':fixed,
         'leading_skew_constant_squared':[8000,1],
         'mean_formula':'-q*binom((q-1)/2,3)*(n)_6/(q)_6; critical weight -1/2',
         'source_sha256':{str(path.relative_to(LAB)):sha256(path.read_bytes()).hexdigest() for path in
              (Path(__file__),source,LAB/'round5/review/review_union_backend.py')},
         'elapsed_seconds':time.perf_counter()-start}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k in ('status','literal_triple_cases','raw_union_equalities','second_moment_union_equalities','elapsed_seconds')},indent=2))


if __name__=='__main__':main()
