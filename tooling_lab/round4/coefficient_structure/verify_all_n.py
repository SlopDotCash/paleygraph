#!/usr/bin/env python3
"""Check every valid n by the general compiler's weighted pair coefficients."""
from fractions import Fraction
import hashlib,json,sys,time
from pathlib import Path
import sympy as sp
from coefficient_polynomial import HERE,ROUND,PATS
from closed_coefficient import contraction_coefficient_degree6,quartic
from contraction_reduction import coefficient,walsh_decomposition

def changes(values):
    segments=[]
    for n,sign in values:
        if not segments or segments[-1]['sign']!=sign:segments.append({'from_n':n,'through_n':n,'sign':sign})
        else:segments[-1]['through_n']=n
    return segments

def main():
    start=time.monotonic();rows=[]
    formula=json.loads((HERE/'symbolic_formula.json').read_text());q,r=sp.symbols('q r')
    raw=sp.sympify(formula['formula_in_q_r'],locals={'q':q,'r':r})
    candidate=-sp.ff(r,4)*sp.ff(q-r-3,4)*quartic(q,r)/(2*sp.ff(q-3,12))
    assert sp.cancel(raw-candidate)==0
    for p in [17,21,25,29,49,101,1297]:
        started=time.monotonic();signs=[];zeros=[]
        for n in range(6,p+1):
            expected=sum((Fraction(u[0]*u[1]*v[0]*v[1]*v[2],64)*
                  (coefficient(p,n,6,u,v,1)-coefficient(p,n,6,u,v,-1)) for u in PATS for v in PATS),Fraction())
            actual=contraction_coefficient_degree6(p,n);assert expected==actual
            sign=(actual>0)-(actual<0);signs.append((n,sign))
            if not actual:zeros.append(n)
        # Check the reducer's full Walsh reconstruction on representative n too.
        extra=sorted({6,7,8,p//2,p-4,p-3,p})
        for n in extra:assert walsh_decomposition(p,n,6)[1]==contraction_coefficient_degree6(p,n)
        quart=sp.Poly(quartic(p,r),r)
        intervals=sp.polys.polytools.intervals(quart,eps=sp.Rational(1,10**10))
        isolations=[{'r_interval':[[int(a.p),int(a.q)],[int(b.p),int(b.q)]],'multiplicity':multiplicity,
                     'approx_n':float((a+b)/2+3)} for (a,b),multiplicity in intervals]
        row={'q':p,'all_n_from':6,'all_n_through':p,'exact_compiler_comparisons':p-5,'full_walsh_reconstructions':extra,
             'integer_zeros_in_valid_domain':zeros,'sign_intervals':changes(signs),
             'quartic_real_root_isolations':isolations,'quartic_degree':int(quart.degree()),
             'c_degree_in_r':8+int(quart.degree()),'seconds':time.monotonic()-started}
        rows.append(row);print(json.dumps({'q':p,'zeros':zeros,'sign_intervals':row['sign_intervals'],'seconds':row['seconds']}),flush=True)
    out={'status':'all exact symbolic and every-cardinality coefficient checks passed','cases':rows,
         'total_exact_cardinality_checks':sum(r['exact_compiler_comparisons'] for r in rows),'elapsed_seconds':time.monotonic()-start,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'closed_formula_sha256':hashlib.sha256((HERE/'closed_coefficient.py').read_bytes()).hexdigest(),
         'compiler_sha256':hashlib.sha256((ROUND/'marked_moments/conditional_moments.py').read_bytes()).hexdigest(),
         'reduction_sha256':hashlib.sha256((ROUND/'marked_moments/contraction_reduction.py').read_bytes()).hexdigest()}
    (HERE/'all_n_validation.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
