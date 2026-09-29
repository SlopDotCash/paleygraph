#!/usr/bin/env python3
"""Exact union coefficients of the degree-six three-mark contraction coefficient."""
from fractions import Fraction
from functools import lru_cache
import hashlib,itertools,json,sys,time
from pathlib import Path
import sympy as sp
HERE=Path(__file__).resolve().parent
ROUND=HERE.parent
sys.dont_write_bytecode=True;sys.path.insert(0,str(ROUND/'marked_moments'))
from conditional_moments import canonical_pair,pair_coefficients
PATS=tuple(itertools.product((-1,1),repeat=3))

@lru_cache(maxsize=None)
def union_coefficients(q):
    values=[Fraction() for k in range(13)]
    for u in PATS:
        for v in PATS:
            sign=u[0]*u[1]*v[0]*v[1]*v[2]
            plus=pair_coefficients(q,6,12,canonical_pair(u,v,1))
            minus=pair_coefficients(q,6,12,canonical_pair(u,v,-1))
            for k in range(13):values[k]+=Fraction(sign*(plus[k]-minus[k]),64)
    return tuple(values)

def symbolic(q,coefs):
    r=sp.Symbol('r');return sp.expand(sum(sp.Rational(a.numerator,a.denominator)*sp.ff(r,k)/sp.ff(q-3,k) for k,a in enumerate(coefs)))

def main():
    start=time.monotonic();rows=[]
    for q in [29,49,101,1297]:
        coefficients=union_coefficients(q);pol=symbolic(q,coefficients);factor=sp.factor(pol)
        rows.append({'q':q,'union_coefficients':[[c.numerator,c.denominator] for c in coefficients],
            'degree_in_r':sp.degree(pol),'factored_in_r':str(factor),'r_definition':'n-3',
            'rational_roots':{str(root):multiplicity for root,multiplicity in sp.roots(pol).items() if root.is_Rational}})
        print(q,factor,flush=True)
    out={'scope':'degree6 and exactly3marks; unweighted union size <=12','cases':rows,
         'elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'compiler_sha256':hashlib.sha256((ROUND/'marked_moments/conditional_moments.py').read_bytes()).hexdigest(),
         'sympy_version':sp.__version__,'python':sys.executable}
    (HERE/'fixed_q_polynomials.json').write_text(json.dumps(out,indent=2,default=int)+'\n')
if __name__=='__main__':main()
