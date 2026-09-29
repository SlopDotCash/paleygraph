#!/usr/bin/env python3
"""Derive the general-q coefficient by an exact common-factor recurrence.
No interpolation: q->q+4 adds one of each nonzero paired column type.
"""
from collections import defaultdict
from fractions import Fraction
import hashlib,json,math,sys,time
from pathlib import Path
import sympy as sp
from coefficient_polynomial import HERE,ROUND,PATS,union_coefficients
from conditional_moments import full_pair_histogram,type_power,multiply,elementary


def base_walsh_polynomial():
    W=defaultdict(int)
    for u in PATS:
        left=elementary(u,6)
        for v in PATS:
            right=elementary(v,6);walshsign=u[0]*u[1]*v[0]*v[1]*v[2]
            marked=[((i,j,0),a*b) for i,a in enumerate(left) for j,b in enumerate(right) if a and b]
            for s in [-1,1]:
                hist=full_pair_histogram(17,s);hist.subtract(zip(u,v));assert min(hist.values())>=0
                pol={(0,0,0):1}
                for (a,b),count in sorted(hist.items()):
                    if count and (a or b):pol=multiply(pol,type_power(a,b,count,6,12),6,12)
                pol=multiply(pol,marked,6,12)
                for key,value in pol.items():W[key]+=s*walshsign*value
    return {key:Fraction(value,64) for key,value in W.items() if value}


def common_factor():
    B={(0,0,0):1}
    for a,b in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        B=multiply(B,[((0,0,0),1),((1,0,1),a),((0,1,1),b),((1,1,1),a*b)],6,12)
    assert B.pop((0,0,0))==1 and all(k>=2 for i,j,k in B)
    return B


def main():
    start=time.monotonic();q,r=sp.symbols('q r');W=base_walsh_polynomial();D=common_factor();T=[];current=W
    for j in range(7):
        row=[current.get((6,6,k),Fraction()) for k in range(13)];T.append(row)
        if j<6:current=multiply(current,list(D.items()),6,12)
    assert tuple(T[0])==union_coefficients(17)
    h=(q-17)/4
    G=[]
    for k in range(13):
        expr=0
        for j in range(7):
            value=T[j][k]
            expr+=sp.Rational(value.numerator,value.denominator)*sp.prod(h-i for i in range(j))/math.factorial(j)
        G.append(sp.factor(expr))
        assert not expr or sp.degree(sp.expand(expr),q)<=k//2
    # Common falling-factorial denominator gives a direct symbolic polynomial.
    denominator=sp.ff(q-3,12)
    numerator=sp.expand(sum(G[k]*sp.ff(r,k)*sp.ff(q-3-k,12-k) for k in range(13)))
    factored=sp.factor(numerator/denominator)
    fixed_checks=[]
    for value in [17,21,25,29,49,101,1297]:
        numerical=union_coefficients(value)
        assert all(sp.Rational(x.numerator,x.denominator)==g.subs(q,value) for x,g in zip(numerical,G))
        fixed_checks.append(value)
    t,u,z=sp.symbols('t u z');B=1+sum(value*t**i*u**j*z**k for (i,j,k),value in D.items())
    out={'status':'general-q formula derived by exact common-factor recurrence, not interpolation',
         'domain':'q>=17 and q=1 mod4; exactly three marks; degree6; r=n-3',
         'q_step':4,'base_q':17,'common_factor_B':str(sp.factor(B)),
         'minimum_v_degree_of_B_minus_one':2,'maximum_binomial_power_needed':6,
         'base_W_nonzero_terms':len(W),'coefficient_degree_bound_in_r':12,
         'binomial_basis_T':[[[a.numerator,a.denominator] for a in row] for row in T],
         'union_coefficients_G_as_polynomials_in_q':[str(g) for g in G],
         'formula_in_q_r':str(factored),'numerator_factored':str(sp.factor(numerator)),
         'common_denominator':str(denominator),'actual_degree_in_r':int(sp.degree(numerator,r)),
         'exact_union_coefficient_replays':fixed_checks,'elapsed_seconds':time.monotonic()-start,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'compiler_sha256':hashlib.sha256((ROUND/'marked_moments/conditional_moments.py').read_bytes()).hexdigest(),
         'sympy_version':sp.__version__,'python':sys.executable}
    (HERE/'symbolic_formula.json').write_text(json.dumps(out,indent=2)+'\n')
    print('B:',sp.expand(B),flush=True)
    print('G:',G,flush=True)
    print('GENERAL:',factored,flush=True)
    print('seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
