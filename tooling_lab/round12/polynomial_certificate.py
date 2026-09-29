#!/usr/bin/env python3
"""Declared affine polynomial-space certificate adapter for both rank-three tools."""
from math import isqrt
from capacitated_pinning import optimize,partition,determinant
from pair_space_preflight import catalogue,cross_key,bound_and_witness,evaluate


def certify(data,method='auto'):
    p,k,s=data['p'],data['k'],data['s'];domain=data['domain'];basis=data['basis'];origin=data.get('origin',[]);n=len(domain)
    if not (type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
            and type(k) is int and 3<=k<=n<=min(p,4096) and type(s) is int and 0<=s<=n):
        raise ValueError('invalid field, degree limit, length or agreement threshold')
    if len(set(domain))!=n or any(type(x) is not int or not 0<=x<p for x in domain):raise ValueError('distinct canonical domain required')
    if len(basis)!=3 or any(len(poly)>k or any(type(a) is not int or not 0<=a<p for a in poly) for poly in [origin,*basis]):raise ValueError('three canonical degree-bounded basis polynomials required')
    columns=[tuple(evaluate(poly,x,p) for poly in basis) for x in domain]
    zeros,groups=partition(p,columns);m=len(groups)
    if m<3:raise ValueError('basis must have evaluation rank three')
    a=groups[0]['direction'];b=groups[1]['direction']
    if not any(determinant(a,b,g['direction'],p) for g in groups[2:]):raise ValueError('basis must have evaluation rank three')
    if method=='auto':method='capacitated_corners' if m<=22 else 'simple_pair_degree'
    if method=='capacitated_corners':proof=optimize(p,columns,s)
    elif method=='simple_pair_degree':
        if zeros or m!=n or n>1024:raise ValueError('pair-degree method requires a simple space of length at most1024')
        ordinary,lines=catalogue(columns,p,cross_key)
        proof={'ordinary_two_point_lines':ordinary,'nontrivial_maximal_lines':[list(L) for L in lines],
               'bound':bound_and_witness(n,s,lines)}
    else:raise ValueError('unknown method')
    return {'schema':'polynomial_rank_three_pinning_v1','method':method,
            'input':{'p':p,'n':n,'k':k,'s':s,'domain':list(domain),'origin':list(origin),'basis':[list(b) for b in basis]},
            'proof':proof,'scope':'Pinning in this declared affine polynomial space only. No discovery, scalar-specific cluster cover, full decoder or prize theorem is asserted.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv) not in (2,3):raise SystemExit('Usage: polynomial_certificate.py input.json [auto|capacitated_corners|simple_pair_degree]')
    print(json.dumps(certify(json.loads(Path(sys.argv[1]).read_text()),sys.argv[2] if len(sys.argv)==3 else 'auto'),separators=(',',':')))
