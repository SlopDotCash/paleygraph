#!/usr/bin/env python3
"""Prime-field rank-two polynomial-cluster certificate interface."""
from math import isqrt
from rank_two_pinning import compile_pinning


def evaluate(coefficients,x,p):
    value=0
    for c in reversed(coefficients):value=(value*x+c)%p
    return value


def certify(data):
    p,k,s=(data[name] for name in ('p','k','s'));domain=data['domain'];basis=data['basis'];origin=data.get('origin',[])
    if not (type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
            and type(k) is int and 2<=k<=len(domain)<=min(p,100_000)
            and type(s) is int and 0<=s<=len(domain)):
        raise ValueError('invalid prime, degree limit, domain length or agreement threshold')
    if len(set(domain))!=len(domain) or any(type(x) is not int or not 0<=x<p for x in domain):
        raise ValueError('domain must contain distinct canonical field elements')
    if len(basis)!=2 or any(len(poly)>k or any(type(c) is not int or not 0<=c<p for c in poly) for poly in [origin,*basis]):
        raise ValueError('two canonical degree-bounded basis polynomials required')
    columns=[(evaluate(basis[0],x,p),evaluate(basis[1],x,p)) for x in domain]
    result=compile_pinning(p,columns,s)
    return {'schema':'rank_two_polynomial_pinning_v1','p':p,'n':len(domain),'k':k,'s':s,
            'domain':list(domain),'origin':list(origin),'basis':[list(b) for b in basis],
            'scope':'Uniform unordered coordinate-pair success within this declared affine polynomial cluster; cluster discovery and coverage are not asserted.',
            'pinning':result}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv)!=2:raise SystemExit('Usage: polynomial_cluster.py input.json')
    print(json.dumps(certify(json.loads(Path(sys.argv[1]).read_text())),indent=2))
