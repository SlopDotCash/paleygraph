#!/usr/bin/env python3
"""Exact independent-sign benchmark, not an actual conference-matrix model."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path

HERE=Path(__file__).resolve().parent


def enc(x):return [x.numerator,x.denominator]


def coefficient(signs,d):
    values=[1]+[0]*d
    for sign in signs:
        for j in range(d,0,-1):values[j]+=sign*values[j-1]
    return values[d]


def main():
    checks=[]
    for n in (6,7,8):
        gram=[[0]*n for _ in range(n)]
        for signs in product((-1,1),repeat=n):
            r=[coefficient(signs[:a]+signs[a+1:],5) for a in range(n)]
            for a in range(n):
                for b in range(n):gram[a][b]+=r[a]*r[b]
        c0=comb(n-2,5) if n>=7 else 0;c1=comb(n-2,4)
        assert all(F(gram[a][b],2**n)==c0+(c1 if a==b else 0) for a in range(n) for b in range(n))
        ratio=F((n-1)*c1*c1,(n*c0+c1)**2+(n-1)*c1*c1)
        assert ratio==F(25,25+(n-1)*(n-5)**2)
        checks.append({'n':n,'sign_vectors':2**n,'gram_entries_checked':n*n})
    path=HERE/'contrast_results.json';data=json.loads(path.read_text());rows=[]
    for r in data['cases'][4:]:
        n=r['n'];reference=F(25,25+(n-1)*(n-5)**2);observed=F(*r['contrast_to_full_squared_norm_fraction'])
        rows.append({'q':r['q'],'n':n,'family':r['family'],'actual_contrast_fraction':enc(observed),
                     'independent_sign_reference_fraction':enc(reference),'actual_over_reference':enc(observed/reference),
                     'actual_over_reference_float':float(observed/reference)})
    out={'status':'passed','scope':'Independent Rademacher deletion-derivative Gram benchmark; not a realized conference sign matrix or a claimed Paley asymptotic.',
         'identity':'E[r_a r_b] = binom(n-2,5) + 1_(a=b) binom(n-2,4)',
         'contrast_fraction':'25 / (25 + (n-1)*(n-5)^2)',
         'reference_asymptotic':'n^3 times this independent-sign fraction tends to25; no transfer to Paley inputs is asserted.',
         'literal_reference_checks':checks,'comparisons':rows,
         'source_sha256':{'reference_model.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'contrast_results.json':sha256(path.read_bytes()).hexdigest()}}
    (HERE/'reference_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if k in ('q','n','family','actual_over_reference_float')} for r in rows]),flush=True)


if __name__=='__main__':main()
