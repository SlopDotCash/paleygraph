#!/usr/bin/env python3
"""Remove the common deletion mode before interpreting insertion covariance."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def enc(x):return [x.numerator,x.denominator]


def analyze(r):
    n=r['n'];C=r['covariance_numerator'];den=r['covariance_denominator'];sums=list(map(sum,C));total=sum(sums)
    centered=[[n*n*C[a][b]-n*(sums[a]+sums[b])+total for b in range(n)] for a in range(n)]
    assert all(sum(row)==0 for row in centered)
    contrast_norm=F(sum(x*x for row in centered for x in row),n**4*den**2)
    full_norm=F(*r['covariance_trace_squared'])
    assert 0<=contrast_norm<=full_norm
    return {'q':r['q'],'n':n,'selected':r['selected'],'family':r.get('family','saved_twin'),
            'contrast_covariance_numerator':centered,'contrast_covariance_denominator':n*n*den,
            'contrast_covariance_trace_squared':enc(contrast_norm),
            'contrast_to_full_squared_norm_fraction':enc(contrast_norm/full_norm) if full_norm else [0,1],
            'total_covariance':enc(F(total,den)),'trace_covariance':enc(F(sum(C[a][a] for a in range(n)),den)),
            'meaning':'P Cov_b(delta(b)) P for P=I-11^T/n; removes the common response shared by all deletion choices.'}


def main():
    source=HERE/'scale_results.json';data=json.loads(source.read_text())
    rows=[analyze(r) for r in data['tiny_witnesses']+data['scale_cases']]
    out={'status':'passed','cases':rows,'source_sha256':{'contrast_analysis.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'scale_results.json':sha256(source.read_bytes()).hexdigest()}}
    (HERE/'contrast_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if k in ('q','n','family','contrast_to_full_squared_norm_fraction')} for r in rows]),flush=True)


if __name__=='__main__':main()
