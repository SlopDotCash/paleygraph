#!/usr/bin/env python3
"""Test labelled insertion correlations on actual one-swap distribution twins."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'round7/local_edits'))
from review_results import literal_signs, elementary, target

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]


def enc(x):return [x.numerator,x.denominator]


def analyze(c,q):
    S=literal_signs(q);n=len(c);outside=sorted(set(range(q))-set(c));m=q-n
    r=[[elementary([S[x][b] for b in c if b!=a],5) for x in range(q)] for a in c]
    V=[[sum(r[a][x]*S[x][b] for x in range(q)) for b in c] for a in range(n)]
    sums=list(map(sum,r));t=list(map(sum,V));u=[V[a][a] for a in range(n)]
    h=[-t[a]-m*u[a] for a in range(n)]
    gram=[[sum(r[a][x]*r[b][x] for x in range(q)) for b in range(n)] for a in range(n)]
    G=[[q*gram[a][b]-sums[a]*sums[b]-sum(V[a][j]*V[b][j] for j in range(n))
        +u[b]*t[a]+u[a]*t[b]+m*u[a]*u[b] for b in range(n)] for a in range(n)]
    covariance=[[F(G[a][b],m)-F(h[a]*h[b],m*m) for b in range(n)] for a in range(n)]
    T=target(S,c,6)
    edits=[[target(S,sorted((set(c)-{a})|{b}),6)-T for b in outside] for a in c]
    for a in range(n):
        assert sum(edits[a])==h[a]
        for b in range(n):assert sum(x*y for x,y in zip(edits[a],edits[b]))==G[a][b]
    row_hist=sorted(sorted(Counter(row).items()) for row in edits)
    column_hist=sorted(sorted(Counter(column).items()) for column in zip(*edits))
    diag=sorted(covariance[a][a] for a in range(n))
    return {'q':q,'n':n,'selected':c,'target':T,'outside':outside,'delta_matrix':edits,
            'derivative_gram':gram,'insertion_cross_second_sum':G,
            'insertion_covariance':[[enc(x) for x in row] for row in covariance],
            'flat_histogram':sorted(Counter(x for row in edits for x in row).items()),
            'unlabelled_row_histograms':row_hist,'unlabelled_column_histograms':column_hist,
            'diagonal_covariance_deck':[enc(x) for x in diag],
            'covariance_trace_squared':enc(sum(covariance[a][b]**2 for a in range(n) for b in range(n))),
            'literal_neighbour_targets':n*m,'exact_cross_moment_entries':n*n}


def main():
    source=LAB/'round7/local_edits/orbit_verification.json'
    twins=json.loads(source.read_text())
    rows=[]
    for case in twins['cases']:
        for pair in case['non_affine_pairs']:
            members=[analyze(p['selected'],case['q']) for p in pair['profiles']]
            comparisons={key:members[0][key]==members[1][key] for key in (
                'flat_histogram','unlabelled_row_histograms','unlabelled_column_histograms',
                'diagonal_covariance_deck','covariance_trace_squared')}
            rows.append({'members':members,'equal_features':comparisons})
    out={'status':'passed','scope':'Actual saved one-swap twins; covariance for edits sharing an insertion, compared with literal targets. Finite diagnostic only.',
         'pairs':rows,'source_sha256':{'round8/shared_insertions/preflight.py':sha256(Path(__file__).read_bytes()).hexdigest(),
                                     'round7/local_edits/review_results.py':sha256((LAB/'round7/local_edits/review_results.py').read_bytes()).hexdigest()},
         'input_sha256':{'round7/local_edits/orbit_verification.json':sha256(source.read_bytes()).hexdigest()}}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([r['equal_features'] for r in rows]),flush=True)


if __name__=='__main__':main()
