#!/usr/bin/env python3
"""Exact covariance of deletion choices sharing a uniform outside insertion."""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
from time import perf_counter

HERE=Path(__file__).resolve().parent


def enc(x):return [x.numerator,x.denominator]


def finish(g):
    q,n=g['q'],g['n'];m=q-n;V=g['internal_contraction'];R=g['derivative_gram']
    assert len(V)==len(R)==n and all(len(row)==n for row in V+R)
    sums=g['insertion_sums'];u=[V[a][a] for a in range(n)];t=list(map(sum,V));h=[-t[a]-m*u[a] for a in range(n)]
    for a in range(n):
        assert R[a][a]==g['insertion_norms_squared'][a]
        for b in range(n):assert R[a][b]==R[b][a]
    G=[[q*R[a][b]-sums[a]*sums[b]-sum(V[a][j]*V[b][j] for j in range(n))
        +u[b]*t[a]+u[a]*t[b]+m*u[a]*u[b] for b in range(n)] for a in range(n)]
    # Keep a single exact denominator m^2, avoiding cancellation in floats.
    C=[[m*G[a][b]-h[a]*h[b] for b in range(n)] for a in range(n)]
    for a in range(n):
        assert C[a][a]>=0
        for b in range(n):assert C[a][b]*C[a][b]<=C[a][a]*C[b][b]
    norm_squared=sum(x*x for row in C for x in row)
    diagonal_squared=sum(C[a][a]**2 for a in range(n))
    drift=F(sum(h),n*m);marginal_var=F(sum(G[a][a] for a in range(n)),n*m)-drift*drift
    return {'q':q,'n':n,'degree':g['degree'],'selected':g['selected'],'target':g['target'],
            'neighbour_count':n*m,'shared_insertions':m,'internal_contraction':V,'derivative_gram':R,
            'insertion_cross_second_sum':G,'insertion_delta_sums':h,'covariance_numerator':C,'covariance_denominator':m*m,
            'covariance_trace_squared':enc(F(norm_squared,m**4)),
            'off_diagonal_squared_energy_fraction':enc(F(norm_squared-diagonal_squared,norm_squared)) if norm_squared else [0,1],
            'marginal_neighbour_mean':enc(g['target']+drift),'marginal_neighbour_variance':enc(marginal_var),
            'scope':'Joint insertion covariance at one actual input; no uniform estimate or target-free predictor.'}


def query(q,c,degree=6):
    c=list(c)
    if type(q)is not int or type(degree)is not int or any(type(x)is not int for x in c):raise ValueError('integer parameters required')
    if not(5<=q<=10000000 and 1<=len(c)<=64 and len(c)<q and 0<=degree<=min(6,len(c))):raise ValueError('prototype bounds exceeded')
    if len(set(c))!=len(c) or any(x<0 or x>=q for x in c):raise ValueError('invalid selected set')
    start=perf_counter()
    run=subprocess.run([str(HERE/'covariance_backend')],input=f'{q} {len(c)} {degree}\n'+' '.join(map(str,c))+'\n',text=True,capture_output=True)
    if run.returncode:raise ValueError(run.stderr.strip())
    g=json.loads(run.stdout);g['selected']=sorted(c)
    assert (g['q'],g['n'],g['degree'])==(q,len(c),degree)
    result=finish(g);result['seconds']=perf_counter()-start
    return result,g
