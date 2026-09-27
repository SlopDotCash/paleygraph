#!/usr/bin/env python3
"""Separate binomial/row-count review of the new derivative Gram data.

Full new matrices through q65537, twelve exact off-diagonal entries on each
larger input, and every Gram row sum through the deletion Euler identity.
This does not claim a complete second large-matrix computation.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random
import sys
from time import perf_counter
import numpy as np

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'round7/local_edits'))
from local_edits import coefficient_table, prime_character


def check_algebra(r):
    q,n=r['q'],r['n'];m=q-n;R=r['derivative_gram'];V=r['internal_contraction']
    C=r['covariance_numerator'];h=r['insertion_delta_sums'];G=r['insertion_cross_second_sum']
    assert r['covariance_denominator']==m*m
    for a in range(n):
        assert h[a]==-sum(V[a])-m*V[a][a]
        for b in range(n):
            assert C[a][b]==m*G[a][b]-h[a]*h[b]
            assert C[a][b]==C[b][a]
    assert F(*r['covariance_trace_squared'])==F(sum(x*x for row in C for x in row),m**4)
    mean=F(sum(h),n*m)+r['target'];var=F(sum(G[a][a] for a in range(n)),n*m)-(mean-r['target'])**2
    assert F(*r['marginal_neighbour_mean'])==mean and F(*r['marginal_neighbour_variance'])==var


def audit(r):
    start=perf_counter();q,n=r['q'],r['n'];c=np.array(r['selected'],dtype=np.int64)
    chi=prime_character(q);table=coefficient_table(n,5)
    bound=comb(n-1,5);chunk=min(8192,((1<<63)-1)//max(1,n*bound*bound))
    rng=random.Random(83100+q);pairs={(0,1),(0,n-1),(n-2,n-1)}
    while len(pairs)<12:
        a,b=sorted(rng.sample(range(n),2));pairs.add((a,b))
    pairs=sorted(pairs);samples=[0]*len(pairs);row_sums=[0]*n;deriv_sums=[0]*n
    full=[[0]*n for _ in range(n)] if q<=65537 else None
    euler_rows=0
    for start_row in range(0,q,chunk):
        x=np.arange(start_row,min(q,start_row+chunk),dtype=np.int64)
        signs=chi[(x[:,None]-c)%q];positive=(signs==1).sum(axis=1);nonzero=(signs!=0).sum(axis=1)
        values=table[nonzero[:,None]-(signs!=0),positive[:,None]-(signs==1)]
        summed=(n-5)*table[nonzero,positive]
        assert np.array_equal(summed,values.sum(axis=1));euler_rows+=len(x)
        totals=(values*summed[:,None]).sum(axis=0)
        for a,total in enumerate(totals):row_sums[a]+=int(total)
        for a,total in enumerate(values.sum(axis=0)):deriv_sums[a]+=int(total)
        if full is not None:
            product=values.T@values
            for a in range(n):
                for b in range(n):full[a][b]+=int(product[a,b])
        else:
            for j,(a,b) in enumerate(pairs):samples[j]+=int((values[:,a]*values[:,b]).sum())
    R=r['derivative_gram'];V=r['internal_contraction'];m=q-n
    assert row_sums==list(map(sum,R))
    if full is not None:assert full==R
    else:
        for value,(a,b) in zip(samples,pairs):assert value==R[a][b]
    # Reconstruct every shared-insertion cross moment from the reviewed ledger.
    for a in range(n):
        for b in range(n):
            u,v=V[a][a],V[b][b]
            expected=q*R[a][b]-deriv_sums[a]*deriv_sums[b]-sum(V[a][j]*V[b][j] for j in range(n))+v*sum(V[a])+u*sum(V[b])+m*u*v
            assert expected==r['insertion_cross_second_sum'][a][b]
    check_algebra(r)
    return {'q':q,'n':n,'family':r['family'],'new_gram_full_entries_checked':n*n if full is not None else 0,
            'new_gram_sampled_off_diagonal_entries_checked':len(pairs) if full is None else 0,
            'sampled_pairs':pairs if full is None else [],'exact_gram_row_sums_checked':n,
            'pointwise_euler_rows_checked':euler_rows,'cross_moment_algebra_entries_checked':n*n,'seconds':perf_counter()-start}


def main():
    data=json.loads((HERE/'scale_results.json').read_text())
    for group in ('source_sha256','input_sha256'):
        for name,value in data[group].items():assert sha256((LAB/name).read_bytes()).hexdigest()==value
    rows=[]
    for r in data['scale_cases']:
        row=audit(r);rows.append(row);print(json.dumps(row),flush=True)
    out={'status':'passed','scope':__doc__,'cases':rows,
         'source_sha256':{'round8/shared_insertions/stream_review.py':sha256(Path(__file__).read_bytes()).hexdigest(),
                          'round7/local_edits/local_edits.py':sha256((LAB/'round7/local_edits/local_edits.py').read_bytes()).hexdigest()},
         'input_sha256':{'round8/shared_insertions/scale_results.json':sha256((HERE/'scale_results.json').read_bytes()).hexdigest()}}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
