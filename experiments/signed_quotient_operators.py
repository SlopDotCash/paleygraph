#!/usr/bin/env python3
"""Exact signed quotient matrices for subgroup and mixed-product periods.

All acceptance checks use integers. Sparse matrix products use int64 only
under an explicit maximum-row-sum bound. This verifies finite identities,
not a uniform bound on the signed spectrum.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json

import numpy as np
import scipy
from scipy.sparse import csr_matrix, identity

from kernel_discriminant import generator, prime

ROOT = Path(__file__).resolve().parents[1]
LIMIT = np.iinfo(np.int64).max


def orbit_data(p, n, g):
    assert prime(p) and n >= 4 and n & (n-1) == 0 and (p-1) % n == 0
    h = [pow(g,j,p) for j in range(n)]
    assert len(set(h)) == n and pow(g,n//2,p) == p-1
    labels = [-1]*p
    signs = [0]*p
    reps = []
    for a in range(1,p):
        if labels[a] < 0:
            index = len(reps)
            reps.append(a)
            for j,u in enumerate(h):
                labels[a*u%p] = index
                signs[a*u%p] = (-1)**j
    assert len(reps) == (p-1)//n
    assert all(signs[x] in (-1,1) for x in range(1,p))
    assert all(labels[g*x%p] == labels[x] and signs[g*x%p] == -signs[x]
               for x in range(1,p))
    return h, reps, labels, signs


def matrices(p, reps, labels, signs, weights):
    rows, cols, signed, unsigned = [], [], [], []
    for i,a in enumerate(reps):
        aa, bb = Counter(), Counter()
        for z,weight in weights.items():
            y = (a-z)%p
            if y:
                j = labels[y]
                aa[j] += weight*signs[y]
                bb[j] += weight
        for j,value in bb.items():
            rows.append(i)
            cols.append(j)
            signed.append(aa[j])
            unsigned.append(value)
    shape = (len(reps),len(reps))
    a = csr_matrix((np.array(signed,dtype=np.int64),(rows,cols)),shape=shape)
    b = csr_matrix((np.array(unsigned,dtype=np.int64),(rows,cols)),shape=shape)
    a.eliminate_zeros()
    assert (a-a.T).nnz == 0 and (b-b.T).nnz == 0
    assert np.all(b.data>0)
    absolute = a.copy()
    absolute.data = np.abs(absolute.data)
    delta = b-absolute
    assert np.all(delta.data>=0)
    return a,b,absolute


def sum_entries(a):
    return sum(map(int,a.data))


def squared_entries(a):
    return sum(int(x)**2 for x in a.data)


def dense_convolve(p, counts, weights):
    following = [0]*p
    for x,count in enumerate(counts):
        if count:
            for z,weight in weights.items():
                following[(x+z)%p] += count*weight
    return following


def matrix_case(p,n,g,kind,depth):
    h,reps,labels,signs = orbit_data(p,n,g)
    k=n//2
    if kind == 'subgroup':
        weights=Counter({x:1 for x in h})
    else:
        kk=[pow(g,2*j,p) for j in range(k)]
        ll=[g*x%p for x in kk]
        weights=Counter((x+y)%p for x in kk for y in ll)
    assert 0 not in weights
    assert all(weights[g*x%p] == value and weights[-x%p] == value for x,value in weights.items())
    d=sum(weights.values())
    energy=sum(x*x for x in weights.values())
    m=len(reps)
    assert d % n == 0 and n*sum(weights.get(x,0) for x in reps) == d
    a,b,absolute=matrices(p,reps,labels,signs,weights)
    assert sum_entries(b)==m*d-d//n
    assert squared_entries(a)==(p*energy-d*d)//n
    assert squared_entries(b)-squared_entries(a)==d*d-2*energy
    loss=sum_entries(b)-sum_entries(absolute)
    assert 0<=2*loss<=d*d-2*energy
    lower_numerator=m*d-d//n-(d*d-2*energy)//2
    assert sum_entries(absolute)>=lower_numerator
    assert int(np.asarray(absolute.sum(axis=1)).max())<=d
    # Exact diagonal sign conjugation, reconstructed from changed field reps.
    switch=[(-1 if i%3==1 else 1) for i in range(m)]
    changed_reps=[g*a0%p if switch[i]<0 else a0 for i,a0 in enumerate(reps)]
    changed_signs=[0]+[switch[labels[x]]*signs[x] for x in range(1,p)]
    changed_a,changed_b,_=matrices(p,changed_reps,labels,changed_signs,weights)
    ss=csr_matrix((np.array(switch,dtype=np.int64),(range(m),range(m))),shape=(m,m))
    assert (changed_a-ss@a@ss).nnz==0 and (changed_b-b).nnz==0
    # Compare matrix traces with independent dense additive convolution.
    power=identity(m,dtype=np.int64,format='csr')
    counts=[1]+[0]*(p-1)
    child_counts=[1]+[0]*(p-1) if kind=='product' else None
    child_weights=Counter({x:1 for x in kk}) if kind=='product' else None
    moments=[]
    for r in range(1,depth+1):
        # Every absolute row sum in A^r is at most d^r, including each
        # partial dot-product sum, so this bounds int64 arithmetic.
        assert d**r<=LIMIT
        power=power@a
        power.eliminate_zeros()
        counts=dense_convolve(p,counts,weights)
        if child_counts is not None:
            child_counts=dense_convolve(p,child_counts,child_weights)
            balanced=sum(child_counts[x]*child_counts[g*x%p] for x in range(p))
            assert counts[0]==balanced
        assert sum(counts)==d**r
        trace=sum(map(int,power.diagonal()))
        numerator=p*counts[0]-d**r
        assert numerator%n==0 and trace==numerator//n
        assert int(np.asarray(abs(power).sum(axis=1)).max())<=d**r
        moments.append({'r':r,'weighted_zero_count':counts[0],
                        'trace':trace,'int64_row_bound':d**r})
    # Enumerate rooted quotient walks of the subgroup graph in two small
    # fields and record their multiplicative endpoint ratio (holonomy).
    holonomy=None
    if p in (17,97):
        holonomy=[]
        count=[1]+[0]*(p-1)
        for r in range(1,5):
            count=dense_convolve(p,count,weights)
            observed={u:sum(count[(a*(1-u))%p] for a in reps) for u in h}
            expected_nonidentity=(d**r-count[0])//n
            assert (d**r-count[0])%n==0
            assert observed[1]==m*count[0]
            assert all(observed[u]==expected_nonidentity for u in h if u!=1)
            signed_trace=sum((-1)**j*observed[u] for j,u in enumerate(h))
            assert signed_trace==moments[r-1]['trace']
            holonomy.append({'r':r,'identity_count':observed[1],
                             'each_nonidentity_count':expected_nonidentity,
                             'signed_trace':signed_trace})
    quartic = n**4 <= 4*p <= 4*n**4
    quartic_bound=None
    if quartic:
        if kind=='subgroup':
            assert n*lower_numerator >= m*(n*n-2)
            quartic_bound={'rho_absolute_at_least':f'{n} - 2/{n}'}
        else:
            assert 4*lower_numerator>=m*(4*k*k-k)
            quartic_bound={'rho_absolute_at_least':f'{k*k} - {k}/4'}
    data={'p':p,'n':n,'generator':g,'kind':kind,'dimension':m,
          'degree':d,'weight_squared_mass':energy,
          'signed_nonzero_entries':a.nnz,'unsigned_nonzero_entries':b.nnz,
          'maximum_absolute_signed_entry':int(abs(a.data).max()),
          'signed_frobenius_squared':squared_entries(a),
          'unsigned_frobenius_squared':squared_entries(b),
          'unsigned_entry_sum':sum_entries(b),'absolute_signed_entry_sum':sum_entries(absolute),
          'entry_cancellation_loss':loss,'loss_upper_bound':(d*d-2*energy)//2,
          'absolute_spectral_lower_bound_numerator':lower_numerator,
          'absolute_spectral_lower_bound_denominator':m,
          'quartic_window':quartic,'quartic_absolute_bound':quartic_bound,
          'gauge_reconstruction_checked':True,
          'independent_child_count_check':kind=='product','moments':moments,
          'holonomy_checks':holonomy}
    if m<=12:
        data['signed_matrix']=a.toarray().tolist()
        data['unsigned_matrix']=b.toarray().tolist()
    print(f'p={p}, n={n}, {kind}: {m}x{m}, {depth} exact traces passed',flush=True)
    return data


def main():
    cases=[]
    for p,n,depth in [(17,8,8),(17,16,8),(97,8,8),(97,32,8),
                       (1049,8,6),(2017,8,6),(17393,16,4)]:
        g=generator(p,n)
        for kind in ('subgroup','product'):
            checked_depth = 6 if n == 32 and kind == 'product' else depth
            cases.append(matrix_case(p,n,g,kind,checked_depth))
    sources=['experiments/signed_quotient_operators.py','experiments/kernel_discriminant.py',
             'experiments/cyclotomic_norm_audit.py']
    output={'status':'Exact finite identities; no uniform signed spectral estimate',
            'numpy_version':np.__version__,'scipy_version':scipy.__version__,
            'cases':cases,'matrix_cases':len(cases),
            'exact_trace_checks':sum(len(x['moments']) for x in cases),
            'exact_holonomy_checks':sum(len(x['holonomy_checks'] or []) for x in cases),
            'independent_balanced_count_checks':sum(len(x['moments']) for x in cases if x['kind']=='product'),
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}}
    (ROOT/'results/signed_quotient_operators.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':
    main()
