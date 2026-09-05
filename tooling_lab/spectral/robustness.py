#!/usr/bin/env python3
"""Post-hoc robustness audit. No thresholds are selected using these results.

Existing pilot/results/witness files are read-only. Affine controls move both
anchors and coordinates, preserving the same graph and Rayleigh witness.
"""
import json,math,hashlib
from fractions import Fraction
from pathlib import Path
import numpy as np
from transport_microscope import Field,HERE,analyze

THRESHOLDS=[(1,3),(1,2),(3,5),(2,3),(3,4)]

def inv(F,a):
    assert a
    if F.degree==1:return pow(int(a),-1,F.p)
    p=F.p;x=int(a)%p;y=int(a)//p;den=pow((x*x-F.nu*y*y)%p,-1,p)
    return x*den%p+p*((-y*den)%p)

def closure(F,T):
    nz=[int(x) for x in T if x]
    if not nz:return [0]
    if F.degree==1:return list(range(F.q))
    p=F.p;a=nz[0];u,v=a%p,a//p
    if any((u*(b//p)-v*(b%p))%p for b in nz):return list(range(F.q))
    return sorted(int(F.mul(k,a)) for k in range(p))

def structure(F,T,anchor_difference=1):
    A=closure(F,T);aset=set(A)
    isfield=1 in aset and len(A)>1 and all(int(F.mul(x,y)) in aset for x in A for y in A) if len(A)<F.q else True
    norm=sorted(int(F.mul(x,inv(F,anchor_difference))) for x in A)
    # Compare exact anchored canonical prime-subfield, not merely cardinality.
    normalized_subfield=(F.degree==2 and norm==list(range(F.p)))
    return {'translation_count':len(T),'additive_closure_size':len(A),
            'proper_subfield_in_current_coordinates':bool(1<len(A)<F.q and isfield),
            'anchor_normalized_prime_subfield_recovered':normalized_subfield}

def correlations(F,C,z):
    full=np.zeros(F.q,dtype=np.int64);full[C]=z
    return [int(z@full[F.add(C,t)]) for t in range(F.q)]

def main():
    threshold_rows=[];affine_rows=[]
    sources=json.loads((HERE/'results.json').read_text())['cases']
    for row in sources:
        w=json.loads((HERE/row['certificate_file']).read_text());F=Field(w['characteristic'],w['degree']);C=np.array(w['C']);z=np.array(w['integer_vector'],dtype=np.int64);n=w['norm_squared']
        corr=correlations(F,C,z);rr=[]
        for num,den in THRESHOLDS:
            T=[t for t,c in enumerate(corr) if den*abs(c)>=num*n]
            rr.append({'threshold':[num,den],**structure(F,T)})
        threshold_rows.append({'q':F.q,'thresholds':rr})
        # Translation plus square scaling transports an actual Paley neighborhood.
        # For extension fields explicitly use a square outside the base subfield.
        if F.degree==2:
            a=next(t for t in range(F.p,F.q) if F.chi[t]==1)
        else:a=next(t for t in range(2,F.q) if F.chi[t]==1)
        b=3;Cp=F.add(F.mul(C,a),b);newcorr=correlations(F,Cp,z)
        assert all(newcorr[int(F.mul(t,a))]==corr[t] for t in range(F.q))
        S=F.chi[F.sub(C[:,None],C[None,:])];Sp=F.chi[F.sub(Cp[:,None],Cp[None,:])]
        assert np.array_equal(S,Sp)
        oldT=[t for t,c in enumerate(corr) if 5*abs(c)>=3*n];newT=[t for t,c in enumerate(newcorr) if 5*abs(c)>=3*n]
        before=structure(F,oldT);after=structure(F,newT,a)
        assert before['anchor_normalized_prime_subfield_recovered']==after['anchor_normalized_prime_subfield_recovered']
        affine_rows.append({'q':F.q,'scale':a,'offset':b,'new_anchors':[b,int(F.add(a,b))],
            'exact_translation_covariance_checks':F.q,'exact_character_entries_checked':len(C)**2,
            'before':before,'after':after})
        print('threshold/affine',F.q,flush=True)
    window_rows=[]
    for F in [Field(101),Field(401),Field(17,2),Field(31,2)]:
        for window in [.01,.025,.05]:
            r=analyze(F,20260905,window,save_vectors=False)
            window_rows.append({k:r[k] for k in ['q','tail_window','tail_rank','radius_H','tail_subspace','tail_witness_rayleigh_H']})
            print('window',F.q,window,r['tail_rank'],r['tail_subspace']['max_abs_overlap'],flush=True)
    out={'status':'post-hoc robustness; no held-out validation or theorem',
         'fixed_witness_thresholds':threshold_rows,'affine_coordinate_controls':affine_rows,
         'tail_window_controls':window_rows,
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'robustness.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
