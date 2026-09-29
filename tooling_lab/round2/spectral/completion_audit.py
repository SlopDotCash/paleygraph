#!/usr/bin/env python3
"""Ablation with the same four translation shifts/signs before and after S3
completion. This removes the top-eigenbasis reselection confound in main logs.
"""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
from fft_transport import HERE,PaleyFFT,joint,complete_symmetry,SEED

def main():
    rows=[]
    for p in [101,401,1009,4001]:
        A=PaleyFFT(p);v0=np.random.default_rng(SEED+p).standard_normal(A.m)
        for sign,which in [(1,'LA'),(-1,'SA')]:
            ev,U=eigsh(A,k=4,which=which,tol=1e-8,ncv=24,v0=v0)
            order=np.argsort(sign*ev)[::-1];U=U[:,order]
            base,_=joint(A,U,SEED+p+sign);Q=complete_symmetry(A,U)
            shifts=base['train_shifts'];signs=base['fixed_alignment_signs']
            K=sum(s*A.overlap_matrix(Q,t) for t,s in zip(shifts,signs))/4
            e,V=np.linalg.eigh(K);v=Q@V[:,-1];q=float(v@(A@v))
            Hq=np.column_stack([A@Q[:,j] for j in range(Q.shape[1])]);G=Q.T@Hq
            min_directional=float(min(np.linalg.eigvalsh(sign*(G+G.T)/2)))
            containment_error=float(np.linalg.norm(U-Q@(Q.T@U)))
            assert e[-1]>=base['joint_train_overlap']-4*containment_error-1e-12
            assert sign*q>=min_directional-1e-10
            rows.append({'p':p,'side':'positive' if sign==1 else 'negative','raw_rank':4,'completed_rank':Q.shape[1],
                'fixed_shifts':shifts,'fixed_signs':signs,'raw_optimum':base['joint_train_overlap'],
                'completed_optimum':float(e[-1]),'gain':float(e[-1]-base['joint_train_overlap']),'raw_span_containment_error':containment_error,
                'completed_witness_directional_rayleigh':sign*q,
                'completed_subspace_min_directional_rayleigh':min_directional})
    output={'status':'same-objective numerical subspace ablation; no asymptotic inference','cases':rows,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'completion_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    for r in rows:print(r['p'],r['side'],r['completed_rank'],r['gain'])
if __name__=='__main__':main()
