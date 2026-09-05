#!/usr/bin/env python3
"""FFT Paley matvec, one-sided Ritz spaces, joint translation witnesses.

Only numerical LOWER witnesses: no global upper eigenvalue certification.
"""
import argparse,hashlib,json,math,struct,sys,time,resource
from pathlib import Path
import numpy as np
import scipy
from scipy.fft import rfft,irfft
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
HERE=Path(__file__).resolve().parent
SEED=20260905

def prime(p):return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))

class PaleyFFT(LinearOperator):
    def __init__(self,p):
        assert prime(p) and p%4==1
        self.p=p;self.chi=np.full(p,-1,dtype=np.int8);self.chi[0]=0
        roots=np.arange(1,(p+1)//2,dtype=np.int64);self.chi[roots*roots%p]=1
        self.C=np.flatnonzero((self.chi==1)&(np.roll(self.chi,1)==1));self.m=len(self.C)
        self.multiplier=self.chi[:p//2+1].astype(float);self.calls=0;self.matvec_seconds=0.
        super().__init__(dtype=np.dtype('float64'),shape=(self.m,self.m))
    def _matvec(self,v):
        start=time.monotonic();v=np.asarray(v).reshape(-1);z=np.zeros(self.p);z[self.C]=v
        answer=irfft(rfft(z)*self.multiplier,n=self.p)[self.C]-v.sum()/self.p
        self.calls+=1;self.matvec_seconds+=time.monotonic()-start
        return answer
    def dense(self):
        assert self.p<=4001,'dense guard prevents accidental large allocation'
        return self.chi[(self.C[:,None]-self.C[None,:])%self.p]/math.sqrt(self.p)-1/self.p
    def overlap_matrix(self,U,t):
        full=np.zeros((self.p,U.shape[1]));full[self.C]=U
        M=U.T@full[(self.C+int(t))%self.p]
        return (M+M.T)/2
    def autocorrelation(self,v):
        z=np.zeros(self.p);z[self.C]=v
        freq=rfft(z);return irfft((freq*freq.conj()).real,n=self.p)

def joint(A,U,seed):
    """Fix shifts/signs from first Ritz vector, then maximize their mean
    signed overlap inside ONE-SIDED computed span. No sign mixing.
    """
    corr=A.autocorrelation(U[:,0]);order=np.argsort(np.abs(corr[1:A.p//2+1]))[::-1]+1
    train=order[:4];holdout=order[4:8]
    signs=np.where(corr[train]>=0,1.,-1.)
    K=sum(float(s)*A.overlap_matrix(U,int(t)) for t,s in zip(train,signs))/len(train)
    vals,vecs=np.linalg.eigh(K);v=U@vecs[:,-1]
    Hv=A@v;ray=float(v@Hv)
    Cnew=A.autocorrelation(v)
    rng=np.random.default_rng(seed);random_shifts=rng.choice(np.arange(1,A.p//2+1),size=4,replace=False)
    randomK=sum(A.overlap_matrix(U,int(t)) for t in random_shifts)/4
    rv,rU=np.linalg.eigh(randomK)
    # Equal-spectrum null: coordinate permutation of the entire one-sided span.
    perm=rng.permutation(A.m);Up=U[perm];cp=A.autocorrelation(Up[:,0]);op=np.argsort(np.abs(cp[1:A.p//2+1]))[::-1]+1;tp=op[:4];sp=np.where(cp[tp]>=0,1.,-1.)
    kp=sum(float(s)*A.overlap_matrix(Up,int(t)) for t,s in zip(tp,sp))/4
    result={'train_shifts':train.tolist(),'fixed_alignment_signs':signs.astype(int).tolist(),
        'heldout_shifts':holdout.tolist(),'seed_train_overlap':float(np.mean(signs*corr[train])),
        'joint_train_overlap':float(vals[-1]),'seed_heldout_abs_overlap':float(np.mean(abs(corr[holdout]))),
        'joint_heldout_abs_overlap':float(np.mean(abs(Cnew[holdout]))),
        'joint_rayleigh_H':ray,'joint_vector_residual':float(np.linalg.norm(Hv-ray*v)),
        'random_fixed_shifts_max_signed_average':float(rv[-1]),
        'permuted_same_spectrum_joint_average':float(np.linalg.eigvalsh(kp)[-1]),
        'method':'fixed four shifts/signs selected from first Ritz vector; one optimizer for their average; four separate heldout shifts'}
    return result,v

def complete_symmetry(A,U):
    """Known S3 action completes symmetry partners omitted by single-start
    Lanczos. Thin Gram SVD discards numerical linear dependencies.
    """
    idx=np.full(A.p,-1,dtype=np.int64);idx[A.C]=np.arange(A.m)
    rp=idx[(1-A.C)%A.p]
    ip=idx[np.array([pow(int(x),-1,A.p) for x in A.C])]
    perms=[np.arange(A.m),rp,ip,rp[ip],ip[rp],rp[ip[rp]]]
    W=np.column_stack([U[g] for g in perms]);G=W.T@W
    ev,V=np.linalg.eigh(G);keep=ev>1e-9*max(ev)
    Q=W@(V[:,keep]/np.sqrt(ev[keep])[None,:]);Q,_=np.linalg.qr(Q)
    return Q

def export_witness(A,v,side):
    if v[np.argmax(abs(v))]<0:v=-v
    z=np.rint(v/np.max(abs(v))*1024).astype('<i4')
    path=HERE/f'witness_{A.p}_{side}.bin'
    pairs=np.column_stack((A.C.astype('<i4'),z)).astype('<i4')
    path.write_bytes(struct.pack('<II',A.p,A.m)+pairs.tobytes())
    norm=int(np.dot(z.astype(np.int64),z.astype(np.int64)));sm=int(z.astype(np.int64).sum())
    zv=z.astype(float)/math.sqrt(norm);quot=float(zv@(A@zv))
    return {'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'format':'little-endian uint32 p,m then m pairs (int32 vertex,int32 coefficient)',
            'norm_squared':norm,'sum_entries':sm,'sum_absolute_entries':int(abs(z.astype(np.int64)).sum()),
            'rounded_rayleigh_H_float':quot,'quantization_max_abs':1024}

def run(p,k,tol):
    start=time.monotonic();A=PaleyFFT(p);rng=np.random.default_rng(SEED+p);v0=rng.standard_normal(A.m)
    result={'p':p,'m':A.m,'k_requested_per_side':k,'eigensolver_tolerance':tol,'dense_validation':None,'sides':[]}
    D=A.dense() if p<=4001 else None
    if D is not None:
        v=rng.standard_normal(A.m);dv=D@v;fv=A@v
        allvals=np.linalg.eigvalsh(D)
        result['dense_validation']={'matvec_max_abs_error':float(np.max(abs(dv-fv))),
              'dense_min_eigenvalue':float(allvals[0]),'dense_max_eigenvalue':float(allvals[-1])}
        assert np.max(abs(dv-fv))<1e-11
    for sign,which,name in [(1,'LA','positive'),(-1,'SA','negative')]:
        side_start=time.monotonic();before=A.calls;converged=True
        try:ev,U=eigsh(A,k=k,which=which,tol=tol,maxiter=1600,ncv=max(24,5*k+1),v0=v0)
        except ArpackNoConvergence as e:
            ev,U=e.eigenvalues,e.eigenvectors;converged=False
            if len(ev)==0:raise
        order=np.argsort(sign*ev)[::-1];ev,U=ev[order],U[:,order]
        # Reorthogonalize and Rayleigh-Ritz within computed span to avoid
        # reliance on exact eigensolver orthogonality/invariance.
        U,_=np.linalg.qr(U);HU=np.column_stack([A@U[:,j] for j in range(U.shape[1])]);G=(U.T@HU+HU.T@U)/2
        eg,V=np.linalg.eigh(G);order=np.argsort(sign*eg)[::-1];eg,V=eg[order],V[:,order];U=U@V;HU=HU@V
        raw_eg=eg.copy();raw_U=U.copy()
        raw_joint,_=joint(A,raw_U,SEED+p+sign)
        U=complete_symmetry(A,U)
        HU=np.column_stack([A@U[:,j] for j in range(U.shape[1])]);G=(U.T@HU+HU.T@U)/2
        eg,V=np.linalg.eigh(G);keep=sign*eg>=min(sign*raw_eg)-1e-6
        eg,V=eg[keep],V[:,keep];order=np.argsort(sign*eg)[::-1];eg,V=eg[order],V[:,order]
        U,HU=U@V,HU@V
        residuals=np.linalg.norm(HU-U*eg[None,:],axis=0)
        jrec,v=joint(A,U,SEED+p+sign)
        # Actual directional Rayleigh is between the compressed eigenvalues.
        directional=sign*jrec['joint_rayleigh_H'];floor=float(min(sign*eg));ceil=float(max(sign*eg))
        assert floor-1e-10<=directional<=ceil+1e-10
        row={'side':name,'solver_reported_convergence':converged,'ritz_eigenvalues':eg.tolist(),'raw_lanczos_eigenvalues':raw_eg.tolist(),'raw_joint':raw_joint,'symmetry_completed_rank':len(eg),
             'residual_norms':residuals.tolist(),'orthogonality_error':float(np.linalg.norm(U.T@U-np.eye(U.shape[1]))),
             'compressed_directional_floor':floor,'compressed_directional_ceiling':ceil,
             'matvec_calls':A.calls-before,'elapsed_seconds':time.monotonic()-side_start,
             'joint':jrec,'witness':export_witness(A,v,name)}
        if D is not None:
            row['dense_eigenvalue_max_abs_error']=float(max(min(abs(allvals-x)) for x in eg))
            row['dense_extremal_error']=float(abs(eg[0]-(allvals[-1] if sign==1 else allvals[0])))
            target=allvals[-k:][::-1] if sign==1 else allvals[:k]
            row['raw_ordered_multiplicity_error']=float(np.max(abs(raw_eg-target)))
            assert row['dense_eigenvalue_max_abs_error']<1e-7 and row['dense_extremal_error']<1e-7
        result['sides'].append(row)
        print('side',p,name,'ritz',eg.tolist(),'residual',max(residuals),'matvec',row['matvec_calls'],'sec',row['elapsed_seconds'],flush=True)
        (HERE/f'results_{p}.json').write_text(json.dumps(result,indent=2)+'\n')
    result.update({'elapsed_seconds':time.monotonic()-start,'total_matvec_calls':A.calls,
        'matvec_seconds':A.matvec_seconds,'process_peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'dense_matrix_bytes_avoided':8*A.m*A.m if D is None else 0,
        'status':'lower witnesses and numerical Ritz estimates only; no spectral upper certificate'})
    (HERE/f'results_{p}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--primes',nargs='+',type=int,default=[101,401,1009,4001,65537]);ap.add_argument('--k',type=int,default=4);ap.add_argument('--tol',type=float,default=1e-8);args=ap.parse_args()
    runtime={'python':sys.executable,'version':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'seed':SEED,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
    for p in args.primes:run(p,args.k,args.tol)
if __name__=='__main__':main()
