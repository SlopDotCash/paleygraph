#!/usr/bin/env python3
"""Edge-conditioned additive transport; deterministic discovery, exact witnesses.

No unproved spectral estimate or novelty theorem is asserted. All mutations and
outputs are confined to this directory. NumPy/SciPy are the only dependencies.
"""
from __future__ import annotations
import argparse, json, math, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

HERE=Path(__file__).resolve().parent

def prime(p):
    return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))

class Field:
    def __init__(self,p,degree=1):
        assert prime(p) and p%2==1 and degree in (1,2)
        self.p,self.degree,self.q=p,degree,p**degree
        self.base=np.array([0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)],dtype=np.int64)
        self.nu=next(x for x in range(2,p) if self.base[x]==-1)
        self.x=np.arange(self.q,dtype=np.int64)
        self.chi=self.base[self.x] if degree==1 else self.base[((self.x%p)**2-self.nu*(self.x//p)**2)%p]
    def sub(self,x,y):
        if self.degree==1:return (x-y)%self.p
        return ((x%self.p-y%self.p)%self.p)+self.p*((x//self.p-y//self.p)%self.p)
    def add(self,x,y):
        if self.degree==1:return (x+y)%self.p
        return ((x%self.p+y%self.p)%self.p)+self.p*((x//self.p+y//self.p)%self.p)
    def mul(self,x,y):
        if self.degree==1:return x*y%self.p
        p=self.p;a=x%p;b=x//p;c=y%p;d=y//p
        return ((a*c+self.nu*b*d)%p)+p*((a*d+b*c)%p)
    def instance(self):
        assert self.q%4==1
        C=self.x[(self.chi==1)&(self.chi[self.sub(self.x,1)]==1)]
        S=self.chi[self.sub(C[:,None],C[None,:])]
        H=S/math.sqrt(self.q)-np.ones(S.shape)/self.q
        return C,S,H

def transport_profile(F,C,U):
    """max |<f,tau_t f>| over unit f in span U, each nonzero t.

U columns orthonormal in R^C, extended by zero. This is the spectral
radius of sym(U^T tau_t U), not ||U^T tau_t U|| (which allows two vectors).
"""
    lookup=np.full(F.q,-1,dtype=np.int64);lookup[C]=np.arange(len(C))
    vals=[];best=(-1,None,None)
    for t in range(1,F.q):
        dest=lookup[F.add(C,t)];ok=dest>=0
        M=U[ok].T@U[dest[ok]];M=(M+M.T)/2
        lam,vec=eigh(M,check_finite=False)
        j=int(np.argmax(np.abs(lam)));score=float(abs(lam[j]))
        vals.append(score)
        if score>best[0]+1e-13:best=(score,t,U@vec[:,j])
    return {'max_abs_overlap':best[0],'maximizing_shift':int(best[1]),
            'median_abs_overlap':float(np.median(vals)),
            'count_ge_0_75':int(np.sum(np.array(vals)>=.75))},best[2]

def certificate(F,C,S,v,label):
    """Round only discovery vector; exact int64 sums independently recomputed
    with Python ints. v^T H v / ||v||² = A/sqrt(q)-B/q.
    """
    if v[np.argmax(np.abs(v))]<0:v=-v
    z=np.rint(v/np.max(np.abs(v))*1024).astype(np.int64)
    n=sum(int(x)**2 for x in z);sm=sum(map(int,z))
    sig=sum(int(z[i])*sum(int(S[i,j])*int(z[j]) for j in range(len(z))) for i in range(len(z)))
    floatq=sig/(math.sqrt(F.q)*n)-sm*sm/(F.q*n)
    # q square permits an exact rational H quotient.
    out={'label':label,'q':F.q,'characteristic':F.p,'degree':F.degree,
         'C':C.tolist(),'integer_vector':z.tolist(),'norm_squared':n,
         'sum_entries':sm,'signed_quadratic_form':sig,
         'rayleigh_H_float':floatq,
         'formula':'signed_quadratic_form/(sqrt(q)*norm_squared) - sum_entries^2/(q*norm_squared)'}
    if F.degree==2:
        from fractions import Fraction
        exact=Fraction(sig,F.p*n)-Fraction(sm*sm,F.q*n)
        out['rayleigh_H_rational']=[exact.numerator,exact.denominator]
    return out

def analyze(F,seed,tail_window=.025,save_vectors=True):
    C,S,H=F.instance();m=len(C)
    ev,U=eigh(H,check_finite=False);k=int(np.argmax(np.abs(ev)))
    radius=float(abs(ev[k]));mask=np.abs(ev)>=radius-tail_window
    # Both signs included so that a largest negative eigenvalue is not missed.
    tail=U[:,mask]
    raw,vraw=transport_profile(F,C,U[:,[k]])
    completed,vcompleted=transport_profile(F,C,tail)
    rng=np.random.default_rng(seed+F.q);perm=rng.permutation(m)
    shuffled,vshuffled=transport_profile(F,C,tail[perm])
    uniform,_=transport_profile(F,C,np.ones((m,1))/math.sqrt(m))
    # Rotation tests distinguish basis invariance from arbitrary eigensolver choices.
    R,_=np.linalg.qr(rng.standard_normal((tail.shape[1],tail.shape[1])))
    rotated,_=transport_profile(F,C,tail@R)
    assert abs(rotated['max_abs_overlap']-completed['max_abs_overlap'])<1e-10
    rec={'field':f'F_{F.p}' if F.degree==1 else f'F_{F.p}^2','q':F.q,'m':m,
         'degree':F.degree,'radius_H':radius,'radius_Z_two_anchor':2*radius/math.sqrt(3),
         'tail_window':tail_window,'tail_rank':int(mask.sum()),
         'top_eigenvector':raw,'tail_subspace':completed,
         'same_spectrum_relabelled_tail':shuffled,'uniform_baseline':uniform,
         'max_basis_rotation_error':abs(rotated['max_abs_overlap']-completed['max_abs_overlap']),
         'spectral_moments_2_4_6':[float(np.mean(ev**j)) for j in [2,4,6]],
         'tail_witness_rayleigh_H':float(vcompleted@H@vcompleted),
         'relabelled_moment_preservation':'exact by simultaneous row/column permutation',
         'subfield_seed':None}
    if F.degree==2:
        v=(C<F.p).astype(float);v/=np.linalg.norm(v)
        pr,_=transport_profile(F,C,v[:,None])
        rec['subfield_seed']={'rayleigh_H':float(v@H@v),'transport':pr,
            'tail_projection_mass':float(np.sum((tail.T@v)**2))}
    if save_vectors:
        cert=certificate(F,C,S,vcompleted,'transport-optimized spectral tail')
        cert['transport_shift']=completed['maximizing_shift']
        zz=np.array(cert['integer_vector'],dtype=np.int64)
        full=np.zeros(F.q,dtype=np.int64);full[C]=zz
        numerator=sum(int(full[x])*int(full[F.add(x,cert['transport_shift'])]) for x in range(F.q))
        cert['translation_inner_product']=numerator
        cert['translation_overlap_exact']=[numerator,cert['norm_squared']]
        rec['certificate_file']=f'witness_{F.q}.json'
        (HERE/rec['certificate_file']).write_text(json.dumps(cert,indent=2)+'\n')
    return rec

def null_counterexample(p=1009,seed=1307):
    """An interval has almost perfect one-step translation; arithmetic alone
    does not imply Paley edge. Also optimize a Paley tail on same dimension.
    """
    F=Field(p);n=math.isqrt(p);C=np.arange(2,2+n,dtype=np.int64)
    S=F.chi[F.sub(C[:,None],C[None,:])]
    H=S/math.sqrt(p)-np.ones((n,n))/p;v=np.ones(n)/math.sqrt(n)
    pr,_=transport_profile(F,C,v[:,None])
    return {'kind':'interval false-positive control','p':p,'size':n,
            'translation':pr,'rayleigh_H':float(v@H@v),
            'radius_H':float(np.max(np.abs(eigh(H,eigvals_only=True))))}

def exact_small_checks():
    checks=0
    for F in [Field(13),Field(17),Field(5,2),Field(7,2)]:
        q=F.q;S=F.chi[F.sub(F.x[:,None],F.x[None,:])]
        assert np.array_equal(S,S.T);assert np.array_equal(S@S,q*np.eye(q,dtype=np.int64)-np.ones((q,q),dtype=np.int64));checks+=2
        C,Sc,H=F.instance();rng=np.random.default_rng(q);perm=np.arange(q);perm[C]=rng.permutation(C)
        Sp=S[perm][:,perm]
        assert np.array_equal(Sp@Sp,S@S);assert np.array_equal(Sp[:,[0,1]],S[:,[0,1]]);checks+=2
        # Conjugacy preserves every polynomial trace, demonstrated exactly at small powers.
        Pc=Sp[np.ix_(C,C)];A=np.eye(len(C),dtype=np.int64);B=A.copy()
        for j in range(1,7):
            A=A@Sc;B=B@Pc;assert np.trace(A)==np.trace(B);checks+=1
    return checks

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--quick',action='store_true');ap.add_argument('--output',default='results.json');args=ap.parse_args()
    fields=[Field(p) for p in ([101,257] if args.quick else [101,257,401,1009,2017,4001])]
    fields += [Field(p,2) for p in ([11,17] if args.quick else [11,17,23,31,43,61])]
    rows=[]
    for F in fields:
        row=analyze(F,20260905);rows.append(row)
        print(F.q,row['tail_rank'],round(row['radius_Z_two_anchor'],6),round(row['top_eigenvector']['max_abs_overlap'],6),round(row['tail_subspace']['max_abs_overlap'],6),round(row['same_spectrum_relabelled_tail']['max_abs_overlap'],6),flush=True)
    out={'status':'experimental mathematical tool; no asymptotic bound or novelty theorem',
         'seed':20260905,'exact_small_checks':exact_small_checks(),'cases':rows,
         'interval_null':null_counterexample(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/args.output).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
