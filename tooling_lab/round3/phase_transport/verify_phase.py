#!/usr/bin/env python3
"""Direct complex sums, Weyl covariance, and Wigner marginals verify the
phase conventions. Nonzero-phase entries remain numerical, not exact cyclotomic certificates.
"""
import cmath,hashlib,json,math
from pathlib import Path
import numpy as np
from scipy.fft import fft,fft2
from ambiguity_lab import HERE,SEED,load,ambiguity

def direct(f,a,b):
    p=len(f);half=pow(2,-1,p)
    return sum(complex(f[x]).conjugate()*complex(f[(x+a)%p])*cmath.exp(2j*math.pi*((b*(x+a*half))%p)/p) for x in range(p))

def displacement(f,a,b):
    p=len(f);x=np.arange(p,dtype=np.int64);half=pow(2,-1,p)
    return np.exp(2j*np.pi*((b*(x+a*half))%p)/p)*f[(x+a)%p]

def wigner(A):
    p=len(A)
    # q first, k second: W(q,k)=p^-1 sum_a f(q+a/2)conj(f(q-a/2)) e_p(-ka).
    return (fft2(A)/(p*p)).T

def main():
    rows=[]
    for p in [101,401,1009]:
      for side in ['positive','negative']:
        C,z,full,n,_=load(p,side);f=full/math.sqrt(n);A=np.load(HERE/f'phase_plane_{p}_{side}.npz')['original']
        rng=np.random.default_rng(SEED+5*p+(side=='negative'));samples=rng.integers(0,p,size=(32,2));err=max(abs(direct(f,int(a),int(b))-A[a,b]) for a,b in samples)
        assert err<1e-10
        W=wigner(A);imag=float(np.max(abs(W.imag)));W=W.real
        marginal_q=float(np.max(abs(W.sum(axis=1)-abs(f)**2)))
        marginal_k=float(np.max(abs(W.sum(axis=0)-abs(fft(f)/math.sqrt(p))**2)))
        purity=abs(float(np.sum(W*W))-1/p)
        assert imag<1e-10 and marginal_q<1e-10 and marginal_k<1e-10 and purity<1e-10
        group_error=0.
        for a,b,c,d in rng.integers(0,p,size=(16,4)):
            lhs=displacement(displacement(f,int(c),int(d)),int(a),int(b))
            phase=np.exp(2j*np.pi*int(((a*d-b*c)*pow(2,-1,p))%p)/p)
            rhs=phase*displacement(f,int((a+c)%p),int((b+d)%p))
            group_error=max(group_error,float(np.max(abs(lhs-rhs))))
        assert group_error<1e-10
        # Fourier covariance tested at p101; all p² points, not sampled.
        fourier_error=None
        if p==101:
            B=ambiguity(fft(f)/math.sqrt(p));x=np.arange(p)
            target=np.array([[A[b,(-a)%p] for b in range(p)] for a in range(p)])
            fourier_error=float(np.max(abs(B-target)));assert fourier_error<1e-10
        rows.append({'p':p,'side':side,'direct_phase_samples':32,'direct_complex_max_error':float(err),
            'Wigner_imaginary_error':imag,'Wigner_position_marginal_error':marginal_q,
            'Wigner_frequency_marginal_error':marginal_k,'Wigner_purity_error':purity,
            'Wigner_l1':float(np.sum(abs(W))),'Weyl_group_samples':16,'Weyl_group_max_error':group_error,
            'Fourier_covariance_full_plane_error':fourier_error})
        print(p,side,'verified',flush=True)
    (HERE/'validation.json').write_text(json.dumps({'status':'phase convention independently checked by direct sums; nonzero phases numerical',
       'cases':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
if __name__=='__main__':main()
