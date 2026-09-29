#!/usr/bin/env python3
"""Known Weyl ambiguity functions applied to saved Paley Rayleigh witnesses.

Finite phase-plane diagnostics only; no conjecture or novelty theorem.
"""
import json,hashlib,math,struct,sys,time
from pathlib import Path
import numpy as np
import scipy
from scipy.fft import fft,ifft
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUT=ROOT/'tooling_lab/round2/spectral'
SEED=20260905


def ambiguity(f):
    """A[a,b]=sum_x conj(f[x]) f[x+a] e_p(b*(x+a/2)).
    f normalized; plus-sign Fourier transform, half computed in F_p.
    """
    p=len(f);assert p<=1009,'p² allocation guard';x=np.arange(p,dtype=np.int64);half=pow(2,-1,p)
    roots=np.exp(2j*np.pi*x/p);out=np.empty((p,p),dtype=np.complex128)
    for a in range(p):
        g=f.conj()*f[(x+a)%p]
        out[a]=p*ifft(g)*roots[(x*a*half)%p]
    return out


def features(A):
    p=len(A);sq=abs(A)**2;nonzero=sq.copy();nonzero[0,0]=0
    a=np.arange(1,p,dtype=np.int64)
    lines=np.array([np.sum(sq[a,(s*a)%p]) for s in range(p)]+[np.sum(sq[0,1:])])
    flat=int(np.argmax(nonzero));peak=list(map(int,np.unravel_index(flat,sq.shape)))
    return {'ordinary_translation_max':float(np.max(abs(A[1:,0]))),
      'all_nonidentity_max':float(np.sqrt(nonzero.max())),'peak_location':peak,
      'both_coordinates_nonzero_max':float(np.max(abs(A[1:,1:]))),
      'moyal_sum_squared':float(sq.sum()),'quartic_concentration':float(np.sum(sq**2)/p),
      'largest_line_energy_excluding_origin':float(max(lines)),
      'largest_line_slope':int(np.argmax(lines)),
      'line_mean_excluding_origin':float(lines.mean()),
      'all_line_energy_sum':float(lines.sum()),
      'phase_plane_rms_excluding_origin':float(np.sqrt((sq.sum()-1)/(p*p-1))),
      'line_energies':lines.tolist()}


def character_matrix(p,C):
    ch=np.full(p,-1,dtype=np.int64);ch[0]=0;x=np.arange(1,(p+1)//2,dtype=np.int64);ch[x*x%p]=1
    return ch[(C[:,None]-C[None,:])%p]/math.sqrt(p)-1/p


def load(p,side):
    path=INPUT/f'witness_{p}_{side}.bin';raw=path.read_bytes();q,m=struct.unpack('<II',raw[:8]);assert q==p
    pair=np.frombuffer(raw[8:],dtype='<i4').reshape(m,2);C=pair[:,0].copy();z=pair[:,1].astype(np.int64)
    full=np.zeros(p,dtype=np.int64);full[C]=z;n=int(z@z)
    return C,z,full,n,hashlib.sha256(raw).hexdigest()


def checks(f,A,C,z=None,full=None,n=None):
    p=len(f);x=np.arange(p,dtype=np.int64);sq=abs(A)**2;norm=float(np.vdot(f,f).real)
    assert abs(norm-1)<1e-12
    herm=float(np.max(abs(A[np.ix_((-x)%p,(-x)%p)]-A.conj())))
    # Row-wise Parseval depends only on |f|, hence is not phase information.
    mass=abs(f)**2;expected=np.array([p*np.sum(mass*mass[(x+a)%p]) for a in range(p)])
    rowerr=float(np.max(abs(sq.sum(axis=1)-expected)))
    assert abs(sq.sum()-p)<1e-9 and herm<1e-10 and rowerr<1e-10
    rec={'moyal_error':float(abs(sq.sum()-p)),'hermitian_phase_error':herm,'row_parseval_error':rowerr}
    if z is not None:
        exact=[sum(int(full[i])*int(full[(i+a)%p]) for i in range(p)) for a in range(p)]
        err=float(max(abs(A[a,0]-exact[a]/n) for a in range(p)))
        assert err<1e-12
        rec.update({'exact_b0_integer_correlations':exact,'norm_squared':n,'exact_b0_max_error':err})
    # One independently computed MUB-basis line purity (slope 3).
    s=3;half=pow(2,-1,p);chirped=f*np.exp(2j*np.pi*((s*half*x*x)%p)/p)
    probs=abs(fft(chirped)/math.sqrt(p))**2
    lhs=float(np.sum(sq[x,(s*x)%p]));rhs=float(p*np.sum(probs**2))
    assert abs(lhs-rhs)<1e-9;rec['MUB_line_parseval_error']=abs(lhs-rhs)
    return rec


def shear_audit(f,A,c=1):
    p=len(f);x=np.arange(p,dtype=np.int64);g=f*np.exp(2j*np.pi*((c*x*x)%p)/p);B=ambiguity(g)
    target=np.empty_like(B)
    for a in range(p):target[a]=A[a,(x+2*c*a)%p]
    err=float(np.max(abs(B-target)));assert err<1e-10
    FA,FB=features(A),features(B)
    return {'chirp_coefficient':c,'complex_shear_max_error':err,
       'ordinary_before':FA['ordinary_translation_max'],'ordinary_after':FB['ordinary_translation_max'],
       'full_plane_max_error':abs(FA['all_nonidentity_max']-FB['all_nonidentity_max']),
       'quartic_invariant_error':abs(FA['quartic_concentration']-FB['quartic_concentration']),
       'max_line_invariant_error':abs(FA['largest_line_energy_excluding_origin']-FB['largest_line_energy_excluding_origin'])},B


def main():
    start=time.monotonic();rows=[];control_checks=[]
    for p in [101,401,1009]:
      x=np.arange(p,dtype=np.int64)
      for side in ['positive','negative']:
        C,z,full,n,sha=load(p,side);f=full/math.sqrt(n);H=character_matrix(p,C)
        A=ambiguity(f);base=features(A);audit=checks(f,A,C,z,full,n);shear,B=shear_audit(f,A)
        source={'p':p,'side':side,'source_witness_sha256':sha,'input_rayleigh_H':float(np.vdot(f[C],H@f[C]).real),
                'baseline':base,'validation':audit,'chirp_shear':shear,'controls':[]}
        rng=np.random.default_rng(SEED+p+(0 if side=='positive' else 1))
        # Same coefficient histogram/support; many independent seeded permutations.
        for seed in range(8):
          g=np.zeros(p,dtype=complex);g[C]=rng.permutation(f[C]);M=ambiguity(g);feat=features(M)
          source['controls'].append({'kind':'shuffled_coefficients','index':seed,'rayleigh_H_same_operator':float(np.vdot(g[C],H@g[C]).real),**feat})
        # Fixed magnitudes, altered phases: row energies must be identical.
        g=abs(f)*np.exp(2j*np.pi*rng.random(p));M=ambiguity(g)
        assert np.max(abs(np.sum(abs(M)**2,axis=1)-np.sum(abs(A)**2,axis=1)))<1e-10
        source['controls'].append({'kind':'random_phases_same_magnitudes','index':0,'rayleigh_H_same_operator':float(np.vdot(g[C],H@g[C]).real),**features(M)})
        # Same-support chirp is a structured phase control, not spectral witness.
        g=np.zeros(p,dtype=complex);g[C]=np.exp(2j*np.pi*((C.astype(np.int64)**2)%p)/p)/math.sqrt(len(C));M=ambiguity(g)
        source['controls'].append({'kind':'uniform_masked_chirp','index':0,'rayleigh_H_same_operator':float(np.vdot(g[C],H@g[C]).real),**features(M)})
        # Store complete complex planes only for original and its chirp shear.
        np.savez_compressed(HERE/f'phase_plane_{p}_{side}.npz',original=A,chirped=B)
        rows.append(source)
        print(p,side,'b0',base['ordinary_translation_max'],'phase',base['all_nonidentity_max'],'line',base['largest_line_energy_excluding_origin'],flush=True)
      # Closed-form full-field controls check normalization and line orientation.
      chirp=np.exp(2j*np.pi*((x*x)%p)/p)/math.sqrt(p);Ac=ambiguity(chirp)
      target=np.zeros((p,p),complex);target[x,(-2*x)%p]=1
      chirperr=float(np.max(abs(Ac-target)));assert chirperr<1e-10
      delta=np.zeros(p);delta[0]=1;Ad=ambiguity(delta);dt=np.zeros((p,p));dt[0,:]=1
      assert np.max(abs(Ad-dt))<1e-12
      control_checks.append({'p':p,'full_chirp_line_max_error':chirperr,'full_chirp_features':features(Ac),'delta_features':features(Ad)})
    out={'status':'known ambiguity/Weyl diagnostics applied to finite Paley witnesses; no inverse theorem or novelty claim',
         'seed':SEED,'cases':rows,'closed_form_controls':control_checks,'elapsed_seconds':time.monotonic()-start,
         'runtime':{'python':sys.executable,'numpy':np.__version__,'scipy':scipy.__version__},
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
