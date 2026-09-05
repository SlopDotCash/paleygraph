#!/usr/bin/env python3
"""Revision: subtract amplitude-only axes and preserve the known S3 action
in null controls before interpreting phase-plane excess as new structure.
"""
import hashlib,json,math
from pathlib import Path
import numpy as np
from ambiguity_lab import HERE,SEED,load,ambiguity,features,character_matrix

def equivariant_permutation(p,C,rng):
    C=list(map(int,C));idx={x:i for i,x in enumerate(C)}
    reflect=lambda x:(1-x)%p;invert=lambda x:pow(x,-1,p)
    unseen=set(C);orbits=[]
    while unseen:
        seed=min(unseen);orbit={seed};todo=[seed]
        while todo:
            x=todo.pop()
            for g in [reflect,invert]:
                y=g(x)
                if y not in orbit:orbit.add(y);todo.append(y)
        unseen-=orbit;orbits.append(sorted(orbit))
    regular=[o for o in orbits if len(o)==6];target=rng.permutation(len(regular));mapping={x:x for x in C}
    for i,o in enumerate(regular):
        pairs={o[0]:regular[target[i]][0]};todo=[o[0]]
        while todo:
            x=todo.pop();y=pairs[x]
            for g in [reflect,invert]:
                xx,yy=g(x),g(y)
                if xx in pairs:assert pairs[xx]==yy
                else:pairs[xx]=yy;todo.append(xx)
        assert set(pairs)==set(o);mapping.update(pairs)
    assert len(set(mapping.values()))==len(C)
    assert all(mapping[reflect(x)]==reflect(mapping[x]) and mapping[invert(x)]==invert(mapping[x]) for x in C)
    return np.array([idx[mapping[x]] for x in C]),len(regular)

def residual_features(A):
    p=len(A);sq=abs(A)**2
    # a=0 is amplitude-only. b=0 is the previously measured real translation axis.
    off=sq[1:,1:]
    full=features(A)
    return {'off_axes_max':float(np.sqrt(off.max())),
        'off_axes_quartic_sum':float(np.sum(off**2)),
        'nonvertical_line_max':float(max(full['line_energies'][:p])),
        'vertical_line':full['line_energies'][p],
        'horizontal_line':full['line_energies'][0]}

def main():
    inp=json.loads((HERE/'results.json').read_text());rows=[]
    for old in inp['cases']:
        p,side=old['p'],old['side'];C,z,full,n,_=load(p,side);f=full/math.sqrt(n);H=character_matrix(p,C)
        A=np.load(HERE/f'phase_plane_{p}_{side}.npz')['original'];base=residual_features(A)
        rng=np.random.default_rng(SEED+2*p+(side=='negative'));controls=[]
        for i in range(8):
            perm,regular=equivariant_permutation(p,C,rng);g=np.zeros(p,complex);g[C]=f[C][perm]
            B=ambiguity(g);controls.append({'index':i,'rayleigh_H_same_operator':float(np.vdot(g[C],H@g[C]).real),
                'regular_orbits_permuted':regular,**residual_features(B)})
        # Chirps preserve all phase-plane magnitude invariants while changing the
        # vector's Rayleigh quotient for the FIXED Paley operator.
        chirp_rays=[];x=np.arange(p,dtype=np.int64)
        for c in [0,1,2,3,4,5,7,11,13,17]:
            g=f*np.exp(2j*np.pi*((c*x*x)%p)/p)
            chirp_rays.append({'c':c,'rayleigh_H_same_operator':float(np.vdot(g[C],H@g[C]).real)})
        rows.append({'p':p,'side':side,'baseline':base,'S3_equivariant_shuffles':controls,
            'chirp_orbit_fixed_operator_rayleighs':chirp_rays,
            'vertical_amplitude_identity_error':abs(base['vertical_line']-(p*np.sum(abs(f)**4)-1))})
        print(p,side,'off4',base['off_axes_quartic_sum'],'null',[round(x['off_axes_quartic_sum'],4) for x in controls],flush=True)
    out={'status':'amplitude-subtracted and symmetry-preserving controls; finite exploratory comparison','cases':rows,
       'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'refinement.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
