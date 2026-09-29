#!/usr/bin/env python3
"""Independent scalar checks on generic complex vectors, not just real witnesses."""
import cmath
from hashlib import sha256
import json
import math
from pathlib import Path
import random
import numpy as np
from ambiguity_lab import ambiguity
from verify_phase import wigner

HERE=Path(__file__).resolve().parent


def root_of_unity(t,p):return cmath.exp(2j*math.pi*(t%p)/p)


def main():
    rng=random.Random(74092026);cases=[]
    for p in (7,13):
        h=pow(2,-1,p)
        for sample in range(3):
            f=[complex(rng.randrange(-5,6),rng.randrange(-5,6)) for _ in range(p)]
            norm=math.sqrt(sum(abs(z)**2 for z in f));f=[z/norm for z in f]
            a=[[sum(f[x].conjugate()*f[(x+u)%p]*root_of_unity(v*(x+u*h),p) for x in range(p))
                for v in range(p)] for u in range(p)]
            actual=ambiguity(np.array(f));err=max(abs(a[u][v]-actual[u,v]) for u in range(p) for v in range(p))
            wd=[[sum(f[(q+t*h)%p]*f[(q-t*h)%p].conjugate()*root_of_unity(-k*t,p) for t in range(p))/p
                 for k in range(p)] for q in range(p)]
            wa=wigner(actual);werr=max(abs(wd[q][k]-wa[q,k]) for q in range(p) for k in range(p))
            ft=[sum(f[x]*root_of_unity(-x*k,p) for x in range(p))/math.sqrt(p) for k in range(p)]
            fa=ambiguity(np.array(ft));ferr=max(abs(fa[u,v]-a[v][(-u)%p]) for u in range(p) for v in range(p))
            assert max(err,werr,ferr)<2e-14
            assert abs(sum(abs(z)**2 for row in a for z in row)-p)<2e-13
            cases.append({'p':p,'sample':sample,'all_phase_points':p*p,
                          'direct_ambiguity_error':err,'direct_Wigner_error':werr,'Fourier_covariance_error':ferr})
    # Check quoted fixed-H chirp effect by direct scalar matrix summation.
    from struct import unpack,iter_unpack
    source=HERE.parents[1]/'round2/spectral/witness_401_negative.bin'
    raw=source.read_bytes();p,m=unpack('<II',raw[:8]);assert p==401
    pairs=list(iter_unpack('<ii',raw[8:]));norm2=sum(z*z for _,z in pairs)
    squares={x*x%p for x in range(1,p)}
    chi=lambda x:0 if x%p==0 else 1 if x%p in squares else -1
    rays=[]
    for c in (0,1):
        f=[z*root_of_unity(c*x*x,p)/math.sqrt(norm2) for x,z in pairs]
        ray=sum(f[i].conjugate()*f[j]*(chi(x-y)/math.sqrt(p)-1/p)
                for i,(x,_) in enumerate(pairs) for j,(y,_) in enumerate(pairs))
        assert abs(ray.imag)<1e-12;rays.append(ray.real)
    assert rays[0]<-.87 and abs(rays[1])<.02
    out={'status':'independent scalar numerical audit passed; no exact nonzero-phase certificate',
         'generic_complex_cases':cases,'fixed_operator_chirp_quotients':rays,
         'source_witness_sha256':sha256(raw).hexdigest(),
         'ambiguity_source_sha256':sha256((HERE/'ambiguity_lab.py').read_bytes()).hexdigest(),
         'Wigner_source_sha256':sha256((HERE/'verify_phase.py').read_bytes()).hexdigest(),
         'reviewer_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'root_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'complete_complex_planes':len(cases),'fixed_H_quotients':rays},indent=2))


if __name__=='__main__':main()
