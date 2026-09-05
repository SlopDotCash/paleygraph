#!/usr/bin/env python3
"""Independent sampled-entry multiplication checks at both large primes.
Dense full matrix is never formed. Explicit row sums cost O(sample_rows*m).
"""
import hashlib,json
from pathlib import Path
import numpy as np
from fft_transport import PaleyFFT,HERE,SEED

def main():
 rows=[]
 for p in [65537,1000033]:
  A=PaleyFFT(p);rng=np.random.default_rng(SEED+17*p)
  # Character array cross-check by Euler criterion, not square enumeration.
  sample=rng.choice(np.arange(1,p),size=256,replace=False)
  assert all(int(A.chi[x])==(1 if pow(int(x),(p-1)//2,p)==1 else -1) for x in sample)
  v=rng.standard_normal(A.m);u=rng.standard_normal(A.m);Av=A@v;Au=A@u
  selected=rng.choice(A.m,size=32,replace=False);err=[]
  for i in selected:
   direct=float(A.chi[(int(A.C[i])-A.C)%p]@v)/np.sqrt(p)-v.sum()/p
   err.append(abs(direct-Av[i]))
  relative_symmetry=abs(float(u@Av-v@Au))/(np.linalg.norm(u)*np.linalg.norm(v))
  assert max(err)<1e-10 and relative_symmetry<1e-12
  rows.append({'p':p,'m':A.m,'Euler_character_checks':256,'direct_rows':32,'direct_character_products':32*A.m,
    'max_row_error':max(err),'relative_bilinear_symmetry_error':relative_symmetry,
    'dense_full_matrix_allocated':False})
  print(rows[-1],flush=True)
 (HERE/'scaled_operator_check.json').write_text(json.dumps({'status':'sampled independent numerical checks, not interval-arithmetic certification','cases':rows,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
if __name__=='__main__':main()
