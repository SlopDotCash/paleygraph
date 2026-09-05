#!/usr/bin/env python3
"""Lower-precision witness enables single-modulus exact ALL translations.
The desired lower Rayleigh quotient is then rechecked, never assumed.
"""
import hashlib,json,struct,subprocess
from pathlib import Path
import numpy as np
from verify_exact import above,HERE
src=HERE/'witness_65537_negative.bin';raw=src.read_bytes();p,m=struct.unpack('<II',raw[:8]);pairs=np.frombuffer(raw[8:],dtype='<i4').reshape(m,2).copy()
pairs[:,1]=np.rint(pairs[:,1]/2).astype('<i4');out=HERE/'witness_65537_negative_512.bin';out.write_bytes(struct.pack('<II',p,m)+pairs.tobytes())
r=json.loads(subprocess.check_output([str(HERE/'exact_ntt'),str(out)],text=True));a,b=866400,1000000
assert above(a,b,p,r['signed_quadratic_form'],r['sum_entries'],r['norm_squared'],-1)
assert 4*a*a>3*b*b
assert r['max_abs_nonzero_translation_numerator']>=0 and r['count_translation_ge_3_5']==1
r.update({'status':'exact finite outlier with no nonzero overlap >=3/5; no asymptotic claim','file':out.name,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
          'strict_directional_rayleigh_lower_bound':[a,b],'rational_bound_exceeds_sqrt3_over2_squared_margin':4*a*a-3*b*b,
          'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
(HERE/'extra_transport_check.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
