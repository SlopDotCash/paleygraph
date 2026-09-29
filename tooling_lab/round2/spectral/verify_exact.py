#!/usr/bin/env python3
"""Run exact NTT replay and attach strict rational directional lower bounds.

No floating point is used to ACCEPT a rational bound. Floating-point values
only propose candidates, then integer squaring with sign cases checks them.
"""
import argparse,json,math,hashlib,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent

def above(a,b,p,A,sm,n,sign):
    # sign*H quotient > a/b iff b*sign*A*sqrt(p) > a*p*n+b*sign*sm².
    X=b*sign*A;Y=a*p*n+b*sign*sm*sm
    if X>=0 and Y<0:return True
    if X<=0 and Y>=0:return False
    if X>0:return X*X*p>Y*Y
    return X*X*p<Y*Y

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--primes',nargs='*',type=int);args=ap.parse_args();rows=[]
    files=[HERE/f'results_{p}.json' for p in args.primes] if args.primes else sorted(HERE.glob('results_*.json'))
    for path in files:
        data=json.loads(path.read_text())
        for side in data['sides']:
            w=side['witness'];wp=HERE/w['file'];assert hashlib.sha256(wp.read_bytes()).hexdigest()==w['sha256']
            rec=json.loads(subprocess.check_output([str(HERE/'exact_ntt'),str(wp)],text=True))
            for k in ['norm_squared','sum_entries','sum_absolute_entries']:assert rec[k]==w[k]
            sign=1 if side['side']=='positive' else -1;b=10**6;a=math.floor(sign*w['rounded_rayleigh_H_float']*b)-1
            assert above(a,b,rec['p'],rec['signed_quadratic_form'],rec['sum_entries'],rec['norm_squared'],sign)
            assert abs(rec['rayleigh_H_display']-w['rounded_rayleigh_H_float'])<1e-10
            rec.update({'side':side['side'],'witness_sha256':w['sha256'],'directional_rayleigh_strict_lower_bound':[a,b],
                        'acceptance':'exact modular convolution and sign-aware integer squaring; numerical search proposes bound only'})
            rows.append(rec);print(rec['p'],side['side'],'exact',rec['signed_quadratic_form'],'lower',a/b,flush=True)
    out={'status':'exact finite lower certificates only','witnesses':rows,
         'cpp_source_sha256':hashlib.sha256((HERE/'exact_ntt.cpp').read_bytes()).hexdigest(),
         'python_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'exact_validation.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
