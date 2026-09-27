#!/usr/bin/env python3
"""Separate full Q CRT readback and blockwise selected-contraction checks."""
from hashlib import sha256
import json
from math import comb,prod
from pathlib import Path
import subprocess
import time
import numpy as np

HERE=Path(__file__).resolve().parent


def coefficient(p,m,d):
    if d<0:return 0
    return sum((-1)**j*comb(m,j)*comb(p,d-j) for j in range(max(0,d-p),min(d,m)+1))


def crt(residues,primes):
    x=0;M=1
    for r,p in zip(residues,primes):x+=M*((r-x)*pow(M,-1,p)%p);M*=p
    return x-M if x>M//2 else x,M


def check(row):
    s=row['statistics'];q=s['q'];c=s['selected'];a=s['deleted'];d=s['degree'];base=sorted(set(c)-set(a));n=len(c)
    start=time.monotonic()
    payload=' '.join(map(str,[q,n,d,*a,*c]))+'\n'
    proc=subprocess.run([str(HERE/'review_contraction')],input=payload,text=True,capture_output=True,check=True)
    review=json.loads(proc.stdout);Q,M=crt(review['residues'],review['primes'])
    bound=q*(q-1)*review['max_abs_g']*review['max_abs_h']
    assert M>2*bound and Q==s['Q']
    # Independently form signs and coefficient-count tables; exact int64 blocks
    # are reduced to Python integers before whole-field accumulation.
    chi=np.full(q,-1,dtype=np.int64);chi[0]=0
    for startx in range(1,(q+1)//2,65536):
        xs=np.arange(startx,min((q+1)//2,startx+65536),dtype=np.int64);chi[xs*xs%q]=1
    table=np.zeros((3,2,len(base)+1),dtype=np.int64)
    for j,degree in enumerate((d,d-1,d-2)):
        for z in range(2):
            for p in range(len(base)+1-z):table[j,z,p]=coefficient(p,len(base)-z-p,degree)
    values={k:0 for k in ('t0','G0','G2','H0','H2')}
    fields={k:[0]*n for k in ('f','v','w')}
    indices=list(np.triu_indices(n));pairs=list(zip(*map(list,indices)))
    full=q<=65537
    if not full:
        off=[(i,j) for i,j in pairs if i!=j];pairs=[off[int(k)] for k in np.linspace(0,len(off)-1,12,dtype=int)]
    gram={ij:0 for ij in pairs}
    for low in range(0,q,4096):
        xs=np.arange(low,min(q,low+4096),dtype=np.int64)
        signs=chi[(xs[:,None]-np.array(c)[None,:])%q]
        rs=signs[:,[i for i,x in enumerate(c) if x in base]]
        z=np.sum(rs==0,axis=1);p=np.sum(rs==1,axis=1)
        target,g,h=(table[j,z,p] for j in range(3))
        for k,vs in [('t0',target),('G0',g),('G2',g*g),('H0',h),('H2',h*h)]:values[k]+=int(vs.sum())
        for k,vs in [('f',g),('v',h),('w',g*h)]:
            block=vs@signs
            for i,v in enumerate(block):fields[k][i]+=int(v)
        for i,j in pairs:gram[i,j]+=int(np.sum(h*signs[:,i]*signs[:,j]))
        for i,x in enumerate(c):
            if low<=x<low+len(xs):assert int(g[x-low])==s['g_selected'][i] and int(h[x-low])==s['h_selected'][i]
    for k,v in values.items():assert s[k]==v,k
    for k,v in fields.items():assert s[k+'_selected']==v,k
    for (i,j),v in gram.items():assert s['K_selected'][i][j]==v
    return {'q':q,'n':n,'family':row['family'],'Q_reconstructed':Q,'Q_bound':bound,'CRT_modulus':M,
            'Q_modular_review':review,'full_K_upper_entries_checked':len(pairs) if full else 0,
            'sampled_K_off_diagonal_entries_checked':len(pairs) if not full else 0,
            'all_selected_fields_checked':5*n,'row_scalar_totals_checked':5,'seconds':time.monotonic()-start}


def main():
    rows=json.loads((HERE/'scale_results.json').read_text())['cases'];checks=[]
    for row in rows:
        r=check(row);checks.append(r);(HERE/'review_partial.json').write_text(json.dumps(checks,indent=2)+'\n')
        print(json.dumps({k:v for k,v in r.items() if k!='Q_modular_review'}),flush=True)
    out={'status':'passed','scope':'Independent full Q reconstruction, full selected vectors and row scalars, full K through65537 and12 off-diagonal K entries per larger input.',
         'cases':checks,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('review_scale.py','review_contraction.cpp','review_contraction')},
         'input_sha256':{'scale_results.json':sha256((HERE/'scale_results.json').read_bytes()).hexdigest()}}
    (HERE/'scale_verification.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'review_partial.json').unlink()


if __name__=='__main__':main()
