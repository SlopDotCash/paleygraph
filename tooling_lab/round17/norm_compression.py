#!/usr/bin/env python3
"""Compress centered coset norms using a checked short kernel relation.

If f(u)=0 mod p and F is supported at u modulo p, then D=fF/p is integral.
Writing Norm(f)=kp and Norm(F)=p^(N-1)tau gives Norm(D)=k*tau.
The digit bound |D_i| <= floor((p-1)||f||_1/(2p)) is exact. This is a
coefficient-height reduction, not a new norm theorem or a cofactor estimate.
"""
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp
from norm_carry import norm,multiply

HERE=Path(__file__).resolve().parent


def kernel_relation(p,N,g):
    if g==2:return [2]+[0]*(N-2)+[1],'specified_2_minus_X_inverse'
    d=min(N,16);u=pow(g,-1,p);rows=[[p]+[0]*(d-1)]
    for j in range(1,d):
        row=[-pow(u,j,p)]+[0]*(d-1);row[j]=1;rows.append(row)
    reduced=sp.Matrix(rows).lll()
    choices=[[int(v) for v in reduced.row(j)] for j in range(d)]
    chosen=min(choices,key=lambda r:(sum(abs(v) for v in r),sum(v*v for v in r),r))
    return chosen+[0]*(N-d),f'bounded_exact_LLL_dimension_{d}'


def main():
    source=HERE/'norm_carry.json';raw=json.loads(source.read_text());cases=[]
    for case in raw['cases']:
        p,N,g=(case[k] for k in ('p','N','g'));f,method=kernel_relation(p,N,g);u=pow(g,-1,p)
        assert any(f) and sum(v*pow(u,j,p) for j,v in enumerate(f))%p==0
        norm_f=norm(f);assert norm_f>0 and norm_f%p==0;k=norm_f//p
        l1=sum(abs(v) for v in f);bound=(p-1)*l1//(2*p);rows=[];seen=set()
        for r in case['records']:
            for a,F,tau in [(r['a'],r['F'],r['norm_defect']),(r['next_a'],r['next_F'],r['next_norm_defect'])]:
                if a in seen:continue
                seen.add(a);numerator=multiply(f,F);assert all(v%p==0 for v in numerator)
                D=[v//p for v in numerator];assert max(map(abs,D))<=bound
                digit_norm=norm(D);assert digit_norm==k*tau
                rows.append({'a':a,'digits':D,'norm_defect':tau,'digit_norm':digit_norm,
                             'original_coefficient_bits':max(abs(v).bit_length() for v in F),
                             'digit_coefficient_bits':max(abs(v).bit_length() for v in D),
                             'original_norm_bits':(tau*p**(N-1)).bit_length(),'digit_norm_bits':digit_norm.bit_length()})
        outcase={'name':case['name'],'p':p,'n':case['n'],'N':N,'g':g,'relation':f,'discovery':method,'relation_l1':l1,
                 'relation_norm':norm_f,'relation_cofactor':k,'digit_height_bound':bound,'records':rows,
                 'maximum_original_coefficient_bits':max(r['original_coefficient_bits'] for r in rows),
                 'maximum_digit_coefficient_bits':max(r['digit_coefficient_bits'] for r in rows),
                 'maximum_original_norm_bits':max(r['original_norm_bits'] for r in rows),
                 'maximum_digit_norm_bits':max(r['digit_norm_bits'] for r in rows)}
        cases.append(outcase);print(json.dumps({k:v for k,v in outcase.items() if k not in ('records','relation')}),flush=True)
    out={'status':'produced','scope':'Exact short-relation digit compression of every distinct coset vector in the carry records. LLL is bounded discovery only; relations and norm identities are checked. No shortest relation, minimal cofactor, runtime gain or uniform shell bound is asserted.',
         'cases':cases,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['norm_compression.py','norm_carry.py']},'input_sha256':{source.name:sha256(source.read_bytes()).hexdigest()}}
    (HERE/'norm_compression.json').write_text(json.dumps(out,separators=(',',':'))+'\n')


if __name__=='__main__':main()
