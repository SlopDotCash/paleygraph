#!/usr/bin/env python3
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify,rational


def main():
    path=HERE/'lattice_information.json';data=json.loads(path.read_text());inputs={path.name:sha256(path.read_bytes()).hexdigest()};cases=[]
    for item in data['cases']:
        source=HERE/item['source'];inputs[item['source']]=sha256(source.read_bytes()).hexdigest();c=json.loads(source.read_text())['certificate']
        certify(c);N,p,f,U=c['N'],c['p'],c['relation'],c['erased']
        M=sp.Matrix([[f[(i-j)%N]*(1 if i>=j else -1) for j in range(N)] for i in range(N)])
        assert len(item['rows'])==len(U);integral=0;congruent=0;carry=0
        for j,record in enumerate(item['rows']):
            v=[rational(x) for x in record['inverse_row']];assert len(v)==N and record['basis_index']==j
            H=[0]*N
            for u,x in zip(U,c['basis'][j]):H[u]=p*x
            # M is invertible by the already certified adjugate identity;
            # checking M*v=p*H_original proves the unique inverse image.
            assert M*sp.Matrix(v)==sp.Matrix(H)
            is_int=all(x.denominator==1 for x in v)
            is_congruent=is_int and all((v[i]-c['scalar_steps'][j]*pow(c['g'],i,p))%p==0 for i in range(N))
            is_carry=c['scalar_steps'][j]==0 and is_int and all(x%p==0 for x in v)
            assert record['integral']==is_int and record['scalar_congruent']==is_congruent and record['integer_carry_direction']==is_carry
            integral+=is_int;congruent+=is_congruent;carry+=is_carry
        assert item['all_kernel_rows_integral_and_congruent']==(integral==congruent==len(U))
        assert item['invisible_directions']==len(c['invisible_directions']) and item['integer_carry_directions']==carry
        cases.append({'source':item['source'],'rows_checked':len(U),'integral_rows':integral,'scalar_congruent_rows':congruent,
                      'integer_carry_directions':carry,'invisible_directions':len(c['invisible_directions'])})
    out={'status':'passed','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['lattice_information_review.py','../round22/conditioned_review.py']}}
    (HERE/'lattice_information_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(cases),flush=True)


if __name__=='__main__':main()
