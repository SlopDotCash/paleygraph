#!/usr/bin/env python3
"""Second iteration: joint structure of ONE certified witness's translations.

The pilot shows high max overlap alone is insufficient (p=101). Here we
retain the actual vector, find its large-overlap translations, and measure
addition/multiplication growth. No optimization is repeated per shift.
"""
import json,math,hashlib
from fractions import Fraction
from pathlib import Path
import numpy as np
from transport_microscope import Field,HERE

def signature(F,C,z):
    n=sum(int(x)**2 for x in z);full=np.zeros(F.q,dtype=np.int64);full[C]=z
    corr=[]
    for t in range(F.q):
        corr.append(int(np.dot(z,full[F.add(C,t)])))
    # Rational threshold 3/5; the comparison uses integer arithmetic.
    T=np.array([t for t,c in enumerate(corr) if 5*abs(c)>=3*n],dtype=np.int64)
    assert len(T)>0 and 0 in T
    sums=np.unique(F.add(T[:,None],T[None,:]));products=np.unique(F.mul(T[:,None],T[None,:]))
    closure=T.copy();steps=0
    while True:
        expanded=np.unique(F.add(closure[:,None],T[None,:]));steps+=1
        if len(expanded)==len(closure):break
        closure=expanded
    nz=closure[closure!=0]
    mult=np.unique(F.mul(closure[:,None],closure[None,:]))
    is_field=1 in closure and np.array_equal(mult,closure) and len(closure)>1
    return {'threshold':[3,5],'large_overlap_translations':T.tolist(),
            'translation_count':len(T),'sumset_size':len(sums),'productset_size':len(products),
            'additive_closure_size':len(closure),'closure_steps':steps,
            'proper_subfield_recovered':bool(is_field and len(closure)<F.q),
            'witness_additive_energy_normalized':str(Fraction(sum(c*c for c in corr),n*n)),
            'fixed_witness_max_abs_overlap':str(Fraction(max(map(abs,corr[1:])),n)),
            'high_overlap_translation_min':str(Fraction(min(abs(corr[t]) for t in T),n))}

def verify(c):
    F=Field(c['characteristic'],c['degree']);C,S,H=F.instance()
    assert C.tolist()==c['C'];z=np.array(c['integer_vector'],dtype=np.int64)
    n=sum(int(x)**2 for x in z);sm=sum(map(int,z));a=int(z@S@z)
    assert n==c['norm_squared'] and sm==c['sum_entries'] and a==c['signed_quadratic_form']
    t=c['transport_shift'];full=np.zeros(F.q,dtype=np.int64);full[C]=z
    corr=sum(int(z[i])*int(full[F.add(int(x),t)]) for i,x in enumerate(C))
    assert corr==c['translation_inner_product']
    return F,C,z

def main():
    records=[]
    for row in json.loads((HERE/'results.json').read_text())['cases']:
        c=json.loads((HERE/row['certificate_file']).read_text());F,C,z=verify(c)
        sig=signature(F,C,z)
        print(F.q,sig['translation_count'],sig['sumset_size'],sig['productset_size'],sig['additive_closure_size'],sig['proper_subfield_recovered'],flush=True)
        records.append({'q':F.q,'signature':sig})
    # Stress control: arithmetic interval has large overlap and small doubling,
    # but the additive closure is the whole prime field; it has no high Paley edge.
    F=Field(1009);C=np.arange(2,33);z=np.ones(len(C),dtype=np.int64)
    output={'status':'iteration after pilot falsified single-score interpretation',
            'exact_witnesses_verified':len(records),'cases':records,
            'interval_signature':signature(F,C,z),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'refinement.json').write_text(json.dumps(output,indent=2)+'\n')
if __name__=='__main__':main()
