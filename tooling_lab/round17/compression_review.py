#!/usr/bin/env python3
"""Check every short relation and digit norm by independent integer resultants."""
from hashlib import sha256
import json
from pathlib import Path
from carry_review import product,result_norm

HERE=Path(__file__).resolve().parent


def main():
    source=HERE/'norm_compression.json';base=HERE/'norm_carry.json';review=HERE/'carry_review.json'
    compressed=json.loads(source.read_text());raw=json.loads(base.read_text());oldreview=json.loads(review.read_text())
    assert oldreview['status']=='passed' and oldreview['input_sha256']['norm_carry.json']==sha256(base.read_bytes()).hexdigest()
    cases=[]
    for c,original in zip(compressed['cases'],raw['cases']):
        assert c['name']==original['name'];p,N,g=(c[k] for k in ('p','N','g'));f=c['relation'];u=pow(g,-1,p)
        assert len(f)==N and any(f) and sum(a*pow(u,j,p) for j,a in enumerate(f))%p==0
        nf=result_norm(f,N);assert nf==c['relation_norm']==p*c['relation_cofactor']
        l1=sum(abs(v) for v in f);bound=(p-1)*l1//(2*p);assert l1==c['relation_l1'] and bound==c['digit_height_bound']
        points={}
        for r in original['records']:
            points[r['a']]=(r['F'],r['norm_defect']);points[r['next_a']]=(r['next_F'],r['next_norm_defect'])
        assert sorted(r['a'] for r in c['records'])==sorted(points)
        for row in c['records']:
            F,tau=points[row['a']];D=row['digits'];assert len(D)==N and all(type(v) is int and abs(v)<=bound for v in D)
            assert product(f,F,N)==[p*v for v in D]
            nd=result_norm(D,N);assert nd==row['digit_norm']==c['relation_cofactor']*tau and row['norm_defect']==tau
            assert row['original_coefficient_bits']==max(abs(v).bit_length() for v in F)
            assert row['digit_coefficient_bits']==max(abs(v).bit_length() for v in D)
            assert row['original_norm_bits']==(tau*p**(N-1)).bit_length() and row['digit_norm_bits']==nd.bit_length()
        for field,record in [('maximum_original_coefficient_bits','original_coefficient_bits'),('maximum_digit_coefficient_bits','digit_coefficient_bits'),('maximum_original_norm_bits','original_norm_bits'),('maximum_digit_norm_bits','digit_norm_bits')]:
            assert c[field]==max(r[record] for r in c['records'])
        record={'case':c['name'],'vectors_checked':len(c['records']),'relation_height':max(map(abs,f)),'relation_l1':l1,'digit_height_bound':bound,
                'original_coefficient_bits':c['maximum_original_coefficient_bits'],'digit_coefficient_bits':c['maximum_digit_coefficient_bits'],
                'original_norm_bits':c['maximum_original_norm_bits'],'digit_norm_bits':c['maximum_digit_norm_bits']};cases.append(record);print(json.dumps(record),flush=True)
    out={'status':'passed','scope':'Each coefficient compression identity and digit-height bound checked by a separate polynomial product implementation; every digit norm and relation norm independently computed by integer resultants. Original large norms are bound to the completed carry review.',
         'cases':cases,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['compression_review.py','carry_review.py']},
         'input_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in [source,base,review]}}
    (HERE/'compression_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
