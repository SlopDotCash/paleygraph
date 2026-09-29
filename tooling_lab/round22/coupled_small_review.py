#!/usr/bin/env python3
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
from coupled_review import setup,target_certificate,rational
from conditioned_review import certify

HERE=Path(__file__).resolve().parent


def main():
    path=HERE/'coupled_small.json'; data=json.loads(path.read_text()); inputs={path.name:sha256(path.read_bytes()).hexdigest()}; cases=[]
    for case in data['cases']:
        source=HERE/case['source']; inputs[source.name]=sha256(source.read_bytes()).hexdigest(); c=json.loads(source.read_text())['certificate']
        prepared=certify(c); words=prepared[2](list(range(c['p']))); K,V,rows,Q=setup(c)
        pairs=0; differences=set()
        for a in range(c['p']):
            for b in range(a+1,c['p']):
                if any(words[a][i]!=words[b][i] for i in K): continue
                values=[]
                for j in V:
                    x=sum(prepared[0][i][j]*(words[b][u]-words[a][u]) for i,u in enumerate(c['erased']))
                    assert x.denominator==1; values.append(int(x))
                assert any(values); pairs+=1
                if next(v for v in values if v)<0: values=[-v for v in values]
                differences.add(tuple(values))
        assert pairs==case['actual_shared_known_pairs'] and len(differences)==case['distinct_signed_actual_differences']
        real=0; separated=0; continuous=0
        for r in case['records']:
            norm=target_certificate(c,r,rows,Q); assert r['separated']==(F(c['p']-1,c['p'])*norm<1)
            separated+=r['separated']
            if r['actual_scalar_pair']:
                a,b=r['actual_scalar_pair']; assert all(words[a][i]==words[b][i] for i in K)
                t=[sum(prepared[0][i][j]*(words[b][u]-words[a][u]) for i,u in enumerate(c['erased'])) for j in V]
                assert t==r['target'] and not r['separated']; real+=1
            if 'normalized_centered_difference' in r:
                v=[rational(x) for x in r['normalized_centered_difference']]
                assert all(abs(x)<=F(c['p']-1,c['p']) for x in v)
                assert all(sum(a*b for a,b in zip(row,v))==0 for row in rows)
                assert all(sum(a*b for a,b in zip(row,v))==x for row,x in zip(Q,r['target']))
                continuous+=1
            if r['separated']:
                mu=[rational(x) for x in r['direction']]
                # Every genuine small difference obeys the proposed separating cut.
                assert all(abs(sum(a*b for a,b in zip(mu,t)))<=F(c['p']-1,c['p'])*norm for t in differences)
        cases.append({'name':case['name'],'all_actual_shared_known_pairs_checked':pairs,
                      'distinct_actual_difference_classes':len(differences),'targets_checked':len(case['records']),
                      'actual_ambiguity_targets_preserved':real,'separators_checked':separated,'continuous_witnesses_checked':continuous})
    out={'status':'passed','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['coupled_small_review.py','coupled_review.py','conditioned_review.py']}}
    (HERE/'coupled_small_review.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(cases),flush=True)


if __name__=='__main__': main()
