#!/usr/bin/env python3
"""Independent complete polynomial/scalar checks of through-code transport."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from copy import deepcopy
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from itertools import product
from pathlib import Path
import json
import random

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'scalar_fibers'/'scalar_fibers.py'


def load(path):
    spec=spec_from_file_location(path.stem,path)
    module=module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_static(F,dom,k,s,word):
    result=set()
    for cs in product(range(F.q),repeat=k):
        cw=tuple(F.eval(cs,x) for x in dom)
        support=tuple(i for i,x in enumerate(cw) if x==word[i])
        if len(support)>=s:result.add((cw,support))
    return result


def full_expansion(candidate,F,record):
    return {(node['scalar'],tuple(node['codeword']),tuple(node['agreement_support']))
            for z in range(F.q) for node in candidate.materialize_scalar(F,record,z)}


def one_pencil(candidate,oracle,OF,F,dom,k,s,u0,u1,pieces):
    result=candidate.certify_pencil(F,dom,k,s,u0,u1,pieces)
    candidate.verify_certificate(F,result)
    expected=oracle.complete_nodes(OF,dom,k,s,u0,u1)
    actual=full_expansion(candidate,F,result)
    assert expected==actual
    assert len(actual)==result['complete_symbolic_node_count']
    static=direct_static(OF,dom,k,s,u1)
    saved={(tuple(x['codeword']),tuple(x['agreement_support']))
           for x in result['static_direction_certificate']['complete_static_list']}
    assert saved==static and len(static)==result['list_size_at_every_nonanchor_scalar']
    counts=Counter(z for z,_,_ in actual)
    zstar=result['anchor']['scalar']
    assert counts[zstar]==1
    assert all(counts[z]==len(static) for z in range(F.q) if z!=zstar)
    assert len(counts)==result['bad_scalar_count']
    return result


def main():
    candidate=load(SOURCE)
    oracle=load(HERE/'review_support_batches.py')
    digest=sha256(SOURCE.read_bytes()).hexdigest()
    F=candidate.PrimeField(3);OF=oracle.OracleField(3);dom=list(range(3))
    static_checks=unavailable=anchor_checks=pencils=nodes=whole=entire=0
    sample=None
    for k in range(1,4):
        for word_tuple in product(range(3),repeat=3):
            word=list(word_tuple)
            pieces=[{'coefficients':[a],'coordinates':[i for i,x in enumerate(word) if x==a]} for a in range(3)]
            cap=sum(min(word.count(a),k-1) for a in set(word))
            for s in range(k,4):
                if s<=cap:
                    try:candidate.certify_static_word(F,dom,k,s,word,pieces)
                    except candidate.UnavailableCertificate:unavailable+=1
                    else:raise AssertionError('accepted uncertifiable root-bound equality')
                    continue
                checked=candidate.certify_static_word(F,dom,k,s,word,pieces)
                actual={(tuple(x['codeword']),tuple(x['agreement_support'])) for x in checked['complete_static_list']}
                assert actual==direct_static(OF,dom,k,s,word)
                static_checks+=1
                for cs in product(range(3),repeat=k):
                    f=[OF.eval(cs,x) for x in dom]
                    for zstar in range(3):
                        u0=[OF.sub(a,OF.mul(zstar,b)) for a,b in zip(f,word)]
                        cert=one_pencil(candidate,oracle,OF,F,dom,k,s,u0,word,pieces)
                        pencils+=1;nodes+=cert['complete_symbolic_node_count']
                        whole+=cert['bad_scalar_count']==3
                        entire+=cert['anchor']['discovery']=='entire_pencil_in_code'
                        if len(cert['nonanchor_fibers'])>=2:sample=cert
    # Exhaust all F3 length-three pencils at k=2. Directly classify every
    # exact codeword intersection, including no intersection and whole lines.
    book={tuple(OF.eval(cs,x) for x in dom) for cs in product(range(3),repeat=2)}
    for u0 in product(range(3),repeat=3):
        for u1 in product(range(3),repeat=3):
            exact=[z for z in range(3) if tuple(OF.add(a,OF.mul(z,b)) for a,b in zip(u0,u1)) in book]
            result=candidate.find_code_anchor(F,dom,2,list(u0),list(u1))
            assert (result is None)==(not exact)
            if exact:
                assert result['scalar'] in exact
                assert (result['discovery']=='entire_pencil_in_code')==(len(exact)==3)
                assert tuple(OF.eval(result['coefficients'],x) for x in dom)==tuple(OF.add(a,OF.mul(result['scalar'],b)) for a,b in zip(u0,u1))
            anchor_checks+=1
    rng=random.Random(493319)
    extension_checks=0
    for q,n,k in ((5,5,2),(9,8,2)):
        F=candidate.PrimeField(q) if q==5 else candidate.Field(3,[1,0,1])
        OF=oracle.OracleField(q);dom=list(range(n))
        for trial in range(20):
            piece_cs=[[rng.randrange(q) for _ in range(k)] for _ in range(2)]
            regions=[[i for i in range(n) if i%2==j] for j in range(2)]
            u1=[OF.eval(piece_cs[i%2],x) for i,x in enumerate(dom)]
            fcs=[rng.randrange(q) for _ in range(k)];zstar=rng.randrange(q)
            u0=[OF.sub(OF.eval(fcs,x),OF.mul(zstar,b)) for x,b in zip(dom,u1)]
            pieces=[{'coefficients':cs,'coordinates':region} for cs,region in zip(piece_cs,regions)]
            for s in (3,n):
                cert=one_pencil(candidate,oracle,OF,F,dom,k,s,u0,u1,pieces)
                pencils+=1;nodes+=cert['complete_symbolic_node_count'];extension_checks+=q==9
                whole+=cert['bad_scalar_count']==q
                entire+=cert['anchor']['discovery']=='entire_pencil_in_code'
                if len(cert['nonanchor_fibers'])>=2:sample=cert
    assert sample is not None
    sample_field=candidate.Field(3,[1,0,1]) if sample['field']['q']==9 else candidate.PrimeField(sample['field']['q'])
    mutations=[('missing_nonanchor_fiber',lambda r:r['nonanchor_fibers'].pop()),
               ('false_scalar_domain',lambda r:r['nonanchor_fibers'][0]['scalar_domain'].__setitem__('excluded',[])),
               ('false_node_count',lambda r:r.__setitem__('complete_symbolic_node_count',0)),
               ('duplicated_partition_coordinate',lambda r:r['static_direction_certificate']['pieces'][0]['coordinates'].append(r['static_direction_certificate']['pieces'][0]['coordinates'][0])),
               ('false_support',lambda r:r['nonanchor_fibers'][0]['agreement_support_for_every_scalar_in_domain'].clear())]
    rejected=[]
    for name,mutation in mutations:
        bad=deepcopy(sample);mutation(bad)
        try:candidate.verify_certificate(sample_field,bad)
        except (AssertionError,candidate.UnavailableCertificate):rejected.append(name)
        else:raise AssertionError('accepted '+name)
    saved_path=SOURCE.parent/'results.json';saved=json.loads(saved_path.read_text())
    for name,saved_hash in saved['source_sha256'].items():
        assert sha256((SOURCE.parent/name).read_bytes()).hexdigest()==saved_hash
    readbacks=[]
    records=[saved['whole_field_extension_example']]+[x['certificate'] for x in saved['inherited_fixtures']]+[x['complete_symbolic_certificate'] for x in saved['large_interleaved']]
    for record in records:
        field=candidate.Field(record['field']['p'],record['field']['modulus']) if record['field']['e']>1 else candidate.PrimeField(record['field']['p'])
        check=candidate.verify_certificate(field,record)
        readbacks.append({'q':field.q,'n':record['n'],**check})
    assert sha256(SOURCE.read_bytes()).hexdigest()==digest
    out={'status':'passed','date':'2026-09-05','complete_polynomial_scalar_comparisons':pencils,
         'complete_nodes_compared':nodes,'whole_field_cases':whole,'entire_pencil_in_code_cases':entire,
         'complete_static_checks':static_checks,'strict_cap_rejections':unavailable,
         'exhaustive_anchor_pencils':anchor_checks,'independent_GF9_pencils':extension_checks,
         'mutations_rejected':rejected,'saved_certificate_readbacks':readbacks,
         'candidate_sha256':digest,'results_sha256':sha256(saved_path.read_bytes()).hexdigest(),
         'reviewer_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'independent_field_oracle_sha256':sha256((HERE/'review_support_batches.py').read_bytes()).hexdigest(),
         'scope':'New tiny oracles use independent field arithmetic and enumerate all polynomials/scalars. Large saved readbacks share candidate helpers and certify supplied pieces; no piece discovery claim.'}
    (HERE/'scalar_fibers_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
