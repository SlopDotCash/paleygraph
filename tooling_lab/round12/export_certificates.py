#!/usr/bin/env python3
"""Export full polynomial certificates and run the separate verifier."""
from copy import deepcopy
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter
from polynomial_certificate import certify
from certificate_verifier import validate,VerificationLimit

HERE=Path(__file__).resolve().parent


def difference(a,b,p,k):return [((a[j] if j<len(a) else 0)-(b[j] if j<len(b) else 0))%p for j in range(k)]


def multiply(a,b,p):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return c


def main():
    actual_path=HERE/'actual_results.json';actual=json.loads(actual_path.read_text());worst=min(actual['rows'],key=lambda r:r['literal_minimum'])
    small_path=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';small=json.loads(small_path.read_text());ids=worst['node_indices']
    a=small['nodes'][ids[0]]['coefficients'];basis=[difference(small['nodes'][i]['coefficients'],a,17,8) for i in ids[1:]]
    hard={'p':17,'k':8,'s':11,'domain':small['domain'],'origin':a,'basis':basis}
    P=[1]
    for x in range(5):P=multiply(P,[(-x)%17,1],17)
    common={'p':17,'k':8,'s':11,'domain':list(range(16)),'basis':[P,[0]+P,[0,0]+P]}
    parallel={'p':17,'k':8,'s':11,'domain':list(range(1,17)),'basis':[[1],[0,0,1],[0,0,0,0,1]]}
    simple={'p':17,'k':8,'s':11,'domain':list(range(1,17)),'basis':[[1],[0,1],[0,0,1,0,0,0,1]]}
    large_path=HERE.parent/'round11/rank_three_actual.json';source=json.loads(large_path.read_text())['large_declared_space'];large={x:source[x] for x in ('p','k','s','domain','origin','basis')}
    cases=[('hard_rank3',hard,'auto',50),('common_roots',common,'auto',20),('parallel_conic',parallel,'auto',120),
           ('simple_corner',simple,'capacitated_corners',147),('simple_pair',simple,'simple_pair_degree',None),
           ('full_length',large,'simple_pair_degree',None)]
    reviews=[];certificates=[];paths=[]
    for name,data,method,expected in cases:
        start=perf_counter();c=certify(data,method);compiled=perf_counter()-start
        start=perf_counter();v=validate(c);reviewed=perf_counter()-start
        if expected is not None:assert v['minimum_interval']==[expected]*2
        if name=='full_length':assert v['minimum_interval']==[11401439,11402277]
        path=HERE/f'{name}.certificate.json';path.write_text(json.dumps(c,separators=(',',':'))+'\n');paths.append(path.name)
        input_path=HERE/f'{name}.input.json';input_path.write_text(json.dumps(data,indent=2)+'\n');paths.append(input_path.name)
        reviews.append({'name':name,'certificate':path.name,'input':input_path.name,'compile_seconds':compiled,'review_seconds':reviewed,**v});certificates.append(c)
        print(json.dumps(reviews[-1]),flush=True)
    invalid=[]
    for label in ('domain','origin_degree','basis_rank','class_missing_coordinate','class_direction','class_triples','occupancy','minimum','probability','line_missing','pair_degree','pair_interval'):
        c=deepcopy(certificates[4] if label.startswith(('line_','pair_')) else certificates[0])
        if label=='domain':c['input']['domain'][0]=c['input']['domain'][1]
        elif label=='origin_degree':c['input']['origin']+=[0]*9
        elif label=='basis_rank':c['input']['basis'][2]=c['input']['basis'][1]
        elif label=='class_missing_coordinate':c['proof']['projective_classes'][0]['coordinates'].pop()
        elif label=='class_direction':c['proof']['projective_classes'][0]['direction']=[0,0,0]
        elif label=='class_triples':c['proof']['independent_class_triples'].pop()
        elif label=='occupancy':c['proof']['class_occupancies'][0]+=1
        elif label=='minimum':c['proof']['minimum_injective_triples']+=1
        elif label=='probability':c['proof']['uniform_triple_success_probability'][0]+=1
        elif label=='line_missing':c['proof']['nontrivial_maximal_lines'].pop()
        elif label=='pair_degree':c['proof']['bound']['dependency_vertex_degrees'][0]+=1
        else:c['proof']['bound']['minimum_injective_triples_interval'][0]+=1
        try:validate(c)
        except (AssertionError,ValueError,IndexError):invalid.append(label)
        else:raise AssertionError(('corruption accepted',label))
    try:validate(certificates[0],max_states=1)
    except VerificationLimit:limited=True
    else:raise AssertionError('a verification budget exhaustion was silently accepted')
    out={'status':'passed','scope':'Standalone polynomial certificates verified without constructor imports. Capacity optima checked by all-integer-occupancy search; full simple input checked with a separate pair grouping and row reduction. Mathematical originality and cluster discovery are not asserted.',
         'cases':reviews,'hard_source_node_indices':ids,'corrupted_certificates_rejected':invalid,'verification_limit_is_not_success':limited,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('export_certificates.py','polynomial_certificate.py','certificate_verifier.py','capacitated_pinning.py','pair_space_preflight.py')},
         'input_sha256':{'actual_results.json':sha256(actual_path.read_bytes()).hexdigest(),
                         '../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(small_path.read_bytes()).hexdigest(),
                         '../round11/rank_three_actual.json':sha256(large_path.read_bytes()).hexdigest()},
         'certificate_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in paths}}
    (HERE/'certificate_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(reviews),'corrupted_certificates_rejected':len(invalid)}),flush=True)


if __name__=='__main__':main()
