#!/usr/bin/env python3
"""Export and independently review actual polynomial-cluster certificates."""
from copy import deepcopy
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from polynomial_cluster import certify
from independent_verifier import validate_certificate,determinant_partition,integer_optimum

HERE=Path(__file__).resolve().parent


def main():
    previous=json.loads((HERE/'actual_results.json').read_text());smallpath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json'
    source=json.loads(smallpath.read_text());largepath=HERE.parent/'round6/pencil_tracks/results.json'
    big=next(r['result'] for r in json.loads(largepath.read_text())['large'] if r['name']=='skew-two')
    nodes=source['nodes'];first=previous['small']['worst_triple_node_indices'][0]
    inputs=[{'p':source['p'],'k':source['k'],'s':source['s'],'domain':source['domain'],
             'origin':nodes[first]['coefficients'],'basis':previous['small']['directions']},
            {'p':big['field']['p'],'k':big['k'],'s':big['s'],'domain':big['domain'],
             'origin':previous['large']['origin'],'basis':previous['large']['directions']}]
    reviews=[];certificates=[]
    for name,data in zip(('hard_f17','large_f65537'),inputs):
        certificate=certify(data);review=validate_certificate(certificate);reviews.append(review);certificates.append(certificate)
        (HERE/f'{name}.input.json').write_text(json.dumps(data,indent=2)+'\n')
        (HERE/f'{name}.certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    # Replay every candidate triple independently using stored codewords, not
    # polynomial normalization or the constructor's greedy occupancy rule.
    values=[];dependent=0;p=source['p']
    for triple in combinations(nodes,3):
        columns=[((b-a)%p,(c-a)%p) for a,b,c in zip(*(t['codeword'] for t in triple))]
        groups,zeros,_=determinant_partition(columns,p)
        if len(groups)<2:dependent+=1;continue
        value,_=integer_optimum([len(g) for g in groups],len(zeros),source['s']);values.append(value)
    assert len(values)==804 and dependent==12 and min(values)==15
    corruptions=0
    for field in ('minimum','probability','class','basis','domain','witness','direction'):
        bad=deepcopy(certificates[0])
        if field=='minimum':bad['pinning']['minimum_injective_pairs']+=1
        elif field=='probability':bad['pinning']['uniform_pair_success_probability'][0]+=1
        elif field=='class':bad['pinning']['projective_classes'][0]['coordinates'].pop()
        elif field=='basis':bad['basis'][1]=bad['basis'][0][:]
        elif field=='domain':bad['domain'][0]=bad['domain'][1]
        elif field=='witness':bad['pinning']['worst_agreement_set'].pop()
        else:bad['pinning']['projective_classes'][0]['direction']=[0,0]
        try:validate_certificate(bad)
        except (AssertionError,ValueError):corruptions+=1
        else:raise AssertionError(('corruption accepted',field))
    out={'status':'passed','scope':'Constructor-free all-pair determinant partition and integer dynamic-program optimum, including the full1024-coordinate certificate. Every804 actual small candidate planes independently replayed.',
         'cases':reviews,'small_rank_two_planes_reviewed':len(values),'small_dependent_triples':dependent,
         'corrupted_certificates_rejected':corruptions,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_certificates.py','polynomial_cluster.py','rank_two_pinning.py','independent_verifier.py')},
         'input_sha256':{'actual_results.json':sha256((HERE/'actual_results.json').read_bytes()).hexdigest(),
                         '../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(smallpath.read_bytes()).hexdigest(),
                         '../round6/pencil_tracks/results.json':sha256(largepath.read_bytes()).hexdigest()},
         'certificate_sha256':{f'{name}.certificate.json':sha256((HERE/f'{name}.certificate.json').read_bytes()).hexdigest() for name in ('hard_f17','large_f65537')}}
    (HERE/'certificate_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('sha256')}),flush=True)


if __name__=='__main__':main()
