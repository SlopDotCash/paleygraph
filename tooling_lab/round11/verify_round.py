#!/usr/bin/env python3
"""Integrate pinning certificates and preserve the frozen preceding snapshot."""
from datetime import datetime,timezone
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import re
from independent_verifier import validate_certificate
from rank_three_verifier import validate

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def bind(path,h):
    assert path.is_file() and digest(path)==h,str(path);BINDINGS.append(str(path.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for item in value:scan(owner,item)
    elif isinstance(value,dict):
        for key,item in value.items():
            if key.endswith('_sha256') and isinstance(item,dict):
                for name,h in item.items():bind(owner.parent/name,h)
            else:scan(owner,item)


def main():
    previous=HERE.parent/'round10';old=json.loads((previous/'manifest.json').read_text())
    for name,record in old['artifacts'].items():bind(previous/name,record['sha256'])
    for path in HERE.glob('*.json'):
        if path.name!='manifest.json':scan(path,json.loads(path.read_text()))
    def read(name):return json.loads((HERE/name).read_text())
    for name in ('certificate_verification.json','rank_two_controls.json','rank_three_results.json','rank_three_review.json','rank_three_scale.json','rank_three_controls.json','rank_three_actual.json','collision_moments.json','plot_metadata.json'):
        assert read(name)['status']=='passed',name
    preflight=read('preflight_results.json');assert preflight['projective_patterns']==3000 and preflight['threshold_checks']==18000 and preflight['agreement_sets_enumerated']==96000
    reviews=read('certificate_verification.json');assert reviews['small_rank_two_planes_reviewed']==804 and reviews['small_dependent_triples']==12 and reviews['corrupted_certificates_rejected']==7
    for name,record in zip(('hard_f17','large_f65537'),reviews['cases']):assert validate_certificate(read(f'{name}.certificate.json'))==record
    assert reviews['cases'][1]['minimum']==83836 and reviews['cases'][1]['full_coordinate_pairs_checked']==523776
    r2=read('rank_two_controls.json');assert len(r2['polynomial_controls'])==8 and len(r2['abstract_scaling_controls'])==2 and len(r2['common_root_extremal_cases'])==38 and len(r2['invalid_inputs_rejected'])==8
    pre=read('rank_three_results.json');assert pre['polynomial_subspaces']==246 and pre['agreement_sets_per_subspace']==4368
    assert [a['ambiguous_fibres'] for a in pre['ablations']]==[1,2,1,0]
    collision=read('rank_three_review.json');assert collision['all_subsets_reviewed']==131072 and collision['all_codewords_reviewed']==9826 and len(collision['controls'])==6
    a,b=collision['rows'];assert a['minimum']==147 and b['minimum']==148
    for key in ('subset_rank_distribution','word_weight_distribution','generalized_hamming_weights'):assert a[key]==b[key]
    assert sum(x[2] for x in a['subset_rank_distribution'])==65536 and sum(count for weight,count in a['word_weight_distribution'])==4913
    moments=read('collision_moments.json');assert moments['moment_equality_1_through_6']==[True,True,False,False,False,False]
    scale=read('rank_three_scale.json');assert len(scale['rows'])==10 and sum(r['status']=='exact' for r in scale['rows'])==4
    total_triples=nodes=0
    for row in scale['rows']:
        checked=validate(read(row['certificate']))
        assert all(row[key]==value for key,value in checked.items())
        assert row['interval'][0]>=comb(row['s']-row['k']+3,3)
        total_triples+=row['all_triples_checked'];nodes+=row['tree_nodes_replayed']
    assert validate(read('rank3_frontier_control.certificate.json'))==scale['frontier_control'] and scale['frontier_control']['status']=='bounded'
    assert len(scale['corrupted_certificates_rejected'])==7
    controls=read('rank_three_controls.json');assert len(controls['boundary_cases'])==70 and controls['literal_subsets']==2760 and len(controls['transformed_cases'])==6 and len(controls['invalid_inputs_rejected'])==8
    actual=read('rank_three_actual.json');assert actual['small_classification']=={'dependent':166,'simple_rank_three':55,'nonsimple_rank_three':2839}
    assert len(actual['small_cases'])==55 and actual['small_agreement_sets_enumerated']==240240 and actual['small_minimum_over_simple_clusters']==120
    for row in actual['small_cases']:
        checked=validate(row['certificate']);assert checked==row['independent_review']
        assert checked['interval']==[row['literal_review']['minimum']]*2
    profile=actual['large_declared_space']['evaluation_profile'];assert profile['rank']==3 and profile['simple'] and profile['projective_classes']==1024
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in HERE.glob('*.md'):
        prose=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('https://','http://','#')):continue
            linked=path.parent/re.sub(r':\d+$','',target.split('#')[0].strip('<>'))
            assert linked.exists() or linked==HERE/'manifest.json',(path,target)
    artifacts={str(path.relative_to(HERE)):{'sha256':digest(path),'bytes':path.stat().st_size} for path in sorted(HERE.rglob('*'))
               if path.is_file() and '__pycache__' not in path.parts and path.name not in ('manifest.json','verification.log')}
    out={'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Declared-space pinning tools, exact finite information obstruction and explicitly bounded rank-three scaling. No cluster coverage, prize proof or historical novelty certification.',
         'previous_round10_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'rank_three_scale_triples_rechecked':total_triples,'rank_three_scale_tree_nodes_replayed':nodes,
         'rank_three_actual_certificates_replayed':55,'rank_two_actual_certificates_replayed':2,
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round11 verified: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings; all {len(old['artifacts'])} round10 artifacts preserved.",flush=True)


if __name__=='__main__':main()
