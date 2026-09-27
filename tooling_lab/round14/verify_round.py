#!/usr/bin/env python3
"""Integrate complete transcript reviews; preserve the preceding frozen round."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re,sys
import numpy as np
from arc_verifier import validate

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def bind(p,h):
    assert p.is_file() and digest(p)==h,str(p);BINDINGS.append(str(p.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for x in value:scan(owner,x)
    elif isinstance(value,dict):
        for k,v in value.items():
            if k.endswith('_sha256') and isinstance(v,dict):
                for name,h in v.items():bind(owner.parent/name,h)
            else:scan(owner,v)


def main():
    previous=HERE.parent/'round13';old=json.loads((previous/'manifest.json').read_text())
    for name,record in old['artifacts'].items():bind(previous/name,record['sha256'])
    def read(name):return json.loads((HERE/name).read_text())
    # Large run files contain no binding dictionaries; their complete bytes
    # are bound by the actual, stress and independent review result files.
    for path in HERE.glob('*.json'):
        if path.name!='manifest.json' and not path.name.endswith('.run.json'):scan(path,read(path.name))
    for name in ('profile_controls.json','zero_results.json','actual_arc_results.json','arc_review.json','arc_controls.json','stress_results.json','filter_comparison.json','higher_dimension_controls.json'):
        assert read(name)['status']=='passed',name
    dimensions=read('dimension_results.json');assert dimensions['status']=='preflight_passed'
    assert [(r['dimension'],r['queries']) for r in dimensions['profiles']]==[(3,2070),(4,7245),(5,26334),(6,112860),(7,452595),(8,2090660)]
    old_controls=read('profile_controls.json')
    assert old_controls['literal_integer_partitions']==7334 and old_controls['profile_queries_checked']==2966
    assert old_controls['exhaustive_ternary_three_by_three_matrices']==19683 and old_controls['finite_field_rank_comparisons']==59049
    actual=read('actual_arc_results.json');assert len(actual['cases'])==7 and sum(len(c['runs']) for c in actual['cases'])==29
    assert [(c['review']['dimension'],c['review']['s'],c['review']['queries_rank_checked']) for c in actual['cases']]==[(3,11,10),(3,410,2070),(4,410,7245),(5,410,26334),(3,192,14685),(4,192,122500),(3,128,37080)]
    assert validate(read('small.certificate.json'))==actual['cases'][0]['review']
    review=read('arc_review.json')
    assert len(review['cases'])==29 and review['queries_independently_replayed']==419998
    assert review['candidate_words_fully_evaluated']==354645 and review['full_word_coordinate_values_checked']==362993184
    for case in actual['cases']:
        if case['name']=='small':continue
        assert all(r['output_count']==2 for r in case['runs'])
        assert case['ambient_johnson_comparison']['below_classical_agreement_threshold']==(case['review']['s'] in (128,192))
    controls=read('arc_controls.json')
    assert controls['total_received_words']==12500 and controls['total_queries_replayed']==81250
    assert len(controls['corrupted_certificates_rejected'])==14 and len(controls['corrupted_transcripts_and_filters_rejected'])==8
    assert len(controls['corrupted_boundary_test_traces_rejected'])==3 and controls['optimized_python_rejected']
    assert controls['rank_deficient_s_set_obstruction']['status']=='all_s_sets_coverage_impossible'
    assert controls['query_budget_control']['status']=='query_budget_exceeded' and controls['bounded_search_control']['status']=='not_found'
    assert validate(controls['common_root_certificate'])==controls['common_root_review']
    boundary=controls['boundary_absence_and_single_occurrence_control'];validate(boundary['certificate'])
    higher=read('higher_dimension_controls.json');assert len(higher['cases'])==4
    for row in higher['cases']:assert validate(row['certificate'])==row['certificate_review']
    stress=read('stress_results.json');assert len(stress['cases'])==3
    assert [r['known_origin_coordinate_tests'] for r in stress['cases']]==[1008]*3
    assert [r['known_origin_accepted'] for r in stress['cases']]==[True,False,True]
    assert [r['known_origin_query_occurrences'] for r in stress['cases']]==[1,1,4]
    fast=read('filter_comparison.json');assert fast['actual_decisions_matched']==354645 and fast['boundary_and_stress_runs_matched']==4
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in HERE.glob('*.md'):
        prose=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http://','https://','#')):continue
            linked=path.parent/target.split('#')[0].strip('<>')
            assert linked.exists() or linked==HERE/'manifest.json',(path,target)
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(HERE.rglob('*'))
               if p.is_file() and '__pycache__' not in p.parts and p.name not in ('manifest.json','verification.log')}
    out={'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'General-dimensional supplied-space arc covers and transcript-based exact filtering, with complete candidate review and below-Johnson boundary stress tests. No cluster discovery, uniform speedup, prize proof or historical novelty certification.',
         'previous_round13_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'small_and_boundary_certificates_rechecked':7,
         'large_review_binding':'Completed independent Cramer replay, complete candidate evaluations and stress checks are bound to all exact sources and inputs; unchanged large numerical reviews are not rerun by integration.',
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,
         'integration_runtime':{'python':sys.version.split()[0],'numpy':np.__version__},'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round14 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round13 artifacts preserved.",flush=True)


if __name__=='__main__':main()
