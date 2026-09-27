#!/usr/bin/env python3
"""Bind completed reviews, recheck small interfaces, preserve the prior snapshot.

The expensive full-length verifier and random-query replay already completed.
This integration checks their exact source/input bindings rather than rerunning
unchanged expensive work. It does not replace those numerical reviews.
"""
from datetime import datetime,timezone
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
from certificate_verifier import validate

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def bind(path,h):
    assert path.is_file() and digest(path)==h,str(path)
    BINDINGS.append(str(path.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for item in value:scan(owner,item)
    elif isinstance(value,dict):
        for key,item in value.items():
            if key.endswith('_sha256') and isinstance(item,dict):
                for name,h in item.items():bind(owner.parent/name,h)
            else:scan(owner,item)


def main():
    previous=HERE.parent/'round11';old=json.loads((previous/'manifest.json').read_text())
    for name,record in old['artifacts'].items():bind(previous/name,record['sha256'])
    def read(name):return json.loads((HERE/name).read_text())
    for path in HERE.glob('*.json'):
        if path.name!='manifest.json':scan(path,read(path.name))
    expected={'preflight_results.json':'toy_passed','actual_results.json':'actual_inputs_passed',
              'pair_space_results.json':'full_input_preflight_passed','pair_space_controls.json':'passed',
              'certificate_verification.json':'passed','pruning_results.json':'passed',
              'pruning_verification.json':'passed','pruning_controls.json':'passed'}
    for name,status in expected.items():assert read(name)['status']==status,name
    pre=read('preflight_results.json')
    assert [pre[k] for k in ('patterns','threshold_checks','literal_subsets','partial_class_witnesses')]==[596,3936,45952,8]
    actual=read('actual_results.json')
    assert [actual[k] for k in ('rank_three_quadruples','dependent_quadruples','unique_labeled_dependency_catalogues','literal_agreement_sets_evaluated','coordinate_triples_checked','previous_simple_cases_matched','minimum_over_all_rank_three_spaces')]==[2894,166,1825,7971600,1620640,55,50]
    reviews=read('certificate_verification.json')
    assert len(reviews['cases'])==6 and len(reviews['corrupted_certificates_rejected'])==12 and reviews['verification_limit_is_not_success']
    for row in reviews['cases'][:5]:
        checked=validate(read(row['certificate']))
        assert all(row[k]==v for k,v in checked.items()),row['name']
    large=reviews['cases'][5]
    assert large['minimum_interval']==[11401439,11402277] and large['all_pairs_checked']==523776
    assert large['all_triple_dependencies_covered']==178433024 and not large['literal_triples_enumerated']
    assert large['success_probability_lower']==[1628777,25490432]
    assert len(read('pair_space_controls.json')['cases'])==73
    runs=read('pruning_results.json');replay=read('pruning_verification.json');controls=read('pruning_controls.json')
    assert len(runs['small_runs'])==17 and len(runs['large_runs'])==2
    assert runs['large_certificate_trial_budget']['trials']==945 and runs['large_universal_baseline_trial_budget']['trials']==1554
    assert replay['queries_replayed']==9234 and replay['candidate_evaluations_recomputed']==2619
    assert replay['small_space_members_enumerated']==4913 and replay['large_global_degree_root_bound_checked']
    assert Fraction(*replay['sum_of_run_failure_bounds'])==Fraction(49,2**46)<Fraction(1,2**40)
    assert controls['all_query_solutions_compared']==9520 and len(controls['corrupted_traces_rejected'])==8
    assert controls['repeated_queries_do_not_inherit_uniform_sampling_guarantee'] and controls['zero_probability_rejected']
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in HERE.glob('*.md'):
        prose=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('https://','http://','#')):continue
            linked=path.parent/target.split('#')[0].strip('<>')
            assert linked.exists() or linked==HERE/'manifest.json',(path,target)
    artifacts={str(path.relative_to(HERE)):{'sha256':digest(path),'bytes':path.stat().st_size}
               for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
               and path.name not in ('manifest.json','verification.log')}
    out={'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Multiplicity-aware and full-length declared-space certificates plus verified actual pruning. No cluster discovery, general prize proof, historical novelty certification or end-to-end speedup claim.',
         'previous_round11_artifacts_preserved':len(old['artifacts']),
         'source_and_input_bindings':len(BINDINGS),'small_certificates_rechecked':5,
         'large_review_and_pruning_replay':'Earlier completed numerical reviews are bound to all exact sources and inputs; unchanged expensive reviews are not rerun by this integration.',
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round12 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round11 artifacts preserved.",flush=True)


if __name__=='__main__':main()
