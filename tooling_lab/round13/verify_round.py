#!/usr/bin/env python3
"""Integrate completed query-cover reviews and preserve the frozen round12."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re
from cover_verifier import validate

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
    previous=HERE.parent/'round12';old=json.loads((previous/'manifest.json').read_text())
    for name,record in old['artifacts'].items():bind(previous/name,record['sha256'])
    def read(name):return json.loads((HERE/name).read_text())
    for path in HERE.glob('*.json'):
        if path.name!='manifest.json':scan(path,read(path.name))
    for name in ('actual_results.json','pruning_review.json','controls.json','family_results.json','plot_metadata.json'):
        assert read(name)['status']=='passed',name
    actual=read('actual_results.json')
    for row in actual['cases']:assert validate(read(row['certificate']))==row['coverage_review']
    a,b=(r['coverage_review'] for r in actual['cases'])
    assert [a['queries_rank_checked'],a['maximum_query_avoiding_set']]==[10,10]
    assert [b['queries_rank_checked'],b['maximum_query_avoiding_set'],b['literal_local_subsets_checked']]==[2070,409,6626]
    assert len(actual['cases'][0]['runs'])==17 and len(actual['cases'][1]['runs'])==2
    replay=read('pruning_review.json')
    assert len(replay['cases'])==19 and replay['queries_replayed']==4310 and replay['candidate_evaluations_recomputed']==3240
    controls=read('controls.json')
    assert controls['literal_global_subsets_checked']==65536 and controls['toy_received_words_exhausted']==3125
    assert controls['toy_queries_independently_replayed']==31250 and controls['toy_space_members']==125
    assert len(controls['corrupted_certificates_rejected'])==15 and len(controls['corrupted_runs_rejected'])==6
    assert validate(controls['common_root_control']['certificate'])==controls['common_root_control']['review']
    family=read('family_results.json')
    assert family['spaces_attempted']==family['spaces_covered']==len(family['rows'])==2894 and family['spaces_not_found']==0
    assert family['literal_agreement_sets_checked']==12640992 and family['queries_rank_checked']==38006 and family['literal_local_subsets_checked']==185216
    for row in family['rows']:assert validate(row['certificate'])==row['review']
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
         'scope':'Deterministic coverage and actual pruning in declared polynomial spaces. Bounded construction and finite controls; no universal cluster discovery, optimum-cover theorem, prize proof or historical novelty certification.',
         'previous_round12_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'standalone_and_common_root_certificates_rechecked':3,'actual_family_certificates_rechecked':2894,
         'execution_review':'Completed Gaussian query replay and complete output oracles are bound to exact sources and inputs; the unchanged execution census is not rerun during integration.',
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,
         'plot_visually_inspected':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round13 verified: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings; all {len(old['artifacts'])} round12 artifacts preserved.",flush=True)


if __name__=='__main__':main()
