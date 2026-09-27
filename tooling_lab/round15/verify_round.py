#!/usr/bin/env python3
"""Freeze round15 evidence after checking bindings and preserving round14."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(path):return sha256(path.read_bytes()).hexdigest()


def bind(path,h):
    assert path.is_file() and digest(path)==h,str(path);BINDINGS.append(str(path.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for v in value:scan(owner,v)
    elif isinstance(value,dict):
        for key,v in value.items():
            if key.endswith('_sha256') and isinstance(v,dict):
                for name,h in v.items():bind(owner.parent/name,h)
            else:scan(owner,v)


def main():
    old=json.loads((HERE.parent/'round14/manifest.json').read_text())
    for name,row in old['artifacts'].items():bind(HERE.parent/'round14'/name,row['sha256'])
    for path in HERE.glob('*.json'):
        if path.name not in ('manifest.json','checkpoint.json'):scan(path,json.loads(path.read_text()))
    root=json.loads((HERE/'root_class_results.json').read_text());search=json.loads((HERE/'capacity_search.json').read_text())
    repair=json.loads((HERE/'exchange_repair.json').read_text());review=json.loads((HERE/'capacity_review.json').read_text())
    controls=json.loads((HERE/'search_controls.json').read_text());parallel=json.loads((HERE/'parallel_controls.json').read_text())
    boundary=json.loads((HERE/'boundary_results.json').read_text())
    assert root['status']=='parallel_preflight_passed' and not root['all_triple_ranks_checked']
    assert search['status']==repair['status']=='found' and review['status']==controls['status']==parallel['status']==boundary['status']=='passed'
    assert len(search['attempts'])==23 and len(repair['trials'])==2 and repair['local_queries_checked']==1682
    assert review['minimum_complete_arc_partition_queries']==136155 and review['feasible_t_u_profiles_checked']==173
    assert review['failed_restart_witnesses_checked']==46 and review['restart_search_queries_reported']==3131565
    assert controls['literal_partitions_checked']==1120 and controls['corrupt_exchange_histories_rejected']==4
    assert parallel['literal_coordinate_set_partitions']==33841 and parallel['threshold_queries']==417
    expected={'queries_independently_replayed':408465,'candidate_words_fully_evaluated':408351,'full_word_coordinate_values_checked':418151424}
    for key,v in expected.items():assert sum(r['review'][key] for r in boundary['cases'])==v
    assert [r['review']['outputs_checked'] for r in boundary['cases']]==[1,0,1]
    assert all(r['known_origin_coordinate_tests']==993 for r in boundary['cases'])
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in [*HERE.glob('*.md'),HERE.parent/'NEXT_ITERATION.md']:
        prose=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http://','https://','#')):continue
            linked=path.parent/target.split('#')[0].strip('<>')
            assert linked.exists() or linked==HERE/'manifest.json',(path,target)
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(HERE.rglob('*'))
               if p.is_file() and '__pycache__' not in p.parts and p.name not in ('manifest.json','verification.log')}
    out={'status':'verified','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Projective capacity obstruction and finite complete-block optimum, checked local exchange repair, complete supplied-space boundary decoding, and controls. Known mathematical ingredients; no historical novelty or prize proof claim.',
         'previous_round14_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'boundary_review_totals':expected,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round15 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round14 artifacts preserved.",flush=True)


if __name__=='__main__':main()
