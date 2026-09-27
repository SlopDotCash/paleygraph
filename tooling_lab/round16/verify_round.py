#!/usr/bin/env python3
"""Freeze the full conditional-decoding round and preserve frozen round15."""
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
    old=json.loads((HERE.parent/'round15/manifest.json').read_text())
    for name,row in old['artifacts'].items():bind(HERE.parent/'round15'/name,row['sha256'])
    reports={}
    for path in HERE.glob('*.json'):
        if path.name in ('checkpoint.json','manifest.json'):continue
        value=json.loads(path.read_text());scan(path,value);reports[path.stem]=value
    for name in ['conditional_review','toy_controls','basis_controls','routing_controls','decoding_controls','decoding_review','guard_review']:
        assert reports[name]['status']=='passed',name
    assert reports['conditional_decoding']['status']==reports['guard_results']['status']=='produced'
    assert reports['toy_controls']['received_words_exhausted']==3125
    assert [r['received_words_exhausted'] for r in reports['basis_controls']['cases']]==[3125,3125]
    assert reports['routing_controls']['supported_cutoffs_tested']==250
    assert sum(r['supported_cutoff_decodings'] for r in reports['decoding_controls']['cases'])==750
    assert reports['decoding_controls']['corrupt_decoding_artifacts_rejected']==4
    original=reports['decoding_review'];guarded=reports['guard_review']
    assert [r['metrics']['query_systems_solved'] for r in original['cases']]==[71113,71113,70892]
    assert [r['metrics']['total_explicit_coordinate_evaluations'] for r in original['cases']]==[452949,450835,1467372]
    assert [r['metrics']['query_systems_solved'] for r in guarded['cases']]==[71520,71520,73554]
    assert [r['metrics']['total_explicit_coordinate_evaluations'] for r in guarded['cases']]==[74419,72345,76276]
    assert all(r['both_complete_supplied_space_lists_matched'] and r['review']['query_replay_performed'] for r in guarded['cases'])
    assert all(r['baseline_complete_supplied_space_list_matched'] and r['review']['query_replay_performed'] for r in original['cases'])
    totals={k:original[k]+guarded[k] for k in ['queries_independently_replayed','candidate_words_fully_evaluated','candidate_coordinate_values_checked']}
    assert totals=={'queries_independently_replayed':429712,'candidate_words_fully_evaluated':429559,'candidate_coordinate_values_checked':413235758}
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
         'scope':'Complete supplied-space conditional decoding, a measured verification regression and its singleton-profile repair. Independent query/full-word review, affine/routing controls and preserved baseline artifacts. No historical novelty, full ambient decoder or prize proof claim.',
         'previous_round15_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'full_scale_review_totals':totals,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round16 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round15 artifacts preserved.",flush=True)


if __name__=='__main__':main()
