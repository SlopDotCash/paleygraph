#!/usr/bin/env python3
"""Bind the projective-capacity preflight and preserve frozen round14."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re

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
    previous=HERE.parent/'round14';old=json.loads((previous/'manifest.json').read_text())
    for name,record in old['artifacts'].items():bind(previous/name,record['sha256'])
    for path in HERE.glob('*.json'):
        if path.name!='checkpoint.json':scan(path,json.loads(path.read_text()))
    root=json.loads((HERE/'root_class_results.json').read_text());controls=json.loads((HERE/'parallel_controls.json').read_text())
    assert root['status']=='parallel_preflight_passed' and not root['all_triple_ranks_checked']
    assert root['input']['n']==1024 and root['input']['s']==96 and root['full_rank_witness']['determinant_mod_p']!=0
    assert sorted(len(g) for g in root['projective_classes'])==[1]*962+[62]
    assert root['unconstrained_profile']['queries']==70070
    best=root['parallel_constrained_profile']['best']
    assert best['queries']==136155 and best['queried_blocks']==33 and best['unqueried_coordinates']==29
    assert controls['status']=='passed' and controls['capacity_patterns']==47
    assert controls['literal_coordinate_set_partitions']==33841 and controls['threshold_queries']==417
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in [*HERE.glob('*.md'),HERE.parent/'NEXT_ITERATION.md']:
        prose=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http://','https://','#')):continue
            linked=path.parent/target.split('#')[0].strip('<>')
            assert linked.exists() or linked==HERE/'checkpoint.json',(path,target)
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(HERE.rglob('*'))
               if p.is_file() and '__pycache__' not in p.parts and p.name not in ('checkpoint.json','checkpoint.log')}
    out={'status':'parallel_preflight_checkpoint_passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Actual-domain polynomial projective-capacity obstruction, relaxed optimal profile and verified class allocation. Higher-rank independence and decoding remain pending; round15 stays active.',
         'previous_round14_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round15 preflight checkpoint: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round14 artifacts preserved.",flush=True)


if __name__=='__main__':main()
