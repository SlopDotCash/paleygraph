#!/usr/bin/env python3
"""Record a bounded preflight checkpoint, preserving the completed round13."""
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
    previous=HERE.parent/'round13';manifest=json.loads((previous/'manifest.json').read_text())
    for name,record in manifest['artifacts'].items():bind(previous/name,record['sha256'])
    for path in HERE.glob('*.json'):
        if path.name!='checkpoint.json':scan(path,json.loads(path.read_text()))
    results=json.loads((HERE/'dimension_results.json').read_text());assert results['status']=='preflight_passed'
    assert [(r['dimension'],r['queries']) for r in results['profiles']]==[(3,2070),(4,7245),(5,26334),(6,112860),(7,452595),(8,2090660)]
    source=json.loads((HERE.parent/'round12/full_length.certificate.json').read_text())['input']
    padded=[b+[0]*(source['k']-len(b)) for b in source['basis']]
    for row in results['cases']:
        assert row['search']['status']=='found' and row['review']['maximum_query_avoiding_set']==409
        artifact=json.loads((HERE/row['artifact']).read_text());data=artifact['input']
        assert artifact['review']==row['review'] and artifact['blocks']==row['search']['blocks']
        assert data['basis'][:3]==padded and data['origin']==source['origin']
        assert all(data[k]==source[k] for k in ('p','n','k','s','domain'))
    controls=json.loads((HERE/'profile_controls.json').read_text())
    assert controls['status']=='passed' and controls['literal_integer_partitions']==7334 and controls['profile_queries_checked']==2966
    assert controls['exhaustive_ternary_three_by_three_matrices']==19683 and controls['finite_field_rank_comparisons']==59049
    zero=json.loads((HERE/'zero_results.json').read_text())
    assert zero['status']=='passed' and zero['naive_profile']['queries']==8 and zero['review']['queries_rank_checked']==10
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
    out={'status':'preflight_checkpoint_passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Source bindings, controlled profile/determinant tests and declared-space enlargement checks. Round14 remains in progress; generic public certificate validation and decoding are pending.',
         'previous_round13_artifacts_preserved':len(manifest['artifacts']),
         'source_and_input_bindings':len(BINDINGS),'python_sources_compiled':len(sources),
         'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round14 preflight checkpoint: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(manifest['artifacts'])} round13 artifacts preserved.",flush=True)


if __name__=='__main__':main()
