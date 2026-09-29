#!/usr/bin/env python3
"""Bind completed round7 outputs and preserve the preceding snapshot."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
LAB=HERE.parent
BINDINGS=[]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def bind(p,value):
    assert p.is_file() and digest(p)==value,str(p)
    BINDINGS.append(str(p.resolve().relative_to(LAB)))


def scan(owner,obj):
    if isinstance(obj,list):
        for x in obj:scan(owner,x)
    elif isinstance(obj,dict):
        for k,v in obj.items():
            if k.endswith('_sha256') and isinstance(v,dict):
                for name,value in v.items():bind(owner.parent/name,value)
            elif k=='source_sha256' and isinstance(v,str):
                assert owner.name=='plot_metadata.json'
                bind(HERE/'plot_results.py',v)
            else:scan(owner,v)


def main():
    previous=json.loads((LAB/'round6/manifest.json').read_text())
    for name,row in previous['artifacts'].items():bind(LAB/'round6'/name,row['sha256'])
    for p in HERE.rglob('*.json'):
        if p.name!='manifest.json':scan(p,json.loads(p.read_text()))
    sources=list(HERE.rglob('*.py'))
    for p in sources:compile(p.read_text(),str(p),'exec')
    def read(name):return json.loads((HERE/'local_edits'/name).read_text())
    for name in ('verification.json','gram_results.json','ablation_results.json','threshold_results.json','orbit_results.json','orbit_verification.json'):
        assert read(name)['status']=='passed',name
    review=read('verification.json');threshold=read('threshold_results.json')
    assert (review['scale_internal_entries_compared'],review['scale_deletion_sums_and_norms_compared'],review['equal_quartic_deck_witness_pairs'])==(7506,206,2)
    assert threshold['finite_distribution_event_bounds_checked']==8298 and threshold['invalid_queries_rejected']==11
    ablation=read('ablation_results.json')['cases']
    assert sum(r['normalized_inputs'] for r in ablation)==9373
    assert sum(r['direct_neighbour_targets'] for r in ablation)==660660
    orbit=read('orbit_verification.json')['cases']
    assert [r['square_affine_orbits'] for r in orbit]==[98,150,190]
    assert all(p['all_neighbour_histograms_equal'] for r in orbit for p in r['non_affine_pairs'])
    broken=[]
    for p in HERE.rglob('*.md'):
        text=re.sub(r'```.*?```|`[^`\n]*`','',p.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if target.startswith(('http://','https://','#')):continue
            linked=p.parent/re.sub(r':\d+$','',target.split('#')[0].strip('<>'))
            if linked!=HERE/'manifest.json' and not linked.exists():broken.append((str(p),target))
    assert not broken,broken
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
               and p!=HERE/'manifest.json' and p!=HERE/'verification.log'}
    out={'created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed',
         'scope':'Integration of completed, separately implemented exact evidence; no human review, proof of either prize, or historical novelty certification.',
         'source_and_input_bindings':len(BINDINGS),'previous_round6_artifacts_preserved':len(previous['artifacts']),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'Round7 verification passed: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings; all {len(previous["artifacts"])} round6 artifacts preserved.',flush=True)


if __name__=='__main__':main()
