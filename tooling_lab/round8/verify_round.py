#!/usr/bin/env python3
"""Integrate completed covariance evidence and preserve round7."""
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


def resolve(owner,name):
    for p in (owner.parent/name,LAB/name,HERE/name):
        if p.is_file():return p
    raise AssertionError((str(owner),name))


def scan(owner,obj):
    if isinstance(obj,list):
        for x in obj:scan(owner,x)
    elif isinstance(obj,dict):
        for key,value in obj.items():
            if key.endswith('_sha256') and isinstance(value,dict):
                for name,expected in value.items():bind(resolve(owner,name),expected)
            elif key=='source_sha256' and isinstance(value,str):
                assert owner.name=='plot_metadata.json';bind(HERE/'plot_results.py',value)
            else:scan(owner,value)


def main():
    previous=json.loads((LAB/'round7/manifest.json').read_text())
    for name,row in previous['artifacts'].items():bind(LAB/'round7'/name,row['sha256'])
    for p in HERE.rglob('*.json'):
        if p.name!='manifest.json':scan(p,json.loads(p.read_text()))
    sources=list(HERE.rglob('*.py'))
    for p in sources:compile(p.read_text(),str(p),'exec')
    def read(name):return json.loads((HERE/'shared_insertions'/name).read_text())
    for name in ('results.json','scale_results.json','verification.json','boundary_verification.json','contrast_results.json','reference_results.json','coupling_verification.json'):
        assert read(name)['status']=='passed',name
    rows=read('verification.json')['cases']
    assert len(rows)==8
    assert sum(r['new_gram_full_entries_checked'] for r in rows)==584
    assert sum(r['new_gram_sampled_off_diagonal_entries_checked'] for r in rows)==48
    assert sum(r['exact_gram_row_sums_checked'] for r in rows)==206
    assert sum(r['pointwise_euler_rows_checked'] for r in rows)==15534568
    boundary=read('boundary_verification.json')
    assert len(boundary['cases'])==25 and sum(r['covariance_entries'] for r in boundary['cases'])==33480
    assert read('coupling_verification.json')['affine_covariance_entries_checked']==392
    assert read('coupling_verification.json')['literal_centered_covariance_entries']==196
    for pair in read('results.json')['pairs']:
        eq=pair['equal_features']
        assert eq['flat_histogram'] and eq['unlabelled_row_histograms'] and eq['diagonal_covariance_deck']
        assert not eq['unlabelled_column_histograms'] and not eq['covariance_trace_squared']
    broken=[]
    for p in HERE.rglob('*.md'):
        prose=re.sub(r'```.*?```|`[^`\n]*`','',p.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http://','https://','#')):continue
            linked=p.parent/re.sub(r':\d+$','',target.split('#')[0].strip('<>'))
            if linked!=HERE/'manifest.json' and not linked.exists():broken.append((str(p),target))
    assert not broken,broken
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size}
               for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts
               and p!=HERE/'manifest.json' and p!=HERE/'verification.log'}
    out={'created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed',
         'scope':'Integration of completed exact arithmetic and stated partial large-matrix readbacks; no human review, prize theorem, or historical originality certification.',
         'source_and_input_bindings':len(BINDINGS),'previous_round7_artifacts_preserved':len(previous['artifacts']),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'Round8 verification passed: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings; all {len(previous["artifacts"])} round7 artifacts preserved.',flush=True)


if __name__=='__main__':main()
