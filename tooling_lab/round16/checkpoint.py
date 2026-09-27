#!/usr/bin/env python3
"""Bind the conditional-cover checkpoint and preserve frozen round15."""
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
    for path in HERE.glob('*.json'):
        if path.name!='checkpoint.json':scan(path,json.loads(path.read_text()))
    compiled=json.loads((HERE/'conditional_results.json').read_text());review=json.loads((HERE/'conditional_review.json').read_text())
    toy=json.loads((HERE/'toy_controls.json').read_text());basis=json.loads((HERE/'basis_controls.json').read_text());routing=json.loads((HERE/'routing_controls.json').read_text())
    assert compiled['status']=='conditional_covers_compiled' and review['status']==toy['status']==basis['status']==routing['status']=='passed'
    assert [r['query_count'] for r in compiled['cases']]==[71113,71113,70892]
    assert [r['review']['queries_rank_checked'] for r in review['cases']]==[71113,71113,70892]
    assert all(r['review']['routing_complete_in_supplied_space'] and r['review']['branches_checked']==2 for r in review['cases'])
    assert toy['received_words_exhausted']==3125 and toy['complete_space_members_per_word']==125 and toy['constituent_queries_checked']==4375
    assert toy['corrupt_conditional_certificates_rejected']==4
    assert [r['received_words_exhausted'] for r in basis['cases']]==[3125,3125]
    assert [r['direction'] for r in basis['cases']]==[[0,1,0],[1,2,3]]
    assert routing['supported_cutoffs_tested']==250 and routing['mixed_light_and_fibre_covers']==25 and routing['multiple_fibre_covers']==100 and routing['corrupt_routing_certificates_rejected']==3
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
    out={'status':'conditional_cover_checkpoint_passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Conditional routing and all residual query certificates independently verified; exhaustive toy decoding and affine-lift controls pass. Full-scale conditional candidate decoding is pending; round16 remains active.',
         'previous_round15_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,
         'large_constituent_queries_independently_checked':sum(r['review']['queries_rank_checked'] for r in review['cases']),
         'exhaustive_received_word_presentation_cases':toy['received_words_exhausted']+sum(r['received_words_exhausted'] for r in basis['cases']),
         'mixed_routing_cover_cases':routing['supported_cutoffs_tested'],'artifacts':artifacts}
    (HERE/'checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round16 checkpoint: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round15 artifacts preserved.",flush=True)


if __name__=='__main__':main()
