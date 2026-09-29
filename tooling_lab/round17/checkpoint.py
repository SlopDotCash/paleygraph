#!/usr/bin/env python3
"""Bind the carry/compression/realizability checkpoint and preserve round16."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import sympy as sp

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
    old=json.loads((HERE.parent/'round16/manifest.json').read_text())
    for name,row in old['artifacts'].items():bind(HERE.parent/'round16'/name,row['sha256'])
    reports={}
    for path in HERE.glob('*.json'):
        if path.name=='checkpoint.json':continue
        value=json.loads(path.read_text());scan(path,value);reports[path.stem]=value
    assert reports['norm_carry']['status']==reports['norm_compression']['status']=='produced'
    assert reports['carry_review']['status']==reports['compression_review']['status']=='passed'
    assert reports['realizability']['status']=='produced' and reports['realizability_review']['status']=='passed'
    realization=reports['realizability_review']
    assert sum(r['positive_vectors_checked'] for r in realization['cases'])==484
    assert sum(r['probes_checked'] for r in realization['cases'])==450
    assert realization['census_words_checked']==realization['census_resultants']==6561
    assert realization['actual_scalar_words_checked']==257 and len(realization['corrupt_artifacts_rejected'])==10
    assert realization['p_divides_cofactor_controls_checked']==3
    assert realization['small_summary']['oracle_actual_norm_set_accepts']==945
    assert realization['small_summary']['oracle_actual_norm_energy_pair_set_accepts']==945
    carry=reports['carry_review']['cases'];compression=reports['compression_review']['cases']
    assert sum(r['steps_checked'] for r in carry)==480 and sum(r['independent_resultants'] for r in carry)==489
    assert sum(r['vectors_checked'] for r in compression)==484
    assert [r['radius_down_norm_up_steps'] for r in carry]==[64,15,13,12,5]
    assert carry[-1]['p']==2013265921 and carry[-1]['n']==128 and carry[-1]['quartic_window']
    assert compression[-1]['digit_height_bound']==6 and compression[-1]['original_coefficient_bits']==30 and compression[-1]['digit_coefficient_bits']==3
    assert compression[-1]['original_norm_bits']==2055 and compression[-1]['digit_norm_bits']==191
    assert reports['norm_compression']['cases'][-1]['relation_cofactor']==9985208709332560769028097
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in [*HERE.glob('*.md'),HERE.parent/'NEXT_ITERATION.md']:
        prose=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('http://','https://','#')):continue
            linked=path.parent/target.split('#')[0].strip('<>');assert linked.exists(),(path,target)
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(HERE.rglob('*'))
               if p.is_file() and '__pycache__' not in p.parts and p.name not in ('checkpoint.json','checkpoint.log')}
    out={'status':'carry_compression_realizability_checkpoint_passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Exact carry cocycles, short-relation digit encodings and inverse membership, with separate resultant and polynomial inverse reviews. Complete small digit-cube census and selected larger probes; no uniform norm or shell theorem is claimed.',
         'previous_round16_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'sympy_version':sp.__version__,'markdown_local_links_passed':True,
         'carry_transitions_checked':480,'carry_resultants_checked':489,'compressed_vector_records_checked':484,
         'realizability_census_checked':6561,'realizability_probes_checked':450,'artifacts':artifacts}
    (HERE/'checkpoint.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round17 checkpoint: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round16 artifacts preserved.",flush=True)


if __name__=='__main__':main()
