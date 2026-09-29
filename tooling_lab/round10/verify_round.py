#!/usr/bin/env python3
"""Integrate exact pair-mean evidence and preserve the preceding snapshot."""
from datetime import datetime,timezone
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(p):return sha256(p.read_bytes()).hexdigest()


def bind(path,h):
    assert path.is_file() and digest(path)==h,str(path);BINDINGS.append(str(path.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for x in value:scan(owner,x)
    elif isinstance(value,dict):
        for k,v in value.items():
            if k.endswith('_sha256') and isinstance(v,dict):
                for name,h in v.items():bind(owner.parent/name,h)
            else:scan(owner,v)


def main():
    old=json.loads((HERE.parent/'round9/manifest.json').read_text())
    for name,r in old['artifacts'].items():bind(HERE.parent/'round9'/name,r['sha256'])
    for p in HERE.glob('*.json'):
        if p.name!='manifest.json':scan(p,json.loads(p.read_text()))
    def read(name):return json.loads((HERE/name).read_text())
    for name in ('boundary_verification.json','scale_results.json','scale_verification.json','projection_toy.json','projection_scale.json','controls_verification.json','shell_mean_verification.json'):
        assert read(name)['status']=='passed',name
    boundary=read('boundary_verification.json');assert len(boundary['cases'])==35 and boundary['all_pair_means_compared']==14855 and boundary['literal_inward_targets_compared']==2229
    scale=read('scale_results.json')['cases'];review=read('scale_verification.json')['cases']
    assert len(scale)==len(review)==8 and sum(len(r['pairs']) for r in scale)==3650
    for r,v in zip(scale,review):
        assert r['q']==v['q'] and r['selected']==v['inward_targets']['selected']
        assert len(r['pairs'])==v['all_pair_means_checked']==len(v['decomposition']['pairs'])
        for pair,other in zip(r['pairs'],v['decomposition']['pairs']):
            assert pair['deleted']==other['deleted'] and pair['target_sum']==other['target_sum']
            assert F(*pair['mean'])==F(pair['target_sum'],r['ordered_insertion_pairs_per_deletion'])
    controls=read('controls_verification.json');assert controls['all_pair_entries_reindexed']==525 and controls['invalid_queries_rejected']==8
    ib=controls['integer_bounds'];assert ib['conservative_accumulator_bound']<ib['signed128_limit'] and ib['inward_target_bound']<ib['signed64_limit']
    toy=read('projection_toy.json')['cases'];assert sum(c['square_affine_representatives'] for c in toy)==438
    assert [r['full_affine_orbits'] for r in toy]==[49,75,95]
    assert [r['mean_deck_fibres'] for r in toy]==[47,74,95]
    assert [r['projected_norm_fibres'] for r in toy]==[41,63,95]
    assert [r['all_projected_means_explained_by_lower_degrees_and_induced_graph'] for r in toy]==[True,True,False]
    assert toy[-1]['projected_norm_fibres_joining_distinct_full_affine_orbits']==0
    assert read('shell_mean_verification.json')['cases'][-1]['distance_two_final_sets']==27498133131022225
    assert read('preflight_results.json')['independent_mean_comparisons']==183
    sources=list(HERE.glob('*.py'))
    for p in sources:compile(p.read_text(),str(p),'exec')
    for p in HERE.glob('*.md'):
        prose=re.sub(r'```.*?```|`[^`\n]*`','',p.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',prose):
            if target.startswith(('https://','http://','#')):continue
            path=p.parent/re.sub(r':\d+$','',target.split('#')[0].strip('<>'))
            assert path.exists() or path==HERE/'manifest.json',(p,target)
    assert not list(HERE.glob('*partial.json'))
    artifacts={str(p.relative_to(HERE)):{'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(HERE.rglob('*'))
               if p.is_file() and '__pycache__' not in p.parts and p.name not in ('manifest.json','verification.log')}
    out={'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Completed query compiler and constructive information audit. Every large pair mean independently checked; no historical originality or prize theorem certification.',
         'previous_round9_artifacts_preserved':len(old['artifacts']),'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f"Round10 verified: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings; all {len(old['artifacts'])} round9 artifacts preserved.",flush=True)


if __name__=='__main__':main()
