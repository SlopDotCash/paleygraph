#!/usr/bin/env python3
"""Bind compressed covers, orbit counts and the integral realizability gap."""
from datetime import datetime,timezone
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import re

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(path):return sha256(path.read_bytes()).hexdigest()


def bind(path,wanted):
    assert path.is_file() and digest(path)==wanted,str(path);BINDINGS.append(str(path.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for v in value:scan(owner,v)
    elif isinstance(value,dict):
        for key,v in value.items():
            if key.endswith('_sha256') and isinstance(v,dict):
                for name,wanted in v.items():bind(owner.parent/name,wanted)
            else:scan(owner,v)


def main():
    preserved={}
    for name,filename in [('round22','manifest.json'),('round21','manifest.json'),('round20','manifest.json'),
                          ('round19','manifest.json'),('round18','manifest.json'),('round17','checkpoint.json'),('round16','manifest.json')]:
        path=HERE.parent/name/filename;data=json.loads(path.read_text())
        for rel,item in data['artifacts'].items():bind(path.parent/rel,item['sha256'])
        preserved[name]=len(data['artifacts'])
    reports={}
    for path in HERE.glob('*.json'):
        if path.name=='manifest.json':continue
        data=json.loads(path.read_text());scan(path,data);reports[path.stem]=data
    for n in ['preflight','cover_summary','small_cover_summary','query_experiments','lattice_information','orbit_summary','gap_certificate','plot_metadata']:
        assert reports[n]['status']=='produced'
    for n in ['preflight_review','cover_review','refinement33_review','refinement40_review','cover36_review',
              'query_review','query36_review','lattice_information_review','orbit_review','gap_review']:
        assert reports[n]['status']=='passed'
    covers=reports['cover_review']['cases']
    for key,v in [('nodes_verified',5748),('cuts_verified',298),('exact_narrowing_steps',24238),
                  ('unresolved_leaves',220),('actual_small_ambiguity_pairs_preserved',550),('continuous_witnesses_checked',83)]:
        assert sum(c[key] for c in covers)==v
    assert covers[0]['universal_unique_completion'] and covers[0]['nodes_verified']==131
    assert not covers[1]['universal_unique_completion'] and covers[1]['unresolved_leaves']==1
    assert covers[2]['universal_unique_completion'] and covers[2]['nodes_verified']==1867
    assert not covers[3]['universal_unique_completion'] and covers[3]['unresolved_integer_points_up_to_sign']==177723199
    assert len(reports['cover_review']['corrupt_coverings_rejected'])==4
    for name,unique,unresolved in [('refinement33_review',True,0),('refinement40_review',False,19)]:
        case=reports[name]['cases'][0]
        assert case['universal_unique_completion']==unique and case['unresolved_leaves']==unresolved
        assert case['added_nodes']==case['added_cuts']==0
    preflight=reports['preflight_review']['cases']
    for key,v in [('queries_replayed',48),('scalar_candidates_checked',225813),('cuts_verified',109),('optimality_certificates',109)]:
        assert sum(c[key] for c in preflight)==v
    query=reports['query_review'];assert query['cyclic_queries_checked']==8 and query['scalar_candidates_checked']==36 and len(query['controls'])==4
    query36=reports['query36_review'];assert query36['erased_coordinates']==36 and query36['scalar']==1234567
    assert reports['erasure36_output']['universal_unique_completion'] and reports['erasure36_output']['count']==1
    assert sum(c['rows_checked'] for c in reports['lattice_information_review']['cases'])==177
    info=reports['lattice_information']['cases']
    for item,d in zip(info[:3],[33,36,40]):
        units=[]
        for row in item['rows']:
            if not row['integer_carry_direction']:continue
            v=[F(*x) for x in row['inverse_row']];support=[i for i,x in enumerate(v) if x]
            assert len(support)==1 and abs(v[support[0]])==2013265921;units.extend(support)
        assert sorted(units)==list(range(d-13))
    orbit=reports['orbit_review']
    assert sum(c['queries_checked'] for c in orbit['small_cases'])==3344
    assert sum(c['direct_scalar_pair_checks'] for c in orbit['small_cases'])==246928
    assert orbit['large_seed_candidates_checked']==27904523 and orbit['large_count']==0
    assert orbit['large_scalar_difference']==345549834 and orbit['fixed_coordinates']==37
    assert reports['gap_review']['full_integer_lattice_checked'] and reports['gap_review']['actual_pair_count']==0
    assert reports['gap_review']['max_digit_difference']==4 and len(reports['gap_review']['corrupt_controls_rejected'])==2
    sources=list(HERE.glob('*.py'))
    for path in sources:compile(path.read_text(),str(path),'exec')
    for path in [*HERE.glob('*.md'),HERE.parent/'NEXT_ITERATION.md']:
        text=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if target.startswith(('http://','https://','#')):continue
            resolved=path.parent/target.split('#')[0].strip('<>')
            assert resolved.exists() or resolved==HERE/'manifest.json',(path,target)
    artifacts={str(path.relative_to(HERE)):{'sha256':digest(path),'bytes':path.stat().st_size}
               for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
               and path.name not in ('manifest.json','verification.log')}
    out={'status':'verified','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Exact compressed difference covers certify unique completion for any36 consecutive cyclic erasures on the saved finite input. Complete orbit counting exhibits an integral40-erasure lattice difference with no actual scalar pair. General40-erasure recovery and prize estimates remain open.',
         'previous_artifacts_preserved':preserved,'source_and_input_bindings':len(BINDINGS),'python_sources_compiled':len(sources),
         'markdown_local_links_passed':True,'cover_nodes_verified':5748,'cover_cuts_verified':298,'cover_narrowing_steps':24238,
         'preflight_queries':48,'small_ambiguity_pairs_preserved':550,'orbit_intersections_checked':3344,
         'large_orbit_candidates_checked':27904523,'integral_gap_checked':True,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'Round23 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; preserved {preserved}.',flush=True)


if __name__=='__main__':main()
