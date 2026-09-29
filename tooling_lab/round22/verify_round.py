#!/usr/bin/env python3
"""Bind the exact conditioning, metric and joint certificates and prior rounds."""
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
BINDINGS=[]


def digest(path): return sha256(path.read_bytes()).hexdigest()


def bind(path,wanted):
    assert path.is_file() and digest(path)==wanted,str(path)
    BINDINGS.append(str(path.resolve()))


def scan(owner,value):
    if isinstance(value,list):
        for v in value: scan(owner,v)
    elif isinstance(value,dict):
        for key,v in value.items():
            if key.endswith('_sha256') and isinstance(v,dict):
                for name,wanted in v.items(): bind(owner.parent/name,wanted)
            else: scan(owner,v)


def main():
    preserved={}
    for name,filename in [('round21','manifest.json'),('round20','manifest.json'),('round19','manifest.json'),
                          ('round18','manifest.json'),('round17','checkpoint.json'),('round16','manifest.json')]:
        path=HERE.parent/name/filename; data=json.loads(path.read_text())
        for rel,item in data['artifacts'].items(): bind(path.parent/rel,item['sha256'])
        preserved[name]=len(data['artifacts'])
    reports={}
    for path in HERE.glob('*.json'):
        if path.name=='manifest.json': continue
        data=json.loads(path.read_text()); scan(path,data); reports[path.stem]=data
    for n in ['conditioned_summary','small_summary','metric_summary','coupled_summary','coupled_small','plot_metadata']:
        assert reports[n]['status']=='produced'
    for n in ['conditioned_review','metric_review','coupled_review','coupled_small_review','query_review','joint_query_review']:
        assert reports[n]['status']=='passed'
    conditioned=reports['conditioned_review']['cases']
    expected={'cuts_verified':138,'optimality_certificates':137,'queries_replayed':1921,
              'full_small_codebook_oracles':1825,'queries_with_nonzero_center_shift':473,'scalar_candidates_checked':11163}
    for key,v in expected.items(): assert sum(c[key] for c in conditioned)==v
    metric=reports['metric_review']['cases']
    for key,v in [('cuts_verified',96),('optimality_certificates',95),('queries_replayed',96),('scalar_candidates_checked',129)]:
        assert sum(c[key] for c in metric)==v
    c=reports['metric_summary']['cases'][0]
    assert c['universal_cap']==512 and c['maximum_observed_box']==64 and c['complete_counts']=={'1':16,'0':16}
    assert all(c['statuses']=={'budget_exceeded':32} for c in reports['metric_summary']['cases'][1:])
    assert all(c['statuses']=={'budget_exceeded':32} for c in reports['conditioned_summary']['cases'][1:])
    joint=reports['coupled_review']
    for key,v in [('initial_nonzero_differences',19682),('sign_representatives',9841),('pair_cuts_checked',156),
                  ('pair_optimality_certificates',156),('remaining_after_pair_cuts',86),('target_separators_checked',86),
                  ('target_optimality_certificates',86),('continuous_feasible_targets',0),('unresolved_targets',0)]:
        assert joint[key]==v
    assert joint['universal_unique_completion']
    small=reports['coupled_small_review']['cases']
    for key,v in [('all_actual_shared_known_pairs_checked',5206),('targets_checked',78),
                  ('actual_ambiguity_targets_preserved',36),('separators_checked',28),('continuous_witnesses_checked',50)]:
        assert sum(c[key] for c in small)==v
    query=reports['query_review']; assert len(query['cyclic_controls'])==64 and len(query['controls'])==8
    assert sum(c['status']=='rejected' for c in query['controls'])==5
    assert query['omitted_center_regression']['anchor_scalar']==20
    assert len(reports['joint_query_review']['controls'])==2
    assert reports['joint_output']['universal_unique_completion'] and reports['joint_output']['count']==1
    assert reports['joint_output']['completions'][0]['scalar']==1234567
    corrupt=(len(reports['conditioned_review']['corrupt_certificates_rejected'])+
             len(reports['metric_review']['corrupt_transformations_rejected'])+
             len(joint['corrupt_controls_rejected'])+1)
    assert corrupt==13
    sources=list(HERE.glob('*.py'))
    for path in sources: compile(path.read_text(),str(path),'exec')
    for path in [*HERE.glob('*.md'),HERE.parent/'NEXT_ITERATION.md']:
        text=re.sub(r'```.*?```|`[^`\n]*`','',path.read_text(),flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if target.startswith(('http://','https://','#')): continue
            resolved=path.parent/target.split('#')[0].strip('<>')
            assert resolved.exists() or resolved==HERE/'manifest.json',(path,target)
    artifacts={str(path.relative_to(HERE)):{'sha256':digest(path),'bytes':path.stat().st_size}
               for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
               and path.name not in ('manifest.json','verification.log')}
    out={'status':'verified','created_utc':datetime.now(timezone.utc).isoformat(),
         'scope':'Exact conditioned intervals, changed metric, complete coupled-difference exclusion and unique completion for any32 consecutive cyclic erasures on the saved finite input. Arbitrary erasures and prize estimates remain open.',
         'previous_artifacts_preserved':preserved,'source_and_input_bindings':len(BINDINGS),
         'python_sources_compiled':len(sources),'markdown_local_links_passed':True,
         'query_replays':2017,'small_codebook_oracles':1825,'coupled_difference_sign_representatives':9841,
         'small_actual_ambiguity_pairs':5206,'cyclic_controls':64,'corrupt_math_and_coverage_controls':corrupt,
         'command_line_rejection_controls':5,'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(f'Round22 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; preserved {preserved}.',flush=True)


if __name__=='__main__': main()
