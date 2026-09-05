#!/usr/bin/env python3
"""Audit the exact-count advance separately from the unproved uniform objective."""
import ast
from datetime import datetime,timezone
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/parallel27_pass_audit_2026_09_05.json'
def H(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def J(p): return json.loads((ROOT/p).read_text())

def main():
    prior=J('results/parallel27_prior_state_2026_09_05.json')
    for p,h in prior['prior_sha256'].items(): assert H(p)==h,('prior input changed',p)
    previous=J('results/parallel26_pass_audit_2026_09_05.json')
    for p,h in previous['artifact_sha256'].items(): assert H(p)==h,('pass26 artifact changed',p)
    assert (ROOT/'sources/manifest.json').read_text()==prior['prior_manifest_text']
    assert len(J('sources/manifest.json'))==36

    inputs=[]
    result_paths=['results/parallel27_six_distinct_2026_09_05.json',
                  'results/parallel27_independent_check_2026_09_05.json']
    for p in result_paths:
        data=J(p)
        assert data['status']=='passed'
        for file,h in data['input_sha256'].items():
            assert H(file)==h,('verification input changed',file)
            inputs.append({'path':file,'sha256':h,'result':p})
    sparse, independent=map(J,result_paths)
    rows=sparse['cases']
    assert len(rows)==17 and len({(r['p'],r['n']) for r in rows})==17
    assert sparse['checks']['independent_six_multiset_enumeration']==42
    assert sparse['checks']['all_distinct_opposite_free_unique_split']==36800
    assert len(independent['direct_cases'])==17
    for r,d in zip(rows,independent['direct_cases']):
        for k in ['p','n','E3','T6','J6','repeated_R6','distinct_balanced_R6','distinct_unbalanced_R6']:
            assert r[k]==d[k],('direct count mismatch',r['p'],k)
        assert r['E3']==sum(r[k] for k in
            ['T6','J6','repeated_R6','distinct_balanced_R6','distinct_unbalanced_R6'])
        for k in ['J6','repeated_R6','distinct_balanced_R6','distinct_unbalanced_R6']:
            assert r[k]>=0
        assert r['distinct_unbalanced_R6']==720*r['n']*r['D6_scaling_orbits']
    larger=[r for r in rows if r['n']>=512]
    assert len(larger)==6
    assert [r['D6_scaling_orbits'] for r in larger]==[5,2,0,8,3,2]
    for r in larger:
        assert r['n']**4//4<=r['p']<=r['n']**4
        assert r['E3']==r['T6']+r['distinct_unbalanced_R6']<15*r['n']**3
    assert any(r['E3']>15*r['n']**3 for r in rows)

    assert independent['mathlib_rev']=='c5ea00351c28e24afc9f0f84379aa41082b1188f'
    axioms=independent['axioms_by_theorem']
    assert set(axioms)=={'PaleyBalancedSplit.'+n for n in
        ['two_equations_force_collision','two_balanced_splits_force_collision','no_two_balanced_splits']}
    assert all(set(v)<={'propext','Classical.choice','Quot.sound'} for v in axioms.values())
    proof='prove2me/Check_paley_balanced_split_collision.lean'
    assert not re.search(r'\b(sorry|admit|axiom|unsafe)\b|import Theorems\.',(ROOT/proof).read_text())
    sources=J('results/parallel27_source_scope_2026_09_05.json')['sources']
    for entry in sources:
        assert H(entry['path'])==entry['sha256']
        assert (ROOT/entry['path']).stat().st_size==entry['bytes']

    central=['README.md','PROVE2ME.md','research/frontier.md','research/source-audit.md',
             'research/checkpoint-2026-09-04.md','sources/manifest.json']
    artifacts={proof}
    for folder in ['research','experiments','results']:
        artifacts.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).glob('parallel27*')
                         if p.is_file() and p!=OUT)
    syntax,links=[],[]
    for p in sorted(artifacts|set(central)):
        value=(ROOT/p).read_text()
        assert not re.search(r'p2m_[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.',value),('credential found',p)
        if p.endswith('.py'): ast.parse(value);syntax.append(p)
        if not p.endswith('.md'): continue
        for line in value.splitlines():
            for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',line):
                if target.startswith(('https://','http://','#','mailto:','codex://')): continue
                if target=='u' and 'RawHyp' in line and p=='research/source-audit.md': continue
                resolved=(ROOT/Path(p).parent/target.split('#')[0]).resolve()
                assert resolved.exists() or resolved==OUT,(p,target)
                links.append({'source':p,'target':target})
    server=J('results/prove2me_projection_server_2026_09_05.json')
    assert server['job_id']=='843b8456-f9ad-4158-809d-8378028cc958' and server['private']
    assert server['publish_job']['status']=='PUBLISHED'
    assert server['theorem_id']=='d025929b-75ed-4c1c-b661-8e491017b153'
    assert server['submission_id']=='27ff60e7-7d6d-48da-8c34-7b123c6c0a3f'
    assert server['verdict']['theorem_id']==server['theorem_id']
    assert server['proof_sha256']==H('prove2me/Sol_paley_projection_survival_count.lean')
    assert server['payload_sha256']==H('prove2me/projection-survival-problem.json')
    verdict=server['verdict']['status']
    assert verdict in ['PENDING','ACCEPTED']
    if verdict=='ACCEPTED':
        assert server['theorem_readback']['status']=='Proved'
        assert server['theorem_readback']['mathlib_rev']==independent['mathlib_rev']
    report={'status':'Passed bounded-result audit; uniform estimate and full goal remain unproved.',
        'audited_at_utc':datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification':'progress',
        'current_turn_classification':'progress: exact count, Lean algebra, and larger independent finite checks',
        'prior_input_hashes_preserved':prior['prior_sha256'],
        'pass26_artifact_hashes_preserved':previous['artifact_sha256'],
        'source_manifest_entries_unchanged':36,'verification_input_checks':inputs,
        'source_checks':sources,'actual_cases':len(rows),'new_larger_cases':len(larger),
        'direct_comparison_cases':len(independent['direct_cases']),'lean_axioms':axioms,
        'syntax_checks':syntax,'local_link_checks':links,
        'artifact_sha256':{p:H(p) for p in sorted(artifacts)},
        'central_file_sha256':{p:H(p) for p in central},
        'hosted_job_id':server['job_id'],'hosted_job_status':'PUBLISHED',
        'hosted_theorem_id':server['theorem_id'],'hosted_submission_id':server['submission_id'],
        'hosted_proof_verdict':verdict,
        'hosted_observation_utc':server['updated_at_utc'],
        'hosted_observation_sha256':H('results/prove2me_projection_server_2026_09_05.json'),
        'worker_state':'All three last reported terminal usage-limit errors; root completed this pass.',
        'goal_status':'active and unachieved',
        'scope_limits':[
            'The computational exponent bounds evaluation cost, not D6 or E3.',
            'Larger cases are selected finite samples, not all quartic primes or an asymptotic theorem.',
            'The algebraic core is Lean-checked; the full combinatorial count is an ordinary proof.',
            'Both implementations have the same author; no separate-author or human review claimed.',
            'No new cancellation exponent, scalar prize bound, Paley equivalence, or full proof.',
        ]}
    OUT.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'passed','prior_pass_artifacts':len(previous['artifact_sha256']),
        'cases':17,'independent_direct_cases':17,'lean_statements':3,
        'artifacts':len(artifacts),'source_entries_unchanged':36,
        'local_links':len(links),'goal_status':report['goal_status']},indent=2))

if __name__=='__main__': main()
