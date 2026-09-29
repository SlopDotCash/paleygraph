#!/usr/bin/env python3
"""Pin the bounded analytic and computational advance without closing the goal."""
import ast
from collections import Counter
from datetime import datetime,timezone
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
OUT='results/parallel28_pass_audit_2026_09_05.json'
def J(p):return json.loads((ROOT/p).read_text())
def H(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def main():
    prior=J('results/parallel28_prior_state_2026_09_05.json')
    for p,h in prior['prior_sha256'].items():assert H(p)==h,p
    assert len(J('sources/manifest.json'))==36
    verification=J('results/parallel28_verification_2026_09_05.json')
    census=J('results/parallel28_complete_remainder_2026_09_05.json')
    for d in [verification,census]:
        assert d['status']=='passed'
        for p,h in d['input_sha256'].items():assert H(p)==h,('changed input',p)
    assert census['total_cases']==len(census['cases'])==28774
    assert census['prior_energy_comparisons']==28753
    assert len({(r['p'],r['n']) for r in census['cases']})==28774
    assert len(verification['direct_cases'])==75
    assert verification['checks']['independent_direct_six_category']==450
    assert verification['exponent_ledger']['best_saving']=='1/72'
    assert F(2849,2880)-F(71,72)==F(1,320)
    for n,expected_count,top in [(4,16,0),(8,95,0),(16,579,0),(32,3705,23040),(64,24379,184320)]:
        rows=[r for r in census['cases'] if r['n']==n]
        assert len(rows)==expected_count
        assert max(r['distinct_unbalanced_R6'] for r in rows)==top
        assert all(n**4//4<=r['p']<=n**4 for r in rows)
        assert all(r['distinct_unbalanced_R6']<=n**3 for r in rows)
        assert all(r['E3']==sum(r[k] for k in ['T6','J6','repeated_R6',
            'distinct_balanced_R6','distinct_unbalanced_R6']) for r in rows)
        assert all(r['distinct_unbalanced_R6']==720*n*r['D6_scaling_orbits'] for r in rows)
        expected={str(k):v for k,v in Counter(r['D6_scaling_orbits'] for r in rows).items()}
        assert expected==census['per_order'][str(n)]['D6_orbit_histogram']
    sources=J('results/parallel28_source_scope_2026_09_05.json')
    for e in sources['sources']:
        assert H(e['path'])==e['sha256']
        assert (ROOT/e['path']).stat().st_size==e['bytes']
        assert (ROOT/e['path']).read_bytes().startswith(b'%PDF')
    for p,h in sources['auxiliary_archive_sha256'].items():assert H(p)==h,p
    server=J('results/prove2me_projection_server_2026_09_05.json')
    assert server['verdict']['status']=='ACCEPTED'
    assert server['theorem_readback']['status']=='Proved'
    assert H('results/prove2me_projection_server_2026_09_05.json')==prior['prior_sha256']['results/prove2me_projection_server_2026_09_05.json']
    central=['README.md','research/frontier.md','research/source-audit.md',
             'research/checkpoint-2026-09-04.md','PROVE2ME.md','sources/manifest.json']
    artifacts={str(p.relative_to(ROOT)) for folder in ['experiments','research','results']
               for p in (ROOT/folder).glob('parallel28*') if p.is_file() and str(p.relative_to(ROOT))!=OUT}
    links=[]
    for p in sorted(artifacts|set(central)):
        value=(ROOT/p).read_text()
        assert not re.search(r'p2m_[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.',value),p
        if p.endswith('.py'):ast.parse(value)
        if not p.endswith('.md'):continue
        for line in value.splitlines():
            for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',line):
                if target.startswith(('http://','https://','#','mailto:','codex://')):continue
                if target=='u' and 'RawHyp' in line and p=='research/source-audit.md':continue
                dest=(ROOT/Path(p).parent/target.split('#')[0]).resolve()
                assert dest.exists() or dest==ROOT/OUT,(p,target)
                links.append({'source':p,'target':target})
    out={'status':'passed bounded-result audit; full goal unproved',
        'audited_at_utc':datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification':'progress',
        'current_turn_classification':'progress: derived stronger subgroup baseline and completed finite D6 census',
        'prior_inputs_preserved':prior['prior_sha256'],
        'artifact_sha256':{p:H(p) for p in sorted(artifacts)},
        'central_file_sha256':{p:H(p) for p in central},
        'source_checks':sources['sources'],'source_manifest_entries_unchanged':36,
        'local_links_checked':links,'actual_prime_order_pairs':28774,
        'independent_direct_fields':75,'uniform_subgroup_power':'71/72',
        'uniform_subgroup_log_power':'1/18','prior_project_power_improvement':'1/320',
        'analytic_proof_kind':'ordinary proof from cited theorems; no separate-author review or Lean verification',
        'hosted_projection_verdict':'ACCEPTED; prior terminal record unchanged',
        'worker_status':'Last reported terminal usage-limit errors; no restart or new worker completion claimed',
        'goal_status':'active and unachieved',
        'scope_limits':['No full two-set Paley bound or prize certificate.',
            'No growing-order cubic D6 bound.',
            'No literature novelty or current-best claim.',
            'Exact finite checks do not prove uniform analytic statements.',
            'No new hosted analytic proof verdict.']}
    (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','artifacts':len(artifacts),'census_cases':28774,
        'direct_fields':75,'local_links':len(links),'goal_status':out['goal_status']}))

if __name__=='__main__':main()
