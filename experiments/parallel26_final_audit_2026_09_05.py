#!/usr/bin/env python3
"""Pin the formalization and source correction without rewriting earlier audits."""
import ast
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/parallel26_pass_audit_2026_09_05.json'
CENTRAL = ['README.md', 'PROVE2ME.md', 'research/frontier.md',
           'research/source-audit.md', 'research/checkpoint-2026-09-04.md',
           'sources/manifest.json']

def H(p): return sha256((ROOT / p).read_bytes()).hexdigest()
def J(p): return json.loads((ROOT / p).read_text())

def main():
    prior = J('results/parallel26_prior_state_2026_09_05.json')
    for p, expected in prior['prior_sha256'].items():
        assert H(p) == expected, ('changed prior input', p)
    old_audits = ['results/parallel25_pass_audit_2026_09_05.json',
                  'results/prove2me_setup_audit_2026_09_05.json']
    allowed_updates = {'PROVE2ME.md': 'Updated current formalization and attribution.',
        'results/prove2me_projection_server_2026_09_05.json':
            'Deliberately refreshed read-only status of the same private job.'}
    preserved, updates = {}, {}
    for audit in old_audits:
        for p, expected in J(audit)['artifact_sha256'].items():
            if H(p) != expected:
                assert p in allowed_updates, ('historical artifact changed', p)
                updates[p] = {'old_sha256': expected, 'new_sha256': H(p),
                              'reason': allowed_updates[p]}
            else:
                preserved[p] = expected
    assert len(J(old_audits[0])['artifact_sha256']) == 20
    assert all(p in preserved for p in J(old_audits[0])['artifact_sha256'])

    # Reconstruct the exact earlier registry from its original ledger.
    pass25prior = J('results/parallel25_prior_state_2026_09_05.json')
    previous_manifest = json.loads(pass25prior['prior_manifest_text'])
    assert sha256(pass25prior['prior_manifest_text'].encode()).hexdigest() == \
        pass25prior['prior_sha256']['sources/manifest.json']
    previous_manifest += [J('results/parallel25_prize_projection_source_2026_09_05.json'),
                          J('results/parallel25_subgroup_source_2026_09_05.json')]
    manifest = J('sources/manifest.json')
    assert len(previous_manifest) == 34 and manifest[:34] == previous_manifest
    new_sources = J('results/parallel26_primary_sources_2026_09_05.json')
    for entry, name in zip(new_sources,
            ['parallel26-arklib-interleaved-code', 'parallel26-arklib-linear-code']):
        entry['name'] = name
    assert manifest[34:] == new_sources and len(manifest) == 36
    inspected_sources = [e['path'] for e in new_sources] + [
        'sources/parallel23-prize-dependencies/ProximityGap_Errors.lean',
        'sources/parallel23-prize-dependencies/ProximityGap_ProximityGenerators.lean',
        'sources/parallel23-prize-dependencies/ListDecodability.lean']
    entries = {e.get('path'): e for e in manifest if e.get('path')}
    source_checks = []
    for p in inspected_sources:
        entry = entries[p]
        assert H(p) == entry['sha256']
        assert (ROOT / p).stat().st_size == entry['bytes']
        source_checks.append({'path': p, 'sha256': H(p), 'bytes': entry['bytes']})
    existing = (ROOT / inspected_sources[2]).read_text()
    assert 'theorem mcaError_interleaved_eq' in existing
    assert 'theorem mcaError_interleaved_le' in existing

    lean = J('results/parallel26_mca_projection_2026_09_05.json')
    proof_path = 'prove2me/Check_paley_mca_projection.lean'
    assert lean['exit_code'] == 0 and lean['status'] == 'passed'
    assert H(proof_path) == lean['proof_sha256']
    assert H('results/parallel26_mca_lean_2026_09_05.log') == lean['log_sha256']
    assert H('experiments/parallel26_lean_check_2026_09_05.py') == lean['checker_sha256']
    assert lean['mathlib_rev'] == 'c5ea00351c28e24afc9f0f84379aa41082b1188f'
    assert lean['toolchain'] == 'leanprover/lean4:v4.30.0'
    assert len(lean['axioms_by_theorem']) == 7
    assert all(set(v) <= {'propext', 'Classical.choice', 'Quot.sound'}
               for v in lean['axioms_by_theorem'].values())
    assert not re.search(r'\b(sorry|admit|axiom|unsafe)\b|import Theorems\.',
                         (ROOT / proof_path).read_text())
    assert (ROOT / proof_path).read_bytes() == \
        (Path.home() / 'prove2me_workspace/Solutions' / Path(proof_path).name).read_bytes()

    server = J('results/prove2me_projection_server_2026_09_05.json')
    assert server['job_id'] == '843b8456-f9ad-4158-809d-8378028cc958'
    assert server['private'] and server['publish_job']['visibility'] == 'private'
    assert server['proof_sha256'] == H('prove2me/Sol_paley_projection_survival_count.lean')
    assert server['payload_sha256'] == H('prove2me/projection-survival-problem.json')
    assert server['publish_job']['status'] == 'PENDING'
    assert not server.get('theorem_id') and not server.get('submission_id') and not server.get('verdict')

    artifacts = {proof_path}
    for folder in ('research', 'experiments', 'results'):
        artifacts.update(str(p.relative_to(ROOT)) for p in (ROOT / folder).glob('parallel26*')
                         if p.is_file() and p != OUT)
    syntax, links, exclusions = [], [], []
    for p in sorted(artifacts | set(CENTRAL)):
        value = (ROOT / p).read_text()
        assert not re.search(r'p2m_[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.', value), ('credential found', p)
        if p.endswith('.py'):
            ast.parse(value); syntax.append(p)
        if p.endswith('.md'):
            for line in value.splitlines():
                for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', line):
                    if target.startswith(('https://', 'http://', '#', 'mailto:', 'codex://')):
                        continue
                    if target == 'u' and 'RawHyp' in line and p == 'research/source-audit.md':
                        exclusions.append({'source': p, 'target': target,
                                           'reason': 'Historical mathematical function notation.'})
                        continue
                    resolved = (ROOT / Path(p).parent / target.split('#')[0]).resolve()
                    assert resolved.exists() or resolved == OUT, (p, target)
                    links.append({'source': p, 'target': target})
    record = {
        'status': 'Passed formalization and source-attribution audit; full goal unproved.',
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification': 'progress',
        'current_turn_classification': 'progress in verification and source accuracy',
        'prior_input_checks': prior['prior_sha256'],
        'prior_artifacts_preserved': preserved, 'documented_mutable_updates': updates,
        'source_manifest_previous_entries_preserved': 34,
        'source_manifest_total_entries': len(manifest), 'source_hash_checks': source_checks,
        'local_lean_principal_statements': lean['axioms_by_theorem'],
        'proof_sha256': lean['proof_sha256'],
        'novelty_correction': 'MCA interleaving invariance already proved in pinned ArkLib; '
                              'independent local verification, not a new mathematical bound.',
        'hosted_incidence_job_status': server['publish_job']['status'],
        'hosted_incidence_job_observed_at_utc': server['updated_at_utc'],
        'hosted_incidence_theorem_id': None, 'hosted_incidence_submission_id': None,
        'hosted_incidence_proof_verdict': None,
        'worker_status': 'All three terminal usage-limit errors; root continued locally.',
        'syntax_checks': syntax, 'local_link_checks': links,
        'mathematical_notation_excluded': exclusions,
        'artifact_sha256': {p: H(p) for p in sorted(artifacts)},
        'central_file_sha256': {p: H(p) for p in CENTRAL},
        'goal_status': 'active and unachieved',
        'scope_limits': [
            'No new uniform scalar MCA/list bound or cancellation exponent.',
            'No full ArkLib dependency build or direct formal identification theorem.',
            'List transfer remains author-reviewed and finitely checked, without Lean verification.',
            'No new coefficient-pigeonhole threshold; already present in earlier project notes.',
            'No Paley-to-prize equivalence or full proof.',
        ],
    }
    OUT.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'prior_artifacts_preserved': len(preserved),
        'documented_updates': list(updates), 'source_files': len(source_checks),
        'manifest_entries': len(manifest), 'local_lean_statements': 7,
        'local_links': len(links), 'goal_status': record['goal_status']}, indent=2))

if __name__ == '__main__':
    main()
