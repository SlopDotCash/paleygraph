#!/usr/bin/env python3
"""Audit saved bounded evidence and provenance; never certify the full conjecture."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def main():
    prior = json.loads((ROOT/'results/parallel29_prior_state_2026_09_05.json').read_text())
    for path, expected in prior['prior_sha256'].items():
        assert digest(path) == expected, ('prior artifact changed', path)
    verified = json.loads((ROOT/'results/parallel29_verification_2026_09_05.json').read_text())
    for path, expected in verified['source_sha256'].items():
        assert digest(path) == expected, ('stale verification', path)
    assert sum(verified['check_counts'].values()) == 3426
    assert verified['incidence']['size_ranges'] == {'N_le_p_squared': 92, 'N_gt_p_squared': 5}
    assert verified['exponents_and_boundaries']['feedback_saving_supremum'] == '1/72'
    assert verified['approximate_fourier']['max_identity_relative_error'] < 1e-9
    source = json.loads((ROOT/'results/parallel29_source_scope_2026_09_05.json').read_text())
    for entry in source['sources']:
        assert digest(entry['path']) == entry['sha256'], entry['path']
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel29-centered-recurrence-2026-09-05.md',
             'research/parallel29-pass-summary-2026-09-05.md']
    links = []
    for path in central+notes:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT/path).read_text()):
            if '://' in target or target.startswith('#'):
                continue
            # Older research prose has RawHyp[...](u), which is math notation.
            # This audit checks local file references, all of which have paths
            # or extensions in the inspected documents.
            if '/' not in target and '.' not in target:
                continue
            dest = ((ROOT/path).parent/target.split('#')[0]).resolve()
            assert dest.exists(), (path, target)
            links.append({'source': path, 'target': str(dest.relative_to(ROOT))})
    artifacts = notes+[
        'experiments/parallel29_verify_2026_09_05.py',
        'experiments/parallel29_final_audit_2026_09_05.py',
        'results/parallel29_prior_state_2026_09_05.json',
        'results/parallel29_verification_2026_09_05.json',
        'results/parallel29_source_scope_2026_09_05.json',
    ]
    result = {
        'status': 'bounded evidence and provenance audit passed; full goal unproved',
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification': 'progress: live Lean diagnosis changed the next execution action',
        'current_turn_classification': 'progress: ordinary centered recurrence and closure of available moment rules',
        'prior_inputs_preserved': prior['prior_sha256'],
        'artifact_sha256': {p: digest(p) for p in artifacts},
        'central_file_sha256': {p: digest(p) for p in central},
        'reused_source_checks': source['sources'],
        'local_links_checked': links,
        'bounded_check_count': sum(verified['check_counts'].values()),
        'analytic_proof_status': 'Written ordinary proof with author audit; no separate-author review or Lean verification.',
        'quantitative_period_status': 'Existing n^(71/72)(log n)^(1/18) bound unchanged.',
        'worker_status': 'All three existing agents inspected this turn; terminal usage-limit errors remain.',
        'new_lean_processes_launched': 0,
        'new_hosted_proof_jobs': 0,
        'goal_status': 'active and unachieved',
        'remaining_obligations': ['Uniform square-root subgroup cancellation',
                                  'Arbitrary two-set classical Paley estimate',
                                  'Official scalar prize bounds and remaining certificate conditions',
                                  'Independent review of the new ordinary argument'],
    }
    dest = ROOT/'results/parallel29_pass_audit_2026_09_05.json'
    dest.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'artifacts': len(artifacts),
                      'checks': result['bounded_check_count'], 'local_links': len(links)}, indent=2))


if __name__ == '__main__':
    main()
