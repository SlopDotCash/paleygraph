#!/usr/bin/env python3
"""Audit the finite counterexample and its provenance, not the full goal."""
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text())


def main():
    prior = read('results/parallel30_prior_state_2026_09_05.json')
    for p, expected in prior['prior_sha256'].items():
        assert digest(p) == expected, ('prior changed', p)
    result = read('results/parallel30_verification_2026_09_05.json')
    source = read('results/parallel30_source_scope_2026_09_05.json')
    for group in [result['source_sha256'], source['local_input_sha256']]:
        for p, expected in group.items():
            assert digest(p) == expected, ('input changed', p)
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    assert sum(result['check_counts'].values()) == 13201
    short = result['short_certificate']
    row = result['independent_direct_counts'][-1]
    assert short['p'] == row['p'] == 215535361 and short['n'] == row['n'] == 128
    assert short['analytic_D6_lower_bound'] == row['triangle_D6'] == 4976640
    assert row['primitive_D6'] == 5990400
    assert row['distinct_unbalanced_R6'] == row['triangle_D6']+row['primitive_D6'] == 10967040
    assert row['primitive_D6'] > 128**3 and short['analytic_D6_lower_bound'] > 128**3
    assert Fraction(row['distinct_unbalanced_R6'], 128**3) == Fraction(5355, 1024)
    assert short['norm_determinant'] == short['p']
    old_direct = read('results/parallel30_direct_witness_2026_09_05.json')
    sparse = read('results/parallel30_sparse_witness_2026_09_05.json')
    for k in ['E3', 'T6', 'J6', 'repeated_R6', 'distinct_balanced_R6', 'distinct_unbalanced_R6']:
        assert old_direct[k] == sparse[k] == row[k], k
    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel30-triangle-remainder-2026-09-05.md',
             'research/parallel30-pass-summary-2026-09-05.md']
    links = []
    for path in central+notes:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT/path).read_text()):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            dest = ((ROOT/path).parent/target.split('#')[0]).resolve()
            assert dest.exists(), (path, target)
            links.append({'source': path, 'target': str(dest.relative_to(ROOT))})
    artifacts = notes+[
        'experiments/parallel30_triangle_search_2026_09_05.py',
        'experiments/parallel30_verify_2026_09_05.py',
        'experiments/parallel30_direct_six.cpp',
        'experiments/parallel30_final_audit_2026_09_05.py',
        'results/parallel30_prior_state_2026_09_05.json',
        'results/parallel30_triangle_search_2026_09_05.json',
        'results/parallel30_sparse_witness_2026_09_05.json',
        'results/parallel30_direct_witness_2026_09_05.json',
        'results/parallel30_verification_2026_09_05.json',
        'results/parallel30_source_scope_2026_09_05.json',
    ]
    audit = {
        'status': 'exact finite counterexample and provenance audit passed; full goal unproved',
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification': 'progress: centered recurrence and closure of moment-bound rules',
        'current_turn_classification': 'progress: arithmetic counterexample and structural triangle bounds',
        'prior_inputs_preserved': prior['prior_sha256'],
        'artifact_sha256': {p: digest(p) for p in artifacts},
        'central_file_sha256': {p: digest(p) for p in central},
        'source_checks': source,
        'local_links_checked': links,
        'bounded_exact_checks': 13201,
        'counterexample': row,
        'short_lower_certificate': short,
        'proof_status': 'Ordinary proof and exact finite arithmetic; independent implementations by same author. No separate-author or Lean verification.',
        'worker_status': 'All three existing agents inspected; terminal usage-limit errors remain.',
        'goal_status': 'active and unachieved',
        'new_lean_jobs': 0,
        'new_hosted_proof_jobs': 0,
        'scope_limits': ['Refutes only the constant-one extrapolation D6<=n^3.',
                         'Does not refute D6=O(n^3), the earlier finite census, or Paley.',
                         'No stronger cancellation exponent or official prize certificate.',
                         'Norm search is exploratory; completeness is not used.'],
    }
    (ROOT/'results/parallel30_pass_audit_2026_09_05.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps({'status': audit['status'], 'artifacts': len(artifacts),
                      'checks': 13201, 'local_links': len(links)}, indent=2))


if __name__ == '__main__':
    main()
