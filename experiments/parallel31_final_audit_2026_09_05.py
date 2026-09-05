#!/usr/bin/env python3
"""Audit the fixed-field certificates and provenance, not the full goal."""
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT / path).read_text())


def main():
    prior = read('results/parallel31_prior_state_2026_09_05.json')
    classification_path = 'results/parallel31_short_multiples_2026_09_05.json'
    classification = read(classification_path)
    certificate = read('results/parallel31_certificate_check_2026_09_05.json')
    source = read('results/parallel31_source_scope_2026_09_05.json')
    for group in [prior['prior_sha256'], classification['input_sha256'],
                  certificate['input_sha256'], source['local_input_sha256']]:
        for path, expected in group.items():
            assert digest(path) == expected, ('input changed', path)
    assert certificate['certificate_sha256'] == digest(classification_path)
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    for item in source['reused_primary_sources']:
        assert digest(item['path']) == item['sha256'], item['path']
    assert sum(classification['check_counts'].values()) == 1148
    assert sum(certificate['check_counts'].values()) == 1026
    assert classification['p'] == classification['norm'] == 215535361
    assert classification['n'] == 128
    assert classification['generator'] == 25525303
    assert classification['normalized_six_sets'] == 714
    orbits = classification['orbits']
    assert len(orbits) == 119
    profile = Counter((r['kind'], r['quotient_l1']) for r in orbits)
    expected = {('triangular', 2): 54, ('primitive', 4): 20,
                ('primitive', 38): 1, ('primitive', 40): 9,
                ('primitive', 42): 17, ('primitive', 44): 18}
    assert profile == expected
    assert {(r['kind'], r['minimum_length']): r['orbits']
            for r in certificate['profile']} == expected
    long = [r for r in orbits if r['kind'] == 'primitive' and r['quotient_l1'] > 4]
    assert len(long) == 45
    assert {r['cancellation_graph']['cycle_rank'] for r in long} == {17, 18, 19, 20}
    for row in orbits:
        length = row['quotient_l1']
        assert sum(map(abs, row['quotient'])) == length
        graph = row['cancellation_graph']
        assert len(graph['edges']) == (3 * length - 6) // 2
        if row['kind'] == 'primitive':
            assert graph['components'] == 1
            assert graph['cycle_rank'] == length // 2 - 2
    previous = read('results/parallel30_verification_2026_09_05.json')
    count = previous['independent_direct_counts'][-1]
    assert len(orbits) * 720 * 128 == count['distinct_unbalanced_R6'] == 10967040
    assert profile[('triangular', 2)] * 720 * 128 == count['triangle_D6']

    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel31-short-multiples-2026-09-05.md',
             'research/parallel31-pass-summary-2026-09-05.md']
    links = []
    for path in central + notes:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT / path).read_text()):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            dest = ((ROOT / path).parent / target.split('#')[0]).resolve()
            assert dest.exists(), (path, target)
            links.append({'source': path, 'target': str(dest.relative_to(ROOT))})
    artifacts = notes + [
        'experiments/parallel31_short_multiples_2026_09_05.py',
        'experiments/parallel31_check_certificate_2026_09_05.py',
        'experiments/parallel31_final_audit_2026_09_05.py',
        'results/parallel31_prior_state_2026_09_05.json',
        classification_path,
        'results/parallel31_certificate_check_2026_09_05.json',
        'results/parallel31_source_scope_2026_09_05.json',
    ]
    audit = {
        'status': 'fixed-field minimum-length certificate and provenance audit passed; full goal unproved',
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification': 'progress: arithmetic counterexample to constant-one D6 bound',
        'current_turn_classification': 'progress: exact minimum triangle lengths for all 119 fixed-field orbits',
        'prior_inputs_preserved': prior['prior_sha256'],
        'artifact_sha256': {p: digest(p) for p in artifacts},
        'central_file_sha256': {p: digest(p) for p in central},
        'source_checks': source,
        'local_links_checked': links,
        'bounded_exact_checks': 2174,
        'orbit_profile': certificate['profile'],
        'long_primitive_orbits': 45,
        'proof_status': 'Ordinary minimum-length proof and exact finite certificates; both implementations by same author. No separate-author or Lean verification.',
        'goal_status': 'active and unachieved',
        'new_lean_jobs': 0,
        'new_hosted_proof_jobs': 0,
        'scope_limits': [
            'One fixed prime and subgroup; no uniform remainder upper bound.',
            'No stronger cancellation exponent or official prize certificate.',
            'The cancellation graphs are not identified with spectral necklace graphs.',
            'No literature novelty or exhaustive higher-energy obstruction claimed.',
        ],
    }
    output = ROOT / 'results/parallel31_pass_audit_2026_09_05.json'
    output.write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({'status': audit['status'], 'artifacts': len(artifacts),
                      'checks': 2174, 'orbits': len(orbits),
                      'local_links': len(links)}, indent=2))


if __name__ == '__main__':
    main()
