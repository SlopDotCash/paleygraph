#!/usr/bin/env python3
"""Provenance and finite-certificate audit for the uniform length argument."""
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT/path).read_text())


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def fraction(value):
    return Fraction(value['numerator'], value['denominator'])


def main():
    prior = read('results/parallel32_prior_state_2026_09_05.json')
    result = read('results/parallel32_verification_2026_09_05.json')
    sources = read('results/parallel32_source_scope_2026_09_05.json')
    for group in [prior['prior_sha256'], result['input_sha256'], sources['local_input_sha256']]:
        for path, expected in group.items():
            assert digest(path) == expected, ('changed input', path)
    assert digest('sources/manifest.json') == sources['main_manifest_sha256']
    assert sum(result['check_counts'].values()) == 2352
    rows = result['trinomial_rows']
    assert len(rows) == 252
    assert {(r['n'], r['b']) for r in rows} == {
        (n, b) for n in [4, 8, 16, 32, 64, 128] for b in range(n)}
    for row in rows:
        n = row['n']
        d = n//2
        epsilon = 1 if n % 3 == 1 else -1
        t = (n-epsilon)//3
        l1 = fraction(row['inverse_l1'])
        l2sq = fraction(row['inverse_l2_squared'])
        assert row['norm'] > 0
        assert l2sq <= 9*t
        assert l1*l1 <= d*l2sq <= 9*d*t <= Fraction(16*n*n, 9)
    assert len(result['conductor_lifts']) == 5
    fixed = result['fixed_field']
    assert fixed['p'] == 215535361 and fixed['n'] == 128
    assert fraction(fixed['inverse_l1']) == Fraction(2878824541, 215535361)
    assert fraction(fixed['six_term_bound']) == 6*fraction(fixed['inverse_l1'])
    assert fixed['integer_six_term_bound'] == fraction(fixed['six_term_bound']).__floor__() == 80
    assert fixed['uniform_six_term_bound'] == 1024
    known = read('results/parallel31_short_multiples_2026_09_05.json')
    assert max(r['quotient_l1'] for r in known['orbits']) == 44
    assert len(known['orbits']) == 119

    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel32-linear-triangle-length-2026-09-05.md',
             'research/parallel32-pass-summary-2026-09-05.md']
    links = []
    for path in central+notes:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT/path).read_text()):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            dest = ((ROOT/path).parent/target.split('#')[0]).resolve()
            assert dest.exists(), (path, target)
            links.append({'source': path, 'target': str(dest.relative_to(ROOT))})
    artifacts = notes+[
        'experiments/parallel32_verify_2026_09_05.py',
        'experiments/parallel32_final_audit_2026_09_05.py',
        'results/parallel32_prior_state_2026_09_05.json',
        'results/parallel32_verification_2026_09_05.json',
        'results/parallel32_source_scope_2026_09_05.json',
    ]
    audit = {
        'status': 'finite certificates and provenance passed; uniform derivation-length proof recorded; full goal unproved',
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification': prior['previous_turn_classification'],
        'current_turn_classification': 'progress: uniform one-orbit derivation length at most 8m for six-term targets',
        'prior_inputs_preserved': prior['prior_sha256'],
        'artifact_sha256': {p: digest(p) for p in artifacts},
        'central_file_sha256': {p: digest(p) for p in central},
        'source_scope': sources, 'local_links_checked': links,
        'bounded_exact_checks': 2352, 'trinomial_cases': 252,
        'conductor_lift_cases': 5, 'prior_orbits_checked': 119,
        'fixed_field': fixed,
        'proof_status': 'Uniform elementary Fourier proof, supported by exact finite calculations. Same-author implementation and review; no Lean or separate-author verification.',
        'worker_status': 'All three existing agents inspected this pass; terminal usage-limit errors remain.',
        'new_lean_jobs': 0, 'new_hosted_proof_jobs': 0,
        'goal_status': 'active and unachieved',
        'scope_limits': ['Target must belong to the integer ideal of the chosen triangle orbit.',
                         'No uniform bound on the number of six-term outputs.',
                         'No stronger period exponent, full Paley proof, or official prize certificate.',
                         'No literature novelty claim.'],
    }
    (ROOT/'results/parallel32_pass_audit_2026_09_05.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps({'status': audit['status'], 'artifacts': len(artifacts),
                      'checks': 2352, 'local_links': len(links)}, indent=2))


if __name__ == '__main__':
    main()
