#!/usr/bin/env python3
"""Read back pass-36 arithmetic, provenance, links, and scope without a build."""
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT/path).read_text())


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def main():
    prior = read('results/parallel36_prior_state_2026_09_05.json')
    resume = read('results/parallel36_resume_state_2026_09_05.json')
    source = read('results/parallel36_source_scope_2026_09_05.json')
    data = read('results/parallel36_verification_2026_09_05.json')
    for mapping in [prior['prior_sha256'], source['local_inputs']]:
        for p, h in mapping.items():
            assert digest(p) == h, p
    for item in source['sources']:
        assert digest(item['path']) == item['sha256']
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    assert digest('experiments/parallel36_verify_2026_09_05.py') == data['source_sha256']
    assert digest('experiments/parallel36_shell_certificate.cpp') == data['cpp_sha256']
    assert digest('results/lean_performance_terminal_2026_09_05_1934.json') == resume['diagnostic_sha256']
    assert data['checks'] == 41574
    assert len(data['full_cpp_comparisons']) == 12
    assert len(data['dyadic_norm_cases']) == 28
    assert len(data['general_symmetric_cases']) == 6
    assert data['dirichlet_inverse_coefficients_checked'] == 4096

    # Check the explicit norm example using ordinary rational elimination,
    # independently of the verifier's fraction-free determinant code.
    p, N, r = 1153, 4, [1, 75, -140, -123]
    matrix = [[F(r[(i-j) % N]*(-1 if i < j else 1)) for j in range(N)] for i in range(N)]
    determinant = F(1)
    for j in range(N):
        pivot = matrix[j][j]
        assert pivot
        determinant *= pivot
        for i in range(j+1, N):
            scale = matrix[i][j]/pivot
            for k in range(j, N):
                matrix[i][k] -= scale*matrix[j][k]
    assert determinant == p**3
    assert sum(x*x for x in r) == 35*p
    assert F(8*(p+1), 24*p) == F(1154, 3459)

    certificates = []
    for entry in data['large_certificates']:
        L = entry['truncation']
        rec = read(f'results/parallel36_shell_6700417_64_{L}_2026_09_05.json')
        p, n, Q, den = rec['p'], rec['n'], rec['weight_scale'], rec['distance_denominator']
        assert p == 6700417 and n == 64 and den == 24*p*p
        assert rec['contains_two'] and rec['mass_checked']
        assert rec['cosets']*n == p-1
        assert rec['weighted_max_coset_index'] == 0
        Dnum, W = rec['max_deviation_numerator'], rec['weight_mass']
        P = F(data['pi_squared_rational_lower'])
        upper_bounds = []
        for key in ['weighted_max_abs_numerator', 'weighted_runner_up_abs_numerator']:
            gain = Dnum*W-int(rec[key])
            assert gain >= 0
            bound = F(36*p*p*Dnum, (p*p+1)*den)-F(2)*P*F(gain, den*Q)+F(n, p-1)
            upper_bounds.append(bound)
        assert upper_bounds == [F(entry['certificate_upper']), F(entry['other_cosets_upper'])]
        lo, hi = map(F, entry['eta_one_interval'])
        assert lo > 43 and upper_bounds[1] < lo <= hi < upper_bounds[0]
        assert entry['maximum_proved_to_occur_on_H']
        if L == 256:
            assert upper_bounds[0] < F(44151, 1000) and upper_bounds[1] < F(40150, 1000)
        else:
            assert L == 4096
            assert upper_bounds[0] < F(43832, 1000) and upper_bounds[1] < F(39838, 1000)
        certificates.append({'truncation': L, 'upper': str(upper_bounds[0]),
                             'other_cosets_upper': str(upper_bounds[1]),
                             'maximum_identified': True})

    old = read('results/period_polynomial_bounds.json')
    old_sources = {'experiments/period_polynomial_bounds.py': old['source_sha256'],
                   **old['dependency_hashes']}
    for p, h in old_sources.items():
        assert digest(p) == h
    assert old['conclusion']['M_equals_period_one']
    assert old['conclusion']['all_other_cosets_absolute_period_less_than'] == 42

    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel36-shell-inversion-2026-09-05.md',
             'research/parallel36-pass-summary-2026-09-05.md']
    links = []
    for p in central+notes:
        text = (ROOT/p).read_text()
        assert 'full' in text.lower() and 'unproved' in text.lower()
        for target in re.findall(r'\]\(([^)]+)\)', text):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            destination = ((ROOT/p).parent/target.split('#')[0]).resolve()
            assert destination.exists(), (p, target)
            links.append({'source': p, 'target': str(destination.relative_to(ROOT))})
    artifacts = notes+[
        'experiments/parallel36_shell_certificate.cpp',
        'experiments/parallel36_verify_2026_09_05.py',
        'experiments/parallel36_final_audit_2026_09_05.py',
        'results/parallel36_prior_state_2026_09_05.json',
        'results/parallel36_resume_state_2026_09_05.json',
        'results/parallel36_source_scope_2026_09_05.json',
        'results/parallel36_verification_2026_09_05.json',
        'results/parallel36_shell_6700417_64_256_2026_09_05.json',
        'results/parallel36_shell_6700417_64_4096_2026_09_05.json',
    ]+[s['path'] for s in source['sources']]
    audit = {
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'finite certificates, source scope and prior-state preservation passed; full goal unproved',
        'previous_goal_turn_classification': resume['previous_goal_turn_classification'],
        'current_turn_classification': 'progress: constant-preserving lattice reformulation and a rational finite certificate without high moments',
        'prior_sha256_preserved': prior['prior_sha256'],
        'artifact_sha256': {p: digest(p) for p in artifacts},
        'central_file_sha256': {p: digest(p) for p in central},
        'source_scope': source, 'old_polynomial_source_hashes_checked': old_sources,
        'large_certificates_recomputed': certificates,
        'alternate_rational_determinant': {'p': 1153, 'norm': str(determinant), 'square_sum': 35*1153},
        'local_links_checked': links, 'bounded_exact_assertions': data['checks'],
        'scope_limits': [
            'The uniform thin-annulus estimate remains unproved.',
            'Pass35 positive moment upper bound remains unproved.',
            'The elementary norm bound is too weak for square-root cancellation.',
            'The existing finite maximum enclosure is already tighter; no improvement to that enclosure.',
            'No improvement to the recorded uniform period exponent.',
            'No classical two-set Paley proof, official prize transfer, or prize proof.',
            'Same-author finite checks do not replace independent review or Lean formalization.'
        ],
        'worker_status': 'All three existing agents freshly inspected: terminal usage-limit errors; root continues.',
        'new_lean_jobs': 0, 'new_hosted_proof_jobs': 0, 'goal_status': 'active and unachieved'
    }
    (ROOT/'results/parallel36_pass_audit_2026_09_05.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps({'status': audit['status'], 'artifacts': len(artifacts),
                      'local_links': len(links), 'exact_assertions': data['checks']}, indent=2))


if __name__ == '__main__':
    main()
