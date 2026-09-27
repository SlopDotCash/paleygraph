#!/usr/bin/env python3
"""Audit centered collision certificates and provenance, not the full goal."""
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from math import isqrt, prod
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads((ROOT/path).read_text())


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def q(value):
    return Fraction(value['numerator'], value['denominator'])


def main():
    prior = read('results/parallel33_prior_state_2026_09_05.json')
    result = read('results/parallel33_verification_2026_09_05.json')
    source = read('results/parallel33_source_scope_2026_09_05.json')
    for group in [prior['prior_sha256'], result['input_sha256'], source['local_input_sha256']]:
        for path, expected in group.items():
            assert digest(path) == expected, ('changed input', path)
    for item in source['sources']:
        assert digest(item['path']) == item['sha256'], item['path']
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    assert sum(result['check_counts'].values()) == 740
    assert len(result['cases']) == 33
    partition_cases = 0
    for row in result['cases']:
        p, n, r = row['p'], row['n'], row['r']
        assert r % 2 == 0 and 4 <= r < p
        alpha = Fraction(n*(p-n), p-1)
        assert q(row['alpha']) == alpha > 0
        full = Fraction(row['energy'])-Fraction(n**r, p)
        fall = prod(n-j for j in range(r)) if r <= n else 0
        distinct = Fraction(row['distinct_zero_count'])-Fraction(fall, p)
        assert full == q(row['centered_energy']) > 0
        assert distinct == q(row['centered_distinct_count'])
        assert sum(z['coefficient']*z['zero_count'] for z in row['partition_rows']) == row['distinct_zero_count']
        assert sum(z['coefficient']*n**len(z['parts']) for z in row['partition_rows']) == fall
        lower = isqrt(alpha.numerator//alpha.denominator)
        epsilon = prod(Fraction(lower+j, lower) for j in range(1, r))-1
        assert epsilon == q(row['rational_epsilon_upper'])
        assert abs(full-distinct) <= epsilon*full
        for item in row['partition_rows']:
            k = len(item['parts'])
            assert sum(item['parts']) == r and all(0 < m < p for m in item['parts'])
            value = Fraction(item['zero_count'])-Fraction(n**k, p)
            assert value == q(item['centered'])
            assert value*value*alpha**(r-k) <= full*full
            partition_cases += 1
    assert partition_cases == 418
    bad = result['excluded_characteristic_example']
    assert bad['p'] == 3 and bad['parts'] == [3, 3]
    assert q(bad['mixed_centered']) == Fraction(8, 3)
    assert q(bad['centered_even_moment']) == Fraction(2, 3)
    for row in result['logarithmic_depth_parameters']:
        n, s = row['n'], row['s']
        assert n == 2**row['log2_n'] and s == 3*row['log2_n']
        lower = isqrt(n-1)
        value = prod(Fraction(lower+j, lower) for j in range(1, 2*s))-1
        assert value == q(row['epsilon_upper'])
        assert (value < 1) == row['two_sided_bound_nontrivial']
    assert Fraction(69, 20)-Fraction(129, 40) == Fraction(9, 40)

    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel33-centered-distinct-moments-2026-09-05.md',
             'research/parallel33-source-check-2026-09-05.md',
             'research/parallel33-pass-summary-2026-09-05.md']
    links = []
    for path in central+notes:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT/path).read_text()):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            dest = ((ROOT/path).parent/target.split('#')[0]).resolve()
            assert dest.exists(), (path, target)
            links.append({'source': path, 'target': str(dest.relative_to(ROOT))})
    artifacts = notes+[
        'experiments/parallel33_verify_2026_09_05.py',
        'experiments/parallel33_final_audit_2026_09_05.py',
        'results/parallel33_prior_state_2026_09_05.json',
        'results/parallel33_verification_2026_09_05.json',
        'results/parallel33_source_scope_2026_09_05.json',
    ]+[item['path'] for item in source['sources'] if '/parallel33-' in item['path']]
    audit = {
        'status': 'centered collision identities, finite bounds, and provenance passed; full goal unproved',
        'audited_at_utc': datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification': prior['previous_turn_classification'],
        'current_turn_classification': 'progress: centered distinct-coordinate reduction at growing depth and repeated-six-word exponent129/40',
        'prior_inputs_preserved': prior['prior_sha256'],
        'artifact_sha256': {p: digest(p) for p in artifacts},
        'central_file_sha256': {p: digest(p) for p in central},
        'source_scope': source, 'local_links_checked': links,
        'bounded_exact_checks': 740, 'set_order_cases': 33,
        'centered_partition_cases': 418,
        'proof_status': 'Ordinary uniform proofs with same-author exact finite implementation and review. Separate-author and Lean verification outstanding.',
        'worker_status': 'All three existing agents inspected this pass; terminal usage-limit errors remain.',
        'new_lean_jobs': 0, 'new_hosted_proof_jobs': 0,
        'goal_status': 'active and unachieved',
        'scope_limits': ['No upper bound on the centered distinct-coordinate aggregate.',
                         'Distinct tuples may have opposite pairs; this is not simply D6.',
                         'No improved full-energy or period exponent, classical Paley proof, or prize certificate.',
                         'No literature novelty claim; logarithmic-depth constants can be large.'],
    }
    (ROOT/'results/parallel33_pass_audit_2026_09_05.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps({'status': audit['status'], 'artifacts': len(artifacts),
                      'checks': 740, 'local_links': len(links)}, indent=2))


if __name__ == '__main__':
    main()
