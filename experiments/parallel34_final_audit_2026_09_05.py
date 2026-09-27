#!/usr/bin/env python3
"""Check pass34 certificates and preserve scope; not a full-goal audit pass."""
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, prod
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def read(p):
    return json.loads((ROOT/p).read_text())


def digest(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()


def q(x):
    return Fraction(x['numerator'], x['denominator'])


def falling(n, r):
    return prod(n-j for j in range(r))


def main():
    prior = read('results/parallel34_prior_state_2026_09_05.json')
    result = read('results/parallel34_verification_2026_09_05.json')
    source = read('results/parallel34_source_scope_2026_09_05.json')
    resume = read('results/parallel34_resume_state_2026_09_05.json')
    for group in [prior['prior_sha256'], result['input_sha256'], source['local_input_sha256']]:
        for p, expected in group.items():
            assert digest(p) == expected, p
    for item in source['sources']:
        assert digest(item['path']) == item['sha256']
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    assert digest(resume['evidence_path']) == resume['evidence_sha256']
    assert sum(result['check_counts'].values()) == 9852
    assert len(result['cases']) == 9
    assert sum(len(row['rows']) for row in result['cases']) == 52
    assert sum(row['subgroup'] for row in result['cases']) == 8

    count_rows = 0
    for case in result['cases']:
        p, n, A = case['p'], case['n'], set(case['elements'])
        assert len(A) == n and 0 not in A and A == {(-a) % p for a in A}
        assert all(a*b % p in A for a in A for b in A) == case['subgroup']
        for row in case['rows']:
            s, r = row['s'], row['r']
            assert r == 2*s and r <= n//2 and r < p
            Q = Fraction(row['distinct'])-Fraction(falling(n, r), p)
            B = Fraction(row['opposite_free'])-Fraction(2**r*falling(n//2, r), p)
            assert q(row['centered_distinct']) == Q
            assert q(row['centered_opposite_free']) == B
            assert 0 <= row['opposite_free'] <= row['distinct'] <= row['energy']
            assert row['gaussian'] == factorial(r)*n**s//(2**s*factorial(s))
            assert Q == sum(factorial(r)//factorial(2*t)*comb(n//2-2*t, s-t)
                            *q(case['rows'][t]['centered_opposite_free']) for t in range(s+1))
            if s >= 2:
                mark = comb(r, 2)
                e, previous = row['energy'], case['rows'][s-1]['energy']
                assert previous**s <= e**(s-1)
                assert n*row['repeated']**2 <= mark**2*e**2
                assert row['repeated'] == e-row['distinct']
                assert row['marked_opposite_upper'] == mark*n*previous
                assert row['opposite_free'] >= e-row['repeated']-row['marked_opposite_upper']
            count_rows += 1
    assert any(c['p'] == 33713 and c['n'] == 16 and
               c['n']**4 <= 4*c['p'] <= 4*c['n']**4 for c in result['cases'])

    # Build the inverse recursively, without using the verifier's closed
    # inverse formula, then compare to its alternate rational expression.
    coefficient_checks = 0
    for item in result['coefficient_parameter_ranges']:
        N, depth = item['N'], item['maximum_s']
        C = [[factorial(2*s)//factorial(2*t)*comb(N-2*t, s-t)
              for t in range(s+1)] for s in range(depth+1)]
        inv = []
        for s in range(depth+1):
            row = []
            for t in range(s+1):
                v = int(s == t)-sum(C[s][k]*inv[k][t] for k in range(t, s))
                j, M = s-t, N-2*t
                c = Fraction(M, M-j)*comb(M-j, j) if j else Fraction(1)
                closed = (-1)**j*factorial(2*s)//factorial(2*t)*c
                assert v == closed
                row.append(v)
                coefficient_checks += 1
            inv.append(row)
    assert coefficient_checks == 2906
    threshold = result['ten_term_threshold']
    assert threshold['n_at_least'] == 2**35
    assert q(threshold['excluded_fraction_upper']) == 45*(Fraction(1, 128)+Fraction(1, 131072)) < Fraction(1, 2)

    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel34-opposite-pair-transform-2026-09-05.md',
             'research/parallel34-pass-summary-2026-09-05.md']
    links = []
    for p in central+notes:
        for target in re.findall(r'\]\(([^)]+)\)', (ROOT/p).read_text()):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            dest = ((ROOT/p).parent/target.split('#')[0]).resolve()
            assert dest.exists(), (p, target)
            links.append(dict(source=p, target=str(dest.relative_to(ROOT))))
    artifacts = notes+[
        'experiments/parallel34_verify_2026_09_05.py',
        'experiments/parallel34_final_audit_2026_09_05.py',
        'results/parallel34_prior_state_2026_09_05.json',
        'results/parallel34_resume_state_2026_09_05.json',
        'results/parallel34_verification_2026_09_05.json',
        'results/parallel34_source_scope_2026_09_05.json',
    ]
    audit = dict(status='opposite-pair certificates and provenance passed; full targets unproved',
                 audited_at_utc=datetime.now(timezone.utc).isoformat(),
                 previous_turn_classification=resume['previous_goal_turn_classification'],
                 current_turn_classification='progress: exact centered opposite-pair transforms, controlled Gaussian hierarchy, and forced opposite-free ten-term relations',
                 prior_inputs_preserved=prior['prior_sha256'],
                 artifact_sha256={p: digest(p) for p in artifacts},
                 central_file_sha256={p: digest(p) for p in central},
                 source_scope=source, local_links_checked=links,
                 bounded_exact_checks=9852, count_rows_rechecked=count_rows,
                 inverse_coefficients_checked_by_recurrence=coefficient_checks,
                 proof_status='Ordinary uniform arguments, same-author exact finite checks and review; separate-author review and Lean verification outstanding.',
                 worker_status=resume['agent_status'], new_lean_jobs=0,
                 new_hosted_proof_jobs=0, goal_status='active and unachieved',
                 scope_limits=['Missing one-sided centered opposite-free hierarchy upper bound.',
                               'No stronger full-energy or period exponent.',
                               'Large-order lower bound conditional on eligible field and set; no subgroup existence theorem.',
                               'No classical Paley proof, full prize reduction, or prize certificate.'])
    (ROOT/'results/parallel34_pass_audit_2026_09_05.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(dict(status=audit['status'], artifacts=len(artifacts),
          verifier_checks=9852, independent_inverse_checks=coefficient_checks,
          count_rows=count_rows, local_links=len(links)),indent=2))


if __name__ == '__main__':
    main()
