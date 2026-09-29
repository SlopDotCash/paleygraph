#!/usr/bin/env python3
"""Pass35 certificate readback with a second root-exclusion check."""
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


def sign(x):
    return (x > 0)-(x < 0)


def radical_sign(coeffs, D, side):
    # Horner reduction in Q[t]/(t^2-D), then choose t=side*sqrt(D).
    a = b = Fraction(0)
    for c in reversed(coeffs):
        a, b = c+D*b, a
    b *= side
    if not a:
        return sign(b)
    if not b or sign(a) == sign(b):
        return sign(a)
    diff = a*a-D*b*b
    return sign(a)*sign(diff)


def main():
    prior = read('results/parallel35_prior_state_2026_09_05.json')
    data = read('results/parallel35_verification_2026_09_05.json')
    source = read('results/parallel35_source_scope_2026_09_05.json')
    for group in [prior['prior_sha256'],data['input_sha256'],source['local_input_sha256']]:
        for p,h in group.items():
            assert digest(p) == h,p
    for item in source['sources']:
        assert digest(item['path']) == item['sha256']
    assert digest('sources/manifest.json') == source['main_manifest_sha256']
    assert sum(data['check_counts'].values()) == 1053
    assert len(data['rational_rows']) == 114 and len(data['finite_field_rows']) == 43
    assert len(data['newton_group_cases']) == 4
    enclosures = derivative_signs = high_sum = 0
    for row in data['rational_rows']:
        N,s = row['N'],row['s']
        n,r = 2*N,2*s
        assert sum(item['count'] for item in row['values']) == N
        assert all(-2 <= q(item['value']) <= 2 for item in row['values'])
        total = sum(q(item['value'])*item['count'] for item in row['values'])
        assert total == q(row['sum'])
        # Grouped binomial products, unlike the per-entry DP in the verifier.
        e = [Fraction(1)]+[Fraction(0)]*r
        for item in row['values']:
            v,m = q(item['value']),item['count']
            e = [sum(e[k-j]*comb(m,j)*v**j for j in range(min(k,m)+1))
                 for k in range(r+1)]
        H,M = factorial(r)*e[r],total**r
        assert H == q(row['H']) and M == q(row['moment'])
        assert H >= -(256*s*n)**s
        assert abs(H) <= 4**s*(M+(64*s*n)**s)
        assert M <= 16**s*H+2*(4096*s*n)**s
        if 'shifted_polynomial' in row:
            assert r <= N//2
            p = [q(c) for c in row['shifted_polynomial']]
            a = Fraction(prod(N-j for j in range(r)),N**r)
            assert p[-1] == a == q(row['normalization'])
            assert Fraction(1,2**r) <= a <= 1
            assert sum(c*total**i for i,c in enumerate(p)) == H
            D = row['root_radius_squared']
            assert D == 64*r*N
            # Every coefficient in p(R+x) and (-1)^r p(-R-x) is
            # nonnegative. Hence there is no real root beyond either
            # endpoint. This uses Taylor signs, not Sturm variation.
            for j in range(r+1):
                derivative = [p[k]*(factorial(k)//factorial(k-j)) for k in range(j,r+1)]
                assert radical_sign(derivative,D,1) >= 0
                assert (-1)**(r-j)*radical_sign(derivative,D,-1) >= 0
                derivative_signs += 2
            enclosures += 1
            high_sum += row['high_sum_branch']
    assert enclosures == 78 and high_sum == 3
    for row in data['finite_field_rows']:
        n,s = row['n'],row['s']
        B,T = q(row['B']),q(row['T'])
        assert B >= -(256*s*n)**s
        assert abs(B) <= 4**s*(T+(64*s*n)**s)
        assert T <= 16**s*B+2*(4096*s*n)**s
    w = data['averaged_polynomial_witness']
    assert w['p'] == 1153 and w['n'] == 8
    assert w['p_times_centered_coefficients'] == [1152,-8,-24,-32,-16]
    assert w['translated_coefficients'] == [1153,0,0,0,-1]
    assert data['ten_term_lower']['constant'] == 90+1280**5

    central = ['README.md','research/frontier.md','research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    notes = ['research/parallel35-single-degree-comparison-2026-09-05.md',
             'research/parallel35-pass-summary-2026-09-05.md']
    links = []
    for p in central+notes:
        for target in re.findall(r'\]\(([^)]+)\)',(ROOT/p).read_text()):
            if '://' in target or target.startswith('#') or ('/' not in target and '.' not in target):
                continue
            dest = ((ROOT/p).parent/target.split('#')[0]).resolve()
            assert dest.exists(),(p,target)
            links.append(dict(source=p,target=str(dest.relative_to(ROOT))))
    artifacts = notes+[
        'experiments/parallel35_verify_2026_09_05.py',
        'experiments/parallel35_final_audit_2026_09_05.py',
        'results/parallel35_prior_state_2026_09_05.json',
        'results/parallel35_verification_2026_09_05.json',
        'results/parallel35_source_scope_2026_09_05.json',
    ]+[item['path'] for item in source['sources']]
    audit = dict(status='single-degree certificates and source scope passed; full targets unproved',
                 audited_at_utc=datetime.now(timezone.utc).isoformat(),
                 previous_turn_classification=prior['previous_turn_classification'],
                 current_turn_classification='progress: unconditional negative-side Gaussian bound and same-degree opposite-free moment criterion',
                 prior_inputs_preserved=prior['prior_sha256'],
                 artifact_sha256={p:digest(p) for p in artifacts},
                 central_file_sha256={p:digest(p) for p in central},
                 source_scope=source,local_links_checked=links,
                 bounded_exact_checks=1053,rational_rows_rechecked=114,
                 finite_field_rows_rechecked=43,
                 alternate_root_enclosures=enclosures,
                 exact_Taylor_derivative_signs=derivative_signs,
                 proof_status='Ordinary uniform application of a primary-source theorem; same-author exact checks and certificate review. Separate-author and Lean verification outstanding.',
                 worker_status='All three prior agents freshly inspected: terminal usage-limit errors. Root continues.',
                 new_lean_jobs=0,new_hosted_proof_jobs=0,goal_status='active and unachieved',
                 scope_limits=['No positive upper bound on B_s.','No improvement to period exponent71/72.',
                               'Degree-ten estimate is a lower bound with a large constant, not two-sided equidistribution.',
                               'No classical Paley proof, full prize reduction, or prize certificate.'])
    (ROOT/'results/parallel35_pass_audit_2026_09_05.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(dict(status=audit['status'],artifacts=len(artifacts),
          checks=1053,root_enclosures=enclosures,Taylor_signs=derivative_signs,
          local_links=len(links)),indent=2))


if __name__ == '__main__':
    main()
