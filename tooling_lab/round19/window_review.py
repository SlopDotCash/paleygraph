#!/usr/bin/env python3
"""Independent lifted-coordinate check of every window and inverse."""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def check(row, N):
    Q = 2**N+1; digits = row['digits']
    assert len(digits) == N and all(type(v) is int and v in (-1, 0, 1) for v in digits)
    F = [sum(v*2**(N-1-j) for j, v in enumerate(digits))]
    for d in digits: F.append(2*F[-1]-Q*d)
    assert F[-1] == -F[0]
    F = F[:-1]; doubled = F+[-v for v in F]; maxima = []
    for s in range(1, N+1):
        raw = [2**s*F[j]-doubled[j+s] for j in range(N)]
        assert all(v % Q == 0 for v in raw)
        maxima.append(max(abs(v//Q) for v in raw))
    assert maxima == row['window_maxima']
    failures = [s for s, v in enumerate(maxima, 1) if v > ((Q-1)*(2**s+1))//(2*Q)]
    assert row['first_failure'] == (min(failures) if failures else None)
    centered = all(2*abs(v) < Q for v in F)
    assert centered == row['realized_over_Q']
    assert all((v-F[0]*pow(2, j, Q)) % Q == 0 for j, v in enumerate(F))
    assert row['odd_support_or_zero'] == (not any(digits) or sum(v*v for v in digits) % 2 == 1)
    assert row['unit_digit_sum_or_zero'] == (not any(digits) or sum(digits) in (-1, 1))
    return {'centered': centered, 'inverse_coefficients': F}


def main():
    source = HERE/'window_hierarchy.json'; data = json.loads(source.read_text()); summaries = []
    for case in data['cases']:
        N = case['N']; assert case['Q'] == 2**N+1 and len(case['records']) == case['cube_words'] == 3**N
        actual_scalars = []
        for digits, row in zip(product((-1, 0, 1), repeat=N), case['records']):
            assert list(digits) == row['digits']; result = check(row, N)
            if result['centered']: actual_scalars.append(result['inverse_coefficients'][0] % case['Q'])
        assert sorted(actual_scalars) == list(range(case['Q']))
        for profile in case['profiles']:
            d = profile['depth']; survivors = [r for r in case['records']
                                               if all(r['window_maxima'][s-1] <= 2**(s-1) for s in range(1, d+1))]
            assert profile['windows_only'] == len(survivors)
            assert profile['with_odd_support'] == sum(r['odd_support_or_zero'] for r in survivors)
            assert profile['with_unit_digit_sum'] == sum(r['unit_digit_sum_or_zero'] for r in survivors)
        summaries.append({'N': N, 'cube_words_checked': len(case['records']),
                          'exact_actual_scalars': len(actual_scalars), 'profiles': case['profiles']})
    controls = []
    for c in data['sharpness_controls']:
        N = c['N']; assert c['Q'] == 2**N+1
        even, odd = check(c['without_parity'], N), check(c['with_parity_and_unit_sum'], N)
        assert not even['centered'] and not odd['centered']
        assert c['without_parity']['first_failure'] == N
        assert c['with_parity_and_unit_sum']['first_failure'] == N//2
        assert even['inverse_coefficients'][1] == -2**(N-1)-1
        controls.append({'N': N, 'no_parity_inverse': even['inverse_coefficients'],
                         'parity_inverse': odd['inverse_coefficients']})
    out = {'status': 'passed', 'scope': 'Every window independently recovered as (2^s F_i-F_(i+s))/Q; all tiny inverse membership decisions and complete scalar coverage checked. Large witness moduli may be composite.',
           'cases': summaries, 'sharpness_controls': controls,
           'source_sha256': {'window_review.py': sha256((HERE/'window_review.py').read_bytes()).hexdigest()},
           'input_sha256': {source.name: sha256(source.read_bytes()).hexdigest()}}
    (HERE/'window_review.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'cube_words': sum(r['cube_words_checked'] for r in summaries),
                      'sharpness_witnesses': 2*len(controls)}), flush=True)


if __name__ == '__main__': main()
