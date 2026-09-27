#!/usr/bin/env python3
"""Measure short-relation windows against exact antiperiodic digit memory."""
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def window_maxima(digits):
    N = len(digits); doubled = list(digits)+[-v for v in digits]
    values = [0]*N; maxima = []
    for s in range(1, N+1):
        values = [2*v+doubled[j+s-1] for j, v in enumerate(values)]
        maxima.append(max(map(abs, values)))
    return maxima


def actual(digits):
    nz = [v for v in digits if v]
    return not nz or len(nz) % 2 == 1 and all(a != b for a, b in zip(nz, nz[1:]))


def row(digits):
    maxima = window_maxima(digits); N = len(digits)
    failures = [s for s, v in enumerate(maxima, 1) if v > 2**(s-1)]
    return {'digits': list(digits), 'window_maxima': maxima,
            'first_failure': min(failures) if failures else None,
            'realized_over_Q': actual(digits),
            'odd_support_or_zero': not any(digits) or sum(abs(v) for v in digits) % 2 == 1,
            'unit_digit_sum_or_zero': not any(digits) or abs(sum(digits)) == 1}


def main():
    cases = []
    for N in (4, 8):
        rows = [row(d) for d in product((-1, 0, 1), repeat=N)]
        profiles = []
        for depth in range(1, N+1):
            remaining = [r for r in rows if r['first_failure'] is None or r['first_failure'] > depth]
            profiles.append({'depth': depth, 'windows_only': len(remaining),
                             'with_odd_support': sum(r['odd_support_or_zero'] for r in remaining),
                             'with_unit_digit_sum': sum(r['unit_digit_sum_or_zero'] for r in remaining)})
        assert profiles[-1]['windows_only'] == 2**N+1
        assert profiles[N//2-1]['with_odd_support'] == 2**N+1
        assert all(r['realized_over_Q'] == (r['first_failure'] is None) for r in rows)
        cases.append({'N': N, 'Q': 2**N+1, 'cube_words': len(rows), 'profiles': profiles, 'records': rows})
        print(json.dumps({k: v for k, v in cases[-1].items() if k != 'records'}), flush=True)
    sharp = []
    for N in (4, 8, 16, 32, 64, 128):
        even = [1, -1]+[0]*(N-2)
        odd = [0]*N; odd[0] = odd[N//2-1] = 1; odd[N//2] = -1
        er, od = row(even), row(odd)
        assert er['first_failure'] == N and not er['realized_over_Q']
        assert od['first_failure'] == N//2 and od['odd_support_or_zero'] and od['unit_digit_sum_or_zero']
        assert not od['realized_over_Q']
        sharp.append({'N': N, 'Q': 2**N+1, 'without_parity': er, 'with_parity_and_unit_sum': od})
    out = {'status': 'produced', 'scope': 'Complete tiny digit-cube depth hierarchy and exact symbolic witness families over Q=2^N+1. The large-Q family does not assert infinitely many prime Fermat moduli or a Paley obstruction.',
           'cases': cases, 'sharpness_controls': sharp,
           'source_sha256': {'window_hierarchy.py': sha256((HERE/'window_hierarchy.py').read_bytes()).hexdigest()}}
    (HERE/'window_hierarchy.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')


if __name__ == '__main__': main()
