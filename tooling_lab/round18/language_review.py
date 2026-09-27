#!/usr/bin/env python3
"""Independent binary-sign transfer matrix, full scalar census, norm review."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from carry_review import result_norm, product as polynomial_product


def binary_profiles(N, k):
    """Use all N binary sign words, whose next endpoint is 1-first.

    S=1-(2^(N-1)+1)*b0 + sum_(j=1)^(N-1) 2^(N-1-j)*bj.
    Support is the number of sign changes, including the twisted endpoint.
    This uses binary digits and fixed position weights, not signed-digit
    paths or their Horner residue update in the producer.
    """
    states = {(b, b, (1-(2**(N-1)+1)*b) % k, 0): 1 for b in (0, 1)}
    for j in range(1, N):
        next_states = defaultdict(int); place = pow(2, N-1-j, k)
        for (first, previous, residue, weight), count in states.items():
            for bit in (0, 1):
                next_states[first, bit, (residue+place*bit) % k, weight+(bit != previous)] += count
        states = dict(next_states)
    rows = [[0]*(N+1) for _ in range(k)]; rows[0][0] = 1
    for (first, previous, residue, weight), count in states.items():
        rows[residue][weight+((1-first) != previous)] += count
    return rows


def recurrence_membership(digits, k):
    N = len(digits); Q = 2**N+1; p = Q//k
    numerator = sum(v*2**(N-1-j) for j, v in enumerate(digits))
    if numerator % k:
        return False, None
    F = [numerator//k]
    for d in digits:
        F.append(2*F[-1]-p*d)
    assert F[-1] == -F[0]
    return all(2*abs(v) < p for v in F), F[:-1]


def direct_scalar_census(p, N):
    """All scalars, including zero, in bounded int64 NumPy batches."""
    assert 2*p < np.iinfo(np.int64).max
    counts = np.zeros(N+1, dtype=np.int64); low = []; coordinates = 0
    for start in range(0, p, 131072):
        a = np.arange(start, min(start+131072, p), dtype=np.int64)
        original = (a+p//2) % p-p//2; x = original.copy(); weight = np.zeros_like(a)
        for j in range(N):
            y = (2*x+p//2) % p-p//2
            delta = 2*x-y
            assert np.all(delta % p == 0)
            d = delta//p; assert np.all(np.abs(d) <= 1)
            weight += d*d; x = y
        assert np.array_equal(x, -original)
        counts += np.bincount(weight, minlength=N+1)
        if N == 32: low.extend(int(v) for v in a[(weight > 0) & (weight <= 5)])
        coordinates += len(a)*N
    assert int(counts.sum()) == p
    return [int(v) for v in counts], low, coordinates


def expected_status(digits, k):
    nz = [v for v in digits if v]
    if not nz: return 'zero'
    if any(x*y > 0 for x, y in zip(nz, nz[1:])): return 'sign_alternation'
    if nz[0] != nz[-1]: return 'twisted_endpoint'
    S = sum(d*2**(len(digits)-1-j) for j, d in enumerate(digits))
    return 'cofactor_congruence' if S % k else 'nonzero'


def main():
    source = HERE/'digit_language.json'; data = json.loads(source.read_text()); cases = []
    largest_low = []
    for c in data['cases']:
        N, p, Q, k = (c[key] for key in ('N', 'p', 'Q', 'k'))
        assert sp.isprime(p) and Q == 2**N+1 == p*k and pow(2, N, p) == p-1
        other = binary_profiles(N, k)
        assert other == c['residue_support_profiles']
        assert all(sum(row) == p for row in other)
        census, low, values = direct_scalar_census(p, N)
        assert census == other[0] and sum(census) == c['actual_section_words'] == p
        weights = [w for w, count in enumerate(census) if w and count]
        assert (min(weights), max(weights)) == (c['minimum_nonzero_digit_support'], c['maximum_digit_support'])
        # Prefix-word totals have a direct support-subset formula.
        assert len(c['layers']) == N
        for j, layer in enumerate(c['layers'], 1):
            assert layer['length'] == j and layer['prefix_words'] == 2**(j+1)-1
        cases.append({'p': p, 'N': N, 'all_residue_profiles_checked': k,
                      'actual_scalar_words_checked': p, 'scalar_coordinate_transitions': values,
                      'actual_support_distribution': census,
                      'minimum_nonzero_support': min(weights), 'maximum_support': max(weights)})
        if N == 32: largest_low = low
        print(json.dumps(cases[-1]), flush=True)
    controls = []
    for c in data['small_cube_controls']:
        N, p = c['N'], c['p']; k = (2**N+1)//p; statuses = Counter(); accepted = 0
        for digits in product((-1, 0, 1), repeat=N):
            member, F = recurrence_membership(digits, k)
            status = expected_status(digits, k); statuses[status] += 1
            assert member == (status in ('zero', 'nonzero'))
            if member:
                assert all((v-F[0]*pow(2, j, p)) % p == 0 for j, v in enumerate(F))
                accepted += 1
        assert dict(sorted(statuses.items())) == c['status_counts']
        assert accepted == p and sum(statuses.values()) == c['cube_words_checked']
        controls.append({'N': N, 'words_checked': sum(statuses.values()), 'accepted': accepted})
    for r in data['saved_g2_records']:
        member, F = recurrence_membership(r['digits'], 641)
        assert member and F[0] % 6700417 == r['a'] and r['result']['scalar_centered'] == F[0]
        assert r['result']['lift_scalar'] == 641*F[0] and r['result']['status'] == 'nonzero'
    for row in data['mutated_probes']:
        member, F = recurrence_membership(row['digits'], 641)
        assert row['status'] == expected_status(row['digits'], 641)
        assert member == (row['status'] in ('zero', 'nonzero'))

    low = data['low_support']['actual_records']; assert sorted(r['a'] for r in low) == largest_low
    assert len(low) == 512 and len(set(largest_low)) == len(low)
    by_weight = defaultdict(Counter); p, N, k = 6700417, 32, 641
    for i, row in enumerate(low):
        digits = row['digits']; member, F = recurrence_membership(digits, k)
        assert member and F[0] == row['scalar_centered'] and F[0] % p == row['a']
        assert sum(v*v for v in digits) == row['weight'] <= 5
        assert sum(v*v for v in F) == row['coefficient_energy']
        assert polynomial_product([2]+[0]*(N-2)+[1], F, N) == [p*v for v in digits]
        nd = result_norm(digits, N)
        assert nd == row['digit_norm'] == k*row['norm_defect']
        by_weight[row['weight']][row['norm_defect']] += 1
        if (i+1) % 128 == 0: print(f'low-support resultants checked {i+1}/{len(low)}', flush=True)
    assert data['low_support']['norm_defect_counts_by_support'] == {
        str(w): [[tau, count] for tau, count in sorted(cs.items())] for w, cs in sorted(by_weight.items())}
    # Boundary and congruence controls each isolate one indispensable condition.
    first = [1]+[0]*31
    endpoint = [0]*32
    for j, digit in zip((0, 1, 8, 31), (1, -1, 1, -1)): endpoint[j] = digit
    repeated = [0]*32
    for j in (22, 24, 31): repeated[j] = 1
    condition_controls = []
    for label, digits, expected in [('cofactor', first, 'cofactor_congruence'),
                                    ('endpoint', endpoint, 'twisted_endpoint'),
                                    ('alternation', repeated, 'sign_alternation')]:
        nz = [v for v in digits if v]
        S = sum(d*2**(31-j) for j, d in enumerate(digits))
        passes = {'cofactor': S % 641 == 0, 'endpoint': len(nz) % 2 == 1,
                  'alternation': all(a != b for a, b in zip(nz, nz[1:]))}
        assert passes == {name: name != label for name in passes}
        assert expected_status(digits, 641) == expected
        assert not recurrence_membership(digits, 641)[0]
        condition_controls.append({'condition_failed': label, 'digits': digits, 'lift_scalar': S,
                                   'conditions_passed': passes})
    out = {'status': 'passed', 'scope': 'Independent binary-sign transfer tables for every residue, direct modular doubling of every actual scalar, complete tiny-cube inverse checks, and independent resultants for the full minimum-support section. No general-prime asymptotic or prize proof.',
           'cases': cases, 'small_controls': controls, 'saved_g2_records_checked': len(data['saved_g2_records']),
           'mutations_checked': len(data['mutated_probes']), 'low_support_resultants': len(low),
           'necessary_condition_controls': condition_controls,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['language_review.py', '../round17/carry_review.py']},
           'input_sha256': {source.name: sha256(source.read_bytes()).hexdigest()}}
    (HERE/'language_review.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'actual_scalars': sum(c['actual_scalar_words_checked'] for c in cases),
                      'low_support_norms': len(low)}), flush=True)


if __name__ == '__main__':
    main()
