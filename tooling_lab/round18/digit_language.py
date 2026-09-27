#!/usr/bin/env python3
"""Compile the generator-2 centered section into digits plus a residue state.

Alternating greedy signed binary expansions are established prior art. This
prototype retains the twisted endpoint and cofactor selection of actual cosets.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import sys
from time import perf_counter
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from norm_carry import norm, multiply


def recognize(digits, k):
    """Exact language membership for Q=2^N+1 and p=Q/k odd."""
    if not digits or type(k) is not int or k < 1 or (2**len(digits)+1) % k:
        raise ValueError('cofactor must be a positive divisor of 2^N+1')
    if any(type(v) is not int or v not in (-1, 0, 1) for v in digits):
        return {'status': 'invalid_digit'}
    nz = [v for v in digits if v]
    if not nz:
        return {'status': 'zero', 'scalar_centered': 0, 'lift_scalar': 0}
    if any(a == b for a, b in zip(nz, nz[1:])):
        return {'status': 'sign_alternation'}
    if len(nz) % 2 == 0:
        return {'status': 'twisted_endpoint'}
    S = sum(v*2**(len(digits)-1-j) for j, v in enumerate(digits))
    if S % k:
        return {'status': 'cofactor_congruence', 'lift_scalar': S}
    return {'status': 'nonzero', 'scalar_centered': S//k, 'lift_scalar': S}


def scalar_digits(a, p, N):
    F = []
    for j in range(N):
        value = a*pow(2, j, p) % p
        F.append(value if 2*value < p else value-p)
    digits = [(2*F[j]-(F[j+1] if j+1 < N else -F[0]))//p for j in range(N)]
    return F, digits


def profiles(N, k):
    """Five sign states, residue mod k, and a digit-support counter.

    State(first,last,residue,weight) records the first and last nonzero
    signs. Zeros preserve the sign state. A further nonzero must flip it.
    Residues are updated most-significant-digit first by r -> 2*r+d.
    """
    dp = {(0, 0, 0, 0): 1}; layers = []
    for position in range(N):
        nxt = defaultdict(int); transitions = 0
        for (first, last, residue, weight), count in dp.items():
            nxt[first, last, 2*residue % k, weight] += count; transitions += 1
            for digit in ((-1, 1) if first == 0 else (-last,)):
                nxt[first or digit, digit, (2*residue+digit) % k, weight+1] += count
                transitions += 1
        dp = dict(nxt)
        layers.append({'length': position+1, 'states': len(dp),
                       'transitions': transitions, 'prefix_words': sum(dp.values())})
    result = [[0]*(N+1) for _ in range(k)]
    for (first, last, residue, weight), count in dp.items():
        if first == last:  # Includes the one empty/zero path.
            result[residue][weight] += count
    assert sum(map(sum, result)) == 2**N+1
    assert [sum(row[w] for row in result) for w in range(N+1)] == [
        1 if w == 0 else 2*comb(N, w) if w % 2 else 0 for w in range(N+1)]
    return result, layers


def low_support_words(N, k, maximum):
    scanned = 0; accepted = []
    for weight in range(1, maximum+1, 2):
        for support in combinations(range(N), weight):
            scanned += 2
            S = sum((-1 if i % 2 else 1)*2**(N-1-j) for i, j in enumerate(support))
            if S % k:
                continue
            for sign in (-1, 1):
                digits = [0]*N
                for i, j in enumerate(support):
                    digits[j] = sign*(-1 if i % 2 else 1)
                accepted.append({'digits': digits, 'weight': weight, 'scalar_centered': sign*S//k})
    return accepted, scanned


def main():
    cases = []
    for N, p in [(4, 17), (8, 257), (16, 65537), (32, 6700417)]:
        Q = 2**N+1; assert Q % p == 0 and sp.isprime(p) and pow(2, N, p) == p-1
        k = Q//p; started = perf_counter(); distribution, layers = profiles(N, k)
        elapsed = perf_counter()-started
        assert all(sum(row) == p for row in distribution)
        positive_weights = [w for w, count in enumerate(distribution[0]) if w and count]
        row = {'N': N, 'n': 2*N, 'p': p, 'Q': Q, 'k': k,
               'residue_support_profiles': distribution, 'layers': layers,
               'actual_section_words': sum(distribution[0]),
               'minimum_nonzero_digit_support': min(positive_weights),
               'maximum_digit_support': max(positive_weights),
               'total_dp_transitions': sum(layer['transitions'] for layer in layers),
               'maximum_dp_states': max(layer['states'] for layer in layers),
               'producer_elapsed_seconds_uncontrolled': elapsed}
        cases.append(row)
        print(json.dumps({key: value for key, value in row.items()
                          if key not in ('residue_support_profiles', 'layers')}), flush=True)

    small_controls = []
    for N, p in [(4, 17), (8, 257)]:
        k = (2**N+1)//p; counts = Counter(); accepted = {}
        for d in product((-1, 0, 1), repeat=N):
            result = recognize(d, k); counts[result['status']] += 1
            if result['status'] in ('zero', 'nonzero'):
                a = result['scalar_centered'] % p
                assert a not in accepted
                F, expected = scalar_digits(a, p, N)
                assert list(d) == expected
                accepted[a] = list(d)
        assert sorted(accepted) == list(range(p))
        small_controls.append({'N': N, 'p': p, 'cube_words_checked': 3**N,
                               'status_counts': dict(sorted(counts.items()))})

    largest = cases[-1]; N, p, k = (largest[key] for key in ('N', 'p', 'k'))
    low, scanned = low_support_words(N, k, 5)
    assert len(low) == sum(largest['residue_support_profiles'][0][1:6])
    tau_by_weight = defaultdict(Counter)
    for row in low:
        a = row['scalar_centered'] % p; F, expected = scalar_digits(a, p, N)
        assert row['digits'] == expected
        nd = norm(row['digits']); assert nd > 0 and nd % k == 0
        row.update({'a': a, 'digit_norm': nd, 'norm_defect': nd//k,
                    'coefficient_energy': sum(v*v for v in F)})
        assert multiply([2]+[0]*(N-2)+[1], F) == [p*v for v in row['digits']]
        tau_by_weight[row['weight']][row['norm_defect']] += 1
    assert len({r['a'] for r in low}) == len(low)

    # Match the two saved g=2 walks and test mutations without presuming all edits fail.
    previous = HERE.parent/'round17/norm_compression.json'
    old = json.loads(previous.read_text()); saved = []; rejected = []
    for case in old['cases']:
        if case['g'] != 2:
            continue
        for r in case['records']:
            result = recognize(r['digits'], k)
            assert result['status'] == 'nonzero' and result['scalar_centered'] % p == r['a']
            saved.append({'case': case['name'], 'a': r['a'], 'digits': r['digits'], 'result': result})
        for r in case['records'][:8]:
            for j in (0, N//2, N-1):
                digits = list(r['digits']); digits[j] = -digits[j] if digits[j] else 1
                rejected.append({'digits': digits, **recognize(digits, k)})
    out = {'status': 'produced', 'scope': 'Exact generator-2 digit language and cofactor-residue transfer counts, including the entire p6700417 arithmetic section. Digit support is not original coefficient energy, and no general-generator or uniform norm bound is inferred.',
           'cases': cases, 'small_cube_controls': small_controls, 'saved_g2_records': saved,
           'mutated_probes': rejected, 'low_support': {'maximum_support': 5,
               'candidate_signed_support_words_examined': scanned,
               'actual_records': low,
               'norm_defect_counts_by_support': {str(w): [[tau, count] for tau, count in sorted(counts.items())]
                                                  for w, counts in sorted(tau_by_weight.items())}},
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['digit_language.py', '../round17/norm_carry.py']},
           'input_sha256': {'../round17/norm_compression.json': sha256(previous.read_bytes()).hexdigest()}}
    (HERE/'digit_language.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print(json.dumps({'low_support_records': len(low), 'signed_support_words_examined': scanned,
                      'norm_defect_counts_by_support': out['low_support']['norm_defect_counts_by_support']}), flush=True)


if __name__ == '__main__':
    main()
