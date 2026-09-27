#!/usr/bin/env python3
"""Independent exact certificate and dense arithmetic replay; no producer import.

The basis spanning proof uses a fraction-free determinant. Scalar candidates
are formed by successive modular set sums, not producer Cartesian iteration.
"""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from math import gcd, prod
import json
from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent


def determinant(matrix):
    A = [row[:] for row in matrix]; n = len(A); sign = 1; previous = 1
    if not n: return 1
    for j in range(n-1):
        if A[j][j] == 0:
            pivot = next(i for i in range(j+1, n) if A[i][j])
            A[j], A[pivot] = A[pivot], A[j]; sign = -sign
        pivot = A[j][j]
        for i in range(j+1, n):
            for k in range(j+1, n):
                numerator = A[i][k]*pivot-A[i][j]*A[j][k]
                assert numerator % previous == 0
                A[i][k] = numerator//previous
            A[i][j] = 0
        previous = pivot
    return sign*A[-1][-1]


def rational(pair):
    assert len(pair) == 2 and all(type(v) is int for v in pair)
    a, b = pair; assert b > 0 and gcd(abs(a), b) == 1
    return Fraction(a, b)


def certify(c):
    p, g, N, f, k = c['p'], c['g'], c['N'], c['relation'], c['k']
    assert sp.isprime(p) and p % 2 and N >= 4 and N & (N-1) == 0 and len(f) == N
    assert pow(g, N, p) == p-1 and sum(f[j]*pow(g, -j, p) for j in range(N)) % p == 0
    matrix = [[f[i-j] if i >= j else -f[N+i-j] for j in range(N)] for i in range(N)]
    product_vector = [sum(row[j]*c['adj'][j] for j in range(N)) for row in matrix]
    assert product_vector == [p*k]+[0]*(N-1) and c['norm_f'] == p*k and k > 0
    w = [c['adj'][0]]+[-c['adj'][N-j] for j in range(1, N)]; assert w == c['projection_weights']
    B = p//2*sum(map(abs, f))//p; assert B == c['digit_bound']
    A = c['erased']; assert A == sorted(set(A)) and all(type(i) is int and 0 <= i < N for i in A)
    d = len(A); divisor = gcd(k, *(w[j] for j in A)); assert divisor == c['kernel_gcd']
    assert c['kernel_index'] == k//divisor
    lift = c['bezout_lift']; assert len(lift) == d and all(type(v) is int for v in lift)
    assert (sum(w[j]*x for j, x in zip(A, lift))-divisor) % k == 0
    basis = c['basis']; assert len(basis) == d and all(len(row) == d and all(type(v) is int for v in row) for row in basis)
    assert abs(determinant(basis)) == k//divisor
    assert len(c['inverse']) == d and all(len(row) == d for row in c['inverse'])
    inv = [[rational(x) for x in row] for row in c['inverse']]
    for i in range(d):
        for j in range(d): assert sum(basis[i][t]*inv[t][j] for t in range(d)) == int(i == j)
    steps = []
    for row in basis:
        dot = sum(w[j]*x for j, x in zip(A, row)); assert dot % k == 0
        steps.append(dot//k % p)
    assert steps == c['scalar_steps']
    radii = [B*sum(abs(inv[i][j]) for i in range(d)) for j in range(d)]
    assert [rational(x) for x in c['coordinate_radii']] == radii
    visible = [j for j in range(d) if steps[j]]
    assert visible == c['visible_directions'] and [j for j in range(d) if not steps[j]] == c['invisible_directions']
    cap = prod((2*radii[j]).numerator//(2*radii[j]).denominator+1 for j in visible)
    assert cap == c['universal_candidate_box_cap']
    assert c['universal_unique_completion'] == all(2*radii[j] < 1 for j in visible)
    limits = {'scalar_product_plus_centering': (p-1)**2+p//2,
              'dense_convolution_absolute_sum': p//2*sum(map(abs, f))}
    assert max(limits.values()) < np.iinfo(np.int64).max
    powers = np.array([pow(g, j, p) for j in range(N)], dtype=np.int64); mat = np.array(matrix, dtype=np.int64)
    def encode_batch(scalars):
        result = []
        for start in range(0, len(scalars), 1024):
            a = np.array(scalars[start:start+1024], dtype=np.int64)
            F = (a[:, None]*powers+p//2) % p-p//2
            values = F @ mat.T; assert np.all(values % p == 0)
            result.extend((values//p).tolist())
        return result
    return inv, radii, encode_batch, limits


def replay(c, known, budget, prepared):
    inv, radii, encoder, _ = prepared; p, N, k, A = c['p'], c['N'], c['k'], c['erased']; B = c['digit_bound']
    assert set(known) == set(range(N))-set(A)
    empty = {'status': 'complete', 'count': 0, 'completions': [], 'candidate_box_size': 0, 'scalar_candidates_checked': 0}
    if any(abs(v) > B for v in known.values()): return {**empty, 'reason': 'known_digit_height'}
    w = c['projection_weights']; known_dot = sum(w[j]*v for j, v in known.items()); rhs = -known_dot % k
    if rhs % c['kernel_gcd']: return {**empty, 'reason': 'syndrome_divisibility'}
    lift = [rhs//c['kernel_gcd']*v for v in c['bezout_lift']]; intervals = []
    for j in range(len(A)):
        origin = sum(lift[i]*inv[i][j] for i in range(len(A)))
        low, high = -origin-radii[j], -origin+radii[j]
        lo = low.numerator//low.denominator
        if Fraction(lo) < low: lo += 1
        intervals.append([lo, high.numerator//high.denominator])
    if any(lo > hi for lo, hi in intervals): return {**empty, 'reason': 'empty_integer_coordinate_interval', 'intervals': intervals}
    size = prod(intervals[j][1]-intervals[j][0]+1 for j in c['visible_directions'])
    if size > budget:
        return {'status': 'budget_exceeded', 'reason': 'candidate_box_too_large', 'candidate_box_size': size,
                'candidate_budget': budget, 'intervals': intervals}
    dot = known_dot+sum(w[j]*v for j, v in zip(A, lift)); assert dot % k == 0
    candidates = {dot//k % p}
    for j in c['visible_directions']:
        lo, hi = intervals[j]; step = c['scalar_steps'][j]
        candidates = {(a+z*step) % p for a in candidates for z in range(lo, hi+1)}
    scalars = sorted(candidates); words = encoder(scalars)
    found = [{'scalar': a, 'digits': word} for a, word in zip(scalars, words) if all(word[j] == v for j, v in known.items())]
    return {'status': 'complete', 'reason': 'exhausted_projected_candidate_box', 'count': len(found), 'completions': found,
            'candidate_box_size': size, 'scalar_candidates_checked': len(candidates), 'intervals': intervals}


def main():
    path = HERE/'erasure_summary.json'; summary = json.loads(path.read_text()); reports = []; sample = None
    inputs = {path.name: sha256(path.read_bytes()).hexdigest()}
    for row in summary['cases']:
        path = HERE/row['artifact']; r = json.loads(path.read_text()); inputs[path.name] = sha256(path.read_bytes()).hexdigest()
        c = r['certificate']; prepared = certify(c); queries = 0; scalar_checks = 0; exhaustive = 0; cyclic = 0
        direct = prepared[2](list(range(c['p']))) if c['p'] < 100 else None
        all_known = set()
        for record in r['records']:
            known = dict(record['known']); expected = replay(c, known, r.get('candidate_budget', 100000), prepared)
            assert expected == record['output']; queries += 1; scalar_checks += expected.get('scalar_candidates_checked', 0)
            if direct is not None:
                true = [a for a, word in enumerate(direct) if all(word[j] == v for j, v in known.items())]
                assert expected['status'] == 'complete' and [x['scalar'] for x in expected['completions']] == true
                exhaustive += 1
            all_known.add(tuple(record['known'][i][1] for i in range(len(record['known']))))
        for record in r.get('initial_budget_controls', []):
            assert replay(c, dict(record['known']), 20000, prepared) == record['output']
        if r['name'].startswith(('p17_basic_', 'p41_general_', 'p97_general_')):
            assert len(all_known) == len(r['records']) == (2*c['digit_bound']+1)**(c['N']-len(c['erased']))
        for record in r.get('cyclic_controls', []):
            s, known, N = record['start'], dict(record['known']), c['N']; d = len(c['erased'])
            rotated = {j: known[(j+s) % N]*(1 if j+s < N else -1) for j in range(d, N)}
            expected = replay(c, rotated, 100000, prepared)
            assert expected['status'] == 'complete'
            scalars = sorted(x['scalar']*pow(c['g'], -s, c['p']) % c['p'] for x in expected['completions'])
            words = prepared[2](scalars)
            expected['completions'] = [{'scalar': a, 'digits': word} for a, word in zip(scalars, words)]
            assert expected == record['output'] and scalars == [record['anchor_scalar']]
            assert all(all(word[j] == v for j, v in known.items()) for word in words)
            cyclic += 1
        if 'bounded_syndrome_alias' in r:
            alias = r['bounded_syndrome_alias']; f = c['relation']; assert alias['digits'] == f and max(map(abs, f)) <= c['digit_bound']
            assert all(f[j] == 0 for j, _ in alias['known'])
            assert alias['lifted_first_coordinate'] == c['p'] > c['p']//2
            expected = replay(c, dict(alias['known']), 100000, prepared)
            assert expected == alias['decoded'] and expected['completions'] == [{'scalar': 0, 'digits': [0]*c['N']}]
            sample = (r, prepared)
        out = {'name': r['name'], 'certificate_verified': True, 'queries_replayed': queries,
               'scalar_candidates_checked': scalar_checks, 'full_small_codebook_oracles': exhaustive,
               'cyclic_controls_checked': cyclic, 'universal_unique_completion': c['universal_unique_completion'],
               'initial_budget_controls_checked': len(r.get('initial_budget_controls', [])),
               'universal_candidate_box_cap': c['universal_candidate_box_cap'], 'overflow_bounds': prepared[3]}
        reports.append(out); print(json.dumps(out), flush=True)
    assert sample; r, prepared = sample; bad = []
    for name, mutate in [
            ('basis_entry', lambda c: c['basis'][0].__setitem__(0, c['basis'][0][0]+1)),
            ('inverse_entry', lambda c: c['inverse'][0][0].__setitem__(0, c['inverse'][0][0][0]+1)),
            ('hidden_scalar_direction', lambda c: c['scalar_steps'].__setitem__(c['visible_directions'][0], 0)),
            ('understated_radius', lambda c: c['coordinate_radii'].__setitem__(c['visible_directions'][0], [0, 1])),
            ('false_bezout', lambda c: c['bezout_lift'].__setitem__(0, c['bezout_lift'][0]+1)),
            ('incorrect_gcd', lambda c: c.__setitem__('kernel_gcd', c['kernel_gcd']+1))]:
        changed = deepcopy(r['certificate']); mutate(changed)
        try: certify(changed)
        except (AssertionError, ZeroDivisionError, StopIteration): bad.append(name)
        else: raise AssertionError('corrupt certificate accepted: '+name)
    out = {'status': 'passed', 'scope': 'Exact basis-spanning, inverse, radius and scalar-step certificates; independent complete candidate-set replay; full small codebook oracles; cyclic transport and a false centered alias. LLL quality is never assumed.',
           'cases': reports, 'corrupt_certificates_rejected': bad,
           'source_sha256': {'erasure_review.py': sha256(Path(__file__).read_bytes()).hexdigest()}, 'input_sha256': inputs,
           'numpy_version': np.__version__}
    (HERE/'erasure_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
