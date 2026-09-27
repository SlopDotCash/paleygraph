#!/usr/bin/env python3
"""Independent batched dense-matrix replay of every reported cross-splice.

The producer uses Horner evaluation, iterative residues, and sparse Python
convolution. This reviewer uses direct powers and NumPy integer matrices,
with explicit int64 bounds before any batch arithmetic.
"""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent


def check_witnesses(r):
    order, N, cut = r['order'], r['N'], r['cut']; M = len(r['scalars'])
    assert sorted(order) == list(range(N)) and 0 < cut < N
    rows = r['cross_acceptance_rows']; selected = r['clique_indices']
    assert len(rows) == M and all(len(row) == M and set(row) <= {'0', '1'} for row in rows)
    assert all(rows[i][i] == '1' for i in range(M))
    assert len(set(selected)) == len(selected) and all(0 <= i < M for i in selected)
    assert r['certified_width_lower_bound'] == len(selected)
    bits = r['pair_witness_selectors']; pos = 0
    assert len(bits) == len(selected)*(len(selected)-1)//2 == r['certified_distinguished_pairs']
    for n, i in enumerate(selected):
        for j in selected[n+1:]:
            if bits[pos] == '0': assert rows[i][j] == '0' and rows[j][j] == '1'
            else: assert bits[pos] == '1' and rows[j][i] == '0' and rows[i][i] == '1'
            pos += 1
    assert r['cross_splices_checked'] == M*M
    assert r['accepted_cross_splices'] == sum(row.count('1') for row in rows)
    prefixes = {tuple(w[j] for j in order[:cut]) for w in r['words_in_original_coordinates']}
    assert len(prefixes) == r['distinct_sampled_prefixes']


def check_dense(r):
    check_witnesses(r)
    N, p, g, f = r['N'], r['p'], r['g'], r['relation']; m = p//2; M = len(r['scalars'])
    assert N >= 4 and N & (N-1) == 0 and sp.isprime(p) and pow(g, N, p) == p-1
    u = pow(g, -1, p); weights = [pow(u, j, p) for j in range(N)]
    assert sum(x*y for x, y in zip(weights, f)) % p == 0
    digit_bound = m*sum(map(abs, f))//p
    limits = {'scalar_times_power_plus_m': (p-1)**2+m,
              'dense_convolution_absolute_sum': m*sum(map(abs, f)),
              'projection_absolute_sum': N*digit_bound*(p-1),
              'projection_scale': (p-1)**2}
    assert all(v < np.iinfo(np.int64).max for v in limits.values())
    matrix = np.array([[f[i-j] if i >= j else -f[N+i-j] for j in range(N)] for i in range(N)], dtype=np.int64)
    powers = np.array([pow(g, j, p) for j in range(N)], dtype=np.int64)
    weights = np.array(weights, dtype=np.int64)
    def reencode(scalars):
        F = (scalars[:, None]*powers+m) % p-m
        numerator = F @ matrix.T; assert np.all(numerator % p == 0)
        return numerator//p
    beta = int((reencode(np.array([1], dtype=np.int64))[0] @ weights) % p)
    assert beta != 0; inv_beta = pow(beta, -1, p)
    words = np.array(r['words_in_original_coordinates'], dtype=np.int64)
    assert words.shape == (M, N) and np.all(np.abs(words) <= digit_bound)
    scalars = r['scalars']; assert len(set(scalars)) == M and all(type(a) is int and 0 <= a < p for a in scalars)
    assert np.array_equal(reencode(np.array(scalars, dtype=np.int64)), words)
    prefix = r['order'][:r['cut']]; flat = ''.join(r['cross_acceptance_rows']); checked = 0
    for start in range(0, M*M, 1024):
        indices = np.arange(start, min(start+1024, M*M)); ii, jj = indices//M, indices % M
        mixed = words[jj].copy(); mixed[:, prefix] = words[ii][:, prefix]
        projected = ((mixed @ weights) % p)*inv_beta % p
        actual = np.all(reencode(projected) == mixed, axis=1)
        expected = np.fromiter((flat[int(i)] == '1' for i in indices), dtype=bool)
        assert np.array_equal(actual, expected); checked += len(indices)
    return {'cross_splices': checked, 'scalar_encodings': M, 'overflow_bounds': limits,
            'calibration_beta': beta, 'distinguished_pairs': r['certified_distinguished_pairs'],
            'certified_width_lower_bound': r['certified_width_lower_bound']}


def main():
    source = HERE/'splice_summary.json'; summary = json.loads(source.read_text()); rows = []; samples = []
    bound = {source.name: sha256(source.read_bytes()).hexdigest()}
    for c in summary['cases']:
        path = HERE/c['artifact']; r = json.loads(path.read_text()); samples.append(r)
        bound[path.name] = sha256(path.read_bytes()).hexdigest()
        result = check_dense(r); rows.append({'name': r['name'], **result}); print(json.dumps(rows[-1]), flush=True)
    bad = []
    for name, kind in [('duplicate_prefix_index', 0), ('false_selector', 1), ('duplicate_read_coordinate', 2), ('corrupt_scalar_word', 3), ('false_acceptance_entry', 4)]:
        r = deepcopy(samples[-1])
        if kind == 0: r['clique_indices'][1] = r['clique_indices'][0]
        if kind == 1:
            i, j = r['clique_indices'][:2]
            rows2 = list(r['cross_acceptance_rows'][i]); rows2[j] = '1'; r['cross_acceptance_rows'][i] = ''.join(rows2)
            r['pair_witness_selectors'] = '0'+r['pair_witness_selectors'][1:]
        if kind == 2: r['order'][1] = r['order'][0]
        if kind == 3: r['words_in_original_coordinates'][0][0] += 1
        if kind == 4:
            # Preserve internal count consistency: arithmetic replay must fail.
            i, j = 0, 1; old = r['cross_acceptance_rows'][i][j]; chars = list(r['cross_acceptance_rows'][i])
            chars[j] = '0' if old == '1' else '1'; r['cross_acceptance_rows'][i] = ''.join(chars)
            r['accepted_cross_splices'] += -1 if old == '1' else 1
            r['clique_indices'] = [0]; r['certified_width_lower_bound'] = 1
            r['pair_witness_selectors'] = ''; r['certified_distinguished_pairs'] = 0
        try: check_dense(r)
        except (AssertionError, IndexError, ValueError): bad.append(name)
        else: raise AssertionError('corruption accepted: '+name)
    out = {'status': 'passed', 'scope': 'Every reported scalar encoding, cross-splice and distinguishing suffix selector replayed using independent direct powers and dense integer arithmetic. No maximum-clique or all-orders claim.',
           'cases': rows, 'corrupt_artifacts_rejected': bad,
           'source_sha256': {'splice_review.py': sha256(Path(__file__).read_bytes()).hexdigest()}, 'input_sha256': bound,
           'numpy_version': np.__version__}
    (HERE/'splice_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
