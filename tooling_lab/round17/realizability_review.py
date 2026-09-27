#!/usr/bin/env python3
"""Separate rational-polynomial inverse / integer-matrix census review."""
from collections import Counter, defaultdict
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as sp
from carry_review import result_norm

HERE = Path(__file__).resolve().parent
X = sp.Symbol('X')


def multiplication_matrix(coefficients):
    N = len(coefficients)
    return [[coefficients[i-j] if i >= j else -coefficients[N+i-j]
             for j in range(N)] for i in range(N)]


def act(matrix, vector):
    return [sum(v*w for v, w in zip(row, vector)) for row in matrix]


def check_certificate(c):
    p, N, g, f = (c[k] for k in ('p', 'N', 'g', 'f'))
    assert type(p) is int and sp.isprime(p) and N >= 4 and N & (N-1) == 0
    assert len(f) == N and all(type(v) is int for v in f) and any(f)
    assert 0 < g < p and pow(g, N, p) == p-1 and (p-1) % (2*N) == 0
    assert sum(v*pow(pow(g, -1, p), j, p) for j, v in enumerate(f)) % p == 0
    nf = result_norm(f, N)
    assert nf == c['norm_f'] == p*c['k'] and nf > 0
    assert c['bound'] == (p-1)*sum(abs(v) for v in f)//(2*p)
    fp = sp.Poly.from_list(f[::-1], X, domain=sp.QQ)
    inverse = sp.invert(fp, sp.Poly(X**N+1, X, domain=sp.QQ))
    adj = [inverse.nth(j)*nf for j in range(N)]
    assert all(v.q == 1 for v in adj) and list(map(int, adj)) == c['adj']
    mf, ma = multiplication_matrix(f), multiplication_matrix(list(map(int, adj)))
    assert act(mf, c['adj']) == [nf]+[0]*(N-1)
    return mf, ma


def check_row(c, matrices, row):
    mf, ma = matrices
    p, k, N, digits = c['p'], c['k'], c['N'], row['digits']
    if len(digits) != N or any(type(v) is not int for v in digits):
        expected = {'status': 'invalid_digits'}
    elif any(abs(v) > c['bound'] for v in digits):
        expected = {'status': 'digit_height'}
    else:
        numerator = act(ma, digits)
        if any(v % k for v in numerator):
            expected = {'status': 'nonintegral_inverse'}
        else:
            F = [v//k for v in numerator]; a = F[0] % p
            canonical = []
            for j in range(N):
                v = a*pow(c['g'], j, p) % p
                canonical.append(v if 2*v < p else v-p)
            centered = all(2*abs(v) < p for v in F)
            scalar = all((v-w) % p == 0 for v, w in zip(F, canonical))
            assert act(mf, F) == [p*v for v in digits]
            status = ('scalar_mismatch' if not scalar else 'outside_centered' if not centered
                      else 'zero' if a == 0 else 'nonzero')
            if status in ('zero', 'nonzero'):
                assert F == canonical
            expected = {'status': status, 'F': F, 'a': a,
                        'centered': centered, 'scalar_congruences': scalar}
    for key, value in expected.items():
        assert row[key] == value, (key, row.get(key), value)
    assert set(row) == {'digits', *expected} or set(row) == {
        'digits', 'digit_norm', 'digit_energy', *expected}
    return expected['status'] in ('zero', 'nonzero')


def main():
    source = HERE/'realizability.json'; prior = HERE/'norm_compression.json'
    data = json.loads(source.read_text()); original = json.loads(prior.read_text())
    assert data['input_sha256'][prior.name] == sha256(prior.read_bytes()).hexdigest()
    reports = []; prepared = []
    for case, base in zip(data['cases'], original['cases']):
        c = case['certificate']; matrices = check_certificate(c); prepared.append(matrices)
        assert (c['p'], c['g'], c['f']) == (base['p'], base['g'], base['relation'])
        assert len(case['positives']) == len(base['records'])
        for row, prior_row in zip(case['positives'], base['records']):
            assert row['digits'] == prior_row['digits'] and row['a'] == prior_row['a']
            assert check_row(c, matrices, row) and row['status'] == 'nonzero'
        for row in case['probes']:
            check_row(c, matrices, row)
        counts = dict(sorted(Counter(r['status'] for r in case['probes']).items()))
        assert counts == case['probe_status_counts']
        reports.append({'case': case['name'], 'positive_vectors_checked': len(case['positives']),
                        'probes_checked': len(case['probes']), 'probe_status_counts': counts})
        print(json.dumps(reports[-1]), flush=True)

    c = data['cases'][0]['certificate']; matrices = prepared[0]
    rows = data['small_census']; assert len(rows) == 3**c['N']
    by_norm = defaultdict(Counter); by_pair = defaultdict(Counter); accepted = []
    for i, (digits, row) in enumerate(zip(product((-1, 0, 1), repeat=c['N']), rows)):
        assert list(digits) == row['digits']
        good = check_row(c, matrices, row)
        nd = result_norm(digits, c['N']); energy = sum(v*v for v in digits)
        assert (nd, energy) == (row['digit_norm'], row['digit_energy'])
        label = 'realized' if good else 'nonrealized'
        by_norm[nd][label] += 1; by_pair[nd, energy][label] += 1
        if good: accepted.append(row)
        if (i+1) % 2000 == 0:
            print(f'census resultants checked {i+1}/{len(rows)}', flush=True)
    # An independent direct scalar enumeration must produce exactly the accepted digit words.
    scalar_words = {}
    for a in range(c['p']):
        F = [((a*pow(c['g'], j, c['p'])+c['p']//2) % c['p'])-c['p']//2
             for j in range(c['N'])]
        raw = act(matrices[0], F); assert all(v % c['p'] == 0 for v in raw)
        scalar_words[a] = [v//c['p'] for v in raw]
    assert {r['a']: r['digits'] for r in accepted} == scalar_words
    summary = data['small_summary']
    assert summary['cube_words'] == len(rows) and summary['actual_words_including_zero'] == len(accepted)
    assert summary['status_counts'] == dict(sorted(Counter(r['status'] for r in rows).items()))
    norms = {r['digit_norm'] for r in accepted}; pairs = {(r['digit_norm'], r['digit_energy']) for r in accepted}
    assert summary['actual_norm_values'] == sorted(norms)
    assert summary['norm_divisible_by_k_accepts'] == sum(r['digit_norm'] % c['k'] == 0 for r in rows)
    assert summary['oracle_actual_norm_set_accepts'] == sum(r['digit_norm'] in norms for r in rows)
    assert summary['oracle_actual_norm_energy_pair_set_accepts'] == sum((r['digit_norm'], r['digit_energy']) in pairs for r in rows)
    assert summary['by_norm'] == [{'norm': k, **dict(v)} for k, v in sorted(by_norm.items())]
    assert summary['by_norm_energy'] == [{'norm': k[0], 'energy': k[1], **dict(v)} for k, v in sorted(by_pair.items())]
    witness = summary['same_norm_and_energy_witness']
    assert check_row(c, matrices, witness['realized'])
    assert not check_row(c, matrices, witness['nonrealized'])
    assert (witness['realized']['digit_norm'], witness['realized']['digit_energy']) == (
        witness['nonrealized']['digit_norm'], witness['nonrealized']['digit_energy'])
    general = data['p_divides_cofactor_control']; gc = general['certificate']; gm = check_certificate(gc)
    assert gc['k'] % gc['p'] == 0
    for row in general['records']: check_row(gc, gm, row)
    assert [r['status'] for r in general['records']] == ['scalar_mismatch', 'zero', 'nonzero']

    # Deliberately corrupt each major certificate/output field; the verifier must reject.
    controls = []
    for field in ('norm_f', 'k', 'adj', 'f', 'bound'):
        corrupted = deepcopy(c)
        if field in ('adj', 'f'): corrupted[field][0] += 1
        else: corrupted[field] += 1
        try: check_certificate(corrupted)
        except AssertionError: controls.append('certificate_'+field)
        else: raise AssertionError('corrupt certificate accepted')
    for field in ('a', 'F', 'status', 'centered', 'scalar_congruences'):
        corrupted = deepcopy(data['cases'][0]['positives'][0])
        if field == 'F': corrupted[field][0] += 1
        elif field == 'status': corrupted[field] = 'zero'
        elif field in ('centered', 'scalar_congruences'): corrupted[field] = False
        else: corrupted[field] += 1
        try: check_row(c, matrices, corrupted)
        except AssertionError: controls.append('output_'+field)
        else: raise AssertionError('corrupt output accepted')
    out = {'status': 'passed', 'scope': 'Rational polynomial Euclidean inverses and integer multiplication matrices, without importing the realizability producer. Every small digit norm independently checked by resultant; accepted set equals the full direct scalar census. Root-run separate implementation, not a human or Lean review.',
           'cases': reports, 'census_words_checked': len(rows), 'census_resultants': len(rows),
           'actual_scalar_words_checked': len(scalar_words), 'small_summary': summary,
           'p_divides_cofactor_controls_checked': len(general['records']),
           'corrupt_artifacts_rejected': controls,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['realizability_review.py', 'carry_review.py']},
           'input_sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in [source, prior]}}
    (HERE/'realizability_review.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'census': len(rows), 'realized': len(accepted),
                      'corruptions_rejected': len(controls)}), flush=True)


if __name__ == '__main__':
    main()
