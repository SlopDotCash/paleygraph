#!/usr/bin/env python3
"""Independent full-inverse oracle for the projection/re-encoding codec."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from carry_review import product, result_norm
from realizability_review import multiplication_matrix, act
X = sp.Symbol('X')


def canonical(a, p, N, g):
    values = [a*pow(g, j, p) % p for j in range(N)]
    return [v if v <= (p-1)//2 else v-p for v in values]


def certificate(c):
    p, N, g, f = (c[key] for key in ('p', 'N', 'g', 'relation'))
    assert type(p) is int and sp.isprime(p) and N >= 4 and N & (N-1) == 0
    assert len(f) == N and all(type(v) is int for v in f) and any(f)
    assert 0 < g < p and pow(g, N, p) == p-1 and c['u'] == pow(g, -1, p)
    u = c['u']; assert sum(v*pow(u, j, p) for j, v in enumerate(f)) % p == 0
    raw = product(f, canonical(1, p, N, g), N); assert all(v % p == 0 for v in raw)
    one = [v//p for v in raw]; assert one == c['calibration_digits']
    beta = sum(v*pow(u, j, p) for j, v in enumerate(one)) % p
    assert beta == c['beta']
    assert c['relation_support'] == sum(v != 0 for v in f)
    assert c['digit_bound'] == (p-1)*sum(map(abs, f))//(2*p)
    nf = result_norm(f, N); assert nf > 0 and nf % p == 0; k = nf//p
    inverse = sp.invert(sp.Poly.from_list(f[::-1], X, domain=sp.QQ), sp.Poly(X**N+1, X, domain=sp.QQ))
    adj = [inverse.nth(j)*nf for j in range(N)]
    assert all(v.q == 1 for v in adj); adj = list(map(int, adj))
    if beta:
        assert c['mode'] == 'field_projection' and 0 < c['inverse_beta'] < p and beta*c['inverse_beta'] % p == 1
        assert 'projection_weights' not in c and 'k' not in c
    else:
        assert c['mode'] == 'single_row_fallback' and c['relation_norm'] == nf and c['k'] == k
        assert c['projection_modulus'] == nf
        expected = [adj[0]]+[-adj[N-j] for j in range(1, N)]
        assert c['projection_weights'] == [v % nf for v in expected]
    return {'adj_matrix': multiplication_matrix(adj), 'f_matrix': multiplication_matrix(f), 'k': k}


def check(c, matrices, record):
    digits, result = record['digits'], record['result']; p, N, k = c['p'], c['N'], matrices['k']
    if len(digits) != N or any(type(v) is not int for v in digits):
        assert result == {'status': 'invalid_digits'}; return False
    numerator = act(matrices['adj_matrix'], digits)
    oracle = False
    if all(v % k == 0 for v in numerator):
        F = [v//k for v in numerator]
        oracle = F == canonical(F[0] % p, p, N, c['g'])
    if any(abs(v) > c['digit_bound'] for v in digits):
        assert result == {'status': 'digit_height'} and not oracle; return False
    if c['beta']:
        root_value = sum(v*pow(c['u'], j, p) for j, v in enumerate(digits)) % p
        a = root_value*pow(c['beta'], -1, p) % p
    else:
        # The full inverse oracle already has this numerator; no producer row-dot used.
        residue = numerator[0] % (k*p)
        if residue % k:
            assert result == {'status': 'projection_nonintegral'} and not oracle; return False
        a = residue//k
    raw = act(matrices['f_matrix'], canonical(a, p, N, c['g']))
    assert all(v % p == 0 for v in raw)
    reen = [v//p for v in raw]; accepted = reen == digits
    assert accepted == oracle
    status = ('accepted_nonzero' if a else 'accepted_zero') if accepted else 'reencoding_mismatch'
    assert result == {'status': status, 'extracted_scalar': a,
                      'coefficient_products': N*sum(v != 0 for v in c['relation']),
                      'projection_terms': N, 'reconstructed_coefficients': N}
    if 'intended_scalar' in record:
        assert accepted and a == record['intended_scalar']
    return accepted


def main():
    source = HERE/'projection_codec.json'; old_source = HERE.parent/'round17/realizability.json'; barrier_source = HERE/'relation_barrier.json'
    data = json.loads(source.read_text()); old = json.loads(old_source.read_text()); barriers = json.loads(barrier_source.read_text())
    cases = []; prepared = []
    assert len(data['cases']) == len(old['cases']) == 5
    for c, prior in zip(data['cases'], old['cases']):
        cert = c['certificate']; matrices = certificate(cert); prepared.append(matrices)
        assert (cert['p'], cert['g'], cert['relation']) == tuple(prior['certificate'][key] for key in ('p', 'g', 'f'))
        counts = {}
        for group in ('positives', 'probes'):
            assert len(c[group]) == len(prior[group]); accepted = 0
            for row, oldrow in zip(c[group], prior[group]):
                assert row['digits'] == oldrow['digits']; good = check(cert, matrices, row)
                assert good == (oldrow['status'] in ('zero', 'nonzero'))
                if good: assert row['result']['extracted_scalar'] == oldrow['a']
                accepted += good
            counts[group] = {'records': len(c[group]), 'accepted': accepted}
        cases.append({'name': c['name'], 'mode': cert['mode'], 'relation_support': cert['relation_support'],
                      'coefficient_products_per_reencoding': cert['N']*cert['relation_support'], 'groups': counts})
        print(json.dumps(cases[-1]), flush=True)
    assert len(data['small_census']) == len(old['small_census']) == 6561
    accepted = 0
    for row, oldrow in zip(data['small_census'], old['small_census']):
        assert row['digits'] == oldrow['digits']; good = check(data['cases'][0]['certificate'], prepared[0], row)
        assert good == (oldrow['status'] in ('zero', 'nonzero')); accepted += good
    assert accepted == 257
    assert len(data['barrier_witnesses']) == len(barriers['cases']) == 5
    for output, prior in zip(data['barrier_witnesses'], barriers['cases']):
        c = output['certificate']; matrices = certificate(c)
        assert (c['p'], c['g'], c['relation']) == (prior['p'], prior['g'], prior['relation'])
        for prefix, digits_key, good in [('positive', 'centered_digits', True), ('negative', 'escaped_digits', False)]:
            assert check(c, matrices, {'digits': prior[digits_key], 'result': output[prefix]}) == good
            assert output[prefix]['extracted_scalar'] == prior['a']
    degeneracy = []
    for control in data['degeneracy_controls']:
        c = control['certificate']; matrices = certificate(c); count = 0
        for row in control['records']: count += check(c, matrices, row)
        degeneracy.append({'name': control['name'], 'mode': c['mode'], 'records_checked': len(control['records']),
                           'accepted': count, 'p_divides_norm_cofactor': matrices['k'] % c['p'] == 0})
    assert [r['mode'] for r in degeneracy] == ['field_projection', 'single_row_fallback', 'single_row_fallback']
    assert all(r['p_divides_norm_cofactor'] for r in degeneracy)
    for row in data['malformed_digit_controls']: check(data['cases'][0]['certificate'], prepared[0], row)
    rejected = []
    base = data['cases'][-1]['certificate']
    for key in ('beta', 'inverse_beta', 'calibration_digits', 'relation', 'digit_bound', 'relation_support'):
        bad = deepcopy(base)
        if key in ('calibration_digits', 'relation'): bad[key][0] += 1
        else: bad[key] += 1
        try: certificate(bad)
        except AssertionError: rejected.append('certificate_'+key)
        else: raise AssertionError('corrupt certificate accepted')
    fallback = data['degeneracy_controls'][1]['certificate']
    for key in ('projection_weights', 'projection_modulus', 'k'):
        bad = deepcopy(fallback)
        if key == 'projection_weights': bad[key][0] += 1
        else: bad[key] += 1
        try: certificate(bad)
        except AssertionError: rejected.append('fallback_'+key)
        else: raise AssertionError('corrupt fallback accepted')
    c = data['cases'][0]['certificate']; matrices = prepared[0]
    for key in ('extracted_scalar', 'coefficient_products', 'status'):
        bad = deepcopy(data['cases'][0]['positives'][0])
        if key == 'status': bad['result'][key] = 'reencoding_mismatch'
        else: bad['result'][key] += 1
        try: check(c, matrices, bad)
        except AssertionError: rejected.append('output_'+key)
        else: raise AssertionError('corrupt output accepted')
    out = {'status': 'passed', 'scope': 'Separate rational polynomial inversion and full multiplication matrices check every membership result; calibration and one-row fallback are independently reconstructed. No producer import. Review is deliberately more expensive than the runtime codec.',
           'cases': cases, 'small_cube_words_checked': 6561, 'small_actual_words': accepted,
           'barrier_pairs_checked': 5, 'degeneracy_controls': degeneracy,
           'malformed_controls_checked': len(data['malformed_digit_controls']), 'corrupt_artifacts_rejected': rejected,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['codec_review.py', '../round17/carry_review.py', '../round17/realizability_review.py']},
           'input_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                            for name in ['projection_codec.json', '../round17/realizability.json', 'relation_barrier.json']}}
    (HERE/'codec_review.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'small_census': 6561, 'accepted': accepted,
                      'degeneracy_records': sum(r['records_checked'] for r in degeneracy),
                      'corruptions_rejected': len(rejected)}), flush=True)


if __name__ == '__main__': main()
