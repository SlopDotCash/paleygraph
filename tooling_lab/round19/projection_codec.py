#!/usr/bin/env python3
"""Exact scalar projection and re-encoding for arbitrary kernel relations.

The fast branch calibrates E(1)(u). Degenerate calibration uses one exact
adjugate row modulo k*p. Neither branch constructs a dense inverse per word
or enumerates cofactor residues. Both require a final exact re-encoding.
"""
from hashlib import sha256
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from realizability import norm_adjugate


def evaluate(digits, root, p):
    value = 0
    for digit in reversed(digits): value = (root*value+digit) % p
    return value


def encode(p, g, f, a):
    N = len(f); F = []; residue = a % p
    for j in range(N):
        F.append(residue if 2*residue < p else residue-p); residue = residue*g % p
    numerator = [0]*N; products = 0
    for j, value in enumerate(f):
        if not value: continue
        for i, coefficient in enumerate(F):
            index = i+j; numerator[index % N] += value*coefficient*(1 if index < N else -1)
            products += 1
    assert all(v % p == 0 for v in numerator)
    return [v//p for v in numerator], products


def prepare(p, g, f):
    N = len(f)
    if not (type(p) is int and sp.isprime(p) and N >= 4 and N & (N-1) == 0 and
            type(g) is int and 0 < g < p and pow(g, N, p) == p-1 and
            all(type(v) is int for v in f) and any(f)):
        raise ValueError('invalid arithmetic input')
    u = pow(g, -1, p)
    if evaluate(f, u, p): raise ValueError('invalid kernel relation')
    one, _ = encode(p, g, f, 1); beta = evaluate(one, u, p)
    c = {'p': p, 'g': g, 'N': N, 'relation': f, 'u': u, 'calibration_digits': one,
         'beta': beta, 'digit_bound': p//2*sum(map(abs, f))//p,
         'relation_support': sum(v != 0 for v in f)}
    if beta:
        c.update({'mode': 'field_projection', 'inverse_beta': pow(beta, -1, p)})
    else:
        nf, adj = norm_adjugate(f); assert nf > 0 and nf % p == 0; k = nf//p
        weights = [adj[0]]+[-adj[N-j] for j in range(1, N)]
        c.update({'mode': 'single_row_fallback', 'relation_norm': nf, 'k': k,
                  'projection_modulus': k*p, 'projection_weights': [v % (k*p) for v in weights]})
    return c


def decode(c, digits):
    N, p = c['N'], c['p']
    if len(digits) != N or any(type(v) is not int for v in digits): return {'status': 'invalid_digits'}
    if any(abs(v) > c['digit_bound'] for v in digits): return {'status': 'digit_height'}
    if c['mode'] == 'field_projection':
        a = evaluate(digits, c['u'], p)*c['inverse_beta'] % p
    else:
        remainder = sum(w*d for w, d in zip(c['projection_weights'], digits)) % c['projection_modulus']
        if remainder % c['k']: return {'status': 'projection_nonintegral'}
        a = remainder//c['k']
    expected, products = encode(p, c['g'], c['relation'], a)
    return {'status': ('accepted_nonzero' if a else 'accepted_zero') if expected == list(digits)
                      else 'reencoding_mismatch', 'extracted_scalar': a,
            'coefficient_products': products, 'projection_terms': N,
            'reconstructed_coefficients': N}


def main():
    source = HERE.parent/'round17/realizability.json'; old = json.loads(source.read_text())
    cases = []; total = 0
    for oldcase in old['cases']:
        base = oldcase['certificate']; c = prepare(base['p'], base['g'], base['f'])
        groups = {}
        for collection in ('positives', 'probes'):
            results = []
            for r in oldcase[collection]:
                result = decode(c, r['digits'])
                assert result['status'].startswith('accepted_') == (r['status'] in ('zero', 'nonzero'))
                if result['status'].startswith('accepted_'): assert result['extracted_scalar'] == r['a']
                results.append({'digits': r['digits'], 'result': result}); total += 1
            groups[collection] = results
        cases.append({'name': oldcase['name'], 'certificate': c, **groups})
        print(json.dumps({'case': oldcase['name'], 'mode': c['mode'], 'support': c['relation_support'],
                          'positive_records': len(groups['positives']), 'probes': len(groups['probes'])}), flush=True)
    small = [{'digits': r['digits'], 'result': decode(cases[0]['certificate'], r['digits'])}
             for r in old['small_census']]
    for r, prior in zip(small, old['small_census']):
        assert r['result']['status'].startswith('accepted_') == (prior['status'] in ('zero', 'nonzero'))

    barriers = HERE/'relation_barrier.json'; adversaries = json.loads(barriers.read_text()); hard = []
    for r in adversaries['cases']:
        c = prepare(r['p'], r['g'], r['relation'])
        yes = decode(c, r['centered_digits']); no = decode(c, r['escaped_digits'])
        assert yes['status'] == 'accepted_nonzero' and no['status'] == 'reencoding_mismatch'
        assert yes['extracted_scalar'] == no['extracted_scalar'] == r['a']
        hard.append({'name': r['name'], 'certificate': c, 'positive': yes, 'negative': no})

    # Test p-dividing cofactors with a good calibration, then two kinds of
    # zero calibration. The fallback must retain exact completeness.
    p, g, f = old['cases'][0]['certificate']['p'], old['cases'][0]['certificate']['g'], old['cases'][0]['certificate']['f']
    N = len(f); squared = [0]*N
    for i, x in enumerate(f):
        for j, y in enumerate(f): squared[(i+j) % N] += x*y*(1 if i+j < N else -1)
    controls = []
    for name, relation in [('scalar_p', [p]+[0]*(N-1)), ('squared_relation', squared), ('scalar_p_squared', [p*p]+[0]*(N-1))]:
        c = prepare(p, g, relation); records = []
        for a in [*range(16), (p+1)//2]:
            D, _ = encode(p, g, relation, a); result = decode(c, D)
            assert result['status'].startswith('accepted_') and result['extracted_scalar'] == a
            records.append({'digits': D, 'result': result, 'intended_scalar': a})
            for j in (0, N//2, N-1):
                changed = D.copy(); changed[j] += 1
                records.append({'digits': changed, 'result': decode(c, changed)})
        controls.append({'name': name, 'certificate': c, 'records': records})
    assert [c['certificate']['mode'] for c in controls] == ['field_projection', 'single_row_fallback', 'single_row_fallback']
    malformed = [[], [True]+[0]*(N-1), [0.5]+[0]*(N-1)]
    invalid = [{'digits': d, 'result': decode(cases[0]['certificate'], d)} for d in malformed]
    assert all(r['result']['status'] == 'invalid_digits' for r in invalid)
    out = {'status': 'produced', 'scope': 'Exact scalar extraction plus full sparse re-encoding, with one-row fallback when the field calibration vanishes. No dense inverse per word or residue-state enumeration; this is a membership algorithm, not a norm estimate.',
           'cases': cases, 'small_census': small, 'barrier_witnesses': hard, 'degeneracy_controls': controls,
           'malformed_digit_controls': invalid,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['projection_codec.py', '../round17/realizability.py', '../round17/norm_carry.py']},
           'input_sha256': {'../round17/realizability.json': sha256(source.read_bytes()).hexdigest(),
                            'relation_barrier.json': sha256(barriers.read_bytes()).hexdigest()}}
    (HERE/'projection_codec.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print(json.dumps({'prior_records': total, 'small_cube_words': len(small), 'barrier_pairs': len(hard),
                      'degeneracy_records': sum(len(c['records']) for c in controls)}), flush=True)


if __name__ == '__main__': main()
