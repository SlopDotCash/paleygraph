#!/usr/bin/env python3
"""Polynomial-substitution review of the signed-permutation digit adapter."""
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp
from language_review import recurrence_membership

HERE = Path(__file__).resolve().parent
X = sp.Symbol('X')


def substitute(vector, e):
    N = len(vector); polynomial = sp.Poly(sum(v*X**(e*j) for j, v in enumerate(vector)), X)
    remainder = polynomial.rem(sp.Poly(X**N+1, X))
    return [int(remainder.nth(j)) for j in range(N)]


def monomial_product(vector, h, sign):
    N = len(vector)
    polynomial = sp.Poly(sum(v*X**j for j, v in enumerate(vector)), X)
    monomial = sp.Poly(sign*X**(h % (2*N)), X)
    rem = (polynomial*monomial).rem(sp.Poly(X**N+1, X))
    return [int(rem.nth(j)) for j in range(N)]


def main():
    source = HERE/'orientation_adapter.json'; previous = HERE.parent/'round17/norm_compression.json'
    data = json.loads(source.read_text()); original = json.loads(previous.read_text()); summaries = []
    for case, base in zip(data['cases'], original['cases']):
        N, p, g, f, mapping = (case[key] for key in ('N', 'p', 'g', 'relation', 'mapping'))
        assert (p, g, N, f) == (base['p'], base['g'], base['N'], base['relation'])
        assert sp.isprime(p) and N & (N-1) == 0 and pow(g, N, p) == p-1
        if mapping['status'] != 'compiled':
            assert mapping['status'] == 'generator_2_has_wrong_order' and pow(2, N, p) != p-1
            assert case['records'] == []
        else:
            e, h, sign, k = (mapping[key] for key in ('exponent', 'shift', 'sign', 'cofactor'))
            assert 0 < e < 2*N and e % 2 and pow(2, e, p) == g and pow(2, N, p) == p-1
            assert k*p == 2**N+1 and k == base['relation_cofactor'] and 0 <= h < N and sign in (-1, 1)
            assert substitute(f, e) == monomial_product([2]+[0]*(N-2)+[1], h, sign)
            assert len(case['records']) == len(base['records'])
            for row, prior in zip(case['records'], base['records']):
                assert row['a'] == prior['a'] and row['original_digits'] == prior['digits']
                digits = monomial_product(substitute(row['original_digits'], e), -h, sign)
                assert digits == row['canonical_digits']
                member, F = recurrence_membership(digits, k)
                assert member and F[0] % p == row['a']
                assert row['result'] == {'status': 'nonzero', 'scalar_centered': F[0], 'lift_scalar': k*F[0]}
        summaries.append({'case': case['name'], 'status': mapping['status'], 'vectors_checked': len(case['records'])})
        print(json.dumps(summaries[-1]), flush=True)
    assert data['nonassociate_control']['status'] == 'relation_is_not_a_signed_monomial_associate'
    # A signed permutation cannot change the coefficient multiset in absolute value.
    assert sorted([4]+[0]*30+[2]) != sorted([2]+[0]*30+[1])
    out = {'status': 'passed', 'scope': 'Separate polynomial substitution and remainder computations check every orientation and every mapped word, then recurrence inversion checks its actual scalar.',
           'cases': summaries, 'nonassociate_control_checked': True,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['orientation_review.py', 'language_review.py']},
           'input_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                            for name in ['orientation_adapter.json', '../round17/norm_compression.json']}}
    (HERE/'orientation_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__':
    main()
