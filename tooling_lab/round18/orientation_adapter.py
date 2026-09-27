#!/usr/bin/env python3
"""Recognize sparse encodings after exact cyclotomic signed permutations."""
from hashlib import sha256
import json
from pathlib import Path
from digit_language import recognize, scalar_digits

HERE = Path(__file__).resolve().parent


def automorphism(vector, exponent):
    N = len(vector); result = [0]*N
    for j, v in enumerate(vector):
        where = j*exponent % (2*N)
        result[where % N] += v if where < N else -v
    return result


def shift(vector, exponent, sign):
    N = len(vector); result = [0]*N
    for j, v in enumerate(vector):
        where = (j+exponent) % (2*N)
        result[where % N] += sign*v if where < N else -sign*v
    return result


def compile_orientation(p, g, f):
    N = len(f)
    if pow(2, N, p) != p-1:
        return {'status': 'generator_2_has_wrong_order'}
    choices = [e for e in range(1, 2*N, 2) if pow(2, e, p) == g]
    if len(choices) != 1:
        return {'status': 'no_unique_generator_exponent'}
    e = choices[0]; transformed = automorphism(f, e); canonical = [2]+[0]*(N-2)+[1]
    for h in range(N):
        for sign in (-1, 1):
            if transformed == shift(canonical, h, sign):
                return {'status': 'compiled', 'exponent': e, 'shift': h, 'sign': sign,
                        'cofactor': (2**N+1)//p}
    return {'status': 'relation_is_not_a_signed_monomial_associate'}


def main():
    source = HERE.parent/'round17/norm_compression.json'; old = json.loads(source.read_text()); cases = []
    for case in old['cases']:
        p, g, N = case['p'], case['g'], case['N']
        mapping = compile_orientation(p, g, case['relation']); rows = []
        if mapping['status'] == 'compiled':
            assert mapping['cofactor'] == case['relation_cofactor']
            for row in case['records']:
                digits = shift(automorphism(row['digits'], mapping['exponent']),
                               -mapping['shift'], mapping['sign'])
                result = recognize(digits, mapping['cofactor'])
                assert result['status'] == 'nonzero' and result['scalar_centered'] % p == row['a']
                F, expected = scalar_digits(row['a'], p, N); assert digits == expected
                rows.append({'a': row['a'], 'original_digits': row['digits'],
                             'canonical_digits': digits, 'result': result})
        cases.append({'name': case['name'], 'p': p, 'g': g, 'N': N,
                      'relation': case['relation'], 'mapping': mapping, 'records': rows})
        print(json.dumps({'case': case['name'], 'mapping': mapping, 'records': len(rows)}), flush=True)
    invalid = compile_orientation(6700417, 2, [4]+[0]*30+[2])
    assert invalid['status'] == 'relation_is_not_a_signed_monomial_associate'
    out = {'status': 'produced', 'scope': 'Signed-permutation adapter for exactly those saved relations that become a monomial associate of 2-X^-1 under X->X^e. Unsupported inputs are explicit; no automatic extension to an arbitrary short relation.',
           'cases': cases, 'nonassociate_control': invalid,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['orientation_adapter.py', 'digit_language.py']},
           'input_sha256': {'../round17/norm_compression.json': sha256(source.read_bytes()).hexdigest()}}
    (HERE/'orientation_adapter.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
