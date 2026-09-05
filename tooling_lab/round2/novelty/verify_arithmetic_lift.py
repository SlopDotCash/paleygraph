#!/usr/bin/env python3
"""Independent standard-library verifier; no discovery-module imports."""
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def mul(a, b, modulus):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) % modulus
             for j in range(len(b[0]))] for i in range(len(a))]


def full_rank_mod_p(a, p):
    a = [[x % p for x in row] for row in a]
    n = len(a)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return False
        a[j], a[pivot] = a[pivot], a[j]
        inv = pow(a[j][j], -1, p)
        for i in range(j+1, n):
            ratio = a[i][j]*inv % p
            a[i] = [(x-ratio*y) % p for x, y in zip(a[i], a[j])]
    return True


def coefficient_matrix(groups, q, modulus):
    cols = []
    for roots in groups:
        coeff = [1]
        for root in roots:
            coeff = [((coeff[k-1] if k else 0)-root*(coeff[k] if k < len(coeff) else 0)) % modulus
                     for k in range(len(coeff)+1)]
        for shift in range(q+1):
            cols.append([0]*shift+coeff+[0]*(q-shift))
    return [list(row) for row in zip(*cols)]


def main():
    path = HERE/'arithmetic_lift_results.json'
    data = json.loads(path.read_text())
    certificates = obstructions = roots_checked = 0
    for profile in data['profiles']:
        p, n, q = profile['p'], profile['n'], profile['cofactor_degree']
        base_roots = profile['precisions'][0]['root_groups']
        previous_roots = None
        previous_modulus = 1
        for case in profile['precisions']:
            e = case['precision']; modulus = p**e
            groups = case['root_groups']
            for gs, base in zip(groups, base_roots):
                for root, a in zip(gs, base):
                    assert root % p == a and pow(root, n, modulus) == 1
                    roots_checked += 1
            if previous_roots:
                assert all(x % previous_modulus == y for gs, hs in zip(groups, previous_roots)
                           for x, y in zip(gs, hs))
            previous_roots, previous_modulus = groups, modulus
            a = coefficient_matrix(groups, q, modulus)
            assert a == case['coefficient_matrix']
            u, v = case['left_transform'], case['right_transform']
            assert full_rank_mod_p(u, p) and full_rank_mod_p(v, p)
            diag = mul(mul(u, a, modulus), v, modulus)
            assert diag == case['diagonal']
            m, nc = len(a), len(a[0])
            assert all(diag[i][j] == 0 for i in range(m) for j in range(nc) if i != j)
            nonzero = [diag[i][i] for i in range(min(m, nc)) if diag[i][i]]
            assert nonzero == [p**v for v in case['pivot_valuations']]
            assert nc-len(nonzero) == case['surviving_initial_dimension']
            certificates += 1
        a = profile['precisions'][1]['coefficient_matrix']
        for obs in profile['first_lift_obstructions']:
            c, l = obs['initial_syzygy'], obs['left_annihilator']
            product = [sum(x*y for x, y in zip(row, c)) for row in a]
            assert all(x % p == 0 for x in product)
            assert all(sum(l[i]*a[i][j] for i in range(len(a))) % p == 0 for j in range(len(c)))
            carry = [x//p % p for x in product]
            assert carry == obs['carry']
            scalar = sum(x*y for x, y in zip(l, carry)) % p
            assert scalar == obs['nonzero_obstruction'] and scalar != 0
            obstructions += 1
    result = {'passed': True, 'modular_diagonal_certificates': certificates,
              'root_checks': roots_checked, 'exact_nonlifting_witnesses': obstructions,
              'results_sha256': sha256(path.read_bytes()).hexdigest(),
              'verifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    HERE.joinpath('arithmetic_lift_validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
