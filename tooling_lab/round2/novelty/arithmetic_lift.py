#!/usr/bin/env python3
"""Domain-preserving arithmetic syzygy lifting; known Hensel/Smith machinery.

New experiment: distinguish lift conventions and classify existing syzygies by
survival modulo p^e while roots remain n-th roots of unity at every precision.
"""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]
sys.path.insert(0, str(LAB / 'proximity'))
from deformation_microscope import F41_GROUPS, domain, coefficient_matrix, kernel, transpose


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def multiply(a, b, modulus):
    return [[sum(x*y for x, y in zip(row, col)) % modulus for col in zip(*b)] for row in a]


def valuation(x, p, e):
    if x == 0:
        return e
    v = 0
    while v < e and x % p == 0:
        x //= p
        v += 1
    return v


def local_diagonal(a, p, e):
    """Return unimodular U,V with U A V diagonal over Z/p^e.

    Minimum-valuation pivoting ensures every remaining entry is divisible
    by the pivot power of p; only units are inverted.
    """
    modulus = p**e
    a = [[x % modulus for x in row] for row in a]
    original = [row[:] for row in a]
    m, n = len(a), len(a[0])
    u, v = identity(m), identity(n)
    pivots = []
    for k in range(min(m, n)):
        val, i, j = min((valuation(a[i][j], p, e), i, j)
                        for i in range(k, m) for j in range(k, n))
        if val == e:
            break
        a[k], a[i] = a[i], a[k]
        u[k], u[i] = u[i], u[k]
        for mat in (a, v):
            for row in mat:
                row[k], row[j] = row[j], row[k]
        power = p**val
        inv = pow(a[k][k] // power, -1, modulus)
        a[k] = [x*inv % modulus for x in a[k]]
        u[k] = [x*inv % modulus for x in u[k]]
        assert a[k][k] == power
        for i in range(k+1, m):
            factor = a[i][k] // power
            a[i] = [(x-factor*y) % modulus for x, y in zip(a[i], a[k])]
            u[i] = [(x-factor*y) % modulus for x, y in zip(u[i], u[k])]
        for j in range(k+1, n):
            factor = a[k][j] // power
            for mat in (a, v):
                for row in mat:
                    row[j] = (row[j]-factor*row[k]) % modulus
        pivots.append(val)
    assert all(a[i][j] == 0 for i in range(m) for j in range(n) if i != j)
    assert multiply(multiply(u, original, modulus), v, modulus) == a
    # Every elementary operation is invertible; check modulo p independently.
    assert not kernel(u, p) and not kernel(v, p)
    survival = n-len(pivots)
    return {'pivot_valuations': pivots, 'precision': e, 'surviving_initial_dimension': survival,
            'diagonal': a, 'left_transform': u, 'right_transform': v}


def root_lift(a, n, p, e):
    """Unique root of X^n-1 modulo p^e reducing to a, for p∤n."""
    assert n % p and pow(a, n, p) == 1
    x, modulus = a % p, p
    for _ in range(1, e):
        next_modulus = modulus*p
        residue = (pow(x, n, next_modulus)-1) % next_modulus
        assert residue % modulus == 0
        digit = -(residue//modulus)*pow(n*pow(x, n-1, p) % p, -1, p) % p
        x += digit*modulus
        modulus = next_modulus
        assert pow(x, n, modulus) == 1
    return x


def first_obstruction(matrix, p):
    """Exact mod-p² obstruction L*(A*c/p), with c a fixed residue lift."""
    right = kernel(matrix, p)
    left = kernel(transpose(matrix), p)
    witnesses = []
    for c in right:
        product_c = [sum(a*b for a, b in zip(row, c)) for row in matrix]
        assert all(x % p == 0 for x in product_c)
        carry = [x//p % p for x in product_c]
        for l in left:
            obstruction = sum(a*b for a, b in zip(l, carry)) % p
            if obstruction:
                witnesses.append({'initial_syzygy': c, 'left_annihilator': l,
                                  'carry': carry, 'nonzero_obstruction': obstruction})
                break
    return witnesses


def profile(name, groups, n, p, q=0, max_e=4):
    modp_matrix = coefficient_matrix(groups, q, p)
    precisions = []
    for e in range(1, max_e+1):
        gs = [[root_lift(a, n, p, e) for a in group] for group in groups]
        matrix = coefficient_matrix(gs, q, p**e)
        cert = local_diagonal(matrix, p, e)
        precisions.append({'root_groups': gs, 'coefficient_matrix': matrix, **cert})
    survival = [r['surviving_initial_dimension'] for r in precisions]
    assert all(a >= b for a, b in zip(survival, survival[1:]))
    naive_matrix = coefficient_matrix(groups, q, p*p)
    naive = local_diagonal(naive_matrix, p, 2)
    return {'name': name, 'p': p, 'n': n, 'cofactor_degree': q,
            'initial_kernel_dimension': len(kernel(modp_matrix, p)),
            'hensel_survival': survival, 'precisions': precisions,
            'first_lift_obstructions': first_obstruction(precisions[1]['coefficient_matrix'], p),
            'naive_residue_survival_mod_p2': naive['surviving_initial_dimension'],
            'naive_roots_on_domain_mod_p2': all(pow(a, n, p*p) == 1 for group in groups for a in group)}


def brute_controls():
    checked = 0
    # Exhaust all 2x2 matrices modulo 4 and compare the projected kernel.
    for vals in product(range(4), repeat=4):
        a = [list(vals[:2]), list(vals[2:])]
        leads = set()
        for x in product(range(4), repeat=2):
            if all(sum(v*w for v, w in zip(row, x)) % 4 == 0 for row in a):
                leads.add(tuple(v % 2 for v in x))
        cert = local_diagonal(a, 2, 2)
        assert len(leads) == 2**cert['surviving_initial_dimension']
        checked += 1
    return checked


def main():
    controls = brute_controls()
    profiles = [profile('F41_linear_middle', F41_GROUPS, 20, 41, q=1)]
    for n, p in ((16, 17), (16, 1009), (24, 97), (48, 193)):
        dom = domain(n, p)
        profiles.append(profile(f'coset_n{n}_p{p}', [dom[i::4] for i in range(3)], n, p))
    antipodal = []
    for p, n in ((17, 16), (41, 20), (97, 24)):
        lifted_minus_one = root_lift(p-1, n, p, 3)
        assert (1+lifted_minus_one) % p**3 == 0
        antipodal.append({'p': p, 'n': n, 'precision': 3,
                          'naive_integer_sum': p, 'domain_preserving_sum_mod_p3': 0,
                          'lifted_minus_one': lifted_minus_one})
    result = {'status': 'finite domain-preserving arithmetic diagnostic; no prize or archimedean bound',
              'known_machinery': 'Hensel lifting, local Smith diagonalization, cokernel lift obstruction',
              'brute_projected_kernel_checks': controls, 'profiles': profiles,
              'lift_convention_controls': antipodal,
              'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    HERE.joinpath('arithmetic_lift_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': controls,
                      'profiles': [{'name': r['name'], 'hensel': r['hensel_survival'],
                                    'naive_p2': r['naive_residue_survival_mod_p2'],
                                    'obstructions': len(r['first_lift_obstructions'])} for r in profiles]}, indent=2))


if __name__ == '__main__':
    main()
