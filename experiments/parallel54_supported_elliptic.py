#!/usr/bin/env python3
"""Exact C-supported witness and Mellin projection/reflection identities.

Cyclotomic polynomial remainders certify character formulas without floating
point arithmetic. Finite tests do not establish the asymptotic spectral edge.
"""
from fractions import Fraction
from pathlib import Path
import json

from parallel53_ambient_elliptic import prime, primitive_root

ROOT = Path(__file__).resolve().parents[1]


def trim(a):
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def divide(a, b):
    a = a[:]
    quotient = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and a != [0]:
        shift = len(a) - len(b)
        assert b[-1] == 1
        lead = a[-1]
        quotient[shift] = lead
        for j, coefficient in enumerate(b):
            a[shift + j] -= lead * coefficient
        trim(a)
    return trim(quotient), trim(a)


def cyclotomic(n):
    polys = {}
    for d in range(1, n + 1):
        if n % d:
            continue
        f = [-1] + [0] * (d - 1) + [1]
        for e, h in polys.items():
            if d % e == 0:
                f, remainder = divide(f, h)
                assert remainder == [0]
        polys[d] = f
    return polys[n]


def character_checks(p, g, chi, C):
    N = p - 1
    q = N // 2
    log = {pow(g, j, p): j for j in range(N)}
    phi = cyclotomic(N)

    def is_zero(coefficients):
        return divide(trim(coefficients), phi)[1] == [0]

    # For gamma = psi restricted to Q:
    # 4 sum_C gamma = 2q delta_gamma,1 + J(psi,chi)
    #                                      + J(psi chi,chi) - 2.
    for a in range(q):
        coefficients = [0] * N
        for x in C:
            coefficients[a * log[x] % N] += 4
        coefficients[0] += 2 - (2 * q if a == 0 else 0)
        for x in range(1, p):
            if x == 1:
                continue
            coefficients[a * log[x] % N] -= chi[(1 - x) % p]
            coefficients[(a + q) * log[x] % N] -= chi[(1 - x) % p]
        assert is_zero(coefficients), ('projection', p, a)

    # 4 sum_C conjugate(alpha(x)) beta(1-x) equals the four
    # complete Jacobi sums from the two extensions of each character.
    for a in range(q):
        for b in range(q):
            coefficients = [0] * N
            for x in C:
                exponent = (-a * log[x] + b * log[(1 - x) % p]) % N
                coefficients[exponent] += 4
            for x in range(2, p):
                lx, ly = log[x], log[(1 - x) % p]
                for s in (0, 1):
                    for t in (0, 1):
                        exponent = ((-a + s * q) * lx + (b + t * q) * ly) % N
                        coefficients[exponent] -= 1
            assert is_zero(coefficients), ('reflection', p, a, b)
    return dict(projection_character_identities=q,
                reflection_character_identities=q * q,
                cyclotomic_degree=len(phi) - 1)


def form(matrix, vector):
    return sum(vector[i] * matrix[i][j] * vector[j]
               for i in range(len(vector)) for j in range(len(vector)))


def verify(p, characters=False):
    assert prime(p) and p % 4 == 1
    g = primitive_root(p)
    q = (p - 1) // 2
    chi = [0] + [1 if pow(x, q, p) == 1 else -1 for x in range(1, p)]
    Q = [x for x in range(1, p) if chi[x] == 1]
    C = [x for x in Q if chi[(1 - x) % p] == 1]
    m = len(C)
    assert m == (p - 5) // 4
    Cset = set(C)
    assert all(2 * (x in Cset) == 1 + chi[(1 - x) % p] - (x == 1) for x in Q)
    index = {x: i for i, x in enumerate(C)}
    R = [index[(1 - x) % p] for x in C]
    I = [index[pow(x, -1, p)] for x in C]
    assert all(R[R[i]] == I[I[i]] == i for i in range(m))
    U = [R[I[i]] for i in range(m)]
    assert all(U[U[U[i]]] == i for i in range(m))
    L = {t: sum(chi[y * (y - 1) * (y - t) % p] for y in range(p)) for t in Q}
    K = [[L[y * pow(x, -1, p) % p] for y in C] for x in C]
    S = [[chi[(x - y) % p] for y in C] for x in C]
    A3 = [[K[i][j] + K[R[i]][R[j]] + K[R[I[i]]][R[I[j]]]
           for j in range(m)] for i in range(m)]
    for i in range(m):
        for j in range(m):
            assert K[i][j] == K[j][i] == K[I[i]][I[j]]
            assert A3[i][j] == A3[I[i]][I[j]] == A3[R[i]][R[j]]
            square = sum(S[i][k] * S[k][j] for k in range(m))
            assert 4 * square == p * (i == j) + A3[i][j] - 6
    pairs = [(i, I[i]) for i in range(m) if i < I[i]]
    for i, ii in pairs:
        for j, jj in pairs:
            lhs = A3[i][j] - A3[ii][j] - A3[i][jj] + A3[ii][jj]
            direct = K[i][j] - K[ii][j] - K[i][jj] + K[ii][jj]
            reflected = K[R[i]][R[j]] - K[R[ii]][R[j]] - K[R[i]][R[jj]] + K[R[ii]][R[jj]]
            assert lhs == direct + 2 * reflected
    out = dict(p=p, generator=g, Q_size=q, C_size=m,
               support_mask_entries=q, full_matrix_entries=m * m,
               odd_basis_entries=len(pairs) ** 2, all_passed=True)
    if characters:
        out.update(character_checks(p, g, chi, C))
    if p == 97:
        positive = [12, 36, 49, 50, 54, 65]
        negative = [2, 3, 9, 33, 62, 89]
        v = [int(x in positive) - int(x in negative) for x in C]
        assert set(positive + negative) <= Cset
        assert all(v[I[i]] == -v[i] for i in range(m))
        norm = sum(x * x for x in v)
        kform = form(K, v)
        reflected = form(K, [v[R[i]] for i in range(m)])
        average3 = form(A3, v)
        assert (norm, kform, reflected, average3) == (12, 876, -356, 164)
        assert 3 * kform > 2 * p * norm
        assert average3 < 2 * p * norm
        out['witness'] = dict(C_in_increasing_order=C, vector=v,
                              positive_support=positive, negative_support=negative,
                              inversion_odd=True, supported_on_C=True,
                              norm_squared=norm, K_quadratic_form=kform,
                              reflected_K_quadratic_form=reflected,
                              three_A_quadratic_form=average3,
                              K_rayleigh=str(Fraction(kform, norm)),
                              A_rayleigh=str(Fraction(average3, 3 * norm)),
                              K_margin_over_two_p_over_three=3 * kform - 2 * p * norm,
                              exact_finite_counterexample_only=True)
    return out


def main():
    primes = [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 109, 113, 193, 241, 257, 641]
    character_primes = {13, 17, 29, 37, 41, 97}
    cases = [verify(p, characters=p in character_primes) for p in primes]
    fields = ['support_mask_entries', 'full_matrix_entries', 'odd_basis_entries',
              'projection_character_identities', 'reflection_character_identities']
    totals = {k: sum(c.get(k, 0) for c in cases) for k in fields}
    out = dict(scope='Exact support/reflection formulas and a C-supported unaveraged finite witness. No asymptotic compressed or averaged norm bound follows from these checks.',
               arithmetic='Standard-library integers, fractions, and exact cyclotomic polynomial remainders. Floating point was used only in a separate exploratory witness search; it is not used in this checker.',
               cases=cases, prime_cases=len(cases), character_prime_cases=len(character_primes),
               totals=totals, all_passed=True)
    (ROOT / 'results/parallel54_supported_elliptic_2026_09_06.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(dict(prime_cases=len(cases), character_prime_cases=len(character_primes), **totals, all_passed=True)))
    print(json.dumps(next(c['witness'] for c in cases if 'witness' in c)))


if __name__ == '__main__':
    main()
