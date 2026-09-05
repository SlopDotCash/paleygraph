#!/usr/bin/env python3
"""Exact checks on actual prime-field two-anchor common neighborhoods."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def chars(p):
    return [0] + [1 if pow(x, (p-1)//2, p) == 1 else -1
                  for x in range(1, p)]


def primitive_root(p):
    n, factors, d = p-1, [], 2
    while d*d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return next(g for g in range(2, p)
                if all(pow(g, (p-1)//q, p) != 1 for q in factors))


def gaussian_jacobi(p, ch):
    g = primitive_root(p)
    eta = [(0, 0)] * p
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    x = 1
    for j in range(p-1):
        eta[x] = units[j % 4]
        assert eta[x][0]**2 - eta[x][1]**2 == ch[x]
        x = x*g % p
    assert x == 1 and len(set(pow(g, j, p) for j in range(p-1))) == p-1
    a = sum(eta[x][0]*ch[(1-x) % p] for x in range(p))
    b = sum(eta[x][1]*ch[(1-x) % p] for x in range(p))
    return g, a, b


def surd_sign(a, b, p):
    """Exact sign of a+b*sqrt(p), for rational a,b and nonsquare p."""
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    if a*a > b*b*p:
        return (a > 0) - (a < 0)
    if a*a < b*b*p:
        return (b > 0) - (b < 0)
    return 0


def quadratic_root_identity():
    reports = []
    for p in (3, 5, 7, 11, 13, 17, 29):
        assert prime(p)
        ch = chars(p)
        cases = 0
        for A in range(p):
            for B in range(1, p):
                left = sum(ch[t*(t*t + A*t + B) % p] for t in range(p))
                right = sum(ch[z*(z*z - 2*A*z + A*A - 4*B) % p]
                            for z in range(p))
                assert left == right
                for z in range(p):
                    direct = sum((t+A+B*pow(t, -1, p)) % p == z
                                 for t in range(1, p))
                    assert direct == 1 + ch[((z-A)**2 - 4*B) % p]
                    COUNTS['quadratic_map_fibres'] += 1
                cases += 1
                COUNTS['quadratic_trace_identities'] += 1
        reports.append({'p': p, 'all_A_nonzero_B_pairs': cases})
    return reports


def actual_neighborhoods():
    reports = []
    for p in range(13, 1001, 4):
        if not prime(p):
            continue
        ch = chars(p)
        g, a, b = gaussian_jacobi(p, ch)
        assert a*a + b*b == p
        f = [ch[x]*ch[(x-1) % p] for x in range(p)]
        assert sum(f) == -1
        Sf = [sum(ch[(x-y) % p]*f[y] for y in range(p)) for x in range(p)]
        U = sum(f[x]*Sf[x] for x in range(p))
        assert U == 2*(a*a-b*b)
        assert -2*p <= U <= 2*p
        C = [x for x in range(p) if ch[x] == ch[(x-1) % p] == 1]
        m = len(C)
        assert m == (p-5)//4 and m > 0
        row_sums = [sum(ch[(x-y) % p] for y in C) for x in C]
        R = sum(row_sums)
        assert 16*R == U - 6*p + 36
        assert 9-2*p <= 4*R <= 9-p
        assert -2-Q(1, p-5) <= Q(R, m) <= -1+Q(4, p-5)
        for x, row_sum in zip(C, row_sums):
            # Polynomial evaluation independent of the f/S convolution.
            Lx = sum(ch[y*(y-1)*(y-x) % p] for y in range(p))
            assert Lx == Sf[x]
            assert 4*row_sum == Lx-6
            COUNTS['actual_row_trace_identities'] += 1
        variance = Q(sum(r*r for r in row_sums), m) - Q(R, m)**2
        trace_variance = Q(sum((Sf[x]-6)**2 for x in C), 16*m) - Q(R, m)**2
        assert variance == trace_variance >= 0
        alpha, beta = Q(3*p+5, 8*p), Q(R, 2*m*p)
        assert alpha == Q(1, 2)*(1-Q(m, p))
        assert surd_sign(alpha, beta, p) > 0
        assert surd_sign(1-alpha, -beta, p) > 0
        assert R < 0
        if p >= 73:
            assert surd_sign(alpha-Q(1, 4), beta, p) >= 0
            assert surd_sign(Q(1, 2)-alpha, -beta, p) > 0
            COUNTS['constant_Fourier_gap_certificates'] += 1
        affine_edges = 0
        if p <= 53:
            for aa in range(p):
                for bb in range(aa+1, p):
                    if ch[(bb-aa) % p] != 1:
                        continue
                    actual = [x for x in range(p)
                              if ch[(x-aa) % p] == ch[(x-bb) % p] == 1]
                    transported = sorted((aa+(bb-aa)*t) % p for t in C)
                    assert actual == transported
                    total = sum(ch[(x-y) % p] for x in actual for y in actual)
                    assert total == R
                    affine_edges += 1
                    COUNTS['actual_affine_edges'] += 1
        reports.append({
            'p': p, 'm': m, 'primitive_root': g, 'Jacobi_sum_Gaussian_integer': [a, b],
            'U': U, 'R_C': R, 'R_C_over_m': str(Q(R, m)),
            'QR_Fourier_mass': {'rational': str(alpha), 'sqrt_p_coefficient': str(beta)},
            'constant_direction_residual_squared': str(variance),
            'affine_edges_checked': affine_edges,
        })
        COUNTS['actual_prime_neighborhoods'] += 1
        COUNTS['exact_Jacobi_norms'] += 1
        COUNTS['double_sum_Jacobi_identities'] += 1
        COUNTS['neighborhood_quadratic_form_identities'] += 1
        COUNTS['constant_direction_variance_identities'] += 1
    return reports


def main():
    fibres = quadratic_root_identity()
    print('Exact quadratic-root identities passed.', flush=True)
    fields = actual_neighborhoods()
    inputs = [
        'research/parallel21-spectral-next-input-2026-09-05.md',
        'experiments/parallel21_spectral_flat_direction_2026_09_05.py',
        'research/parallel17-prime-uncertainty-2026-09-05.md',
        'research/parallel13-principal-budget-2026-09-05.md',
    ]
    result = {
        'status': 'All finite exact checks passed; no all-vector Fourier bound or spectral edge is proved.',
        'arithmetic': 'Python integers, Gaussian integers, Fraction, exact signs in Q(sqrt(p)).',
        'input_sha256': {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs},
        'counts': dict(COUNTS),
        'quadratic_root_fields': fibres,
        'actual_prime_fields': fields,
        'limitations': [
            'The uniform theorem concerns the constant vector on two-anchor common neighborhoods only.',
            'The verified residual identity supplies no asymptotic estimate for its right side.',
            'Passing to more anchors or arbitrary supported vectors is not justified.',
            'No determinant bound for the full matrix or logarithmic depth range is improved.',
            'No arbitrary projection, square-field example or abstract sign kernel is used as a counterexample.',
            'This lane has no separate human or formal proof review.',
        ],
    }
    output = ROOT/'results/parallel21_spectral_flat_direction_2026_09_05.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'counts': dict(COUNTS)}, indent=2), flush=True)


if __name__ == '__main__':
    main()
