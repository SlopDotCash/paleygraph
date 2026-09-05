#!/usr/bin/env python3
"""Exact actual-field checks; no independent modular-form computation."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def field(p):
    assert prime(p) and p % 4 == 1
    ch = [0] + [1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1, p)]
    # Direct polynomial trace, independent of the neighborhood row sums.
    L = [sum(ch[y*(y-1)*(y-t) % p] for y in range(p)) for t in range(p)]
    C = [t for t in range(2, p) if ch[t] == ch[t-1] == 1]
    m = len(C)
    assert m == (p-5)//4 and m > 0
    assert L[0] == L[1] == -1
    T = sum(x*x for x in L[2:])
    W = sum(L[s*s % p]**2 for s in range(2, p-1))
    Z = sum(L[t]**2 for t in C)
    A = sum(ch[t]*L[t]**2 for t in range(2, p))
    A1 = sum(ch[t-1]*L[t]**2 for t in range(2, p))
    A01 = sum(ch[t]*ch[t-1]*L[t]**2 for t in range(2, p))
    assert A == A1 == A01
    assert W == T+A and 4*Z == T+3*A == 3*W-2*T
    assert T == p*p-2*p-3
    COUNTS['prime_fields'] += 1
    COUNTS['restricted_second_moment_identities'] += 1
    for t in range(2, p):
        assert L[(1-t) % p] == L[t]
        assert L[pow(t, -1, p)] == ch[t]*L[t]
        assert sum(s*s % p == t for s in range(2, p-1)) == 1+ch[t]
        COUNTS['parameter_symmetries_and_square_fibres'] += 1
    # U from a double sum; check the bounded mean with the previous
    # independently computed Gaussian-integer Jacobi data where available.
    f = [ch[t]*ch[(t-1) % p] for t in range(p)]
    U = sum(f[t]*L[t] for t in range(p))
    assert abs(U) <= 2*p
    ell = F(sum(L[t] for t in C), m)
    assert 4*sum(L[t] for t in C) == U+6
    rows = [sum(ch[(x-y) % p] for y in C) for x in C]
    for t, row in zip(C, rows):
        assert 4*row == L[t]-6
        COUNTS['actual_neighborhood_rows'] += 1
    R = sum(rows)
    variance = F(sum(x*x for x in rows), m)-F(R, m)**2
    assert variance == (F(Z, m)-ell**2)/16 >= 0
    b2 = variance/(4*p)
    assert b2 == (F(Z, m)-ell**2)/(64*p)
    # Algebraic label only: this is not an independent Hecke computation.
    a8 = p*(p-3)-4-W
    assert A == -a8-p-1
    assert 4*Z == p*p-5*p-6-3*a8
    exact = F(1, 64)-F(6+3*a8, 64*p*(p-5))-F((U+6)**2, 64*p*(p-5)**2)
    assert b2 == exact
    # Multiply Bu by sqrt(m). Its entries are a+d_x*sqrt(p),
    # where a=(1-m/p)/2 and d_x=row/(2p). Subtracting q leaves
    # only the nonconstant d_x term, exactly in Q(sqrt(p)).
    aa = (1-F(m, p))/2
    ds = [F(row, 2*p) for row in rows]
    dbar = sum(ds)/m
    assert F(p, m)*sum((d-dbar)**2 for d in ds) == b2
    # Squared norm of Bu = rational part + sqrt(p) part.
    norm2_r = aa*aa + F(p, m)*sum(d*d for d in ds)
    norm2_s = 2*aa*sum(ds)/m
    assert norm2_r == aa*aa+p*dbar*dbar+b2
    assert norm2_s == 2*aa*dbar
    COUNTS['exact_projection_variances'] += 1
    # The removed block is (c_x+c_y)/(2m sqrt(p)), c_x=row-R/m.
    # Its squared Frobenius norm is 2*b^2 and its action on 1 is c/(2sqrt(p)).
    if p <= 101:
        cs = [F(row)-F(R,m) for row in rows]
        assert sum(cs) == 0
        frob = F(0)
        for c in cs:
            assert sum(c+d for d in cs) == m*c
            for d in cs:
                frob += (c+d)**2/(4*m*m*p)
                COUNTS['exact_deleted_block_entries'] += 1
        assert frob == 2*b2
    return {'p':p,'m':m,'T':T,'W':W,'Z':Z,'common_twisted_sum':A,
            'a8_inferred_from_imported_formula_not_independently_computed':a8,
            'U':U,'ell':str(ell),'b_squared':str(b2),
            'Z_over_mp':str(F(Z,m*p)),
            'q':{'rational':str(aa),'sqrt_p_coefficient':str(dbar)},
            'u_B_squared_u':{'rational':str(norm2_r),'sqrt_p_coefficient':str(norm2_s)}}


def main():
    reports = [field(p) for p in range(13,1001,4) if prime(p)]
    prior_path = 'results/parallel21_spectral_flat_direction_2026_09_05.json'
    prior = json.loads((ROOT/prior_path).read_text())
    old = {f['p']:f for f in prior['actual_prime_fields']}
    for f in reports:
        z = old[f['p']]
        a,b = z['Jacobi_sum_Gaussian_integer']
        assert f['U'] == z['U'] == 2*(a*a-b*b)
        assert F(f['b_squared'])*4*f['p'] == F(z['constant_direction_residual_squared'])
        COUNTS['prior_independent_Jacobi_and_variance_comparisons'] += 1
    paths = ['research/parallel23-spectral-coupling-2026-09-05.md',
             'experiments/parallel23_spectral_coupling_2026_09_05.py',
             'sources/grove-hypergeometric-moments-2026.html',
             'results/parallel23_spectral_source_2026_09_05.json',
             'research/parallel21-spectral-next-input-2026-09-05.md',prior_path]
    result = {'status':'All exact finite identities passed; asymptotic coupling limit relies on the cited existing theorem.',
              'arithmetic':'Python integers and Fraction; Q(sqrt(p)) pairs. No floating point.',
              'input_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
              'counts':dict(COUNTS),'prime_fields':reports,
              'limitations':['No independent Hecke trace computation or check of the published trace formula.',
                             'Finite checks do not establish the asymptotic second-moment theorem.',
                             'No all-vector upper estimate, spectral edge, classical Paley or prize proof.']}
    out=ROOT/'results/parallel23_spectral_coupling_2026_09_05.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(out),'counts':dict(COUNTS)},indent=2))


if __name__ == '__main__':
    main()
