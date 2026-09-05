#!/usr/bin/env python3
"""Exact actual-field checks; published cohomology is not verified by this run."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def check(ok, name):
    assert ok, name
    COUNTS[name] += 1


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def surd_mul(x, y, p):
    return (x[0]*y[0]+p*x[1]*y[1], x[0]*y[1]+x[1]*y[0])


def field(p):
    chi = [0] + [1 if pow(t, (p-1)//2, p) == 1 else -1 for t in range(1, p)]
    C = [t for t in range(2, p) if chi[t] == chi[t-1] == 1]
    m = len(C)
    check(prime(p) and p % 4 == 1 and m == (p-5)//4, 'actual_prime_neighborhood')
    L = [sum(chi[y*(y-1)*(y-t) % p] for y in range(p)) for t in range(p)]
    # Full S D_0 action, independently of the C-compressed recurrence.
    k = [sum(chi[(t-y) % p]*chi[y]*L[y] for y in range(p)) for t in range(p)]
    f0 = [chi[t]*chi[(t-1) % p] for t in range(p)]
    UJ = sum(f0[t]*L[t] for t in range(p))
    check(L[0] == L[1] == -1 and k[0] == 1 and k[1] == UJ,
          'full_field_boundary_values')
    check(abs(UJ) <= 2*p, 'Jacobi_mean_bound_finite_check')
    check(sum(k[t]**2 for t in range(1, p)) == p**3-2*(p*p+p+1),
          'exact_rank_three_Mellin_norm')
    if p <= 53:
        inv = {x: pow(x, -1, p) for x in range(1, p)}
        f = [chi[(1-t) % p] for t in range(p)]
        k2 = [0]+[sum(f[a]*f[t*inv[a] % p] for a in range(1, p)) for t in range(1, p)]
        conv3 = [0]+[sum(k2[a]*f[t*inv[a] % p] for a in range(1, p)) for t in range(1, p)]
        check(k2[1:] == L[1:] and conv3[1:] == k[1:], 'independent_multiplicative_convolution')
    pulls = [[k[t] for t in C], [k[(1-t) % p] for t in C],
             [k[(1-pow(t, -1, p)) % p] for t in C]]
    K = [sum(q[i] for q in pulls) for i in range(m)]
    LC = [L[t] for t in C]
    SC = [[chi[(x-y) % p] for y in C] for x in C]
    rows = [sum(row) for row in SC]
    SL = [sum(a*b for a, b in zip(row, LC)) for row in SC]
    for i, t in enumerate(C):
        check(4*rows[i] == L[t]-6, 'actual_neighborhood_row')
        check(4*SL[i] == p+6+K[i], 'actual_rank_three_recurrence')
        check(L[pow(t, -1, p)] == L[t] and L[(1-t) % p] == L[t]
              and k[pow(t, -1, p)] == k[t], 'actual_S3_parameter_symmetries')
    # All finite analytic inequalities are checked after squaring; no float.
    D = range(2, p)
    for e in range(2):
        for f in range(2):
            weights = {t: chi[t]**e*chi[t-1]**f for t in D}
            mixed = sum(weights[t]*L[t]*k[t] for t in D)
            single = sum(weights[t]*k[t] for t in D)
            check(mixed*mixed <= 36*p**4, 'mixed_rank_two_three_twist_bound')
            check(single*single <= 9*p**3, 'single_rank_three_twist_bound')
            if e or f:
                z2 = sum(weights[t]*L[t]**2 for t in D)
                z3 = sum(weights[t]*k[t]**2 for t in D)
                check(z2*z2 <= 16*p**3, 'nontrivial_rank_two_diagonal_twist')
                check(z3*z3 <= 81*p**5, 'nontrivial_rank_three_diagonal_twist')
            all_pulls = [[k[t] for t in D], [k[(1-t) % p] for t in D],
                         [k[(1-pow(t, -1, p)) % p] for t in D]]
            for i in range(3):
                for j in range(i+1, 3):
                    corr = sum(weights[t]*all_pulls[i][t-2]*all_pulls[j][t-2] for t in D)
                    check(corr*corr <= 81*p**5, 'distinct_reflection_pullback_twist')
    ell = F(sum(LC), m)
    check(ell == F(UJ+6, p-5) and abs(ell) <= 4, 'exact_bounded_mean')
    g = [F(x)-ell for x in LC]
    d = sum(x*x for x in g)
    s0 = sum(pulls[0])
    s = sum(a*b for a, b in zip(g, pulls[0]))
    check(sum(K) == 3*s0 and sum(a*b for a, b in zip(g, K)) == 3*s,
          'S3_symmetrized_mean_and_pairing')
    b2 = d/(64*m*p)
    q = ((1-F(m, p))/2, (ell-6)/(8*p))
    report = {'p': p, 'm': m, 'd': str(d), 'ell': str(ell), 's0': s0,
              's': str(s), 'b1_squared': str(b2), 'q': [str(x) for x in q],
              'K_squared_norm': sum(x*x for x in K)}
    if d == 0:
        COUNTS['zero_first_residual_fields'] += 1
        report['next_direction'] = 'undefined: d=0; not divided by zero'
        return report
    Sg = [sum(a*b for a, b in zip(row, g)) for row in SC]
    aS = sum(a*b for a, b in zip(g, Sg))/d
    check(aS == (3*s/d-ell)/4, 'exact_next_diagonal_coefficient')
    check(sum(Sg) == d/4, 'exact_uniform_coupling')
    h = [F(K[i])-F(3*s0, m)-3*s/d*g[i] for i in range(m)]
    h2 = sum(x*x for x in h)
    check(sum(h) == sum(a*b for a, b in zip(g, h)) == 0, 'exact_second_residual_orthogonality')
    check(h2 == sum(x*x for x in K)-F(9*s0*s0, m)-9*s*s/d,
          'exact_second_residual_norm_ledger')
    z = [Sg[i]-d/(4*m)-aS*g[i] for i in range(m)]
    check(all(4*z[i] == h[i] for i in range(m)), 'direct_compressed_second_residual')
    c2 = h2/(64*p*d)
    check(c2 == sum(x*x for x in z)/(4*p*d) >= 0, 'exact_beta_two_squared')
    a = (F(1, 2), aS/(2*p))
    # ||Bv||² = 1/4 + ||Sg||²/(4pd) + aS/(2sqrt(p)).
    a_sq = surd_mul(a, a, p)
    direct_v_norm2 = (F(1, 4)+sum(x*x for x in Sg)/(4*p*d), aS/(2*p))
    check(direct_v_norm2 == (a_sq[0]+b2+c2, a_sq[1]), 'retained_leakage_Gram_diagonal')
    q_sq = surd_mul(q, q, p)
    check(q_sq[0]+b2 == ((1-F(m, p))/2)**2+F(sum(x*x for x in rows), 4*m*p),
          'actual_uniform_B_squared_norm')
    if p >= 1024:
        check(d >= F(p*p, 8), 'explicit_positive_variance_threshold')
        check(abs(aS) <= 40, 'explicit_diagonal_error_constant')
    report.update({'next_direction': 'defined', 'a1': [str(x) for x in a],
                   'a1_S_rayleigh': str(aS), 'beta2_squared': str(c2),
                   'h_squared_norm': str(h2),
                   'Bv_squared_norm': [str(x) for x in direct_v_norm2],
                   'Bu_squared_norm': [str(q_sq[0]+b2), str(q_sq[1])]})
    return report


def main():
    reports = [field(p) for p in range(13, 1101, 4) if prime(p)]
    inputs = ['research/parallel24-spectral-operator-2026-09-05.md',
              'experiments/parallel24_spectral_operator_2026_09_05.py',
              'research/parallel23-spectral-coupling-2026-09-05.md',
              'research/parallel21-spectral-next-input-2026-09-05.md',
              'research/parallel6-seeded-kernels-2026-09-04.md',
              'sources/katz-g2-hypergeometric.pdf',
              'sources/katz-gauss-kloosterman-monodromy.pdf',
              'sources/katz-g2-hypergeometric.txt',
              'sources/katz-gauss-kloosterman-monodromy.txt']
    result = {'status': 'PASS', 'arithmetic': 'Integers, Fraction, exact Q(sqrt(p)) pairs.',
              'scope': 'Actual K2 coefficients and retained leakage; no all-vector spectral edge.',
              'counts': dict(COUNTS), 'fields': reports,
              'source_urls': ['https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf',
                              'https://web.math.princeton.edu/~nmk/Katz-GKM.pdf'],
              'input_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
              'limitations': ['Finite checks do not prove the imported cohomological results.',
                              'No inference of asymptotic limits from finite trends.',
                              'The next h-S-h input and full operator remain unbounded here.']}
    out = ROOT/'results/parallel24_spectral_operator_2026_09_05.json'
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'counts': dict(COUNTS), 'output': str(out)}, indent=2))


if __name__ == '__main__':
    main()
