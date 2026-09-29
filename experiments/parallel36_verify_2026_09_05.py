#!/usr/bin/env python3
"""Independent rational checks of the pass-36 formulas and C++ certificates.

No Lean run, external package, or floating-point acceptance criterion is used.
The universal analytic argument remains an ordinary mathematical proof.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
BINARY = Path('/private/tmp/paley-parallel36-shell')
CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def atan_bounds(q, terms=32):
    s = sum((F((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(terms)), F(0))
    e = F(1, (2*terms+1)*q**(2*terms+1))
    return (s, s+e) if terms % 2 == 0 else (s-e, s)


a5, b5 = atan_bounds(5)
a239, b239 = atan_bounds(239)
PI_LO, PI_HI = 16*a5-4*b239, 16*b5-4*a239
PI2_FLOOR = F(9869604401089, 10**12)


@lru_cache(maxsize=None)
def cosine_bounds(p, r):
    r = min(r % p, p-r % p)
    # Interval evaluation of the Taylor polynomial, with an absolute remainder.
    xl, xu = 2*PI_LO*r/p, 2*PI_HI*r/p
    lo = hi = F(0)
    for j in range(21):
        l, u = xl**(2*j)/math.factorial(2*j), xu**(2*j)/math.factorial(2*j)
        if j % 2:
            lo -= u
            hi -= l
        else:
            lo += l
            hi += u
    error = xu**42/math.factorial(42)
    return lo-error, hi+error


def mu_sieve(limit):
    mu = [1]*(limit+1)
    primes = [True]*(limit+1)
    for q in range(2, limit+1):
        if primes[q]:
            for j in range(q, limit+1, q):
                mu[j] *= -1
                if j > q:
                    primes[j] = False
            for j in range(q*q, limit+1, q*q):
                mu[j] = 0
    return mu


MU = mu_sieve(4096)


def b_formula(k):
    t, m = 0, k
    while m % 2 == 0:
        t, m = t+1, m//2
    return F(-MU[m], m*m*(2**(t+1) if t else 1))


def trunc_weight(k, scale):
    b = b_formula(k)
    return int(b*scale)  # Fraction -> int rounds toward zero.


def certificate_upper(record, maximum=None):
    p, n = record['p'], record['n']
    D = F(record['max_deviation_numerator'], record['distance_denominator'])
    maximum = int(record['weighted_max_abs_numerator']) if maximum is None else maximum
    benefit = record['max_deviation_numerator']*record['weight_mass']-maximum
    check(benefit >= 0, 'nonnegative fixed-weight benefit')
    R = F(benefit, record['distance_denominator']*record['weight_scale'])
    return 36*F(p*p, p*p+1)*D-2*PI2_FLOOR*R+F(n, p-1)


def direct_rows(p, A):
    return [12*sum(min(a*h % p, p-a*h % p)**2 for h in A)
            -len(A)*p*(p+1) for a in range(1, p)]


def abs_interval(lo, hi):
    return (0 if lo <= 0 <= hi else min(abs(lo), abs(hi)), max(abs(lo), abs(hi)))


def det_bareiss(matrix):
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a)-1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k+1, len(a)) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        v = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                num = a[i][j]*v-a[i][k]*a[k][j]
                check(num % previous == 0, 'Bareiss exact division')
                a[i][j] = num//previous
        for i in range(k+1, len(a)):
            a[i][k] = 0
        previous = v
    return sign*a[-1][-1]


def norm_case(p, n, g, a):
    N = n//2
    check(n & (n-1) == 0 and pow(g, N, p) == p-1, 'dyadic generator')
    r = [((a*pow(g, j, p)+p//2) % p)-p//2 for j in range(N)]
    values = [sum(r[j]*pow(pow(g, k, p), j, p) for j in range(N)) % p
              for k in range(1, n, 2)]
    check(sum(v != 0 for v in values) == 1, 'exact rank-one evaluation')
    check(values[-1] == a*N % p, 'exceptional root is g inverse')
    matrix = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            k = (i-j) % N
            matrix[i][j] = r[k]*(-1 if i < j else 1)
    norm = det_bareiss(matrix)
    S = sum(x*x for x in r)
    check(norm != 0 and norm % p**(N-1) == 0, 'nonzero norm divisibility')
    check(S**N >= norm**2 >= p**(2*N-2), 'Parseval AMGM norm bound')
    check(S % p == 0, 'squared-distance congruence')
    return {'p': p, 'n': n, 'generator': g, 'frequency': a,
            'centered_power_coordinates': r, 'square_sum': S,
            'norm': norm, 'norm_over_p_to_N_minus_1': str(F(norm, p**(N-1)))}


def main():
    started = time.perf_counter()
    check(PI_LO > 3 and PI_HI < F(22, 7), 'Machin interval')
    check(PI2_FLOOR < PI_LO**2, 'explicit rational pi-squared lower bound')
    # Recover the inverse recursively, independently of its claimed closed form.
    inverse = [F(0)]*4097
    inverse[1] = -1
    for k in range(2, 4097):
        total = F(0)
        for d in range(2, math.isqrt(k)+1):
            if k % d == 0:
                total += F((-1)**d, d*d)*inverse[k//d]
                if d*d != k:
                    e = k//d
                    total += F((-1)**e, e*e)*inverse[d]
        total += F((-1)**k, k*k)*inverse[1]
        inverse[k] = total
    for k in range(1, 4097):
        check(inverse[k] == b_formula(k), 'Dirichlet inverse formula')

    # General symmetric sets, including non-subgroups and all nonzero residues.
    sets = [(3, [1, 2]), (7, [1, 2, 5, 6]), (13, [1, 3, 4, 9, 10, 12]),
            (17, [1, 2, 3, 14, 15, 16]), (29, [1, 4, 7, 22, 25, 28]),
            (17, list(range(1, 17)))]
    analytic_cases = []
    for p, A in sets:
        check(len(set(A)) == len(A) and set(A) == {p-a for a in A}, 'symmetric distinct nonzero set')
        ds = direct_rows(p, A)
        check(sum(ds) == 0, 'exact nonprincipal mean')
        D = F(max(map(abs, ds)), 24*p*p)
        mf_lo, mf_hi, m_hi = 0, 0, 0
        for a in range(1, p):
            intervals = [cosine_bounds(p, a*h) for h in A]
            lo, hi = sum(v[0] for v in intervals), sum(v[1] for v in intervals)
            m_hi = max(m_hi, abs_interval(lo, hi)[1])
            lo, hi = lo+F(len(A), p-1), hi+F(len(A), p-1)
            l, u = abs_interval(lo, hi)
            mf_lo, mf_hi = max(mf_lo, l), max(mf_hi, u)
        if D:
            invariant_two = {2*a % p for a in A} == set(A)
            coefficient = 24 if invariant_two else 12
            check(mf_lo >= coefficient*F(p*p, p*p-1)*D, 'forward norm bound with rigorous cosine intervals')
            check(mf_hi <= 36*F(p*p, p*p+1)*D, 'inverse norm bound with rigorous cosine intervals')
        analytic_cases.append({'p': p, 'A': A, 'D': str(D)})

    # Recompute full small-field C++ output using direct cosets and rational weights.
    cpp_cases = []
    for p, n, L in [(3, 2, 16), (5, 2, 16), (5, 4, 16), (7, 6, 16),
                    (17, 4, 32), (17, 8, 32), (17, 16, 32), (41, 8, 64),
                    (97, 16, 128), (257, 16, 256), (1153, 8, 128), (33713, 16, 256)]:
        out = subprocess.run([str(BINARY), str(p), str(n), str(L)], check=True,
                             capture_output=True, text=True, timeout=15)
        rec = json.loads(out.stdout)
        g, h, q = rec['primitive_generator'], rec['subgroup_generator'], rec['cosets']
        H = [pow(h, j, p) for j in range(n)]
        check(len(set(H)) == n and pow(h, n, p) == 1, 'independent subgroup order')
        representatives = [pow(g, j, p) for j in range(q)]
        label = {}
        cs = []
        for j, a in enumerate(representatives):
            coset = [a*x % p for x in H]
            for x in coset:
                check(x not in label, 'disjoint cosets')
                label[x] = j
            cs.append(sum(min(x, p-x)**2 for x in coset))
        check(len(label) == p-1, 'all nonzero elements covered')
        ds = [12*c-n*p*(p+1) for c in cs]
        weights = [(k, trunc_weight(k, rec['weight_scale'])) for k in range(1, L+1) if k % p]
        ws = [sum(w*ds[label[k*a % p]] for k, w in weights) for a in representatives]
        check(rec['max_deviation_numerator'] == max(map(abs, ds)), 'full D independently recomputed')
        check(rec['minimum_all_H_square_sum'] == min(cs) and rec['maximum_all_H_square_sum'] == max(cs), 'distance extrema')
        check(rec['weight_mass'] == sum(abs(w) for _, w in weights), 'coefficient mass')
        check(int(rec['weighted_max_abs_numerator']) == max(map(abs, ws)), 'full weighted maximum')
        check(int(rec['benefit_numerator']) == max(map(abs, ds))*rec['weight_mass']-max(map(abs, ws)), 'benefit')
        cpp_cases.append({'p': p, 'n': n, 'truncation': L, 'cosets': q, 'upper': str(certificate_upper(rec))})

    norms = []
    for p, n, g in [(5, 4, 2), (17, 8, 2), (17, 16, 3), (41, 8, 3),
                    (97, 16, 8), (257, 16, 2), (1153, 8, 75)]:
        # Choose an exact-order generator rather than assuming the listed candidate works.
        if pow(g, n//2, p) != p-1:
            g = next(x for x in range(2, p) if pow(x, n//2, p) == p-1)
        for a in [1, 2, 3, p//2]:
            norms.append(norm_case(p, n, g, a))

    large = []
    for L in [256, 4096]:
        out = subprocess.run([str(BINARY), '6700417', '64', str(L)], check=True,
                             capture_output=True, text=True, timeout=15)
        rec = json.loads(out.stdout)
        result_path = ROOT/f'results/parallel36_shell_6700417_64_{L}_2026_09_05.json'
        result_path.write_text(json.dumps(rec, indent=2)+'\n')
        p, n, g, h = rec['p'], rec['n'], rec['primitive_generator'], rec['subgroup_generator']
        H = [pow(h, j, p) for j in range(n)]
        check(set(H) == {pow(2, j, p) for j in range(n)}, 'large subgroup is generated by two')
        def d_at(a):
            return 12*sum(min(a*x % p, p-a*x % p)**2 for x in H)-n*p*(p+1)
        for key in ['minimum_coset_index', 'maximum_coset_index']:
            a = pow(g, rec[key], p)
            c = (d_at(a)+n*p*(p+1))//12
            expected = rec['minimum_all_H_square_sum' if key.startswith('minimum') else 'maximum_all_H_square_sum']
            check(c == expected, 'large distance extremizer by direct residues')
        weights = [(k, trunc_weight(k, rec['weight_scale'])) for k in range(1, L+1)]
        check(sum(abs(w) for _, w in weights) == rec['weight_mass'], 'large weight mass')
        for key, value in [('weighted_max_coset_index', 'weighted_max_abs_numerator'),
                           ('weighted_runner_up_coset_index', 'weighted_runner_up_abs_numerator')]:
            a = pow(g, rec[key], p)
            direct = sum(w*d_at(k*a % p) for k, w in weights if w)
            check(abs(direct) == int(rec[value]), 'large weighted extremizer by direct residues')
        upper = certificate_upper(rec)
        other_upper = certificate_upper(rec, int(rec['weighted_runner_up_abs_numerator']))
        eta_intervals = [cosine_bounds(p, pow(2, j, p)) for j in range(32)]
        eta_lo, eta_hi = 2*sum(v[0] for v in eta_intervals), 2*sum(v[1] for v in eta_intervals)
        check(eta_lo > 43 and eta_hi < upper, 'rigorous finite period enclosure')
        check(upper**2 < 1970, 'strict improvement on prior finite moment certificate')
        grid = 10**12
        rounded_lo = F((eta_lo*grid).__floor__(), grid)
        rounded_hi = F((eta_hi*grid).__ceil__(), grid)
        maximum_identified = rec['weighted_max_coset_index'] == 0 and other_upper < eta_lo
        large.append({'p': p, 'n': n, 'truncation': L, 'certificate_upper': str(upper),
                      'upper_float_for_display_only': float(upper),
                      'other_cosets_upper': str(other_upper),
                      'other_cosets_upper_float_for_display_only': float(other_upper),
                      'eta_one_interval': [str(rounded_lo), str(rounded_hi)],
                      'maximum_proved_to_occur_on_H': maximum_identified,
                      'cpp_elapsed_seconds': rec['elapsed_seconds']})

    out = {'created_at_utc': datetime.now(timezone.utc).isoformat(), 'checks': CHECKS,
           'elapsed_seconds': time.perf_counter()-started,
           'dirichlet_inverse_coefficients_checked': 4096,
           'general_symmetric_cases': analytic_cases, 'full_cpp_comparisons': cpp_cases,
           'dyadic_norm_cases': norms, 'large_certificates': large,
           'pi_squared_rational_lower': str(PI2_FLOOR),
           'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'cpp_sha256': hashlib.sha256((ROOT/'experiments/parallel36_shell_certificate.cpp').read_bytes()).hexdigest(),
           'scope': 'Exact finite verification by the same author using independent computational paths. The analytic argument is not Lean-verified or independently reviewed; the uniform bound and full conjectures remain unproved.'}
    (ROOT/'results/parallel36_verification_2026_09_05.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'checks': CHECKS, 'elapsed_seconds': out['elapsed_seconds'], 'large_certificates': large}, indent=2))


if __name__ == '__main__':
    main()
