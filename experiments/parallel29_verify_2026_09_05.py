#!/usr/bin/env python3
"""Bounded checks of the centered recurrence; not an analytic or Lean proof.

One process, standard library only. Integer/rational identities are exact.
Complex Fourier evaluations are explicitly approximate diagnostics.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import cmath
import json
import math
import random
import time

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def check(ok, label):
    assert ok, label
    COUNTS[label] += 1


def subgroup(p, n):
    # Independent small-field construction by testing roots, not a generator.
    h = [x for x in range(1, p) if pow(x, n, p) == 1]
    check(len(h) == n, 'subgroup_order')
    check(all(a*b % p in h for a in h for b in h), 'subgroup_closure')
    return h


def cosets(p, h):
    left = set(range(1, p))
    out = []
    while left:
        c = {min(left)*a % p for a in h}
        left -= c
        out.append(sorted(c))
    return out


def conv(f, g):
    p = len(f)
    out = [0]*p
    for x, a in enumerate(f):
        for y, b in enumerate(g):
            out[(x+y) % p] += a*b
    return out


def norm2(f):
    return sum(a*a for a in f)


def e0(f):
    return norm2(conv(f, f)) - F(sum(f)**4, len(f))


def incidence_checks():
    ranges = Counter()
    max_ratio_sq = F(0)
    for p in [3, 5, 7]:
        for n in range(1, p):
            if (p-1) % n:
                continue
            h = subgroup(p, n)
            orbits = cosets(p, h)
            for mask in range(1, 1 << len(orbits)):
                q = sorted(x for i, c in enumerate(orbits)
                           if mask >> i & 1 for x in c)
                m = len(q)
                points = list(product(q, q, h))
                planes = list(product(h, q, q))  # h,u,v in x+h*y-v*z=u
                total = len(points)
                incident = lambda pt, pl: (pt[0]+pl[0]*pt[1]-pl[2]*pt[2]-pl[1]) % p == 0
                actual = sum(incident(pt, pl) for pt in points for pl in planes)
                f = [int(x in q) for x in range(p)]
                energy = norm2(conv(f, f))
                check(actual == n*n*energy, 'exact_incidence_encoding')
                check(len(set(planes)) == total == n*m*m, 'balanced_distinct_planes')
                directions = Counter((a, v) for a, u, v in planes)
                # Integer-scaled variance avoids rounding M/p.
                field_counts = [sum(incident(pt, pl) for pl in planes)
                                for pt in product(range(p), repeat=3)]
                variance_scaled = sum((p*c-total)**2 for c in field_counts)
                variance_formula = p*p*(p*p*total-p*sum(v*v for v in directions.values()))
                check(variance_scaled == variance_formula, 'exact_plane_variance')
                check((p*actual-total*total)**2 <= p**4*total**2,
                      'centered_incidence_cauchy_bound')
                t = F(energy)-F(m**4, p)
                check(t >= 0, 'centered_set_energy_nonnegative')
                ratio_sq = t*t*n/F(m**6)
                max_ratio_sq = max(max_ratio_sq, ratio_sq)
                # Merely a finite K=1 check, not the unknown incidence constant.
                check(ratio_sq <= 1, 'selected_sets_K_one_sanity_only')
                ranges['N_le_p_squared' if total <= p*p else 'N_gt_p_squared'] += 1
    return {'size_ranges': dict(ranges), 'largest_squared_centered_set_ratio': str(max_ratio_sq)}


def weight_and_convolution_checks():
    rng = random.Random(290905)
    data = []
    for p, n in [(5, 4), (7, 3), (13, 4), (17, 4), (17, 16), (29, 7), (31, 5)]:
        h = subgroup(p, n)
        orbits = cosets(p, h)
        for _ in range(12):
            g = [0]*p
            for c in orbits:
                w = rng.randrange(-4, 5)
                for x in c:
                    g[x] = w
            a, b = sum(abs(v) for v in g), norm2(g)
            check(e0(g)**2*n <= (10368*a*a*b)**2,
                  'signed_off_origin_K_one_sanity_only')
            g[0] = rng.randrange(-8, 9)
            a, b = sum(abs(v) for v in g), norm2(g)
            residual = max(F(0), e0(g)-8*g[0]**4)
            check(residual**2*n <= (82944*a*a*b)**2,
                  'signed_with_origin_K_one_sanity_only')
        indicator = [int(x in h) for x in range(p)]
        r = [1]+[0]*(p-1)
        rs, ts, fs = {}, {}, {}
        for s in range(1, 13):
            r = conv(r, indicator)
            rs[s] = r
            ts[s] = F(norm2(r))-F(n**(2*s), p)
            fs[s] = [F(v)-F(n**s, p) for v in r]
            check(sum(r) == n**s, 'exact_convolution_mass')
            check(norm2(fs[s]) == ts[s], 'exact_centered_second_moment')
            check(sum(fs[s]) == 0, 'exact_centered_zero_mass')
            check(all(fs[s][a*x % p] == fs[s][x] for a in h for x in range(p)),
                  'centered_H_invariance')
            check(abs(fs[s][0]) <= n**(s-1), 'centered_origin_pointwise')
            check(fs[s][0]**4 <= n**(2*s-2)*ts[s], 'centered_origin_absorption')
            check(sum(abs(v) for v in fs[s]) <= 2*n**s, 'centered_L1_bound')
        for s in range(1, 7):
            check(conv(fs[s], fs[s]) == fs[2*s], 'exact_centered_convolution_identity')
            check(e0(fs[s]) == ts[2*s], 'exact_centered_fourth_moment')
            check(ts[2*s]**2*n <= (331784*n**(2*s)*ts[s])**2,
                  'centered_recurrence_K_one_sanity_only')
        data.append({'p': p, 'n': n, 'T_s': {str(s): str(t) for s, t in ts.items()}})
    return data


def fourier_checks():
    # Independent complex check of centering in the mixed-moment identity.
    worst_relative_error = 0.0
    for p, n in [(5, 4), (7, 3), (13, 4), (17, 8), (31, 5)]:
        h = subgroup(p, n)
        exp = [cmath.exp(2j*math.pi*t/p) for t in range(p)]
        eta = [sum(exp[a*x % p] for x in h) for a in range(p)]
        indicator = [int(x in h) for x in range(p)]
        r = [1]+[0]*(p-1)
        rows, energies = {}, {}
        for order in range(1, 5):
            r = conv(r, indicator)
            rows[order] = r
            energies[order] = float(F(norm2(r))-F(n**(2*order), p))
        for k, ell in product(range(2, 5), range(1, 5)):
            aa = [sum(abs(eta[t])**k*exp[-t*x % p] for t in range(p))/p for x in range(p)]
            ac = [v-n**k/p for v in aa]
            rc = [v-n**ell/p for v in rows[ell]]
            for a in [1, p-1]:
                lhs = sum(rows[ell][y]*abs(eta[a*y % p])**k for y in range(p))
                bilinear = sum(ac[x]*rc[y]*exp[a*x*y % p] for x in range(p) for y in range(p))
                rhs = bilinear+n**ell*aa[0]+n**k*rows[ell][0]-n**(k+ell)/p
                err = abs(lhs-rhs)/max(1.0, abs(lhs), abs(rhs))
                worst_relative_error = max(worst_relative_error, err)
                check(err < 1e-9, 'approximate_centered_mixed_identity')
                gate = 2/n+math.sqrt(p*energies[k]*energies[ell])/n**(k+ell)
                check((abs(eta[a])/n)**(k*ell) <= gate+1e-9,
                      'approximate_centered_mixed_inequality')
    return {'relative_tolerance': 1e-9, 'max_identity_relative_error': worst_relative_error,
            'scope': 'Floating-point diagnostics, not certified inequalities.'}


def exponent_and_boundary_checks():
    check(6**4 == 1296, 'layer_constant')
    check(8*1296 == 10368, 'sign_split_constant')
    check(8*10368 == 82944, 'origin_split_constant')
    check(4*82944+8 == 331784, 'recurrence_constant')
    mixed = [F(s, 36*2**s) for s in range(25)]
    order_two = [F(10*j-9, 240*2**j) for j in range(25)]
    direct = [F(j-2, 12*2**j) for j in range(25)]
    check(max(mixed) == F(1, 72), 'dyadic_pair_maximum')
    check(max(order_two) == F(11, 960), 'order_two_pair_maximum')
    check(max(direct) == F(1, 96), 'direct_coset_maximum')
    knots = {1: F(1), 2: F(31, 20), 3: F(2), 6: F(5, 2),
             12: F(3), 24: F(7, 2)}

    def envelope(s):
        s = F(s)
        if s >= 24:
            return s/36+F(17, 6)
        for lo, hi in zip(list(knots), list(knots)[1:]):
            if lo <= s <= hi:
                return knots[lo]+(knots[hi]-knots[lo])*(s-lo)/(hi-lo)
        raise ValueError(s)

    candidates = []
    orders = [2, 3, 6, 12, 24, None]
    for k, ell in product(orders, repeat=2):
        if k is None and ell is None:
            saving = F(0)
        elif k is None or ell is None:
            saving = F(1, 72*(k or ell))
        else:
            saving = (envelope(k)+envelope(ell)-4)/(2*k*ell)
        candidates.append({'k': k, 'l': ell, 'saving': str(saving)})
    check(max(F(c['saving']) for c in candidates) == F(1, 72),
          'feedback_mixed_endpoint_maximum')
    # These are rational spot checks supporting the interval proofs in the note.
    for s in [F(i, 2) for i in range(2, 121)]:
        check(envelope(2*s) >= envelope(s)+F(1, 2), 'feedback_doubling_closure')
        check(envelope(s)-s/36 <= F(17, 6), 'feedback_intercept_bound')
        check((envelope(s)-3)/(2*s) < F(1, 72), 'feedback_direct_saving_below_limit')
        for t in [F(1), F(3, 2), F(3), F(12), F(25)]:
            check(envelope(s+t) >= envelope(s)+t/36, 'feedback_amplitude_closure')
            check(envelope(s+t) >= envelope(s)+envelope(t)-3, 'feedback_coset_product_closure')
    # Prove the recurrence exponent by checking the closed-form sum identity.
    for s, j in product(range(1, 9), range(13)):
        check(sum((2*s*2**i-F(1, 2) for i in range(j)), F(0))
              == 2*s*(2**j-1)-F(j, 2), 'iteration_exponent_identity')
    origin_ratios = {}
    subfield_ratios = {}
    for p in [3, 5, 7, 11, 17, 31, 101]:
        f = [F(int(x == 0))-F(1, p) for x in range(p)]
        check(e0(f) == 1-F(1, p), 'origin_counterexample_energy')
        a, b = sum(abs(v) for v in f), norm2(f)
        ratio_sq = e0(f)**2*(p-1)/(a*a*b)**2
        check(ratio_sq == F(p-1)/(16*(1-F(1, p))**4), 'origin_counterexample_ratio')
        origin_ratios[str(p)] = str(ratio_sq)
        q = [0]+[1]*(p-1)
        m = p-1
        check(norm2(conv(q, q)) == m**3-m**2+m, 'subfield_energy_formula')
        factor = 1-F(1, m)+F(1, m*m)-F(m, p**6)
        subfield_ratios[str(p)] = str(m*factor*factor)
    return {'best_mixed_saving': str(max(mixed)), 'best_direct_saving': str(max(direct)),
            'feedback_endpoint_candidates': candidates,
            'feedback_saving_supremum': '1/72',
            'origin_counterexample_squared_ratios': origin_ratios,
            'subfield_squared_ratios': subfield_ratios,
            'scope': 'Finite arithmetic checks; infinite tail comparisons and divergence proved in prose.'}


def main():
    started = time.monotonic()
    result = {
        'incidence': incidence_checks(),
        'centered_convolutions': weight_and_convolution_checks(),
        'approximate_fourier': fourier_checks(),
        'exponents_and_boundaries': exponent_and_boundary_checks(),
    }
    result.update({
        'status': 'bounded checks passed; analytic theorem not formally verified; full goal unproved',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'elapsed_seconds': time.monotonic()-started,
        'check_counts': dict(COUNTS),
        'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in [
            'experiments/parallel29_verify_2026_09_05.py',
            'research/parallel29-centered-recurrence-2026-09-05.md',
            'sources/parallel28-moment-recurrence/shkredov-published-2019.pdf',
        ]},
        'limitations': ['K=1 checks apply only to the selected finite examples.',
                        'No new Lean process or hosted proof job was launched.',
                        'No independent reviewer completed this pass.'],
    })
    dest = ROOT/'results/parallel29_verification_2026_09_05.json'
    dest.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'path': str(dest), 'checks': sum(COUNTS.values()),
                      'elapsed_seconds': result['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
