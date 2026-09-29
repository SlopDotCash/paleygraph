#!/usr/bin/env python3
"""Exact checks for the bounded pass23 prize-bridge investigation."""
from collections import Counter
from itertools import combinations, product
from fractions import Fraction
from math import comb, isqrt
from pathlib import Path
from hashlib import sha256
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def check(test, key):
    assert test, key
    COUNTS[key] += 1


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def domain(p, n):
    assert prime(p) and (p-1) % n == 0
    for z in range(1, p):
        vals = sorted({pow(z, j, p) for j in range(n)})
        if len(vals) == n and pow(z, n, p) == 1:
            return vals
    raise AssertionError('subgroup generator')


def ev(c, x, p):
    y = 0
    for v in reversed(c):
        y = (y*x+v) % p
    return y


def rootpoly(T, p):
    q = [1]
    for x in T:
        new = [0]*(len(q)+1)
        for j, v in enumerate(q):
            new[j] = (new[j]-x*v) % p
            new[j+1] = (new[j+1]+v) % p
        q = new
    return q


def rem(c, q, p):
    vals = list(c)
    s = len(q)-1
    for i in range(len(vals)-1, s-1, -1):
        v = vals[i]
        for j, w in enumerate(q):
            vals[i-s+j] = (vals[i-s+j]-v*w) % p
    return tuple((vals+[0]*s)[:s])


def newton(c, p):
    ans = []
    for j in range(1, len(c)+1):
        ans.append((-j*c[j-1]-sum(c[i-1]*ans[j-i-1] for i in range(1, j))) % p)
    return tuple(ans)


def coefficient_lists():
    reports = []
    for p, n in ((13, 12), (17, 8)):
        D = domain(p, n)
        k = 2
        code = [(c, tuple(ev(c, x, p) for x in D)) for c in product(range(p), repeat=k)]
        for a in (1, 2):
            s = k+a
            hist = Counter()
            for T in combinations(D, s):
                q = rootpoly(T, p)
                c = tuple(q[s-j] for j in range(1, a+1))
                moments = tuple(sum(pow(x, j, p) for x in T) % p for j in range(1, a+1))
                check(newton(c, p) == moments, 'Newton_subset_identities')
                hist[c] += 1
            targets = {(0,)*a, (1,)*a, tuple(range(1, a+1))}
            for c in targets:
                W = [0]*(s+1)
                W[s] = 1
                for j, v in enumerate(c, 1):
                    W[s-j] = v
                received = tuple(ev(W, x, p) for x in D)
                direct = sum(sum(u == v for u, v in zip(received, vals)) >= s for _, vals in code)
                check(direct == hist[c], 'monic_list_classifications')
            reports.append({'p': p, 'n': n, 'k': k, 'a': a, 'subsets': comb(n, s)})
    p, D, s = 17, list(range(1, 17)), 8
    hist = Counter((sum(T) % p, sum(x*x for x in T) % p) for T in combinations(D, s))
    check(hist[(0, 0)] == 54, 'false_linear_bound_fibre')
    check(p*p*54-comb(16, 8) == 2736 > p*p-1, 'false_linear_bound_deviation')
    for a, U in ((1, 1), (2, 5)):
        p, D, s = 13, list(range(1, 13)), 6
        hist = Counter(tuple(sum(pow(x, j, p) for x in T) % p for j in range(1, a+1))
                       for T in combinations(D, s))
        for target in product(range(p), repeat=a):
            check(abs(p**a*hist[target]-comb(12, 6)) <= (p**a-1)*comb(U+5, 6),
                  'finite_multivariate_count_bounds')
    return reports


def extension_centers():
    p = 5
    def add(x, y): return ((x[0]+y[0]) % p, (x[1]+y[1]) % p)
    def mul(x, y): return ((x[0]*y[0]+2*x[1]*y[1]) % p, (x[0]*y[1]+x[1]*y[0]) % p)
    def evaluate(c, x):
        v = (0, 0)
        for z in reversed(c): v = add(mul(v, (x, 0)), z)
        return v
    theta = (0, 1)
    # The first center has a non-base-field top coefficient; the second only
    # has an extension-valued lower coefficient, which the code can absorb.
    fixtures = [([(0, 0), theta, (1, 0)], 0), ([theta, (0, 0), (1, 0)], 2)]
    reports = []
    for W, expected in fixtures:
        words = [evaluate(W, x) for x in range(1, 5)]
        constants = [c for c in product(range(p), repeat=2) if sum(c == y for y in words) >= 2]
        check(len(constants) == expected, 'F25_extension_center_lists')
        reports.append({'center': W, 'constant_codewords': constants})
    return reports


def remainder_tests():
    reports = []
    for p, n, k in ((5, 4, 1), (7, 6, 2)):
        D = domain(p, n)
        code = [(c, tuple(ev(c, x, p) for x in D)) for c in product(range(p), repeat=k)]
        for width in (1, 2):
            words = [tuple(tuple((seed*(j+1)**2+i*(j+2)) % p for j in range(n))
                           for i in range(width)) for seed in (0, 1, 2)]
            for s in (k+1, k+2):
                Ts = list(combinations(range(n), s))
                qs = [rootpoly([D[i] for i in T], p) for T in Ts]
                restricted_codes = [{tuple(vals[i] for i in T) for _, vals in code} for T in Ts]
                rems = [[tuple(rem(c, q, p) for c in W) for q in qs] for W in words]
                values = [tuple(tuple(ev(c, x, p) for x in D) for c in W) for W in words]
                for wi, W in enumerate(words):
                    predicted = {tuple(row[:k] for row in rr) for rr in rems[wi]
                                 if all(not any(row[k:]) for row in rr)}
                    direct = set()
                    for chosen in product(code, repeat=width):
                        agreements = sum(all(chosen[r][1][x] == values[wi][r][x]
                                             for r in range(width)) for x in range(n))
                        if agreements >= s:
                            direct.add(tuple(item[0] for item in chosen))
                    check(predicted == direct, 'distinct_remainder_list_equalities')
                for wi in range(len(words)):
                    for zi in range(len(words)):
                        predicted = set()
                        for rrW, rrZ in zip(rems[wi], rems[zi]):
                            a = tuple(v for row in rrW for v in row[k:])
                            b = tuple(v for row in rrZ for v in row[k:])
                            pos = next((i for i, v in enumerate(b) if v), None)
                            if pos is not None:
                                gamma = -a[pos]*pow(b[pos], -1, p) % p
                                if all((u+gamma*v) % p == 0 for u, v in zip(a, b)):
                                    predicted.add(gamma)
                        direct = set()
                        for gamma in range(p):
                            for T, allowed in zip(Ts, restricted_codes):
                                fold = all(tuple((values[wi][r][x]+gamma*values[zi][r][x]) % p for x in T)
                                           in allowed for r in range(width))
                                joint = all(tuple(values[ii][r][x] for x in T) in allowed
                                            for ii in (wi, zi) for r in range(width))
                                if fold and not joint:
                                    direct.add(gamma)
                        check(predicted == direct, 'nontrivial_remainder_MCA_equalities')
                reports.append({'p': p, 'n': n, 'k': k, 'rows': width, 's': s, 'word_fixtures': len(words)})
    return reports


def box_collision():
    p, n, a, R = 17, 4, 3, 8
    D = domain(p, n)
    assert R**n < p**a
    boxes = {}
    for index, coeff in enumerate(product(range(p), repeat=a), 1):
        vals = tuple(ev((0,)+coeff, x, p) for x in D)
        box = tuple(R*v//p for v in vals)
        if box in boxes:
            prev, prev_vals = boxes[box]
            P = tuple((u-v) % p for u, v in zip(coeff, prev))
            differences = tuple(u-v for u, v in zip(vals, prev_vals))
            check(any(P) and any(differences), 'nonzero_actual_box_polynomial')
            check(all(R*abs(d) < p for d in differences), 'exact_phase_arc_bounds')
            check(all(ev((0,)+P, x, p) == d % p for x, d in zip(D, differences)),
                  'box_polynomial_evaluation_identity')
            return {'p': p, 'D': D, 'a': a, 'R': R, 'tested_vectors': index,
                    'zero_constant_polynomial_coefficients': [0]+list(P),
                    'small_integer_evaluations': differences,
                    'real_sum_lower_bound': str(Fraction(n)*(1-Fraction(20, R*R)))}
        boxes[box] = (coeff, vals)
    raise AssertionError('pigeonhole collision not found')


def official():
    folder = ROOT/'sources/official-prize-2026-09-04'
    manifest = json.loads((folder/'manifest.json').read_text())
    files = manifest['files'] + [f for dep in manifest['dependencies'] for f in dep['files']]
    for f in files:
        check(sha256((folder/f['path']).read_bytes()).hexdigest() == f['sha256'], 'official_archived_source_hashes')
    p, n, k, a, R = 2130706433, 262144, 131072, 26215, 8
    check(prime(p) and (p-1) % n == 0, 'official_base_prime_and_domain_divisibility')
    check(p > 2**30 and 30*a-3*n == 18 > 0, 'official_box_cardinality_certificate')
    check(k+a == 157287 and n-k-a == 104857, 'official_degree_and_radius')
    check(Fraction(1)-Fraction(20, R*R) == Fraction(11, 16), 'official_real_phase_lower_fraction')
    budget = p**6//2**128
    check(budget == 274980728111395087, 'official_combined_integer_budget')
    check(budget*2**128 <= p**6 < (budget+1)*2**128, 'official_budget_floor')
    return {'p': p, 'q': p**6, 'n': n, 'k': k, 'a': a, 's': k+a,
            'radius': [n-k-a, n], 'R': R, 'real_phase_sum_strict_lower': 180224,
            'combined_MCA_and_list_budget': budget,
            'contract_commit': manifest['commit'], 'current_status_claim': False}


def main():
    result = {'status': 'Exact checks passed; the full prize bridge remains unproved.',
              'coefficient_list_fixtures': coefficient_lists(), 'extension_centers': extension_centers(),
              'general_remainder_fixtures': remainder_tests(), 'actual_box_collision': box_collision(),
              'pinned_official_parameters': official()}
    inputs = ['research/parallel23-prize-bridge-2026-09-05.md',
              'experiments/parallel23_prize_bridge_2026_09_05.py',
              'research/official-profile-and-trace.md', 'research/subgroup-target.md',
              'research/subset-sums-and-lists.md', 'research/prize-reduction-audit.md',
              'sources/official-prize-2026-09-04/manifest.json']
    result['input_sha256'] = {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs}
    result['counts'] = dict(COUNTS)
    result['limitations'] = ['No production-sized list or MCA maximum was evaluated.',
                            'The polynomial-phase obstruction does not contradict linear subgroup cancellation.',
                            'The obstruction does not prove a list or soundness lower bound.',
                            'No current prize status or formal verification is claimed.']
    output = ROOT/'results/parallel23_prize_bridge_2026_09_05.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'output': str(output), 'counts': dict(COUNTS),
                      'box': result['actual_box_collision']}, indent=2))


if __name__ == '__main__':
    main()
