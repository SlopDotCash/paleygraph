#!/usr/bin/env python3
"""Exact finite tests of triangle-generated D6 and its explicit counterexample.

No asymptotic assertion, probabilistic primality certificate, or Lean proof.
The direct executable is built from the separately enumerating C++ source.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import isqrt
from pathlib import Path
import importlib.util
import json
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
CHECKS = Counter()


def check(ok, label):
    assert ok, label
    CHECKS[label] += 1


def balanced(word, p):
    for tail in combinations(range(1, 6), 2):
        left = {0, *tail}
        a = b = 1
        for i, x in enumerate(word):
            if i in left:
                a = a*x % p
            else:
                b = b*x % p
        if (a+b) % p == 0:
            return True
    return False


def classify(a, b, p):
    word = tuple(a)+tuple(b)
    if len(set(word)) != 6:
        return 'overlap'
    if any(-x % p in word for x in word):
        return 'opposite'
    if balanced(word, p):
        return 'balanced'
    return 'good'


def group(p, n):
    check(p > 3 and all(p % d for d in range(2, isqrt(p)+1)), 'trial_division_prime')
    check(n >= 4 and n & (n-1) == 0 and (p-1) % n == 0, 'dyadic_splitting_order')
    g = next(pow(a, (p-1)//n, p) for a in range(2, p)
             if pow(pow(a, (p-1)//n, p), n//2, p) != 1)
    h = sorted({pow(g, j, p) for j in range(n)})
    check(len(h) == n and pow(g, n, p) == 1, 'exact_subgroup_order')
    return g, h


def inspect(p, n):
    g, h = group(p, n)
    hs = set(h)
    kappa = sum((1-x) % p in hs for x in h)
    tau = kappa-3*int(2 in hs)
    triangles = sorted({tuple(sorted((x, y, (-x-y) % p)))
                        for x in h for y in h if (-x-y) % p in hs
                        and len({x, y, (-x-y) % p}) == 3})
    check(len(triangles)*6 == n*tau and tau % 6 == 0, 'exact_triangle_count')
    check(all(sum(t) % p == 0 and len(t) == 3 for t in triangles), 'actual_zero_triangles')
    canonical = lambda t: min(tuple(sorted(s*x % p for x in t)) for s in h)
    representatives = sorted({canonical(t) for t in triangles})
    check(len(representatives)*6 == tau, 'triangle_orbit_count')
    for a in representatives:
        check(len({tuple(sorted(s*x % p for x in a)) for s in h}) == n, 'free_triangle_scaling')
    relative = Counter()
    diagonal = []
    for a, b in product(representatives, repeat=2):
        counts = Counter(classify(a, [t*x % p for x in b], p) for t in h)
        check(counts['good'] >= max(0, n-28), 'relative_pair_lower_bound')
        relative.update(counts)
        if a == b:
            ratios = {x*pow(y, -1, p) % p for x in a for y in a}
            bad = ratios | {-x % p for x in ratios} | {-x*x % p for x in ratios}
            actual_bad = {t for t in h if classify(a, [t*x % p for x in a], p) != 'good'}
            check(bad == actual_bad, 'exact_diagonal_forbidden_ratios')
            check(len(bad) <= 20, 'diagonal_bad_ratio_bound')
            diagonal.append({'representative': a, 'bad_count': len(bad), 'good_count': n-len(bad)})
    union_count = Counter()
    for a, b in product(triangles, repeat=2):
        if classify(a, b, p) == 'good':
            union_count[tuple(sorted(a+b))] += 1
    check(all(c == 2 for c in union_count.values()), 'unique_zero_triple_partition')
    for word in union_count:
        zeros = sum(sum(word[i] for i in (0, *tail)) % p == 0
                    for tail in combinations(range(1, 6), 2))
        check(zeros == 1, 'direct_zero_partition_count')
    triangle_d6 = 720*len(union_count)
    check(triangle_d6 == 360*n*relative['good'], 'relative_and_absolute_pair_counts_agree')
    check(10*n*max(0, n-28)*tau*tau <= triangle_d6 <= 10*n*n*tau*tau,
          'general_triangle_D6_bounds')
    if representatives:
        check(triangle_d6 >= 360*n*max(0, n-20), 'single_orbit_D6_lower_bound')
    return {'p': p, 'n': n, 'generator': g, 'kappa': kappa, 'tau': tau,
            'triangles': len(triangles), 'triangle_scaling_orbits': len(representatives),
            'relative_pair_categories': dict(relative), 'diagonal': diagonal,
            'triangle_D6': triangle_d6, 'triangle_D6_scaling_orbits': len(union_count)//n}


def witness():
    p, n, g = 215535361, 128, 25525303
    check(n**4//4 <= p <= n**4, 'witness_quartic_window')
    check(pow(g, 64, p) == p-1 and pow(g, 128, p) == 1, 'witness_generator_orders')
    a = [1, g, pow(g, 19, p)]
    check(a == [1, 25525303, 190010057] and sum(a) == p, 'witness_zero_triangle')
    exponents = {0, 1, 19}
    ratio_exponents = {(x-y) % n for x in exponents for y in exponents}
    bad = ratio_exponents | {(e+64) % n for e in ratio_exponents} | {(2*e+64) % n for e in ratio_exponents}
    check(len(bad) == 20, 'witness_twenty_forbidden_scales')
    good = sorted(set(range(n))-bad)
    for e in good:
        check(classify(a, [pow(g, e, p)*x % p for x in a], p) == 'good',
              'witness_each_good_scale')
    analytic_lower = 360*n*(n-20)
    check(analytic_lower > n**3, 'analytic_counterexample_without_full_enumeration')
    # Independent integer determinant for Norm(1+X+X^19) in degree64.
    spec = importlib.util.spec_from_file_location('determinant', ROOT/'experiments/cyclotomic_norm_audit.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    d = n//2
    matrix = [[0]*d for _ in range(d)]
    for j in range(d):
        for e in [0, 1, 19]:
            q, r = divmod(e+j, d)
            matrix[r][j] += (-1)**q
    determinant = abs(module.determinant(matrix))
    check(determinant == p, 'independent_cyclotomic_norm_determinant')
    return {'p': p, 'n': n, 'generator': g, 'triangle': a,
            'ratio_exponents': sorted(ratio_exponents), 'bad_scale_exponents': sorted(bad),
            'good_scale_exponents': good, 'analytic_D6_lower_bound': analytic_lower,
            'n_cubed': n**3, 'norm_determinant': determinant}


def main():
    started = time.monotonic()
    cases = [inspect(p, n) for p, n in [(5, 4), (17, 16), (97, 32), (257, 64), (215535361, 128)]]
    cert = witness()
    with tempfile.TemporaryDirectory(prefix='paley-pass30-') as tmp:
        exe = str(Path(tmp)/'direct')
        subprocess.run(['/usr/bin/clang++', '-O2', '-std=c++17', '-Wall', '-Wextra',
                        str(ROOT/'experiments/parallel30_direct_six.cpp'), '-o', exe], check=True)
        run = subprocess.run([exe], input=''.join(f"{r['p']} {r['n']}\n" for r in cases),
                             capture_output=True, text=True, check=True)
    direct = [json.loads(line) for line in run.stdout.splitlines()]
    for row, independent in zip(cases, direct, strict=True):
        check(row['triangle_D6'] == independent['triangle_D6'], 'independent_Cpp_triangle_D6')
    old = json.loads((ROOT/'results/parallel30_sparse_witness_2026_09_05.json').read_text())
    for key in ['E3', 'T6', 'J6', 'repeated_R6', 'distinct_balanced_R6', 'distinct_unbalanced_R6']:
        check(old[key] == direct[-1][key], 'independent_Cpp_full_remainder')
    check(direct[-1]['primitive_D6'] > 128**3, 'primitive_remainder_also_exceeds_n_cubed')
    paths = ['experiments/parallel30_verify_2026_09_05.py', 'experiments/parallel30_direct_six.cpp',
             'experiments/cyclotomic_norm_audit.py', 'results/parallel30_sparse_witness_2026_09_05.json']
    result = {'status': 'exact finite checks passed; no uniform upper bound or full conjecture proof',
              'cases': cases, 'short_certificate': cert, 'independent_direct_counts': direct,
              'witness_D6_ratio': str(Fraction(direct[-1]['distinct_unbalanced_R6'], 128**3)),
              'check_counts': dict(CHECKS), 'elapsed_seconds': time.monotonic()-started,
              'source_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}
    (ROOT/'results/parallel30_verification_2026_09_05.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': sum(CHECKS.values()), 'elapsed_seconds': result['elapsed_seconds'],
                      'witness': direct[-1], 'short_certificate': cert}, indent=2))


if __name__ == '__main__':
    main()
