#!/usr/bin/env python3
"""Bounded exploratory norm search for a quartic subgroup with a zero triangle.

Factorization is discovery only. Every reported field prime is separately
verified by exhaustive trial division and its triangle is checked directly.
No completeness of the prime search or arbitrary-size primality claim.
"""
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt, prod
from pathlib import Path
import importlib.util
import json
import random
import time

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('prior_norm', ROOT/'experiments/parallel25_subgroup_growing_orders_2026_09_05.py')
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)
RNG = random.Random(300905)


def probable_prime(v):
    if v < 2:
        return False
    for b in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if v % b == 0:
            return v == b
    d, s = v-1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        x = pow(b, d, v)
        if x in [1, v-1]:
            continue
        for _ in range(s-1):
            x = x*x % v
            if x == v-1:
                break
        else:
            return False
    return True


def factor(v):
    if v == 1:
        return []
    if probable_prime(v):
        return [v]
    for small in [2, 3, 5, 7, 11, 13]:
        if v % small == 0:
            return [small]+factor(v//small)
    for attempt in range(12):
        c = RNG.randrange(1, v)
        x = y = RNG.randrange(2, v)
        for _ in range(100000):
            x = (x*x+c) % v
            y = (y*y+c) % v
            y = (y*y+c) % v
            d = gcd(abs(x-y), v)
            if d == v:
                break
            if d > 1:
                return factor(d)+factor(v//d)
    raise RuntimeError(('bounded factor search exhausted', v))


def main():
    started = time.monotonic()
    n = 128
    rows, found = [], {}
    cache = {}
    for m in [4, 8, 16, 32, 64, 128]:
        for b in range(2, m):
            if b in [m//2, m//2+1]:
                continue
            coefficients = [0]*(m//2)
            for e in [0, 1, b]:
                coefficients[e % (m//2)] += 1 if e < m//2 else -1
            norm = MOD.norm_recursive(coefficients)
            assert norm > 0
            if norm not in cache:
                cache[norm] = factor(norm)
            factors = cache[norm]
            assert prod(factors) == norm
            rows.append({'m': m, 'b': b, 'norm': norm, 'discovery_factors': dict(Counter(factors))})
            for p in factors:
                if p in found or not (p % n == 1 and n**4//4 <= p <= n**4):
                    continue
                assert all(p % d for d in range(2, isqrt(p)+1))
                g = next(pow(a, (p-1)//n, p) for a in range(2, p)
                         if pow(pow(a, (p-1)//n, p), n//2, p) != 1)
                h = sorted({pow(g, j, p) for j in range(n)})
                assert len(h) == n and pow(g, n, p) == 1
                pairs = [(x, (1-x) % p) for x in h if (1-x) % p in h]
                tau = len(pairs)-3*int(2 in h)
                assert tau > 0 and tau % 6 == 0
                triangle = next([1, -x % p, -y % p] for x, y in pairs
                                if len({1, -x % p, -y % p}) == 3)
                assert sum(triangle) % p == 0 and all(x in h for x in triangle)
                found[p] = {'p': p, 'n': n, 'generator': g, 'triangle': triangle,
                            'kappa': len(pairs), 'tau': tau,
                            'predicted_D6_lower_bound': 10*n*(n-28)*tau*tau,
                            'n_cubed': n**3, 'discovery_m': m, 'discovery_b': b}
                print(json.dumps(found[p]), flush=True)
    result = {'n': n, 'scope': 'Exploratory candidate search; no completeness claim.',
              'rows': rows, 'found': list(found.values()), 'elapsed_seconds': time.monotonic()-started,
              'norm_routine_sha256': sha256((ROOT/'experiments/parallel25_subgroup_growing_orders_2026_09_05.py').read_bytes()).hexdigest()}
    (ROOT/'results/parallel30_triangle_search_2026_09_05.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'found_count': len(found), 'norm_rows': len(rows),
                      'distinct_norms': len(cache), 'elapsed_seconds': result['elapsed_seconds']}))


if __name__ == '__main__':
    main()
