#!/usr/bin/env python3
"""Exact elementary covering lower bounds; no claim about arbitrary decoders."""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    results_path = HERE / 'results.json'
    raw = json.loads(results_path.read_text())
    finite = []
    for x in raw['results']:
        n, k, s = x['n'], x['k'], x['s']
        volumes = []
        for m, num in x['candidate_support_sizes'].items():
            m = int(m)
            volume = sum(math.comb(m, i) * math.comb(n - m, s - i)
                         for i in range(k, min(m, s) + 1) if 0 <= s - i <= n - m)
            volumes.extend([volume] * num)
        total = 0
        for i, volume in enumerate(sorted(volumes, reverse=True)):
            total += volume
            if total >= math.comb(n, s):
                break
        M = x['maximal_joint_support_size']
        full, remainder = divmod(n, M)
        cap = full * min(M, k - 1) + min(remainder, k - 1)
        finite.append({'configuration': x['configuration'], 'sample': x['sample'],
                       'catalog_aware_volume_lower_bound': i + 1,
                       'exact_integer_partition_cap_minimum': cap,
                       'strict_partition_cap_fails': cap >= s,
                       'caveat': 'Size-only cap lower bound; actual attainment is separately certified by cover_barriers.'})
    R, alpha = .5, 11 / 16
    entropy = lambda t: -t * math.log2(t) - (1 - t) * math.log2(1 - t)
    exponent = entropy(R) - alpha * entropy(R / alpha)
    scale = []
    for n in [16, 32, 64, 128, 256, 512, 1024]:
        k, s = n // 2, 11 * n // 16
        num, den = math.comb(n, k), math.comb(s, k)
        scale.append({'n': n, 'k': k, 's': s, 'exact_integer_lower_bound': (num + den - 1) // den,
                      'log2_counting_ratio': math.log2(num) - math.log2(den)})
    output = {'finite_bounds': finite,
              'pure_k_support_fixed_rate': {'R': R, 'alpha': alpha,
                  'entropy_exponent_bits_per_coordinate': exponent,
                  'formula': 'H2(R)-alpha*H2(R/alpha)>0 for 0<R<alpha<1', 'finite_instances': scale},
              'scope': 'Pure k-support lower bound is conditional on that structural regime; no saved hard-pencil sequence at growing n is claimed.',
              'source_bindings': [{'path': p.name, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                                  for p in [Path(__file__).resolve(), results_path]]}
    (HERE / 'counting_bounds.json').write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
