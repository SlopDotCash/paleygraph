#!/usr/bin/env python3
"""Exact finite witnesses for the missing published incidence hypothesis.

This checks image cardinalities and weighted counting identities. It neither
refutes an asymptotic energy estimate nor proves a Paley conjecture.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import json
from math import isqrt
from pathlib import Path


def subgroup(p, generator, order):
    assert p > 2 and all(p % q for q in range(2, isqrt(p) + 1))
    h = {pow(generator, j, p) for j in range(order)}
    assert len(h) == order and pow(generator, order, p) == 1
    assert all(x * y % p in h for x in h for y in h)
    return h


def analyze(p, h, s):
    t = len(h)
    assert 0 not in s and s and all(x * z % p in s for x in s for z in h)

    @lru_cache(None)
    def coset(x):
        return min(x * z % p for z in h)

    c = {coset(x) for x in s}
    assert len(s) == t * len(c)
    image_fibers = Counter(
        (x * pow(z, -1, p) % p, y * pow(z, -1, p) % p)
        for x in s for y in s for z in s
    )
    multiplicities = Counter(
        (coset(a * pow(c0, -1, p) % p), coset(b * pow(c0, -1, p) % p))
        for a in c for b in c for c0 in c
    )
    rho = {
        (u, v): sum((u * x + v * y) % p == 1 for x in h for y in h)
        for u, v in multiplicities
    }
    direct_incidence = sum((x + y) % p in s for x in s for y in s)
    weighted_incidence = t * sum(m * rho[uv] for uv, m in multiplicities.items())
    unweighted_incidence = t * sum(rho.values())
    assert direct_incidence == weighted_incidence
    assert sum(multiplicities.values()) == len(c) ** 3
    assert len(image_fibers) == t * t * len(multiplicities)
    assert sum(image_fibers.values()) == len(s) ** 3
    assert all(m == t * multiplicities[coset(u), coset(v)]
               for (u, v), m in image_fibers.items())
    assert multiplicities[coset(1), coset(1)] == len(c)
    required = len(s) ** 3 // t
    holds = len(image_fibers) == required
    assert holds == all(m == 1 for m in multiplicities.values())
    assert holds == (len(c) == 1)  # All three sets are the same here.
    layer_sizes = [sum(m >= j for m in multiplicities.values())
                   for j in range(1, max(multiplicities.values()) + 1)]
    assert sum(layer_sizes) == sum(multiplicities.values())
    return {
        'p': p, 'h': sorted(h), 's': sorted(s), 't': t,
        'coset_representatives': sorted(c),
        'image_cardinality': len(image_fibers),
        'required_image_cardinality': required,
        'published_extra_hypothesis_holds': holds,
        'raw_fiber_histogram': dict(sorted(Counter(image_fibers.values()).items())),
        'quotient_multiplicity_histogram': dict(sorted(Counter(multiplicities.values()).items())),
        'quotient_layer_sizes': layer_sizes,
        'direct_incidence': direct_incidence,
        'weighted_incidence': weighted_incidence,
        'unweighted_incidence': unweighted_incidence,
    }


def main():
    p = 1153
    h = subgroup(p, 75, 8)
    r = Counter((x - y) % p for x in h for y in h)
    energy = sum(v * v for v in r.values())
    d = Fraction(energy, 16 * len(h) ** 2)
    levels = {}
    for x, v in r.items():
        if x == 0 or v <= d:
            continue
        i = 1
        while v > 2 ** i * d:
            i += 1
        levels.setdefault(i, set()).add(x)
    assert energy == 168 and d == Fraction(21, 128)
    examples = []
    for i, s in sorted(levels.items()):
        result = analyze(p, h, s)
        result.update(case='old_proof_dyadic_level', level=i,
                      difference_multiplicities=sorted({r[x] for x in s}),
                      energy=energy, cutoff=str(d))
        examples.append(result)
    assert examples[0]['published_extra_hypothesis_holds']
    assert examples[1]['image_cardinality'] == 1600
    assert examples[1]['required_image_cardinality'] == 1728
    assert not examples[1]['published_extra_hypothesis_holds']

    # A nested subgroup also witnesses an actual loss in the incidence count
    # when the normalization multiplicity is omitted.
    p = 97
    h = subgroup(p, pow(5, 12, p), 8)
    s = subgroup(p, pow(5, 4, p), 24)
    result = analyze(p, h, s)
    result['case'] = 'nested_subgroups'
    assert result['quotient_multiplicity_histogram'] == {3: 9}
    assert result['weighted_incidence'] == 3 * result['unweighted_incidence']
    assert result['unweighted_incidence'] > 0
    examples.append(result)
    assert Fraction(22, 9) < Fraction(49, 20) < Fraction(32, 13)
    assert Fraction(32, 13) - Fraction(49, 20) == Fraction(3, 260)
    report = {
        'scope': 'exact finite counting; missing-hypothesis witness only',
        'full_conjecture_proved': False,
        'asymptotic_22_over_9_refuted': False,
        'examples': examples,
    }
    root = Path(__file__).resolve().parents[1]
    target = root / 'results/parallel38_coset_multiplicity_2026_09_06.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    for ex in examples:
        print(ex['case'], 'p=', ex['p'], '|S|=', len(ex['s']),
              'image=', ex['image_cardinality'],
              'required=', ex['required_image_cardinality'],
              'I=', ex['direct_incidence'], 'unweighted=', ex['unweighted_incidence'])
    print(target)


if __name__ == '__main__':
    main()
