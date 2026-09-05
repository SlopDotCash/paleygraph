#!/usr/bin/env python3
"""Independent finite codeword/cover oracles for compressed_support.py."""
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
LANE = HERE.parent/'proximity'
sys.path.insert(0, str(LANE))
import compressed_support as tool


def arithmetic(q):
    if q == 9:
        # Independent GF9 = F3[a]/(a²+1), encoded x+3y.
        def add(x, y):
            return (x % 3+y % 3) % 3+3*((x//3+y//3) % 3)
        def mul(x, y):
            a, b, c, d = x % 3, x//3, y % 3, y//3
            return (a*c-b*d) % 3+3*((a*d+b*c) % 3)
        return add, mul
    return lambda a, b: (a+b) % q, lambda a, b: a*b % q


def expected_nodes(q, dom, k, s, u0, u1):
    add, mul = arithmetic(q)
    code = []
    for coeff in product(range(q), repeat=k):
        values = []
        for x in dom:
            y = 0
            for c in reversed(coeff):
                y = add(mul(y, x), c)
            values.append(y)
        code.append(tuple(values))
    answer = set()
    for z in range(q):
        word = [add(a, mul(z, b)) for a, b in zip(u0, u1)]
        for cw in code:
            support = tuple(i for i, (a, b) in enumerate(zip(word, cw)) if a == b)
            if len(support) >= s:
                answer.add((z, cw, support))
    return answer


def expand_record(record):
    q = record['field']['q']
    add, mul = arithmetic(q)
    out = {(r['scalar'], tuple(r['codeword']), tuple(r['agreement_support'])) for r in record['nodes']}
    for track in record['all_field_track_certificates']:
        a, b = track['intercept'], track['slope']
        for z in range(q):
            cw = tuple(add(x, mul(z, y)) for x, y in zip(a, b))
            word = tuple(add(x, mul(z, y)) for x, y in zip(record['u0'], record['u1']))
            support = tuple(i for i, (x, y) in enumerate(zip(cw, word)) if x == y)
            assert len(support) >= record['s']
            out.add((z, cw, support))
    return out


def cover_checks():
    cases = supports = 0
    for n in range(1, 10):
        for k in range(1, n+1):
            for s in range(k, n+1):
                for mode in ('partition', 'anchor', 'full_k'):
                    plan = tool.cover_plan(n, k, s, mode)
                    # Check the universal sufficient condition independently.
                    occupied = set()
                    for block in plan['blocks']:
                        assert not occupied.intersection(block)
                        occupied.update(block)
                    assert n-len(occupied)+sum(min(len(b), k-1) for b in plan['blocks']) < s
                    for support in combinations(range(n), s):
                        assert any(len(set(support).intersection(block)) >= k for block in plan['blocks'])
                        supports += 1
                    cases += 1
    return {'plans': cases, 'supports': supports}


def run_one(q, dom, k, s, u0, u1):
    F = tool.Field(3, [1, 0, 1]) if q == 9 else tool.PrimeField(q)
    expected = expected_nodes(q, dom, k, s, u0, u1)
    checks = symbolic = 0
    for mode in ('partition', 'anchor', 'full_k'):
        for limit in (0, q):
            result = tool.run_census(F, dom, k, s, u0, u1, mode, materialize_limit=limit)
            actual = expand_record(result)
            assert actual == expected, (q, k, s, mode, limit, len(actual), len(expected))
            expected_scalars = {r[0] for r in expected}
            assert result['finite_bad_scalar_count'] == len(expected_scalars)
            if result['node_ledger_materialized']:
                explicit = {(r['scalar'], tuple(r['codeword']), tuple(r['agreement_support'])) for r in result['nodes']}
                assert explicit == expected
            else:
                assert result['whole_field_correlated'] and result['bad_scalars'] is None
                symbolic += 1
            checks += 1
    return checks, symbolic


def main():
    covers = cover_checks()
    counts = {'distinct_pencils': 0, 'census_oracle_equalities': 0, 'symbolic_all_field_cases': 0}
    words = list(product(range(3), repeat=3))
    for u0 in words:
        for u1 in words:
            a, b = run_one(3, [0, 1, 2], 1, 2, u0, u1)
            counts['distinct_pencils'] += 1
            counts['census_oracle_equalities'] += a
            counts['symbolic_all_field_cases'] += b
    rng = random.Random(20260905)
    for q, k, s, samples in ((5, 2, 3, 60), (5, 2, 2, 12), (9, 2, 3, 20), (9, 2, 2, 6)):
        dom = [0, 1, 2, 3]
        for _ in range(samples):
            u0, u1 = [[rng.randrange(q) for _ in dom] for _ in range(2)]
            a, b = run_one(q, dom, k, s, u0, u1)
            counts['distinct_pencils'] += 1
            counts['census_oracle_equalities'] += a
            counts['symbolic_all_field_cases'] += b
    result = {'date': '2026-09-05', 'passed': True, 'cover_checks': covers, **counts,
              'scope': 'Complete tiny polynomial/scalar oracles, all three cover modes and both symbolic/materialized branches; no production theorem',
              'reviewed_source_sha256': sha256((LANE/'compressed_support.py').read_bytes()).hexdigest(),
              'reviewer_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    HERE.joinpath('compressed_review.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
