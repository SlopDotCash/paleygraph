#!/usr/bin/env python3
"""Independent certificate readback: no constructor or scalar-bucket reuse.

The checker evaluates polynomials by powers, uses literal Python-set
intersections, and scans every field scalar and selected track. It does not
trust the candidate catalog, greedy run, stored covering witness, or node list.
"""
import copy
import hashlib
import json
import random
import sys
import time
from collections import Counter
from itertools import combinations, product
from math import comb, isqrt
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]


def values(coefficients, domain, p):
    return tuple(sum(a * pow(x, j, p) for j, a in enumerate(coefficients)) % p for x in domain)


def verify(record):
    p, n, k, s = [record[key] for key in ['p', 'n', 'k', 's']]
    assert all(type(a) is int for a in [p, n, k, s])
    assert 2 <= p <= 5000 and all(p % d for d in range(2, isqrt(p) + 1))
    assert 1 <= k <= s <= n and comb(n, s) <= 200000
    domain, u0, u1 = [record[key] for key in ['domain', 'u0', 'u1']]
    assert len(domain) == len(u0) == len(u1) == n and len(set(domain)) == n
    assert all(type(a) is int and 0 <= a < p for a in domain + u0 + u1)
    assert record['event_definition'] == 'plain at-least-s agreement with degree-less-than-k polynomials'
    supports, evals, unique_tracks = [], [], set()
    for tr in record['tracks']:
        aa, bb = tr['intercept_coefficients'], tr['slope_coefficients']
        assert len(aa) == len(bb) == k
        assert all(type(a) is int and 0 <= a < p for a in aa + bb)
        key = tuple(aa), tuple(bb)
        assert key not in unique_tracks
        unique_tracks.add(key)
        av, bv = values(aa, domain, p), values(bb, domain, p)
        support = [i for i in range(n) if av[i] == u0[i] and bv[i] == u1[i]]
        assert support == tr['joint_support'] and len(support) >= k
        assert len(tr['base']) == k and len(set(tr['base'])) == k
        assert set(tr['base']).issubset(support)
        supports.append(set(support))
        evals.append((av, bv))
    multiplicities, actual_witnesses = Counter(), []
    for aa in combinations(range(n), s):
        hits = [j for j, support in enumerate(supports) if len(set(aa).intersection(support)) >= k]
        assert hits, ('uncovered', aa)
        multiplicities[len(hits)] += 1
        actual_witnesses.append(hits[0])
    assert record['coverage_count'] == comb(n, s)
    assert record['lexicographic_s_set_track_witnesses'] == actual_witnesses
    assert {int(a): b for a, b in record['coverage_multiplicity_histogram'].items()} == dict(multiplicities)
    actual, representations = {}, {}
    for z in range(p):
        received = tuple((a + z * b) % p for a, b in zip(u0, u1))
        for j, (av, bv) in enumerate(evals):
            cw = tuple((a + z * b) % p for a, b in zip(av, bv))
            support = tuple(i for i, (a, b) in enumerate(zip(cw, received)) if a == b)
            if len(support) >= s:
                aa, bb = record['tracks'][j]['intercept_coefficients'], record['tracks'][j]['slope_coefficients']
                coeff = tuple((a + z * b) % p for a, b in zip(aa, bb))
                key = z, cw
                actual[key] = coeff, support
                representations.setdefault(key, []).append(j)
    exported = {}
    for r in record['nodes']:
        key = r['scalar'], tuple(r['codeword'])
        assert key not in exported
        exported[key] = tuple(r['coefficients']), tuple(r['agreement_support'])
        assert r['represented_by_tracks'] == representations.get(key)
    assert exported == actual and record['node_count'] == len(actual)
    assert record['whole_field_track_ids'] == [j for j, support in enumerate(supports) if len(support) >= s]
    return {'tracks': len(evals), 'covered_s_sets': len(actual_witnesses), 'nodes': len(actual),
            'scalar_track_pairs_checked': p * len(evals),
            'coverage_minimum_multiplicity': min(multiplicities),
            'whole_field_tracks': len(record['whole_field_track_ids'])}


def literal_codebook(record):
    p, k, s = [record[key] for key in ['p', 'k', 's']]
    output = set()
    for h in product(range(p), repeat=k):
        cw = values(h, record['domain'], p)
        for z in range(p):
            support = tuple(i for i, (a, b, c) in enumerate(zip(cw, record['u0'], record['u1']))
                            if a == (b + z * c) % p)
            if len(support) >= s:
                output.add((z, cw, support))
    return output


def main():
    started = time.perf_counter()
    outputs = []
    bindings = []
    for path in sorted(HERE.glob('*.certificate.json')):
        record = json.loads(path.read_text())
        outputs.append({'file': path.name, **verify(record)})
        bindings.append({'path': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    # Small arbitrary words and boundary s=k, k=1; compare all q^k codewords.
    from overlap_cover import track_from_base, compile_certificate
    rng = random.Random(611816)
    tiny_count = tiny_nodes = whole_count = 0
    for p, n, k, s in [(3, 3, 1, 2), (3, 3, 2, 2), (3, 3, 2, 3),
                        (5, 4, 1, 3), (5, 4, 2, 3), (5, 4, 3, 4)]:
        for _ in range(12):
            domain = list(range(n))
            u0, u1 = [[rng.randrange(p) for _ in domain] for _ in range(2)]
            tracks = {}
            for base in combinations(range(n), k):
                tr = track_from_base(p, domain, k, u0, u1, base)
                key = tuple(tr['intercept_coefficients']), tuple(tr['slope_coefficients'])
                tracks[key] = tr
            record = compile_certificate(p, domain, k, s, u0, u1, list(tracks.values()))
            stats = verify(record)
            expected = literal_codebook(record)
            exported = {(r['scalar'], tuple(r['codeword']), tuple(r['agreement_support'])) for r in record['nodes']}
            assert exported == expected
            tiny_count += 1
            tiny_nodes += len(expected)
            whole_count += stats['whole_field_tracks']
    # Deterministic whole-field pencil control.
    p, dom, k, s, u0, u1 = 5, [0, 1, 2, 3], 2, 3, [1, 3, 0, 2], [2, 3, 4, 0]
    tr = track_from_base(p, dom, k, u0, u1, [0, 1])
    rec = compile_certificate(p, dom, k, s, u0, u1, [tr])
    assert verify(rec)['whole_field_tracks'] == 1 and rec['node_count'] == p
    assert len(literal_codebook(rec)) == p
    tiny_count += 1
    tiny_nodes += p
    whole_count += 1
    original = json.loads((HERE / 'coset_medium_characteristic_sample0.certificate.json').read_text())
    bads = []
    r = copy.deepcopy(original); r['tracks'][0]['joint_support'].append(15); bads.append(('false_support', r))
    r = copy.deepcopy(original); r['tracks'][0]['intercept_coefficients'][0] ^= 1; bads.append(('false_polynomial', r))
    r = copy.deepcopy(original); r['tracks'].pop(); bads.append(('missing_required_track', r))
    r = copy.deepcopy(original); r['nodes'].pop(); bads.append(('missing_node', r))
    r = copy.deepcopy(original); r['nodes'].append(r['nodes'][0]); bads.append(('duplicate_node', r))
    r = copy.deepcopy(original); r['domain'][0] += r['p']; bads.append(('noncanonical_domain', r))
    rejected = []
    for name, r in bads:
        try:
            verify(r)
        except (AssertionError, ValueError):
            rejected.append(name)
        else:
            raise AssertionError('Accepted corruption: ' + name)
    bindings.extend({'path': name, 'sha256': hashlib.sha256((HERE / name).read_bytes()).hexdigest()}
                    for name in ['overlap_cover.py', 'run_overlap.py', 'verify_overlap.py'])
    result = {'status': 'passed', 'independent_readbacks': outputs, 'tiny_full_codebook_cases': tiny_count,
              'tiny_nodes': tiny_nodes, 'tiny_whole_field_tracks': whole_count,
              'corrupted_certificates_rejected': rejected, 'source_bindings': bindings,
              'scope': 'All field scalars for every selected track, all s-sets, separate full-codebook tiny controls.',
              'elapsed_seconds': time.perf_counter() - started}
    (HERE / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
