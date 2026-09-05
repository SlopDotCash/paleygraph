#!/usr/bin/env python3
"""Hidden-witness fixtures, complete tiny oracles and honest failure controls."""
from copy import deepcopy
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import random
import sys
import time
sys.dont_write_bytecode = True
from piece_discovery import *
from scalar_fibers import materialize_scalar, interpolate
from compressed_support import domain


def book(F, dom, k):
    return [(cs, tuple(F.sum(F.mul(c, F.power(x, j)) for j, c in enumerate(cs)) for x in dom))
            for cs in product(range(F.q), repeat=k)]


def static_oracle(codebook, word, s):
    return {(cw, tuple(i for i, v in enumerate(cw) if v == word[i]))
            for _, cw in codebook if sum(v == word[i] for i, v in enumerate(cw)) >= s}


def check_tiny(F, dom, k, s, u0, word, codebook, result):
    if result['complete_static_list_certified']:
        got = {(tuple(x['codeword']), tuple(x['agreement_support']))
               for x in result['static_certificate']['complete_static_list']}
        assert got == static_oracle(codebook, word, s)
    if result['complete_scalar_fibers_certified']:
        actual = {(z, cw, support) for z in range(F.q)
                  for cw, support in static_oracle(codebook, [F.add(a, F.mul(z, b)) for a, b in zip(u0, word)], s)}
        got = {(z, tuple(node['codeword']), tuple(node['agreement_support'])) for z in range(F.q)
               for node in materialize_scalar(F, result['scalar_certificate'], z)}
        assert got == actual
        return len(actual)
    return 0


def summary(result):
    return {'status': result['status'], 'seconds': result['seconds'],
            'polynomial_cover_size_lower_bound': result['polynomial_cover_size_lower_bound'],
            'minimum_polynomial_cover_size_certified': result['minimum_polynomial_cover_size_certified'],
            'complete_static_list_certified': result['complete_static_list_certified'],
            'complete_scalar_fibers_certified': result['complete_scalar_fibers_certified'],
            'anchor_found': result['anchor'] is not None,
            'attempts': [{'r': a['piece_bound'], 'status': a['status'],
                          'rank': a.get('linear_certificate', {}).get('rank'),
                          'nullity': len(a.get('linear_certificate', {}).get('kernel_basis', [])),
                          'factor_status': a.get('factor_recovery', {}).get('status')}
                         for a in result['attempts']]}


def chunk_baseline(F, dom, word, k, s):
    """Every k arbitrary points interpolate, but this rarely certifies the list."""
    begin = time.perf_counter()
    grouped = {}
    for start in range(0, len(dom), k):
        indices = list(range(start, min(start + k, len(dom))))
        cs = interpolate(F, [dom[i] for i in indices], [word[i] for i in indices])
        cs += [0] * (k - len(cs))
        assert all(evaluate(F, cs, dom[i]) == word[i] for i in indices)
        grouped.setdefault(tuple(cs), []).extend(indices)
    cap = sum(min(len(indices), k - 1) for indices in grouped.values())
    return {'method': 'consecutive groups of at most k points, ordinary interpolation',
            'piece_count': len(grouped), 'root_count_cap': cap,
            'strict_list_certificate_available': s > cap, 'seconds': time.perf_counter() - begin}


def hidden_fixture(F, dom, k, r, seed):
    rng = random.Random(seed)
    polys = [[rng.randrange(F.q) for _ in range(k)] for _ in range(r)]
    labels = [i % r for i in range(len(dom))]
    rng.shuffle(labels)
    word = [evaluate(F, polys[labels[i]], x) for i, x in enumerate(dom)]
    f = [rng.randrange(F.q) for _ in range(k)]
    zstar = rng.randrange(F.q)
    u0 = [F.sub(evaluate(F, f, x), F.mul(zstar, y)) for x, y in zip(dom, word)]
    return u0, word, polys, labels


def main():
    start = time.perf_counter()
    output = {'schema': 'piece-discovery-experiments/v1', 'tiny': [], 'inherited': [],
              'hidden_random': [], 'generic': [], 'hard_cosets': [], 'large': [], 'negative_controls': []}
    tiny_counts = {'instances': 0, 'complete_static': 0, 'complete_scalar': 0,
                   'nodes_compared': 0, 'discovery_failures': 0, 'whole_field': 0}
    for F, n, k, s in [(PrimeField(3), 3, 1, 2), (PrimeField(5), 5, 2, 3),
                        (Field(3, [1, 0, 1]), 8, 2, 3)]:
        dom = list(range(n)); codebook = book(F, dom, k)
        inputs = []
        if F.q == 3:
            for word in product(range(F.q), repeat=n):
                for f in range(F.q):
                    zstar = (sum(word) + f) % F.q
                    u0 = [F.sub(f, F.mul(zstar, y)) for y in word]
                    inputs.append((u0, list(word)))
        else:
            inputs = [hidden_fixture(F, dom, k, 2, 51231 + seed)[:2] for seed in range(32)]
            rng = random.Random(5529 + F.q)
            for _ in range(16):
                word = [rng.randrange(F.q) for _ in dom]
                u0 = [F.sub(2, y) for y in word]
                inputs.append((u0, word))
        outcomes = {}
        for u0, word in inputs:
            result = discover(F, dom, k, s, word, u0)
            verify_export(F, result)
            nodes = check_tiny(F, dom, k, s, u0, word, codebook, result)
            tiny_counts['instances'] += 1
            tiny_counts['complete_static'] += result['complete_static_list_certified']
            tiny_counts['complete_scalar'] += result['complete_scalar_fibers_certified']
            tiny_counts['nodes_compared'] += nodes
            tiny_counts['discovery_failures'] += not result['complete_static_list_certified']
            if result['complete_scalar_fibers_certified']:
                tiny_counts['whole_field'] += result['scalar_certificate']['bad_scalar_count'] == F.q
            outcomes[result['status']] = outcomes.get(result['status'], 0) + 1
        output['tiny'].append({'field': field_record(F), 'n': n, 'k': k, 's': s,
                               'instances': len(inputs), 'all_successful_lists_and_all_scalar_nodes_compared': True,
                               'outcomes': outcomes})
    output['tiny_totals'] = tiny_counts
    print('tiny', tiny_counts, flush=True)

    old = json.loads((LAB / 'round4/support_batches/results.json').read_text())
    for item in old['cost_experiments']:
        d = item['result']; F = PrimeField(d['field']['p'])
        result = discover(F, d['domain'], d['k'], d['s'], d['u1'], d['u0'])
        verify_export(F, result)
        output['inherited'].append({'name': item['name'], 'summary': summary(result), 'result': result,
                                    'chunk_baseline': chunk_baseline(F, d['domain'], d['u1'], d['k'], d['s']),
                                    'planting_not_passed_to_discovery': True})
        print(item['name'], summary(result), flush=True)
    for seed in range(8):
        F = PrimeField(257); dom = list(range(192)); k, s = 8, 80
        u0, word, polys, labels = hidden_fixture(F, dom, k, 3, 80431 + seed)
        result = discover(F, dom, k, s, word, u0)
        assert result['complete_scalar_fibers_certified']
        assert {tuple(p['coefficients']) for p in result['discovered_pieces']} == {tuple(p) for p in polys}
        verify_export(F, result)
        output['hidden_random'].append({'seed': 80431 + seed, 'summary': summary(result),
                                       'exact_planted_polynomials_recovered_after_blind_run': True,
                                       'planting_not_passed_to_discovery': True})
    print('hidden random passed', flush=True)

    for seed in range(4):
        F = PrimeField(257); dom = list(range(64)); rng = random.Random(7811 + seed)
        word = [rng.randrange(F.q) for _ in dom]
        result = discover(F, dom, 4, 20, word)
        verify_export(F, result)
        output['generic'].append({'seed': 7811 + seed, 'summary': summary(result), 'result': result})
    hard = json.loads((LAB / 'proximity/stack_results.json').read_text())['records']
    for d in hard:
        F = PrimeField(d['p']); dom = domain(d['n'], d['p'])
        result = discover(F, dom, d['k'], d['s'], d['u1'], d['u0'])
        verify_export(F, result)
        output['hard_cosets'].append({'configuration': d['configuration'], 'sample': d['sample'],
                                     'summary': summary(result), 'result': result})
    print('generic and hard controls', len(output['hard_cosets']), flush=True)

    old_scalar = json.loads((LAB / 'round4/scalar_fibers/results.json').read_text())
    for item in old_scalar['large_interleaved']:
        d = item['complete_symbolic_certificate']; F = PrimeField(d['field']['p'])
        result = discover(F, d['domain'], d['k'], d['s'], d['u1'], d['u0'])
        assert result['complete_scalar_fibers_certified']
        assert result['static_certificate']['complete_static_list'] == d['static_direction_certificate']['complete_static_list']
        assert result['scalar_certificate']['complete_symbolic_node_count'] == d['complete_symbolic_node_count']
        verify_export(F, result)
        output['large'].append({'name': f"interleaved-n1024-k64-{item['parts']}-pieces", 'summary': summary(result),
                                'chunk_baseline': chunk_baseline(F, d['domain'], d['u1'], d['k'], d['s']),
                                'result': result, 'planting_not_passed_to_discovery': True,
                                'exact_equality_to_frozen_supplied_witness_list': True})
        print('large', item['parts'], summary(result), flush=True)

    # A known cover can fail the strict list threshold; success is not fabricated.
    F = PrimeField(17); dom = list(range(12)); k = 3
    u0, word, _, _ = hidden_fixture(F, dom, k, 2, 50591)
    cap_result = discover(F, dom, k, 4, word, u0, max_pieces=2)
    assert cap_result['status'] == 'piece_cover_discovered' and not cap_result['complete_static_list_certified']
    output['negative_controls'].append({'name': 'cover-found-list-cap-not-strict', 'result': cap_result})
    no_anchor_u0 = [F.mul(x, x) for x in dom]
    # k2 means X² is outside the code; this direction cannot cancel it globally.
    _, word, _, _ = hidden_fixture(F, dom, 2, 2, 50131)
    result = discover(F, dom, 2, 5, word, no_anchor_u0)
    assert result['complete_static_list_certified'] and result['anchor'] is None
    output['negative_controls'].append({'name': 'static-discovery-without-code-anchor', 'result': result})
    # Distinct linear branches whose difference vanishes at the only allowed center.
    word = [0 if i % 2 else x for i, x in enumerate(dom)]
    result = discover(F, dom, 2, 5, word, center_budget=1, max_pieces=2)
    assert not result['complete_static_list_certified']
    repaired = discover(F, dom, 2, 5, word, center_budget=8, max_pieces=2)
    assert repaired['complete_static_list_certified']
    output['negative_controls'].append({'name': 'single-colliding-center-fails-expanded-center-budget-recovers',
                                       'initial': result, 'iteration': repaired})
    for item in output['negative_controls']:
        for key in ['result', 'initial', 'iteration']:
            if key in item:
                verify_export(F, item[key])

    base = output['large'][0]['result']
    bad = deepcopy(base); bad['attempts'][-1]['linear_certificate']['particular'][0] ^= 1
    corruptions = [bad]
    bad = deepcopy(base); bad['attempts'][0]['linear_certificate']['inconsistency_witness']['rhs_residual'] = 0
    corruptions.append(bad)
    bad = deepcopy(base); bad['discovered_pieces'][0]['coordinates'].pop(); corruptions.append(bad)
    rejected = 0
    for bad in corruptions:
        try:
            verify_export(PrimeField(base['field']['p']), bad)
        except (AssertionError, ValueError):
            rejected += 1
    assert rejected == len(corruptions)
    output['tampered_exports_rejected'] = rejected
    output['source_sha256'] = {name: sha256((HERE / name).read_bytes()).hexdigest()
                               for name in ['piece_discovery.py', 'run_experiments.py']}
    output['elapsed_seconds'] = time.perf_counter() - start
    (HERE / 'results.json').write_text(json.dumps(output, indent=2) + '\n')
    print('complete', output['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
