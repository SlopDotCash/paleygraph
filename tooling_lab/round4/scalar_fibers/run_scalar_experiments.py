#!/usr/bin/env python3
"""Exact tiny oracles, inherited failure fixtures and interleaved scale cases."""
from copy import deepcopy
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import random
import sys
import time

sys.dont_write_bytecode = True
from scalar_fibers import *

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]
sys.path.insert(0, str(LAB / 'round4' / 'support_batches'))
from run_experiments import piecewise
from compressed_support import domain


def codebook(F, dom, k):
    # Independent direct monomial evaluation, not the production Horner helper.
    return [(list(cs), [F.sum(F.mul(c, F.power(x, j)) for j, c in enumerate(cs)) for x in dom])
            for cs in product(range(F.q), repeat=k)]


def oracle(F, dom, k, s, u0, u1, book):
    nodes = set()
    for z in range(F.q):
        word = [F.add(a, F.mul(z, b)) for a, b in zip(u0, u1)]
        for _, cw in book:
            support = tuple(i for i in range(len(dom)) if cw[i] == word[i])
            if len(support) >= s:
                nodes.add((z, tuple(cw), support))
    return nodes


def expanded(F, certificate):
    return {(node['scalar'], tuple(node['codeword']), tuple(node['agreement_support']))
            for z in range(F.q) for node in materialize_scalar(F, certificate, z)}


def from_direction_pieces(F, dom, k, cs_list, labels, fcs, zstar):
    direction = [evaluate(F, cs_list[labels[i]], x) for i, x in enumerate(dom)]
    anchor = [evaluate(F, fcs, x) for x in dom]
    u0 = [F.sub(f, F.mul(zstar, d)) for f, d in zip(anchor, direction)]
    pieces = [{'coefficients': cs, 'coordinates': [i for i, label in enumerate(labels) if label == j]}
              for j, cs in enumerate(cs_list)]
    return u0, direction, pieces


def inherited_pieces(F, planting):
    target, labels = planting['target_polynomial'], planting['labels']
    return [{'coefficients': [F.sub(f, a) for f, a in zip(target, cs)],
             'coordinates': [i for i, label in enumerate(labels) if label == j]}
            for j, cs in enumerate(planting['intercept_polynomials'])]


def main():
    start = time.perf_counter()
    output = {'schema': 'scalar-fiber-experiments/v1', 'tiny_oracles': [],
              'inherited_fixtures': [], 'large_interleaved': [], 'negative_controls': []}
    tiny_count = static_checks = anchor_checks = node_count = 0
    whole_field = 0
    # Every length-three direction over F3, every constant codeword anchor and
    # every designated anchor scalar. Constant pieces certify every static list.
    F = PrimeField(3); dom = list(range(3)); book = codebook(F, dom, 1)
    for word in product(range(3), repeat=3):
        pieces = [{'coefficients': [c], 'coordinates': [i for i, x in enumerate(word) if x == c]}
                  for c in range(3)]
        for f in range(3):
            for zstar in range(3):
                u0 = [F.sub(f, F.mul(zstar, x)) for x in word]
                cert = certify_pencil(F, dom, 1, 2, u0, list(word), pieces)
                actual = oracle(F, dom, 1, 2, u0, list(word), book)
                assert expanded(F, cert) == actual
                verify_certificate(F, cert)
                tiny_count += 1; node_count += len(actual)
                whole_field += cert['bad_scalar_count'] == F.q
    output['tiny_oracles'].append({'field': 3, 'n': 3, 'k': 1, 's': 2,
        'complete_pencil_oracles': tiny_count, 'all_direction_words_and_constant_anchors_and_anchor_scalars': True})
    saved_whole = None
    for F, n, trials in [(PrimeField(5), 5, 32), (Field(3, [1, 0, 1]), 8, 20)]:
        dom = list(range(n)); k, s = 2, 3
        book = codebook(F, dom, k)
        for seed in range(trials):
            rng = random.Random(19307 + seed + F.q)
            cs = [[rng.randrange(F.q) for _ in range(k)] for _ in range(2)]
            labels = [i % 2 for i in range(n)]
            fcs = [rng.randrange(F.q) for _ in range(k)]
            zstar = rng.randrange(F.q)
            u0, u1, pieces = from_direction_pieces(F, dom, k, cs, labels, fcs, zstar)
            cert = certify_pencil(F, dom, k, s, u0, u1, pieces)
            actual = oracle(F, dom, k, s, u0, u1, book)
            assert expanded(F, cert) == actual
            verify_certificate(F, cert)
            static_actual = {(tuple(cw), tuple(i for i in range(n) if cw[i] == u1[i]))
                             for _, cw in book if sum(a == b for a, b in zip(cw, u1)) >= s}
            static_cert = {(tuple(h['codeword']), tuple(h['agreement_support']))
                           for h in cert['static_direction_certificate']['complete_static_list']}
            assert static_actual == static_cert
            anchor_cs = interpolate(F, dom[:k], cert['anchor']['codeword'][:k])
            assert [evaluate(F, anchor_cs, x) for x in dom] == cert['anchor']['codeword']
            tiny_count += 1; static_checks += 1; anchor_checks += 1; node_count += len(actual)
            whole_field += cert['bad_scalar_count'] == F.q
            if F.q == 9 and cert['list_size_at_every_nonanchor_scalar'] == 2:
                saved_whole = cert
        output['tiny_oracles'].append({'field': F.q, 'n': n, 'k': k, 's': s,
            'complete_pencil_oracles': trials, 'seed_start': 19307 + F.q,
            'independent_static_list_equalities': trials})
    assert saved_whole is not None
    output['whole_field_extension_example'] = saved_whole
    output['tiny_totals'] = {'complete_pencil_oracles': tiny_count,
        'static_list_oracles': static_checks, 'anchor_interpolant_checks': anchor_checks,
        'exact_nodes_compared': node_count, 'whole_field_pencils': whole_field}
    inherited = json.loads((LAB / 'round4' / 'support_batches' / 'results.json').read_text())
    F = PrimeField(65537)
    for case in inherited['cost_experiments']:
        if case['planting']['layout'] == 'generic':
            continue
        r = case['result']; pieces = inherited_pieces(F, case['planting'])
        cert = certify_pencil(F, r['domain'], r['k'], r['s'], r['u0'], r['u1'], pieces)
        assert cert['bad_scalar_set'] == {'kind': 'singleton', 'scalar': 1}
        assert materialize_scalar(F, cert, 1) == r['nodes']
        verify_certificate(F, cert)
        output['inherited_fixtures'].append({'name': case['name'],
            'previous_track_interpolations': r['stats']['tracks_interpolated'],
            'exact_all_nodes_equal': True, 'certificate': cert})
        print(case['name'], 'previous tracks', r['stats']['tracks_interpolated'],
              'piece checks', cert['static_direction_certificate']['distinct_piece_polynomials'],
              'seconds', cert['elapsed_seconds'], flush=True)
    n, k, s = 1024, 64, 410
    dom = domain(n, F.q)
    for parts in [3, 2]:
        u0, u1, planting = piecewise(F, dom, k, s, 'interleaved', parts=parts)
        cert = certify_pencil(F, dom, k, s, u0, u1, inherited_pieces(F, planting))
        verify_certificate(F, cert)
        expected = 0 if parts == 3 else 2
        assert cert['list_size_at_every_nonanchor_scalar'] == expected
        # A few evaluated fibers check transport implementation; the algebraic
        # certificate, not these samples, establishes every scalar's behavior.
        for z in [0, 1, 2, 65536]:
            for node in materialize_scalar(F, cert, z):
                word = [F.add(a, F.mul(z, b)) for a, b in zip(u0, u1)]
                support = [i for i in range(n) if word[i] == node['codeword'][i]]
                assert support == node['agreement_support']
        output['large_interleaved'].append({'parts': parts, 'planting': planting,
            'complete_symbolic_certificate': cert, 'transport_spotcheck_scalars': [0, 1, 2, 65536]})
        print('interleaved1024', parts, 'static list', expected, 'symbolic nodes',
              cert['complete_symbolic_node_count'], 'seconds', cert['elapsed_seconds'], flush=True)
    # Strictness at the piece root bound is necessary: x agrees twice with
    # [0,0,3,3], yet is neither of its two constant pieces.
    F = PrimeField(7); dom = list(range(4)); word = [0, 0, 3, 3]
    pieces = [{'coefficients': [0], 'coordinates': [0, 1]},
              {'coefficients': [3], 'coordinates': [2, 3]}]
    try:
        certify_static_word(F, dom, 2, 2, word, pieces)
        raise AssertionError('accepted threshold equality')
    except UnavailableCertificate:
        pass
    support = [i for i, x in enumerate(dom) if x == word[i]]
    assert support == [0, 3]
    output['negative_controls'].append({'name': 'strict_root_bound_needed', 'rejected': True,
        'field': 7, 'n': 4, 'k': 2, 's': 2, 'word': word,
        'outside_piece_polynomial': [0, 1], 'agreement_support': support})
    assert find_code_anchor(F, dom, 2, [x*x % 7 for x in dom], [1]*4) is None
    output['negative_controls'].append({'name': 'no_exact_anchor', 'detected': True})
    F = PrimeField(5); dom = list(range(5)); book = codebook(F, dom, 2)
    anchor_nodes = [cw for _, cw in book if sum(x == 0 for x in cw) >= 1]
    assert len(anchor_nodes) == 21
    try:
        certify_pencil(F, dom, 2, 1, [0]*5, [1]*5,
                        [{'coefficients': [1], 'coordinates': list(range(5))}])
        raise AssertionError('accepted s<k')
    except UnavailableCertificate:
        pass
    output['negative_controls'].append({'name': 'anchor_s_below_k', 'rejected': True,
        'field': 5, 'n': 5, 'k': 2, 's': 1, 'actual_anchor_list_size': 21})
    original = saved_whole
    def reject(name, mutate):
        bad = deepcopy(original); mutate(bad)
        try:
            verify_certificate(Field(3, [1, 0, 1]), bad)
        except (AssertionError, UnavailableCertificate):
            return name
        raise AssertionError('accepted corrupted certificate ' + name)
    output['tampering_rejected'] = [
        reject('missing_fiber', lambda r: r['nonanchor_fibers'].pop()),
        reject('false_static_cap', lambda r: r['static_direction_certificate'].__setitem__('outside_piece_agreement_cap', -1)),
        reject('false_anchor_polynomial', lambda r: r['anchor']['coefficients'].__setitem__(0, (r['anchor']['coefficients'][0] + 1) % 9)),
        reject('false_symbolic_support', lambda r: r['nonanchor_fibers'][0]['agreement_support_for_every_scalar_in_domain'].pop()),
        reject('missing_piece_coordinate', lambda r: r['static_direction_certificate']['pieces'][0]['coordinates'].pop()),
        reject('noncanonical_coordinate', lambda r: r['domain'].__setitem__(0, 9)),
    ]
    output['elapsed_seconds'] = round(time.perf_counter() - start, 6)
    output['source_sha256'] = {name: sha256((HERE / name).read_bytes()).hexdigest()
                               for name in ['scalar_fibers.py', 'run_scalar_experiments.py']}
    (HERE / 'results.json').write_text(json.dumps(output, indent=2) + '\n')
    print('tiny', output['tiny_totals'], 'saved', output['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
