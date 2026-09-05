#!/usr/bin/env python3
"""Whole-field oracle checks, no-anchor discovery and partition failures."""
from copy import deepcopy
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import random
import sys
import time
sys.dont_write_bytecode = True
from pencil_tracks import *
from scalar_fibers import interpolate
from compressed_support import domain


def codebook(F, dom, k):
    return [tuple(F.sum(F.mul(c, F.power(x, j)) for j, c in enumerate(cs)) for x in dom)
            for cs in product(range(F.q), repeat=k)]


def node_set(nodes):
    return {(node['scalar'], tuple(node['codeword']), tuple(node['agreement_support'])) for node in nodes}


def oracle(F, dom, s, u0, u1, book):
    out = set()
    for z in range(F.q):
        word = [F.add(a, F.mul(z, b)) for a, b in zip(u0, u1)]
        for cw in book:
            supp = tuple(i for i, y in enumerate(cw) if y == word[i])
            if len(supp) >= s:
                out.add((z, cw, supp))
    return out


def check_all_scalars(F, cert, book=None):
    nodes = set()
    # Direct equality and support calculations at every scalar; when book is
    # supplied this also compares every degree<k codeword independently.
    for z in range(F.q):
        got = materialize_scalar(F, cert, z)
        word = [F.add(a, F.mul(z, b)) for a, b in zip(cert['u0'], cert['u1'])]
        direct = set()
        for t in cert['tracks']:
            cs = [F.add(a, F.mul(z, b)) for a, b in zip(t['intercept_coefficients'], t['slope_coefficients'])]
            cw = tuple(F.sum(F.mul(c, F.power(x, j)) for j, c in enumerate(cs)) for x in cert['domain'])
            supp = tuple(i for i, y in enumerate(cw) if y == word[i])
            if len(supp) >= cert['s']:
                direct.add((z, cw, supp))
        assert node_set(got) == direct
        nodes.update(direct)
    assert len(nodes) == cert['complete_symbolic_node_count']
    assert len({z for z, _, _ in nodes}) == cert['bad_scalar_count']
    if book is not None:
        assert nodes == oracle(F, cert['domain'], cert['s'], cert['u0'], cert['u1'], book)
    return len(nodes)


def chunk_tracks(F, dom, k, u0, u1):
    tracks = []
    for start in range(0, len(dom), k):
        indices = list(range(start, min(start + k, len(dom))))
        a = interpolate(F, [dom[i] for i in indices], [u0[i] for i in indices])
        b = interpolate(F, [dom[i] for i in indices], [u1[i] for i in indices])
        tracks.append({'intercept_coefficients': a, 'slope_coefficients': b, 'coordinates': indices})
    return tracks


def linear_combination(F, a, b, c):
    return [F.add(x, F.mul(c, y)) for x, y in zip(a, b)]


def planted(F, dom, k, parts, seed, geometry='skew', grid=False):
    rng = random.Random(seed)
    f, g, h = [[rng.randrange(F.q) for _ in range(k)] for _ in range(3)]
    if geometry == 'collision':
        assert parts == 3
        aa = [linear_combination(F, f, g, c) for c in [0, 1, 3]]
        bb = [linear_combination(F, h, g, c) for c in [0, 1, 2]]
    elif parts == 2 and not grid:
        assert k >= 2
        g[1] = 1  # Difference g(X)+z cannot vanish identically at a scalar.
        aa = [f, linear_combination(F, f, g, 1)]
        one = [1] + [0] * (k - 1)
        bb = [h, linear_combination(F, h, one, 1)]
    else:
        aa = [[rng.randrange(F.q) for _ in range(k)] for _ in range(parts)]
        bb = [[rng.randrange(F.q) for _ in range(k)] for _ in range(parts)]
    labels = [(i % parts, (i // parts) % parts if grid else i % parts) for i in range(len(dom))]
    # Labels are interleaved in the domain, never passed to the discovery tool.
    u0 = [evaluate(F, aa[a], x) for (a, _), x in zip(labels, dom)]
    u1 = [evaluate(F, bb[b], x) for (_, b), x in zip(labels, dom)]
    witness = [{'intercept_coefficients': a, 'slope_coefficients': b,
                'coordinates': [i for i, pair in enumerate(labels) if pair == (j, ell)]}
               for j, a in enumerate(aa) for ell, b in enumerate(bb)
               if any(pair == (j, ell) for pair in labels)]
    assert find_code_anchor(F, dom, k, u0, u1) is None
    return u0, u1, witness


def summary(d):
    out = {'status': d['status'], 'seconds': d['seconds'], 'code_anchor': d['code_anchor'],
           'mode': d['mode'], 'raw_tracks': len(d.get('raw_intersection_tracks', [])),
           'retained_tracks': len(d.get('verified_tracks', [])), 'root_count_cap': d.get('joint_root_count_cap')}
    if d['complete']:
        c = d['certificate']
        out.update({'generic_list_size': c['generic_list_size'], 'exceptional_scalars': len(c['exceptional_scalars']),
                    'collision_events': len(c['track_collisions']), 'bad_scalars': c['bad_scalar_count'],
                    'symbolic_nodes': c['complete_symbolic_node_count']})
    return out


def main():
    start = time.perf_counter()
    output = {'schema': 'pencil-track-experiments/v1', 'tiny': [], 'blind_tiny': [],
              'medium': [], 'large': [], 'iterations': [], 'negative_controls': [], 'anchor_comparison': []}
    totals = {'complete_polynomial_scalar_oracles': 0, 'nodes': 0, 'no_anchor_cases': 0,
              'whole_field_cases': 0, 'collision_cases': 0}
    saved = None
    for F, n, k, s, count in [(PrimeField(3), 3, 1, 2, 729), (PrimeField(5), 5, 2, 4, 64),
                              (Field(3, [1, 0, 1]), 8, 2, 5, 40)]:
        dom = list(range(n)); book = codebook(F, dom, k)
        rng = random.Random(33811 + F.q)
        inputs = product(range(F.q), repeat=2 * n) if F.q == 3 else (
            [rng.randrange(F.q) for _ in range(2 * n)] for _ in range(count))
        local_nodes = 0
        for values in inputs:
            u0, u1 = list(values[:n]), list(values[n:])
            cert = compile_certificate(F, dom, k, s, u0, u1, chunk_tracks(F, dom, k, u0, u1))
            verify_certificate(F, cert)
            nodes = check_all_scalars(F, cert, book)
            local_nodes += nodes; totals['nodes'] += nodes
            totals['complete_polynomial_scalar_oracles'] += 1
            totals['no_anchor_cases'] += cert['code_anchor'] is None
            totals['whole_field_cases'] += cert['bad_scalar_count'] == F.q
            totals['collision_cases'] += bool(cert['track_collisions'])
            if F.q == 9 and saved is None:
                saved = cert
        output['tiny'].append({'field': field_record(F), 'n': n, 'k': k, 's': s,
                               'cases': count, 'nodes': local_nodes, 'all_polynomials_and_scalars_compared': True})
    output['tiny_totals'] = totals; output['saved_extension_certificate'] = saved
    print('tiny', totals, flush=True)

    for F, n, trials in [(PrimeField(5), 5, 20), (Field(3, [1, 0, 1]), 8, 20)]:
        dom = list(range(n)); k, s = 2, 3; book = codebook(F, dom, k)
        successes = failures = nodes = 0
        for seed in range(trials):
            u0, u1, _ = planted(F, dom, k, 2, 90191 + seed)
            d = discover_pencil(F, dom, k, s, u0, u1)
            verify_discovery(F, d)
            if d['complete']:
                successes += 1; nodes += check_all_scalars(F, d['certificate'], book)
            else:
                failures += 1
        output['blind_tiny'].append({'field': field_record(F), 'n': n, 'k': k, 's': s,
                                     'trials': trials, 'complete': successes, 'incomplete': failures,
                                     'nodes_compared': nodes, 'all_polynomials_and_scalars_compared_on_success': True})
    print('blind tiny', output['blind_tiny'], flush=True)

    F = PrimeField(257); dom = list(range(192)); k, s = 8, 80
    for name, parts, geometry in [('skew-two', 2, 'skew'), ('pair-colliding-three', 3, 'collision')]:
        u0, u1, _ = planted(F, dom, k, parts, 44801, geometry)
        d = discover_pencil(F, dom, k, s, u0, u1)
        assert d['complete'] and d['code_anchor'] is None
        verify_discovery(F, d)
        count = check_all_scalars(F, d['certificate'])
        output['medium'].append({'name': name, 'summary': summary(d), 'result': d,
                                 'all_scalars_direct_track_enumeration_nodes': count,
                                 'all_codewords_enumerated': False, 'hidden_labels_not_passed_to_discovery': True})
        print('medium', name, summary(d), flush=True)

    # Failure-driven iteration: ambiguous points create singleton intersection
    # fragments; full joint supports allow verified greedy reassignment.
    F = PrimeField(17); dom = list(range(12)); k, s = 2, 3
    labels = [int(i % 2 == 0) for i in range(12)]
    u0 = [x if label else 0 for x, label in zip(dom, labels)]
    u1 = labels
    raw = discover_pencil(F, dom, k, s, u0, u1, reassign=False)
    fixed = discover_pencil(F, dom, k, s, u0, u1, reassign=True)
    assert not raw['complete'] and fixed['complete']
    assert len(raw['verified_tracks']) == 3 and len(fixed['verified_tracks']) == 2
    verify_discovery(F, raw); verify_discovery(F, fixed)
    check_all_scalars(F, fixed['certificate'], codebook(F, dom, k))
    output['iterations'].append({'name': 'ambiguous-coordinate-fragment-reassignment', 'initial': raw, 'iteration': fixed})
    print('greedy iteration', summary(raw), summary(fixed), flush=True)

    # The same joint cover can be easier to find in one coordinate basis.
    F = PrimeField(257); dom = list(range(192)); k, s = 8, 80
    u0, u1, _ = planted(F, dom, k, 2, 47711, grid=True)
    good = discover_pencil(F, dom, k, s, u0, u1)
    bad = discover_pencil(F, dom, k, s, u0, u1, slice_scalars=(0, 1))
    assert good['complete'] and not bad['complete']
    verify_discovery(F, good); verify_discovery(F, bad)
    output['iterations'].append({'name': 'four-joint-regions-two-marginals-versus-four-piece-slice',
                                 'initial': bad, 'iteration': good})
    print('basis iteration', summary(bad), summary(good), flush=True)

    F = PrimeField(65537); dom = domain(1024, F.p); k, s = 64, 410
    for name, parts, geometry in [('skew-two', 2, 'skew'), ('pair-colliding-three', 3, 'collision')]:
        u0, u1, planting = planted(F, dom, k, parts, 118921, geometry)
        d = discover_pencil(F, dom, k, s, u0, u1)
        assert d['complete'] and d['code_anchor'] is None
        verify_discovery(F, d)
        c = d['certificate']
        # Compare selected scalar materializations to the supplied track witness,
        # which is revealed only after blind discovery. Completeness is certified.
        reference = compile_certificate(F, dom, k, s, u0, u1, planting)
        assert reference['complete_symbolic_node_count'] == c['complete_symbolic_node_count']
        assert reference['bad_scalar_count'] == c['bad_scalar_count']
        zs = sorted({0, 1, 2, 3, 8, 137, F.q - 1} | {x['scalar'] for x in c['track_collisions']})
        for z in zs:
            assert node_set(materialize_scalar(F, c, z)) == node_set(materialize_scalar(F, reference, z))
        output['large'].append({'name': name, 'summary': summary(d), 'result': d,
                               'spotchecked_scalars': zs, 'hidden_labels_not_passed_to_discovery': True,
                               'reference_symbolic_count_equality': True, 'field_scalars_enumerated': False})
        print('large', name, summary(d), flush=True)

    F = PrimeField(257); dom = list(range(192)); k, s = 8, 50
    u0, u1, _ = planted(F, dom, k, 3, 1291, grid=True)
    d = discover_pencil(F, dom, k, s, u0, u1)
    assert d['status'] == 'joint_cover_found_certificate_unavailable'
    verify_discovery(F, d)
    output['negative_controls'].append({'name': 'nine-joint-regions-cap-failure', 'result': d})
    for seed in range(3):
        F = PrimeField(257); dom = list(range(64)); rng = random.Random(7711 + seed)
        u0, u1 = [[rng.randrange(F.q) for _ in dom] for _ in range(2)]
        d = discover_pencil(F, dom, 4, 20, u0, u1)
        assert not d['complete']; verify_discovery(F, d)
        output['negative_controls'].append({'name': f'generic-no-small-marginal-cover-{seed}', 'result': d})

    old = json.loads((LAB / 'round4/scalar_fibers/results.json').read_text())
    inherited = old['inherited_fixtures'][1]
    # The saved certificate is used only as a reference, not to supply pieces.
    cert0 = inherited['certificate']
    F = PrimeField(cert0['field']['p'])
    d = discover_pencil(F, cert0['domain'], cert0['k'], cert0['s'], cert0['u0'], cert0['u1'])
    assert d['complete']; verify_discovery(F, d)
    assert d['certificate']['complete_symbolic_node_count'] == cert0['complete_symbolic_node_count']
    import scalar_fibers as old_route
    for z in range(F.q):
        assert node_set(materialize_scalar(F, d['certificate'], z)) == node_set(old_route.materialize_scalar(F, cert0, z))
    output['anchor_comparison'].append({'source_name': inherited.get('name'), 'result': d,
                                        'all_field_scalar_equality_to_frozen_anchor_route': True})

    base = output['large'][1]['result']['certificate']
    corruptions = []
    bad = deepcopy(base); bad['tracks'][0]['always_coordinates'].pop(); corruptions.append(bad)
    bad = deepcopy(base); bad['generic_scalar_domain']['excluded'].pop(); corruptions.append(bad)
    bad = deepcopy(base); bad['track_collisions'][0]['scalar'] ^= 1; corruptions.append(bad)
    bad = deepcopy(base); bad['exceptional_scalars'][0]['nodes'][0]['equal_track_ids'].pop(); corruptions.append(bad)
    bad = deepcopy(base); bad['complete_symbolic_node_count'] += 1; corruptions.append(bad)
    rejected = 0
    for bad in corruptions:
        try:
            verify_certificate(PrimeField(base['field']['p']), bad)
        except (AssertionError, ValueError):
            rejected += 1
    assert rejected == len(corruptions)
    output['corruptions_rejected'] = rejected
    output['elapsed_seconds'] = time.perf_counter() - start
    output['source_sha256'] = {name: sha256((HERE / name).read_bytes()).hexdigest()
                               for name in ['pencil_tracks.py', 'run_experiments.py']}
    (HERE / 'results.json').write_text(json.dumps(output, indent=2) + '\n')
    print('complete', output['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
