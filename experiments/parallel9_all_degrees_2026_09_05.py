#!/usr/bin/env python3
"""Exact rank bookkeeping and raw character sums; no numerical proof of weights."""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import random

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def dim(rep):
    return sum(length*n for (sign, length), n in rep.items())


def blocks(rep, sign):
    return sum(n for (s, length), n in rep.items() if s == sign)


def ones(rep, sign):
    return rep.get((sign, 1), 0)


def twist(rep):
    return Counter({(1-sign, length): n for (sign, length), n in rep.items()})


def quotient(rep):
    out = Counter()
    for (sign, length), n in rep.items():
        length -= sign == 0
        if length:
            out[sign, length] += n
    return out


def reconstruct(quot, size):
    out = Counter()
    for (sign, length), n in quot.items():
        out[sign, length+(sign == 0)] += n
    fill = size-dim(out)
    assert fill >= 0
    if fill:
        out[0, 1] += fill
    return out


def mc(local, a):
    finite = list(range(a))+['y']
    ambient = sum(dim(local[v])-blocks(local[v], 0) for v in finite)
    size = ambient-blocks(local['infinity'], 1)
    out = {v: reconstruct(twist(quotient(local[v])), size) for v in finite}
    out['infinity'] = quotient(twist(reconstruct(local['infinity'], ambient)))
    assert all(dim(rep) == size for rep in out.values())
    return out, size


def mask_twist(local, mask, a):
    out = dict(local)
    selected = [v for v in range(a) if mask >> v & 1]
    for v in selected:
        out[v] = twist(local[v])
    if len(selected) % 2:
        out['infinity'] = twist(local['infinity'])
    return out, selected


def initial(a):
    return {**{v: Counter({(0, 1): 1}) for v in range(a)},
            'y': Counter({(1, 1): 1}), 'infinity': Counter({(1, 1): 1})}


def packed(local):
    return {str(v): [[sign, length, n] for (sign, length), n in sorted(rep.items())]
            for v, rep in local.items()}


def rank_case(a, depth):
    totals = dict(transitions=0, inverse_checks=0, local_block_identities=0,
                  new_finite_invariant_dimensions=0)
    levels = {j: dict(count=0, min_rank=None, max_rank=0, min_jump=None)
              for j in range(1, depth+1)}
    witness = None
    def visit(local, word, previous_jump):
        nonlocal witness
        if len(word) == depth:
            return
        old = dim(local['y'])
        for mask in range(1, 1 << a):
            g, selected = mask_twist(local, mask, a)
            out, new = mc(g, a)
            jump = new-old
            step = len(word)+1
            assert new >= step+1 and new <= a*old+1
            if step == 1:
                assert new == len(selected)+(len(selected) % 2)
                assert jump >= 1
            else:
                fixed_t = {v: blocks(local[v], 0)-blocks(local[v], 1) for v in range(a)}
                infinity_t = blocks(local['infinity'], 1)-blocks(local['infinity'], 0)
                assert min(list(fixed_t.values())+[infinity_t]) >= previous_jump
                predicted = -previous_jump+sum(fixed_t[v] for v in selected)
                if len(selected) % 2:
                    predicted += infinity_t
                assert jump == predicted and jump >= previous_jump
            assert dim(quotient(out['y'])) == 1
            expected_y = Counter({(1, 1): 1, (0, 1): new-1}) if step % 2 == 0 else Counter({(0, 2): 1})
            if step % 2:
                if new > 2:
                    expected_y[0, 1] = new-2
            assert out['y'] == expected_y
            inverse, inverse_rank = mc(out, a)
            assert inverse_rank == old and inverse == g
            totals['inverse_checks'] += 1
            for v in range(a):
                assert blocks(out[v], 0) == jump+blocks(g[v], 0)
                assert blocks(out[v], 1) == blocks(g[v], 0)-ones(g[v], 0)
                assert blocks(out[v], 0)-blocks(out[v], 1) == jump+ones(g[v], 0)
                totals['local_block_identities'] += 3
            gi, oi = g['infinity'], out['infinity']
            assert blocks(oi, 1) == jump+blocks(gi, 1)
            assert blocks(oi, 0) == blocks(gi, 1)-ones(gi, 1)
            assert blocks(oi, 1)-blocks(oi, 0) == jump+ones(gi, 1)
            totals['local_block_identities'] += 3
            punctual = sum(blocks(g[v], 0) for v in selected)
            infinity = blocks(g['infinity'], 1)
            assert new+punctual+infinity <= (2*a+2)*old
            totals['new_finite_invariant_dimensions'] += punctual
            if a == 2 and word+(mask,) == (1, 2, 1):
                assert punctual == 1 and infinity == 0 and new == 6
                witness = dict(sequence=['A', 'B', 'A'], before=packed(local),
                               twisted=packed(g), after=packed(out),
                               new_invariant_anchor=0, new_invariant_dimension=1,
                               infinity_invariant_dimension=0, output_rank=new)
            level = levels[step]
            level['count'] += 1
            level['min_rank'] = new if level['min_rank'] is None else min(level['min_rank'], new)
            level['max_rank'] = max(level['max_rank'], new)
            level['min_jump'] = jump if level['min_jump'] is None else min(level['min_jump'], jump)
            totals['transitions'] += 1
            visit(out, word+(mask,), jump)
    visit(initial(a), (), 0)
    for j, level in levels.items():
        assert level['count'] == ((1 << a)-1)**j
    return dict(a=a, depth=depth, levels=levels, counts=totals, finite_boundary_witness=witness)


def chi_table(p):
    return [0]+[1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1, p)]


def field_case(p, anchors, seed):
    a = len(anchors)
    chi = chi_table(p)
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    masks = {}
    for mask in range(1, 1 << a):
        masks[mask] = np.array([int(np.prod([chi[(x-v) % p] for j, v in enumerate(anchors)
                                                     if mask >> j & 1])) for x in range(p)], dtype=object)
    rng = random.Random(seed)
    endpoint_checks = trace_sign_checks = literal_chains = zero_masks = small_orders = 0
    original_necklaces = 0
    nonzero_values = []
    # Trace orders one and two retain all diagonal zeros.
    for m, dm in masks.items():
        assert sum((dm[:, None]*s).diagonal()) == 0
        small_orders += 1
        for n, dn in masks.items():
            value = int(sum(((dm[:, None]*s)@(dn[:, None]*s)).diagonal()))
            assert value == int(sum(dm)*sum(dn)-dm@dn)
            assert abs(value) <= a*a*p
            small_orders += 1
    words = [tuple(rng.randrange(1, 1 << a) for _ in range(k))
             for k in range(3, 11) for _ in range(20)]
    if p == 5:
        words += list(product(range(1, 1 << a), repeat=3))
    for word in words:
        k = len(word)
        actual_matrix = np.eye(p, dtype=object)
        for mask in word:
            actual_matrix = actual_matrix@(masks[mask][:, None]*s)
        actual = int(sum(actual_matrix.diagonal()))
        # Build the full path from its right endpoint toward the left.
        path = s.copy()
        trace_function = -s.copy()
        for step, mask in enumerate(reversed(word[1:-1]), 1):
            masked = masks[mask][:, None]*trace_function
            for j, x in enumerate(anchors):
                if mask >> j & 1:
                    assert np.all(masked[x] == 0)
                    zero_masks += 1
            path = s@(masks[mask][:, None]*path)
            trace_function = -s@masked
            assert np.array_equal(trace_function, (-1)**(step+1)*path)
            trace_sign_checks += p*p
        assert np.array_equal(trace_function, (-1)**(k-1)*path)
        endpoint = int(sum((np.outer(masks[word[0]], masks[word[-1]])*s*path).flat))
        assert endpoint == actual
        endpoint_checks += 1
        coefficient = 3*a*(2*a+2)**(k-2)
        assert actual*actual <= coefficient**2*p**(k+1)
        if actual and len(nonzero_values) < 8:
            nonzero_values.append(dict(word=list(word), value=actual))
        if p == 5 and k == 3:
            literal = 0
            for xs in product(range(p), repeat=k):
                term = 1
                for i, mask in enumerate(word):
                    term *= int(masks[mask][xs[i]])*chi[(xs[i]-xs[(i+1) % k]) % p]
                literal += term
            assert literal == actual
            original_necklaces += 1
    if p == 5:
        y = next(x for x in range(p) if x not in anchors)
        for length in (1, 2):
            for word in product(range(1, 1 << a), repeat=length):
                path = s.copy()
                for mask in word:
                    path = s@(masks[mask][:, None]*path)
                for x in range(p):
                    value = 0
                    for zs in product(range(p), repeat=length):
                        vertices = (x,)+zs+(y,)
                        term = 1
                        for i in range(length+1):
                            term *= chi[(vertices[i]-vertices[i+1]) % p]
                        for i, mask in enumerate(reversed(word), 1):
                            term *= int(masks[mask][vertices[i]])
                        value += term
                    assert value == path[x, y]
                    literal_chains += 1
    return dict(p=p, anchors=list(anchors), endpoint_identities=endpoint_checks,
                trace_sign_entries=trace_sign_checks, literal_chain_entries=literal_chains,
                zero_mask_rows=zero_masks, orders_one_two=small_orders,
                original_literal_necklaces=original_necklaces, examples=nonzero_values)


def main():
    ranks = []
    for a, depth in ((1, 40), (2, 8), (3, 5), (4, 4), (5, 3)):
        case = rank_case(a, depth)
        ranks.append(case)
        print(json.dumps(dict(a=a, depth=depth, counts=case['counts'])), flush=True)
    fields = []
    for seed, (p, anchors) in enumerate(((5, (0, 1)), (5, (0, 1, 2)),
                                       (13, (0, 2, 5)), (17, (1, 3, 7, 10)),
                                       (29, (0, 1, 4)))):
        fields.append(field_case(p, anchors, seed))
        print(json.dumps(dict(p=p, anchors=anchors, status='exact raw-sum checks passed')), flush=True)
    inputs = [
        'research/parallel9-all-degrees-2026-09-05.md',
        'experiments/parallel9_all_degrees_2026_09_05.py',
        'sources/katz-rigid-local-systems.pdf',
        'sources/bbd-faisceaux-pervers.pdf',
        'sources/katz-gauss-kloosterman-monodromy.pdf',
        'sources/kunisky-2303.16475v1.html',
    ]
    out = ROOT/'results/parallel9_all_degrees_2026_09_05.json'
    result = dict(
        status='All finite exact checks passed; uniform theorem and weight assertions require the written source-dependent proof.',
        scope='Individual necklaces at arbitrary fixed degree and length; no full spectral-aggregate or Paley conclusion.',
        arithmetic='Integer Jordan-block multiplicities and original character sums; no floating-point acceptance.',
        limitations='Small-field upper-bound checks with exponential constants can be vacuous. They are not evidence for the asymptotic theorem. Rank, zero-mask and normalization checks are the substantive finite checks.',
        rank_cases=ranks, field_cases=fields,
        rank_totals={key: sum(c['counts'][key] for c in ranks) for key in ranks[0]['counts']},
        field_totals={key: sum(c[key] for c in fields) for key in (
            'endpoint_identities', 'trace_sign_entries', 'literal_chain_entries',
            'zero_mask_rows', 'orders_one_two', 'original_literal_necklaces')},
        input_sha256={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(written=str(out), rank_totals=result['rank_totals'],
                          field_totals=result['field_totals'])), flush=True)


if __name__ == '__main__':
    main()
