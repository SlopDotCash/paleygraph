#!/usr/bin/env python3
"""Exact abstract cap minimum and maximum-agreement support certificates."""
from itertools import combinations
from math import comb
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]
sys.path.insert(0, str(LAB / 'round4/scalar_fibers'))
from scalar_fibers import PrimeField, Field, evaluate, interpolate, canonical, field_record


def minimum_cap(n, k, M):
    assert n >= 1 and 1 <= k <= n and 1 <= M <= n
    full, rest = divmod(n, M)
    t = k - 1
    parts = [M] * full + ([rest] if rest else []) if full <= 1000 else None
    return {'minimum': full * min(M, t) + min(rest, t),
            'attaining_integer_partition': parts, 'block_size_upper_bound': M,
            'compressed_attaining_partition': {'full_blocks': full, 'full_block_size': M, 'remainder': rest},
            'scope': 'integer partitions only; the minimizing sizes need not be realized by polynomial supports'}


def cap_barrier(n, k, s, M, proof_kind):
    assert 1 <= k <= s <= n
    if proof_kind not in ['exact_maximum', 'certified_upper_bound']:
        raise ValueError('A candidate maximum is only a lower bound and cannot certify this obstruction.')
    result = minimum_cap(n, k, M)
    return {'n': n, 'k': k, 's': s, 'M': M, 'M_evidence': proof_kind, **result,
            'status': 'entire_root_count_cap_family_blocked' if result['minimum'] >= s else 'abstract_size_bound_does_not_block',
            'meaning': 'Every verified polynomial or affine-track partition has cap at least this minimum. Only a minimum >=s rules out the family; the opposite inequality does not certify an achievable partition.'}


def enumerate_bases(F, dom, k, words):
    """Tiny generic-field oracle: exact scalar maxima and their joint maximum."""
    n = len(dom)
    assert 1 <= k <= n and len(set(dom)) == n and words
    canonical(F, dom)
    for w in words:
        assert len(w) == n; canonical(F, w)
    maxima = [0] * (len(words) + 1)
    witnesses = [None] * len(maxima)
    ledger = []
    for base in combinations(range(n), k):
        cs = [interpolate(F, [dom[i] for i in base], [w[i] for i in base]) for w in words]
        masks = [sum(1 << i for i, x in enumerate(dom) if evaluate(F, h, x) == w[i]) for h, w in zip(cs, words)]
        joint = (1 << n) - 1
        for mask in masks:
            joint &= mask
        masks.append(joint)
        ledger.append({'base': list(base), 'support_masks': masks})
        for j, mask in enumerate(masks):
            size = mask.bit_count()
            if size > maxima[j]:
                maxima[j] = size
                witnesses[j] = {'base': list(base), 'coefficients': cs if j == len(words) else [cs[j]],
                                'support': [i for i in range(n) if mask >> i & 1]}
    assert len(ledger) == comb(n, k) and min(maxima) >= k
    return {'field': field_record(F), 'domain': list(dom), 'k': k, 'words': words,
            'maxima': maxima, 'witnesses': witnesses, 'basis_count': len(ledger), 'ledger': ledger,
            'completeness': 'Every polynomial vector with >=k joint matches interpolates on one enumerated k-base; every word vector has an interpolant on each base.'}


def polynomial_partition_bound(F, dom, k, u0, u1, tracks):
    """Exact/upper maximum-joint bound from a supplied complete track partition.

    This does not require the scalar-list threshold inequality.
    """
    sys.path.insert(0, str(LAB / 'round6/pencil_tracks'))
    from pencil_tracks import normalize_tracks
    tracks = normalize_tracks(F, dom, k, u0, u1, tracks)
    cap = sum(min(len(t['coordinates']), k - 1) for t in tracks)
    supports = []
    for t in tracks:
        supports.append([i for i, x in enumerate(dom) if evaluate(F, t['intercept_coefficients'], x) == u0[i]
                         and evaluate(F, t['slope_coefficients'], x) == u1[i]])
    observed = max(map(len, supports))
    upper = max(observed, cap)
    return {'proof_kind': 'exact_maximum' if observed >= cap else 'certified_upper_bound',
            'M': upper, 'candidate_maximum': observed, 'outside_list_cap': cap,
            'maximal_joint_supports': supports, 'tracks': tracks,
            'proof': 'A polynomial pair distinct from every listed pair differs in at least one component on each region and has at most k-1 common roots there. The listed pairs are checked on the full domain.'}


def attain_from_maximum(F, dom, k, words, witness):
    """When M>=n/2 and its complement has <=k points, realize the cap minimum."""
    n = len(dom)
    cs = witness['coefficients']
    assert len(cs) == len(words) and all(len(h) <= k for h in cs)
    support = [i for i, x in enumerate(dom) if all(evaluate(F, h, x) == w[i] for h, w in zip(cs, words))]
    assert support == witness['support']
    M = len(support)
    complement = sorted(set(range(n)) - set(support))
    if 2 * M < n or len(complement) > k:
        return None
    pieces = [{'coefficients': cs, 'coordinates': support}]
    if complement:
        other = [interpolate(F, [dom[i] for i in complement], [w[i] for i in complement]) for w in words]
        pieces.append({'coefficients': other, 'coordinates': complement})
    for piece in pieces:
        for i in piece['coordinates']:
            assert all(evaluate(F, h, dom[i]) == w[i] for h, w in zip(piece['coefficients'], words))
    cap = sum(min(len(piece['coordinates']), k - 1) for piece in pieces)
    assert cap == minimum_cap(n, k, M)['minimum']
    return {'pieces': pieces, 'cap': cap, 'attains_abstract_minimum': True,
            'requires_upper_bound_evidence_for_optimality': True,
            'method': 'maximum support plus polynomial interpolation on its at-most-k-point complement'}
