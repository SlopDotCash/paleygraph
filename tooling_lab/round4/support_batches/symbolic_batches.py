#!/usr/bin/env python3
"""Certified same-track batching of a complete interpolation-base cover.

ZDD family operations and polynomial interpolation are established methods.
This local instrument links exact symbolic subtraction to checkable track witnesses.
"""
from collections import Counter
from functools import lru_cache
from math import comb
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
ROUND2 = Path(__file__).resolve().parents[2] / 'round2' / 'proximity'
sys.path.insert(0, str(ROUND2))
from compressed_support import (PrimeField, Field, cover_plan, verify_cover,
                                interpolation_matrix, decode, node_set)


class FamilyZDD:
    """Canonical zero-suppressed set-family DAG, terminals 0=∅ and 1={∅}."""
    def __init__(self):
        self.nodes = [None, None]
        self.counts = [0, 1]
        self.supports = [0, 0]
        self.min_sizes = [float('inf'), 0]
        self.unique = {}
        self.visits = 0

    def make(self, variable, low, high):
        if high == 0:
            return low
        key = (variable, low, high)
        if key not in self.unique:
            self.unique[key] = len(self.nodes)
            self.nodes.append(key)
            self.counts.append(self.counts[low] + self.counts[high])
            self.supports.append((1 << variable) | self.supports[low] | self.supports[high])
            self.min_sizes.append(min(self.min_sizes[low], 1 + self.min_sizes[high]))
        return self.unique[key]

    def choose(self, coordinates, k):
        coordinates = tuple(sorted(coordinates))
        @lru_cache(None)
        def recur(position, remaining):
            if remaining == 0:
                return 1
            if remaining > len(coordinates) - position:
                return 0
            return self.make(coordinates[position], recur(position + 1, remaining),
                             recur(position + 1, remaining - 1))
        return recur(0, k)

    def first(self, root):
        assert root != 0
        answer = []
        while root > 1:
            variable, _, high = self.nodes[root]
            answer.append(variable)
            root = high
        assert root == 1
        return tuple(answer)

    def subtract_subsets(self, root, allowed_mask):
        """Remove all represented subsets of allowed without enumerating them."""
        @lru_cache(None)
        def recur(node):
            self.visits += 1
            if node == 0 or not (self.supports[node] & ~allowed_mask):
                return 0
            # Measured iteration: if too few allowed coordinates remain, this
            # entire branch is disjoint from the family being removed.
            if (self.supports[node] & allowed_mask).bit_count() < self.min_sizes[node]:
                return node
            variable, low, high = self.nodes[node]
            if allowed_mask & (1 << variable):
                return self.make(variable, recur(low), recur(high))
            return self.make(variable, recur(low), high)
        return recur(root)

    def materialize(self, root):
        if root == 0:
            return set()
        if root == 1:
            return {()}
        variable, low, high = self.nodes[root]
        return self.materialize(low) | {(variable,) + tail for tail in self.materialize(high)}


def track_coordinates(F, A, B, u0, u1):
    common, buckets = [], {}
    for i, (a, b, c, d) in enumerate(zip(A, B, u0, u1)):
        d0, d1 = F.sub(c, a), F.sub(d, b)
        if d1 == 0:
            if d0 == 0:
                common.append(i)
        else:
            z = F.mul(F.sub(0, d0), F.inv(d1))
            buckets.setdefault(z, []).append(i)
    return common, buckets


def run_census(F, dom, k, s, u0, u1, mode='partition', materialize_limit=1000):
    started = time.perf_counter()
    n = len(dom)
    assert len(set(dom)) == n and len(u0) == len(u1) == n
    plan = cover_plan(n, k, s, mode)
    verify_cover(plan)
    zdd = FamilyZDD()
    roots = [zdd.choose(block, k) for block in plan['blocks']]
    block_sets = list(map(set, plan['blocks']))
    initial_nodes = len(zdd.nodes)
    assert sum(zdd.counts[root] for root in roots) == plan['basis_count']
    nodes, tracks, batches, all_field, finite = {}, set(), [], [], []
    stats = Counter()
    while any(roots):
        root = next(root for root in roots if root)
        base = zdd.first(root)
        assert len(base) == k
        matrix = interpolation_matrix(F, [dom[i] for i in base], dom)
        A = decode(F, matrix, [u0[i] for i in base])
        B = decode(F, matrix, [u1[i] for i in base])
        assert (A, B) not in tracks
        tracks.add((A, B))
        common, buckets = track_coordinates(F, A, B, u0, u1)
        common_set = set(common)
        mask = sum(1 << i for i in common)
        expected = sum(comb(len(block & common_set), k)
                       for block in block_sets if len(block & common_set) >= k)
        before = sum(zdd.counts[root] for root in roots)
        roots = [zdd.subtract_subsets(root, mask) for root in roots]
        removed = before - sum(zdd.counts[root] for root in roots)
        # Distinct degree<k tracks cannot share k common coordinates. Thus no
        # selected base inside this common set was removed by an earlier track.
        assert removed == expected >= 1
        stats['tracks_interpolated'] += 1
        stats['bases_certified'] += removed
        record = {'first_base': list(base), 'base_multiplicity': removed,
                  'intercept': list(A), 'slope': list(B), 'common_coordinates': common}
        batches.append(record)
        if len(common) >= s:
            record['all_field'] = True
            all_field.append(record)
            chosen = range(F.q) if F.q <= materialize_limit else []
        else:
            chosen = [z for z, indices in buckets.items() if len(indices) + len(common) >= s]
            if chosen:
                record['all_field'] = False
                record['qualifying_scalar_buckets'] = {str(z): buckets[z] for z in sorted(chosen)}
                finite.append(record)
        for z in chosen:
            cw = tuple(F.add(a, F.mul(z, b)) for a, b in zip(A, B))
            support = sorted(common + buckets.get(z, []))
            assert len(support) >= s
            nodes[(z, cw)] = {'scalar': z, 'codeword': list(cw), 'agreement_support': support}
    assert stats['bases_certified'] == plan['basis_count']
    stats['interpolations_avoided'] = plan['basis_count'] - stats['tracks_interpolated']
    stats['zdd_nodes_initial'] = initial_nodes
    stats['zdd_nodes_allocated'] = len(zdd.nodes)
    stats['zdd_subtraction_node_visits'] = zdd.visits
    node_list = [node for _, node in sorted(nodes.items())]
    whole = bool(all_field)
    return {'field': {'p': F.p, 'e': F.e, 'q': F.q}, 'n': n, 'k': k, 's': s,
            'domain': dom, 'u0': u0, 'u1': u1, 'cover_certificate': plan,
            'cover_verified': True, 'completeness': 'complete exact selected-cover census; all skipped bases carry same-track certificates',
            'batch_certificates': batches, 'finite_track_certificates': finite,
            'all_field_track_certificates': all_field,
            'finite_bad_scalar_count': F.q if whole else len({x['scalar'] for x in node_list}),
            'whole_field_correlated': whole,
            'bad_scalars': list(range(F.q)) if whole and F.q <= materialize_limit else
                (None if whole else sorted({x['scalar'] for x in node_list})),
            'nodes': node_list, 'node_ledger_materialized': not whole or F.q <= materialize_limit,
            'stats': dict(stats), 'elapsed_seconds': round(time.perf_counter() - started, 6)}


def verify_batches(F, record):
    """Certificate readback; does not use the ZDD or enumerate covered k-bases."""
    plan = record['cover_certificate']
    verify_cover(plan)
    k, n = record['k'], record['n']
    dom, u0, u1 = record['domain'], record['u0'], record['u1']
    assert len(dom) == n and len(set(dom)) == n
    tracks, total = set(), 0
    for batch in record['batch_certificates']:
        base = batch['first_base']
        assert len(base) == len(set(base)) == k
        assert any(set(base) <= set(block) for block in plan['blocks'])
        matrix = interpolation_matrix(F, [dom[i] for i in base], dom)
        A, B = tuple(batch['intercept']), tuple(batch['slope'])
        assert A == decode(F, matrix, [u0[i] for i in base])
        assert B == decode(F, matrix, [u1[i] for i in base])
        assert (A, B) not in tracks
        tracks.add((A, B))
        common = {i for i in range(n) if A[i] == u0[i] and B[i] == u1[i]}
        assert sorted(common) == batch['common_coordinates']
        size = sum(comb(len(common & set(block)), k) for block in plan['blocks']
                   if len(common & set(block)) >= k)
        assert size == batch['base_multiplicity'] >= 1
        total += size
    assert total == plan['basis_count']
    return {'distinct_tracks_checked': len(tracks), 'covered_bases_certified': total,
            'all_skips_certified_without_basis_enumeration': True}
