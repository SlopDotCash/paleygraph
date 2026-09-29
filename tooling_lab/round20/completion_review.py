#!/usr/bin/env python3
"""Separate implementation: dense arithmetic and binary endpoint languages.

Does not import the producer, the diagram library, or the projection codec.
Induction on independently reconstructed binary suffixes verifies the large
symbolic languages. Exact graph reduction and reachability certify minimal
nonempty residual widths for these fixed layered readers.
"""
from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def validate_graph(g):
    N = g['dimension']; layers = g['layers']; counts = g['counts']
    assert len(layers) == len(counts) == N+1 and layers[-1] == [[]] and counts[-1] == [1]
    assert g['root'] == 0 and len(layers[0]) == 1
    reached = {0}
    for d in range(N):
        assert reached == set(range(len(layers[d])))
        following = set(); signatures = set()
        assert len(layers[d]) == len(counts[d])
        for s, node in enumerate(layers[d]):
            sig = tuple(map(tuple, node)); assert sig and sig == tuple(sorted(sig))
            assert len({x for x, _ in sig}) == len(sig) and sig not in signatures
            signatures.add(sig)
            assert all(type(x) is int and type(t) is int and 0 <= t < len(layers[d+1]) for x, t in sig)
            assert counts[d][s] == sum(counts[d+1][t] for _, t in sig)
            following.update(t for _, t in sig)
        reached = following
    assert reached == {0}
    assert g['widths'] == list(map(len, layers)) and g['word_count'] == counts[0][0]
    assert g['nodes'] == sum(map(len, layers))
    assert g['edges'] == sum(len(node) for layer in layers[:-1] for node in layer)


def direct_codebook(c, order):
    p, g, f = c['p'], c['g'], c['relation']; N = len(f); m = p//2
    assert p*p < np.iinfo(np.int64).max and m*sum(map(abs, f)) < np.iinfo(np.int64).max
    matrix = np.array([[f[i-j] if i >= j else -f[N+i-j] for j in range(N)] for i in range(N)], dtype=np.int64)
    powers = np.array([pow(g, j, p) for j in range(N)], dtype=np.int64)
    F = (np.arange(p, dtype=np.int64)[:, None]*powers+m) % p-m
    numerators = F @ matrix.T; assert np.all(numerators % p == 0)
    words = numerators//p
    return words[:, order]


def all_paths(g):
    def walk(d, s, prefix):
        if d == g['dimension']:
            yield tuple(prefix); return
        for x, t in g['layers'][d][s]: yield from walk(d+1, t, prefix+[x])
    return list(walk(0, g['root'], []))


def binary_counter(N, k, fixed=None):
    fixed = fixed or {}
    @lru_cache(None)
    def count(d, start, bit, residue):
        if d == N: return int(bit == 1-start and residue == 0)
        total = 0
        for nxt in (0, 1):
            digit = nxt-bit
            if d in fixed and digit != fixed[d]: continue
            r = (residue+pow(2, N-1-d, k)*digit) % k
            total += count(d+1, start, nxt, r)
        return total
    return count


def verify_binary_language(record):
    g = record['graph']; N = g['dimension']; k = record['config']['k']; count = binary_counter(N, k)
    initial = ((0, 0, 0), (1, 1, 0)); pending = {(0, initial, True)}; pairs = 0
    for d in range(N+1):
        following = set()
        for node, states, zero in pending:
            pairs += 1
            assert g['counts'][d][node] == sum(count(d, *s) for s in states)+int(zero)
            if d == N: continue
            expected = {}
            for digit in (-1, 0, 1):
                after = []
                for start, bit, residue in states:
                    nxt = bit+digit
                    if nxt not in (0, 1): continue
                    r = (residue+pow(2, N-1-d, k)*digit) % k
                    state = (start, nxt, r)
                    if count(d+1, *state): after.append(state)
                z = zero and digit == 0
                if after or z: expected[digit] = (tuple(sorted(after)), z)
            edges = dict(g['layers'][d][node]); assert set(edges) == set(expected)
            for digit, (states2, z) in expected.items(): following.add((edges[digit], states2, z))
        pending = following
    return pairs


def verify_queries(record, codebook=None):
    N = record['graph']['dimension']; checked_words = 0
    for q in record['queries']:
        fixed = dict(q['fixed']); words = q['words']
        assert len(words) == min(12, q['count']) and q['truncated'] == (len(words) < q['count'])
        assert len({tuple(w) for w in words}) == len(words) and words == sorted(words)
        if codebook is not None:
            mask = np.ones(len(codebook), dtype=bool)
            for i, v in fixed.items(): mask &= codebook[:, i] == v
            expected = sorted(map(tuple, codebook[mask].tolist()))
            assert q['count'] == len(expected) and list(map(tuple, words)) == expected[:12]
        else:
            k = record['config']['k']; count = binary_counter(N, k, fixed)
            expected_count = count(0, 0, 0, 0)+count(0, 1, 1, 0)+int(all(v == 0 for v in fixed.values()))
            assert q['count'] == expected_count
            # Independent greedy binary completion counts certify the listed
            # lexicographic first words, including queries with no completion.
            expected_words = []
            def visit(prefix, states, zero):
                if len(expected_words) == 12: return
                d = len(prefix)
                if d == N: expected_words.append(prefix); return
                for digit in (-1, 0, 1):
                    if d in fixed and fixed[d] != digit: continue
                    after = []
                    for start, bit, residue in states:
                        nxt = bit+digit
                        if nxt not in (0, 1): continue
                        state = (start, nxt, (residue+pow(2, N-1-d, k)*digit) % k)
                        if count(d+1, *state): after.append(state)
                    z = zero and digit == 0
                    size = sum(count(d+1, *s) for s in after)+int(z)
                    if size: visit(prefix+[digit], after, z)
                    if len(expected_words) == 12: return
            if expected_count: visit([], [(0, 0, 0), (1, 1, 0)], all(v == 0 for v in fixed.values()))
            assert expected_words == words
        checked_words += len(words)
    return checked_words


def main():
    summary_path = HERE/'completion_summary.json'; summary = json.loads(summary_path.read_text())
    reports = []; bound = {summary_path.name: sha256(summary_path.read_bytes()).hexdigest()}; small = None
    for c in summary['cases']:
        path = HERE/c['artifact']; r = json.loads(path.read_text()); bound[path.name] = sha256(path.read_bytes()).hexdigest()
        validate_graph(r['graph']); p = r['config']['p']; N = r['graph']['dimension']
        book = None; paths_checked = 0; grammar_pairs = 0
        if r['complete_census']:
            book = direct_codebook(r['config'], r['order'])
            paths = all_paths(r['graph']); actual = set(map(tuple, book.tolist()))
            assert len(actual) == p == len(paths) and set(paths) == actual
            paths_checked = len(paths)
        else:
            grammar_pairs = verify_binary_language(r)
        query_words = verify_queries(r, book)
        row = {'name': r['name'], 'direct_scalar_words': paths_checked, 'binary_language_pairs': grammar_pairs,
               'query_counts': len(r['queries']), 'listed_completions': query_words, 'minimal_widths': r['graph']['widths']}
        reports.append(row); print(json.dumps(row), flush=True)
        if small is None: small = r
    corruptions = []
    for name, mutate in [
            ('root_count', lambda g: g['counts'][0].__setitem__(0, g['counts'][0][0]+1)),
            ('edge_target', lambda g: g['layers'][0][0][0].__setitem__(1, 99999)),
            ('duplicate_unreachable_residual', lambda g: g['layers'][2].append(deepcopy(g['layers'][2][0]))),
            ('repeated_edge_label', lambda g: g['layers'][0][0].append(deepcopy(g['layers'][0][0][0]))),
            ('false_width', lambda g: g['widths'].__setitem__(2, g['widths'][2]-1))]:
        changed = deepcopy(small['graph']); mutate(changed)
        try: validate_graph(changed)
        except (AssertionError, IndexError): corruptions.append(name)
        else: raise AssertionError('corruption accepted: '+name)
    out = {'status': 'passed', 'scope': 'Full direct dense scalar-codebook equality for eight finite cases; independent binary-endpoint residual languages for five symbolic cases; all reported completion counts and first twelve words checked. Minimality is within fixed layered read order only.',
           'cases': reports, 'corrupt_artifacts_rejected': corruptions,
           'source_sha256': {'completion_review.py': sha256(Path(__file__).read_bytes()).hexdigest()}, 'input_sha256': bound,
           'numpy_version': np.__version__}
    (HERE/'completion_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
