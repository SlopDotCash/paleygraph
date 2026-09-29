#!/usr/bin/env python3
"""Exact layered residual languages. No skipped levels; reject state implicit.

This is standard bottom-up decision-diagram minimization, not new automata
theory. Construction from a codebook still costs a complete scalar census.
"""
from functools import lru_cache


def reduce_layers(raw, accept):
    """raw[d][state] maps symbols to next states; discard empty residuals."""
    N = len(raw); layers = [None]*(N+1); layers[N] = [[]]
    counts = [None]*(N+1); counts[N] = [1]
    names = {s: 0 for s in accept}
    for d in range(N-1, -1, -1):
        intern = {}; renaming = {}; nodes = []; sizes = []
        for state, edges in raw[d].items():
            sig = tuple(sorted((x, names[t]) for x, t in edges.items() if t in names))
            if not sig: continue
            if sig not in intern:
                intern[sig] = len(nodes); nodes.append([list(e) for e in sig])
                sizes.append(sum(counts[d+1][t] for _, t in sig))
            renaming[state] = intern[sig]
        layers[d] = nodes; counts[d] = sizes; names = renaming
    assert len(layers[0]) == 1
    return {'dimension': N, 'root': 0, 'layers': layers, 'counts': counts,
            'widths': list(map(len, layers)), 'word_count': counts[0][0],
            'nodes': sum(map(len, layers)),
            'edges': sum(len(node) for layer in layers[:-1] for node in layer)}


def from_words(words):
    words = sorted(set(map(tuple, words))); N = len(words[0])
    assert all(len(w) == N for w in words)
    raw = [{} for _ in range(N)]
    for w in words:
        for d in range(N):
            raw[d].setdefault(w[:d], {})[w[d]] = w[:d+1]
    result = reduce_layers(raw, set(words))
    result['construction'] = 'complete_scalar_codebook'
    result['raw_states'] = sum(len(layer) for layer in raw)+len(words)
    return result


def generator_two(N, k):
    """All centered Q=2^N+1 encodings with S(D) divisible by k.

    State is (first nonzero sign, last nonzero sign, weighted residue).
    Nonzero signs alternate and their number is odd. Zero is included.
    """
    assert N >= 4 and N & (N-1) == 0 and k > 0 and ((1 << N)+1) % k == 0
    states = {(0, 0, 0)}; raw = []; transition_count = 0
    for d in range(N):
        layer = {}; following = set(); weight = pow(2, N-1-d, k)
        for state in sorted(states):
            first, last, r = state; edges = {}
            for digit in (-1, 0, 1):
                if digit and digit == last: continue
                nxt = (first or digit, digit or last, (r+weight*digit) % k)
                edges[digit] = nxt; following.add(nxt); transition_count += 1
            layer[state] = edges
        raw.append(layer); states = following
    accept = {s for s in states if s[2] == 0 and s[0] == s[1]}
    result = reduce_layers(raw, accept)
    result.update({'construction': 'alternating_sign_residue_grammar', 'cofactor': k,
                   'raw_states': sum(len(x) for x in raw)+len(states),
                   'raw_transitions': transition_count})
    assert result['word_count'] == ((1 << N)+1)//k
    return result


def complete(diagram, fixed, limit=0):
    """Count all completions of arbitrary coordinate assignments in read order.

    Return up to limit lexicographic complete words, with a truncation flag.
    This never labels a truncated list as a complete list.
    """
    N = diagram['dimension']; layers = diagram['layers']
    if type(limit) is not int or limit < 0: raise ValueError('invalid limit')
    if any(type(i) is not int or i < 0 or i >= N or type(v) is not int
           for i, v in fixed.items()): raise ValueError('invalid fixed coordinates')
    visits = 0
    @lru_cache(None)
    def count(d, s):
        nonlocal visits
        visits += 1
        if d == N: return 1
        return sum(count(d+1, t) for x, t in layers[d][s]
                   if d not in fixed or x == fixed[d])
    total = count(0, diagram['root']); words = []
    def walk(d, s, prefix):
        if len(words) >= limit: return
        if d == N: words.append(prefix); return
        for x, t in layers[d][s]:
            if (d not in fixed or x == fixed[d]) and count(d+1, t):
                walk(d+1, t, prefix+[x])
                if len(words) >= limit: return
    if limit: walk(0, diagram['root'], [])
    return {'count': total, 'words': words, 'truncated': len(words) < total,
            'count_states_visited': visits}


def accepted(diagram, word):
    if len(word) != diagram['dimension']: return False
    s = diagram['root']
    for d, x in enumerate(word):
        edges = dict(diagram['layers'][d][s])
        if x not in edges: return False
        s = edges[x]
    return True
