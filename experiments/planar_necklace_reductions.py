#!/usr/bin/env python3
"""Exact affine/projective graph checks and binary length-six necklaces.

All arithmetic is integer. Graph partition functions are summed directly,
without substituting the formulas under test. Finite checks supplement the
ordinary proofs in research/planar-necklace-reductions.md.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

import numpy as np

from paley_exact import character
from localized_necklace_identities import convolution_moments, paired_trace

ROOT = Path(__file__).resolve().parents[1]


def edge(u, v):
    return tuple(sorted((u, v)))


def frame_and_dual(word):
    k = len(word)
    assert k >= 3 and set(word) == {'0', '1'}
    edges = [edge(i, (i + 1) % k) for i in range(k)]
    edges += [edge(i, k + int(word[i])) for i in range(k)]
    faces = []
    for label in ('0', '1'):
        positions = [i for i, x in enumerate(word) if x == label]
        for start, end in zip(positions, positions[1:] + positions[:1]):
            arc = [start]
            while len(arc) == 1 or arc[-1] != end:
                arc.append((arc[-1] + 1) % k)
            if label == '1':
                arc.reverse()
            faces.append([k + int(label)] + arc)
    incidence = {e: [] for e in edges}
    for f, cycle in enumerate(faces):
        for u, v in zip(cycle, cycle[1:] + cycle[:1]):
            incidence[edge(u, v)].append((f, u, v))
    assert len(faces) == k and k + 2 - len(edges) + len(faces) == 2
    dual = []
    for pair in incidence.values():
        assert len(pair) == 2
        (f, u, v), (g, x, y) = pair
        assert (u, v) == (y, x)
        dual.append(edge(f, g))
    return edges, dual, faces


def delete(vertices, edges, removed):
    removed = set(removed)
    return [v for v in vertices if v not in removed], [e for e in edges if not removed.intersection(e)]


def graph_sum(p, chi, vertices, edges, projective=False, fixed=None):
    """Literal recursive sum, pruning only zero edge weights.

    The affine version with no prescribed values fixes one vertex to 0
    and multiplies by p, using translation invariance alone.
    Infinity is represented by p. Repeated edges are retained.
    """
    vertices = list(vertices)
    fixed = dict(fixed or {})
    assert all(v in vertices for v in fixed)
    if any(u == v for u, v in edges):
        return 0
    if not vertices:
        return 1
    if not projective and not fixed:
        v = max(vertices, key=lambda x: sum(x in e for e in edges))
        return p * graph_sum(p, chi, vertices, edges, fixed={v: 0})

    def weight(x, y):
        if x == p or y == p:
            return 0 if x == y else 1
        return chi[(x - y) % p]

    assigned = dict(fixed)
    factor = 1
    for u, v in edges:
        if u in assigned and v in assigned:
            factor *= weight(assigned[u], assigned[v])
    if not factor:
        return 0
    remaining = [v for v in vertices if v not in assigned]
    # Select an elimination ordering using adjacency alone, before enumeration.
    order = []
    known = set(assigned)
    while remaining:
        v = max(remaining, key=lambda x: (sum((x == a and b in known) or (x == b and a in known)
                                               for a, b in edges), sum(x in e for e in edges), -x))
        order.append(v)
        known.add(v)
        remaining.remove(v)
    backward = []
    known = set(assigned)
    for v in order:
        backward.append([b if a == v else a for a, b in edges
                         if (a == v and b in known) or (b == v and a in known)])
        known.add(v)
    domain = range(p + int(projective))

    def visit(index):
        if index == len(order):
            return 1
        v = order[index]
        total = 0
        for value in domain:
            sign = 1
            for u in backward[index]:
                sign *= weight(value, assigned[u])
                if sign == 0:
                    break
            if sign:
                assigned[v] = value
                total += sign * visit(index + 1)
        assigned.pop(v, None)
        return total

    return factor * visit(0)


def canonical(word):
    words = [word, word[::-1], ''.join('1' if x == '0' else '0' for x in word)]
    words.append(words[-1][::-1])
    return min(w[i:] + w[:i] for w in words for i in range(len(w)))


def formula(word, p, t):
    k = len(word)
    minority = min(word.count('0'), word.count('1'))
    if minority == 0:
        return (p - 1) * t[k]
    if minority == 1:
        return -t[k]
    if minority == 2:
        mark = '0' if word.count('0') == 2 else '1'
        places = [i for i, x in enumerate(word) if x == mark]
        ell = places[1] - places[0]
        r = k - ell
        s, d = min(ell, r), abs(ell - r)
        return p*t[ell]*t[r] - p**s*t[d] - t[k] + 2*(-1)**k*sum(p**j for j in range(1, s))
    assert k == 6 and minority == 3
    name = canonical(word)
    return {'000111': p*t[3]**2-t[6],
            '001011': p*(t[5]+2*t[3])-t[6],
            '010101': p*(p-5)*t[4]+p*p*(p-2)-t[6]}[name]


def matrix_case(p):
    chi = character(p)
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    matrices = [np.array([chi[(x-b) % p] for x in range(p)], dtype=object)[:, None]*s for b in (0, 1)]
    powers = [[np.eye(p, dtype=object)] for _ in range(2)]
    t = convolution_moments(p, chi, 12)
    for b in (0, 1):
        for _ in range(12):
            powers[b].append(powers[b][-1] @ matrices[b])
    blocks = 0
    for a in range(1, 12):
        for b in range(1, 13-a):
            observed = paired_trace(powers[0][a], powers[1][b])
            assert observed == p*t[a]*t[b]-t[a+b]
            assert observed**2 <= (a+b)**4*p**(a+b)
            blocks += 1
    values = {}
    def descend(word, power):
        if len(word) == 6:
            value = sum(int(x) for x in power.diagonal())
            assert value == formula(word, p, t), (p, word)
            assert value**2 <= 36**2*p**7
            values[word] = value
            return
        for b in (0, 1):
            descend(word+str(b), power @ matrices[b])
    descend('', np.eye(p, dtype=object))
    # Affine averaging of powers, evaluated directly over all p anchors.
    averaging = 0
    if p <= 17:
        sums = [np.zeros((p, p), dtype=object) for _ in range(6)]
        for b in range(p):
            cb = np.array([chi[(x-b) % p] for x in range(p)], dtype=object)[:, None]*s
            power = np.eye(p, dtype=object)
            for j in range(6):
                power = power @ cb
                sums[j] += power
        base = p*np.eye(p, dtype=object)-np.ones((p, p), dtype=object)
        for j in range(1, 7):
            assert np.array_equal(sums[j-1], t[j]*base)
            averaging += 1
    return {'p': p, 'two_block_checks': blocks, 'binary_six_checks': len(values),
            'affine_power_averaging_checks': averaging, 'binary_six_values': values,
            'traces_t0_through_t12': t}


def graph_case(p):
    chi = character(p)
    t = convolution_moments(p, chi, 6)
    # Original planar frames and their face-by-face constructed duals.
    dual_results = []
    for word in ('001011', '010101'):
        original, dual, faces = frame_and_dual(word)
        z = graph_sum(p, chi, range(8), original)
        zd = graph_sum(p, chi, range(6), dual)
        assert z == p*zd
        assert z == p*(p-1)*formula(word, p, t)+p*(p-1)*t[6]
        dual_results.append({'word': word, 'faces': faces, 'dual_edges': dual,
                             'original_partition': z, 'dual_partition': zd})
    h = [edge(a,b) for a,b in [(0,1),(0,2),(1,2),(3,4),(3,5),(4,5),
                               (0,5),(1,5),(1,3),(2,3),(2,4),(2,5)]]
    assert Counter(h) == Counter(frame_and_dual('001011')[1])
    completed = h + [edge(6,v) for v in (0,2,4,5)]
    assert all(sum(v in e for e in completed) % 2 == 0 for v in range(7))
    assert all(edge(2,v) in completed for v in range(7) if v != 2)
    fixed_new = graph_sum(p, chi, range(7), completed, projective=True, fixed={6:p})
    fixed_universal = graph_sum(p, chi, range(7), completed, projective=True, fixed={2:p})
    assert fixed_new == fixed_universal == p*(p-1)*t[5]
    zh = graph_sum(p, chi, range(6), h)
    deletions = [graph_sum(p, chi, *delete(range(6), h, [v])) for v in (1,3)]
    assert deletions == [-p*(p-1)*t[3]]*2
    assert zh+sum(deletions) == fixed_new
    # Octahedron: three disjoint nonedge pairs, all other pairs are edges.
    octa = [edge(a,b) for a,b in combinations(range(6),2) if a//2 != b//2]
    z_affine = graph_sum(p, chi, range(6), octa)
    z_projective = graph_sum(p, chi, range(6), octa, projective=True)
    elliptic = [sum(chi[u]*chi[(u-1)%p]*chi[(u-a)%p] for u in range(p)) for a in range(p)]
    weighted = sum(chi[a]*elliptic[a]**2 for a in range(1,p))
    assert weighted == t[4]
    assert z_projective == p*(p*p-1)*(weighted+p)
    assert z_projective-z_affine == 6*p*(p-1)*t[4]+3*p*p*(p-1)
    assert z_affine == dual_results[1]['dual_partition']
    return {'p':p, 'duality_checks':dual_results, 'completed_condition_at_new_vertex':fixed_new,
            'completed_condition_at_universal_vertex':fixed_universal,
            'deleted_partitions':deletions, 'octahedron_affine':z_affine,
            'octahedron_projective':z_projective, 'elliptic_twisted_second_moment':weighted}


def main():
    balanced = {canonical(''.join(w)) for w in product('01',repeat=6) if w.count('0')==3}
    assert balanced == {'000111','001011','010101'}
    # Verify oriented embeddings for every nonconstant binary length-six word,
    # including bridges (and hence loops in the dual) when an anchor occurs once.
    embeddings = 0
    for letters in product('01', repeat=6):
        word = ''.join(letters)
        if len(set(word)) == 2:
            original, dual, faces = frame_and_dual(word)
            assert len(original) == len(dual) == 12 and len(faces) == 6
            embeddings += 1
    assert embeddings == 62
    matrices=[]
    for p in (5,13,17,29,41,61,97):
        matrices.append(matrix_case(p))
        print(json.dumps({'p':p,'matrix_checks':'passed'}),flush=True)
    graphs=[]
    for p in (5,13):
        graphs.append(graph_case(p))
        print(json.dumps({'p':p,'literal_graph_checks':'passed'}),flush=True)
    sources=['experiments/planar_necklace_reductions.py','experiments/localized_necklace_identities.py',
             'experiments/paley_exact.py','sources/lu-zheng-zheng-1305.3405v3.html']
    result={'status':'Exact checks passed; no full Paley/localization theorem',
            'arithmetic':'Python integers and NumPy object matrices',
            'balanced_binary_bracelets':sorted(balanced),'matrix_cases':matrices,'graph_cases':graphs,
            'oriented_embedding_checks':embeddings,
            'totals':{key:sum(c[key] for c in matrices) for key in
                      ('two_block_checks','binary_six_checks','affine_power_averaging_checks')},
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in sources}}
    (ROOT/'results/planar_necklace_reductions.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['totals']),flush=True)


if __name__=='__main__':
    main()
