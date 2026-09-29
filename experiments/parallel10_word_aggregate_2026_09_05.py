#!/usr/bin/env python3
"""Exact finite bookkeeping and original-path Gram certificates, not a sheaf prover."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import isqrt, gcd, lcm
from pathlib import Path
import importlib.util
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('pass9', ROOT/'experiments/parallel9_all_degrees_2026_09_05.py')
P9 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P9)


def vertices(a):
    return list(range(a))+['y', 'infinity']


def freeze(local, a):
    return tuple(tuple((s, k, n) for (s, k), n in sorted(local[v].items()) if n)
                 for v in vertices(a))


def thaw(key, a):
    return {v: Counter({(s, k): n for s, k, n in row}) for v, row in zip(vertices(a), key)}


def rank(local):
    return P9.dim(local['y'])


def recover(local, a):
    reverse = []
    while rank(local) > 1:
        old_rank = rank(local)
        g, _ = P9.mc(local, a)
        differences = [P9.blocks(g[v], 0)-P9.blocks(g[v], 1) for v in range(a)]
        assert all(differences)
        mask = sum(1 << v for v, t in enumerate(differences) if t < 0)
        assert mask
        local, _ = P9.mask_twist(g, mask, a)
        assert rank(local) < old_rank
        assert all(P9.blocks(local[v], 0) > P9.blocks(local[v], 1) for v in range(a))
        reverse.append(mask)
    assert local == P9.initial(a)
    return tuple(reversed(reverse))


def kummer(v, a):
    local = {u: Counter({(int(u == v), 1): 1}) for u in vertices(a)}
    local['infinity'] = Counter({(1, 1): 1})
    return freeze(local, a)


def constant(a):
    return freeze({v: Counter({(0, 1): 1}) for v in vertices(a)}, a)


def inventory_step(inventory, mask, a):
    """Numerical local-type inventory. Equal types are not assumed isomorphic."""
    out = Counter()
    selected = {v for v in range(a) if mask >> v & 1}
    for key, multiplicity in inventory.items():
        if key[0] == 'point':
            v = key[1]
            if v not in selected:
                out[kummer(v, a)] += multiplicity
            continue
        local = thaw(key, a)
        g, _ = P9.mask_twist(local, mask, a)
        for v in selected:
            out[kummer(v, a)] += multiplicity*P9.blocks(g[v], 0)
        if rank(g) == 1:
            ramified = [v for v in vertices(a)[:-1] if P9.blocks(g[v], 1)]
            if not ramified:
                continue  # Constant middle extension has zero compact convolution.
            if len(ramified) == 1:
                out['point', ramified[0]] += multiplicity
                out[constant(a)] += multiplicity
                continue
        middle, _ = P9.mc(g, a)
        out[freeze(middle, a)] += multiplicity
        out[constant(a)] += multiplicity*P9.blocks(g['infinity'], 1)
    return +out


def inventory_invariants(inventory, a):
    total_rank = 0
    conductor = dict.fromkeys(vertices(a)[:-1], 0)
    for key, multiplicity in inventory.items():
        if key[0] == 'point':
            conductor[key[1]] += multiplicity
        else:
            local = thaw(key, a)
            r = rank(local)
            total_rank += multiplicity*r
            for v in conductor:
                conductor[v] += multiplicity*(r-P9.blocks(local[v], 0))
    return total_rank, conductor


def recurrent_levels(a, depth):
    w = r = q = 1
    z = 0
    out = []
    for m in range(depth+1):
        out.append(dict(depth=m, words=w, rank_sum=r, rank_square_sum=q, conductor_square_sum=z))
        if a == 1:
            w, r, q, z = 1, m+2, (m+2)**2, (m+1)**2
        else:
            h, n = 2**(a-1), 2**a-1
            a2 = 2**(a-2)*(a*a+3*a-1)-1
            w, r, q, z = n*w, (h*(a+1)-1)*r+h*w, a2*q+h*(a+2)*r+(h//2)*(z+w), a*h*q+(h-1)*z
    return out


def word_audit(a, depth, inventory_depth):
    levels = [dict(depth=m, words=0, rank_sum=0, rank_square_sum=0, conductor_square_sum=0,
                   principal_rank_square_sum=0, error_rank_square_sum=0) for m in range(depth+1)]
    signatures = {}
    recovery_checks = inventory_checks = 0
    def visit(word, local, conductors, inventory):
        nonlocal recovery_checks, inventory_checks
        m, d = len(word), rank(local)
        r = 1+sum(conductors)
        assert d <= r
        key = freeze(local, a)
        assert key not in signatures
        signatures[key] = word
        assert recover(local, a) == word
        recovery_checks += 1
        row = levels[m]
        for k, val in dict(words=1, rank_sum=r, rank_square_sum=r*r,
                           conductor_square_sum=sum(c*c for c in conductors),
                           principal_rank_square_sum=d*d, error_rank_square_sum=(r-d)**2).items():
            row[k] += val
        if inventory is not None:
            ri, ci = inventory_invariants(inventory, a)
            assert ri == r and ci['y'] == 1
            assert [ci[v] for v in range(a)] == conductors
            assert inventory[key] >= 1
            inventory_checks += 1
        if m == depth:
            return
        for mask in range(1, 1 << a):
            g, _ = P9.mask_twist(local, mask, a)
            nxt, _ = P9.mc(g, a)
            cs = [r if mask >> v & 1 else conductors[v] for v in range(a)]
            inv = inventory_step(inventory, mask, a) if m < inventory_depth else None
            visit(word+(mask,), nxt, cs, inv)
    start = P9.initial(a)
    visit((), start, [0]*a, Counter({freeze(start, a): 1}))
    expected = recurrent_levels(a, depth)
    for actual, prescribed in zip(levels, expected):
        for k, val in prescribed.items():
            assert actual[k] == val, (a, actual, prescribed)
    if a == 2:
        assert [r['rank_square_sum'] for r in levels[:4]] == [1, 17, 199, 2001]
    return dict(a=a, depth=depth, inventory_depth=inventory_depth, recovery_checks=recovery_checks,
                distinct_signatures=len(signatures), inventory_checks=inventory_checks, levels=levels)


def sqrt_upper(x, scale=10**9):
    x = F(x)
    assert x >= 0
    v = isqrt(x.numerator*scale*scale//x.denominator)
    if F(v*v, scale*scale) < x:
        v += 1
    out = F(v, scale)
    assert out*out >= x
    return out


def positive_definite_integer(matrix):
    """Bareiss leading-minor test with exact division and positive pivots."""
    b = [[int(x) for x in row] for row in matrix]
    previous = 1
    for k in range(len(b)):
        pivot = b[k][k]
        assert pivot > 0, ('nonpositive leading minor', k, pivot)
        for i in range(k+1, len(b)):
            for j in range(k+1, len(b)):
                numerator = pivot*b[i][j]-b[i][k]*b[k][j]
                value, rem = divmod(numerator, previous)
                assert rem == 0
                b[i][j] = value
        previous = pivot
    return previous.bit_length()


def gram_certificate(paths, rows, p, m, center, bound):
    count = len(paths)
    x = np.stack([v[rows] for v in paths], axis=1)
    max_entry = int(np.abs(x).max()) if x.size else 0
    if len(rows)*max_entry**2 >= 2**63:
        x = x.astype(object)
    numerator = x.T@x
    denom = p**(m+1)
    residual = [[F(int(numerator[i,j]), denom)-(center if i == j else 0)
                 for j in range(count)] for i in range(count)]
    # All denominators divide this positive integer.
    common = denom*center.denominator*bound.denominator
    bits = []
    for sign in (-1, 1):
        matrix = []
        for i in range(count):
            row = []
            for j in range(count):
                value = common*((bound if i == j else 0)+sign*residual[i][j])
                assert value.denominator == 1
                row.append(value.numerator)
            matrix.append(row)
        bits.append(positive_definite_integer(matrix))
    return dict(dimension=count, points=len(rows), center=str(center), upper_error=str(bound),
                relative_error_less_than_one=bound<center, determinant_bit_lengths=bits)


def field_audit(p, anchors, y, m, fresh, level):
    assert p % 4 == 1 and all(p % d for d in range(2, isqrt(p)+1))
    a = len(anchors)
    chi = np.array(P9.chi_table(p), dtype=np.int64)
    masks = {mask: np.array([int(np.prod([chi[(x-v) % p] for i,v in enumerate(anchors)
                                         if mask >> i & 1])) for x in range(p)], dtype=np.int64)
             for mask in range(1, 1 << a)}
    paths = {(): chi[(np.arange(p)-y) % p]}
    convolution_checks = literal_entries = mask_rows = 0
    for step in range(1, m+1):
        nxt = {}
        for word, path in paths.items():
            for mask, dm in masks.items():
                masked = dm*path
                assert p*int(np.abs(masked).max()) < 2**62
                convolution = np.convolve(chi, masked)
                result = convolution[:p].copy()
                result[:p-1] += convolution[p:]
                nxt[word+(mask,)] = result
                for i,v in enumerate(anchors):
                    if mask >> i & 1:
                        assert masked[v] == 0
                        mask_rows += 1
                for x in sorted({0, y, p//2, p-1}):
                    expected = sum(int(chi[(x-t) % p])*int(masked[t]) for t in range(p))
                    assert result[x] == expected
                    literal_entries += 1
                convolution_checks += 1
        paths = nxt
    assert len(paths) == (2**a-1)**m
    if p == 13 and m == 2:
        for word, path in paths.items():
            for x in range(p):
                value = 0
                for z1,z2 in product(range(p), repeat=2):
                    value += int(chi[(x-z1)%p])*int(masks[word[1]][z1])*int(chi[(z1-z2)%p])*int(masks[word[0]][z2])*int(chi[(z2-y)%p])
                assert value == path[x]
                literal_entries += 1
    dw, lw = level['principal_rank_square_sum'], level['error_rank_square_sum']
    sqrtinv = sqrt_upper(F(1,p))
    eps, gam = (a*dw+1)*sqrtinv, sqrt_upper(F(lw,p))
    bound = eps+2*gam*sqrt_upper(1+eps)+gam*gam
    rows = [x for x in range(p) if x not in set(anchors)|{y}]
    values = list(paths.values())
    certificates = [dict(kind='full_open', **gram_certificate(values, rows, p, m, F(1), bound))]
    for size in range(1,len(fresh)+1):
        beta = F(1,2**size)
        e = ((a+F(size,2))*dw+beta)*sqrtinv+F(size*dw,2*p)
        boundcell = e+2*gam*sqrt_upper(beta+e)+gam*gam
        for signs in product((-1,1), repeat=size):
            chosen = fresh[:size]
            cell = [x for x in rows if x not in chosen and all(chi[(x-v)%p] == s for v,s in zip(chosen,signs))]
            certificates.append(dict(kind='cell', fresh=chosen, signs=signs,
                                      **gram_certificate(values,cell,p,m,beta,boundcell)))
    return dict(p=p, anchors=anchors, y=y, depth=m, words=len(paths), convolution_checks=convolution_checks,
                independent_literal_entries=literal_entries, zero_mask_rows=mask_rows,
                certificates=certificates)


def null_vector_certificate():
    p, y, m = 13, 2, 3
    chi = P9.chi_table(p)
    s = np.array([[chi[(x-z) % p] for z in range(p)] for x in range(p)], dtype=np.int64)
    masks = {mask: np.array([int(np.prod([chi[(x-v)%p] for v in (0,1) if mask >> v & 1]))
                             for x in range(p)],dtype=np.int64) for mask in (1,2,3)}
    paths = {():s[:,y]}
    for step in range(m):
        paths = {word+(mask,):s@(dm*v) for word,v in paths.items() for mask,dm in masks.items()}
    words = list(paths)
    points = [x for x in range(p) if x not in (0,1,y)]
    original = [[int(paths[w][x]) for w in words] for x in points]
    mat = [[F(x) for x in row] for row in original]
    pivots = []
    for col in range(len(words)):
        row = next((i for i in range(len(pivots),len(points)) if mat[i][col]), None)
        if row is None:
            continue
        cur = len(pivots)
        mat[cur],mat[row] = mat[row],mat[cur]
        value = mat[cur][col]
        mat[cur] = [x/value for x in mat[cur]]
        for i in range(len(points)):
            if i != cur:
                value = mat[i][col]
                mat[i] = [x-value*y for x,y in zip(mat[i],mat[cur])]
        pivots.append(col)
        if len(pivots) == len(points):
            break
    free = next(i for i in range(len(words)) if i not in pivots)
    z = [F(0)]*len(words)
    z[free] = F(1)
    for row,pivot in enumerate(pivots):
        z[pivot] = -mat[row][free]
    den = lcm(*(v.denominator for v in z))
    ints = [int(v*den) for v in z]
    divisor = gcd(*ints)
    ints = [v//divisor for v in ints]
    residual = [sum(x*v for x,v in zip(row,ints)) for row in original]
    assert any(ints) and all(v == 0 for v in residual)
    assert len(words)>len(points)
    return dict(p=p,anchors=[0,1],y=y,depth=m,rows=len(points),columns=len(words),
                exact_rank=len(pivots),points=points,
                coefficients=[dict(word=list(w),coefficient=v) for w,v in zip(words,ints) if v],
                coefficient_squared_norm=sum(v*v for v in ints),integer_residual=residual,
                meaning='The actual raw-path evaluation map has a nonzero integer kernel vector; a relative isometry error below one is impossible for all coefficients.')


def main():
    audits=[]
    for a,depth,idepth in [(1,30,10),(2,6,4),(3,4,3),(4,3,2)]:
        audit=word_audit(a,depth,idepth)
        audits.append(audit)
        print(json.dumps({k:audit[k] for k in ('a','depth','recovery_checks','inventory_checks')}),flush=True)
    lookup={case['a']:case for case in audits}
    fields=[]
    for p,anchors,y,m,fresh in [(13,[0,1],2,2,[3,4]),(101,[0,1],2,3,[3,4]),
                                (257,[0,2,5],1,2,[3,4]),(10009,[0,1],2,1,[3,4]),
                                (10009,[0],1,4,[2,3])]:
        case=field_audit(p,anchors,y,m,fresh,lookup[len(anchors)]['levels'][m])
        fields.append(case)
        print(json.dumps(dict(p=p,depth=m,words=case['words'],certificates=len(case['certificates']),
                              nonvacuous=sum(c['relative_error_less_than_one'] for c in case['certificates']))),flush=True)
    inputs=['research/parallel10-word-aggregate-2026-09-05.md','experiments/parallel10_word_aggregate_2026_09_05.py',
            'research/parallel9-all-degrees-2026-09-05.md','experiments/parallel9_all_degrees_2026_09_05.py',
            'sources/katz-rigid-local-systems.pdf','sources/bbd-faisceaux-pervers.pdf',
            'sources/katz-gauss-kloosterman-monodromy.pdf']
    result=dict(status='All exact finite checks passed. Uniform geometric and weight claims require the source-dependent written proof.',
                scope='All-word quadratic aggregates on the open curve and fresh adjacency cells; not the full cyclic spectral aggregate.',
                arithmetic='Integer local block counts, original integer circular convolutions, rational upper square roots, and fraction-free positive-definiteness certificates.',
                limitations='Some bounds are vacuous at small primes. Independent literal path values test normalization; no numerical test proves sheaf irreducibility or purity.',
                word_audits=audits,field_audits=fields,null_vector_certificate=null_vector_certificate(),
                totals=dict(recovery_checks=sum(x['recovery_checks'] for x in audits),
                            distinct_signatures=sum(x['distinct_signatures'] for x in audits),
                            inventory_checks=sum(x['inventory_checks'] for x in audits),
                            convolutions=sum(x['convolution_checks'] for x in fields),
                            literal_entries=sum(x['independent_literal_entries'] for x in fields),
                            zero_mask_rows=sum(x['zero_mask_rows'] for x in fields),
                            gram_or_cell_bounds=sum(len(x['certificates']) for x in fields),
                            positive_definiteness_certificates=2*sum(len(x['certificates']) for x in fields),
                            bounds_with_relative_error_below_one=sum(c['relative_error_less_than_one'] for x in fields for c in x['certificates'])),
                input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    out=ROOT/'results/parallel10_word_aggregate_2026_09_05.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(written=str(out),totals=result['totals'])),flush=True)


if __name__=='__main__':
    main()
