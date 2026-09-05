#!/usr/bin/env python3
"""Exact finite checks of seeded kernels; the uniform proof is in the note."""
from hashlib import sha256
from itertools import combinations_with_replacement, product
from pathlib import Path
import json
import numpy as np

from parallel_necklace_2026_09_04 import matrices
from parallel3_necklace_2026_09_04 import convolution_kernels
from parallel5_anchor_aggregate_2026_09_04 import psd_certificate

ROOT = Path(__file__).resolve().parents[1]


def one_field(p, R):
    chi, diag, cs = matrices(p)
    conv = convolution_kernels(p, chi, R)
    seeds = [1, 2] if p == 5 else [1, 2, 3]
    free = [t for t in range(1, p) if t not in seeds]
    supports = [(), (free[0],), (free[1],), tuple(free[:2])]
    index = [(b, r) for b in seeds for r in range(1, R+1)]
    kernels = {(b, r): [0] + [int(conv[r][t*pow(b, -1, p) % p])
                              for t in range(1, p)] for b, r in index}
    ratios = [[y*pow(x, -1, p) % p for y in range(1, p)] for x in range(1, p)]
    I = np.eye(p-1, dtype=object)
    J = np.zeros((p-1, p-1), dtype=object)
    for t in range(1, p):
        J[t-1, pow(t, -1, p)-1] = int(chi[t])
    assert np.array_equal(J@J, I)
    assert np.trace(J) == 2
    assert ((p-1)+int(np.trace(J)))//2 == (p+1)//2
    for r in range(1, R+1):
        assert np.array_equal(J@np.array(conv[r][1:], dtype=object), conv[r][1:])
    point_checks = 0
    for b, r in index:
        for t in range(1, p):
            rank = r-1 if t == b else r
            assert kernels[b, r][t]**2 <= rank**2*p**(r-1)
            point_checks += 1
    S = np.array([[int(chi[(x-y) % p]) for y in range(p)] for x in range(p)], dtype=object)
    v = np.array(chi[1:], dtype=object)
    full_seed = np.array([[int(chi[(1-t*pow(b, -1, p)) % p])
                           for b in range(1, p)] for t in range(1, p)], dtype=object)
    assert np.array_equal(full_seed, S[1:, 1:]*v[None, :])
    assert np.array_equal(full_seed@full_seed.T, p*I-np.ones((p-1, p-1), dtype=object)-np.outer(v, v))
    assert sum(v) == 0
    odd_example = None
    if p == 13:
        Q = [t for t in range(p) if chi[t] == chi[(t-1) % p] == 1]
        assert Q == [4, 10]
        sub = S[np.ix_(Q, Q)]
        assert sub.tolist() == [[0, -1], [-1, 0]]
        odd = np.array([1, -1], dtype=object)
        assert np.array_equal(sub@odd, odd)
        assert all(conv[r][4] == conv[r][10] for r in range(1, R+1))
        assert [kernels[3, 1][t] for t in Q] == [1, -1]
        odd_example = dict(Q=Q, restricted_S=sub.tolist(), new_seed=3,
                           old_sector_eigenvalue=-1, missing_sector_eigenvalue=1)
    sums, certs, aggregate, cells = [], [], [], []
    for support in supports:
        m = len(support)
        w = [1]*p
        for t in range(1, p):
            for a in support:
                w[t] *= int(chi[(t-a) % p])
        for (b, r), (c, s) in combinations_with_replacement(index, 2):
            f = [w[t]*kernels[b, r][t]*kernels[c, s][t] for t in range(p)]
            corrected = not support and b == c and r == s
            if corrected:
                f = [0] + [f[t]+p**(r-1)*((t == b)-1) for t in range(1, p)]
                d = 2*r-2
            else:
                d = m*r*s+r+s-int(b == c)
            bound = d*d*p**(r+s-1)
            for twist in ['trivial', 'quadratic']:
                total = sum(f[t]*(1 if twist == 'trivial' else int(chi[t])) for t in range(1, p))
                assert total*total <= bound, (p, support, b, r, c, s, twist, total, bound)
                sums.append(dict(anchors=support, left=[b, r], right=[c, s],
                                 corrected_diagonal=corrected, d=d, twist=twist, value=total))
            if p <= 17 or (r == s and r in {1, R}):
                op = np.array([[f[t] for t in row] for row in ratios], dtype=object)
                assert np.array_equal(op@op.T, op.T@op)
                pivots = psd_certificate(bound*I-op@op.T)
                certs.append(dict(anchors=support, left=[b, r], right=[c, s],
                                  corrected_diagonal=corrected, d=d, scaled_schur_pivots=pivots))
        for L in sorted({1, len(seeds)}):
            for depth in sorted({1, R}):
                ix = [(b, r) for b in seeds[:L] for r in range(1, depth+1)]
                U = depth*(depth+1)*(2*depth+1)//6
                C = L*depth*(3*depth+1)//2-depth-1
                coefficient_cases = [[(1, 0)]*len(ix),
                                     [((-1)**j, 0) for j in range(len(ix))],
                                     [(j+1, (-1)**j) for j in range(len(ix))]]
                for coeff in coefficient_cases:
                    norm = sum((x*x+y*y)*p**(r-1) for (b, r), (x, y) in zip(ix, coeff))
                    values = [0]*p
                    for t in range(1, p):
                        re = sum(x*kernels[b, r][t] for (b, r), (x, y) in zip(ix, coeff))
                        im = sum(y*kernels[b, r][t] for (b, r), (x, y) in zip(ix, coeff))
                        values[t] = re*re+im*im
                    for twist in ['trivial', 'quadratic']:
                        total = sum(w[t]*values[t]*(1 if twist == 'trivial' else int(chi[t]))
                                    for t in range(1, p))
                        if m:
                            assert total*total <= ((m+2)*L*U)**2*p*norm*norm
                        else:
                            deviation = abs(total-(p*norm if twist == 'trivial' else 0))
                            residual = max(deviation-2*norm, 0)
                            assert residual*residual <= C*C*p*norm*norm
                        aggregate.append(dict(anchors=support, L=L, R=depth, coefficients=coeff,
                                              twist=twist, value=total, coefficient_norm=norm))
                    if not m:
                        continue
                    for signs in product([-1, 1], repeat=m):
                        for zero_sign in [None, -1, 1]:
                            e = int(zero_sign is not None)
                            denominator = 2**(m+e)
                            cell = [t for t in range(1, p) if t not in support and
                                    all(chi[(t-a) % p] == eps for a, eps in zip(support, signs)) and
                                    (zero_sign is None or chi[t] == zero_sign)]
                            actual = sum(values[t] for t in cell)
                            numerator = [0]*p
                            for t in range(1, p):
                                value = 1 if zero_sign is None else 1+zero_sign*int(chi[t])
                                for a, eps in zip(support, signs):
                                    value *= 1+eps*int(chi[(t-a) % p])
                                numerator[t] = value
                            soft = sum(numerator[t]*values[t] for t in range(1, p))
                            border = sum(numerator[t]*values[t] for t in support)
                            assert denominator*actual == soft-border
                            assert all(0 <= numerator[t] <= denominator//2 for t in support)
                            expansion = 0
                            for mask in range(2**m):
                                for use_zero in range(1+e):
                                    subtotal = 0
                                    for t in range(1, p):
                                        value = values[t]
                                        if use_zero:
                                            value *= zero_sign*int(chi[t])
                                        for j, (a, eps) in enumerate(zip(support, signs)):
                                            if mask >> j & 1:
                                                value *= eps*int(chi[(t-a) % p])
                                        subtotal += value
                                    expansion += subtotal
                            assert expansion == soft
                            alpha = 2**e*(C+(m*2**(m-1)+2*(2**m-1))*L*U)
                            beta = 2**e*(2+m*2**(m-1)*L*U)
                            residual = max(abs(denominator*actual-p*norm)-beta*norm, 0)
                            assert residual*residual <= alpha*alpha*p*norm*norm
                            cells.append(dict(anchors=support, signs=signs, zero_sign=zero_sign,
                                              L=L, R=depth, coefficients=coeff, cardinality=len(cell),
                                              cell_mass=actual, boundary_numerator=border,
                                              coefficient_norm=norm, cleared_error=[alpha, beta]))
    b, c = seeds[:2]
    collision = sum(kernels[b, 1][t]*kernels[c, 1][t]*int(chi[(t-b) % p])*int(chi[(t-c) % p])
                    for t in range(1, p))
    assert collision == int(chi[b*c % p])*(p-3)
    return dict(p=p, maximum_rank=R, seeds=seeds, inversion_checks=R,
                pointwise_checks=point_checks, full_seed_gram_entries=(p-1)**2,
                missing_sector_example=odd_example, pair_sums=sums,
                mellin_operator_certificates=certs, quadratic_combinations=aggregate,
                adjacency_cells=cells, seed_anchor_collision=dict(seeds=[b, c], value=collision))


def main():
    cases = []
    for p, R in [(5, 6), (13, 3), (17, 3), (29, 3), (41, 3)]:
        cases.append(one_field(p, R))
        print(f'passed p={p}, maximum rank={R}', flush=True)
    inputs = ['research/parallel6-seeded-kernels-2026-09-04.md',
              'experiments/parallel6_seeded_kernels_2026_09_04.py',
              'research/parallel5-anchor-aggregate-2026-09-04.md',
              'experiments/parallel5_anchor_aggregate_2026_09_04.py',
              'research/parallel4-kernel-aggregate-2026-09-04.md',
              'research/parallel3-necklace-2026-09-04.md',
              'experiments/parallel_necklace_2026_09_04.py',
              'experiments/parallel3_necklace_2026_09_04.py',
              'sources/katz-g2-hypergeometric.pdf',
              'sources/katz-gauss-kloosterman-monodromy.pdf']
    result = dict(status='Finite exact checks passed; uniform seeded correlation proof is in the note.',
                  arithmetic='Integer kernels, exact fraction-free PSD certificates for every Mellin mode, integer-cleared cell inequalities.',
                  limitations='The old span misses an explicit largest eigenvector. The larger family has a limited isometry range; no full Paley spectral bound follows.',
                  cases=cases,
                  input_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs})
    dest = ROOT/'results/parallel6_seeded_kernels_2026_09_04.json'
    dest.write_text(json.dumps(result, indent=2)+'\n')
    print(dest)


if __name__ == '__main__':
    main()
