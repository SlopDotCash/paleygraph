#!/usr/bin/env python3
"""Exact anchor-weight and adjacency-cell checks; finite evidence, not the proof."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import numpy as np

from parallel_necklace_2026_09_04 import matrices
from parallel3_necklace_2026_09_04 import convolution_kernels

ROOT = Path(__file__).resolve().parents[1]


def psd_certificate(matrix):
    """Fraction-free symmetric elimination, retaining exact zero pivots."""
    a = [[int(x) for x in row] for row in matrix]
    n = len(a)
    assert all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    previous = 1
    pivots = []
    for j in range(n):
        pivot = a[j][j]
        assert pivot >= 0, (j, pivot)
        pivots.append(pivot)
        if pivot == 0:
            assert all(a[j][i] == 0 for i in range(j + 1, n))
            continue
        for i in range(j + 1, n):
            for k in range(i, n):
                numerator = pivot * a[i][k] - a[i][j] * a[j][k]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                a[i][k] = a[k][i] = quotient
        previous = pivot
    return pivots


def coeff_vectors(depth):
    return [[(1, 0)] * depth,
            [((-1) ** r, 0) for r in range(depth)],
            [(r + 1, (-1) ** r) for r in range(depth)]]


def one_field(p, R):
    chi, diag, cs = matrices(p)
    conv = convolution_kernels(p, chi, R)
    S = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    kernels = [None, S]
    for r in range(2, R+1):
        kernels.append(kernels[-1] @ cs['A'])
    ratios = [[y * pow(x, -1, p) % p for y in range(1, p)] for x in range(1, p)]
    kernel_checks = point_checks = 0
    for r in range(1, R+1):
        for x in range(1, p):
            for y in range(1, p):
                assert kernels[r][x, y] == chi[x] * conv[r][ratios[x-1][y-1]]
                kernel_checks += 1
        for t in range(1, p):
            rank_bound = r-1 if t == 1 else r
            assert conv[r][t] ** 2 <= rank_bound ** 2 * p ** (r-1)
            point_checks += 1
    supports = [(1,), (2,), (1, 2), (2, 3), (1, 2, 3), (2, 3, 4)]
    certificates = []
    sums = []
    aggregate_checks = []
    cell_checks = []
    for support in supports:
        m = len(support)
        w = [1] * p
        for t in range(p):
            for anchor in support:
                w[t] *= chi[(t-anchor) % p]
        for r in range(1, R+1):
            for s in range(r, R+1):
                d = m*r*s + (0 if 1 in support else r+s-1)
                bound_squared = d*d*p**(r+s-1)
                entry_kernel = [int(w[t] * conv[r][t] * conv[s][t]) for t in range(p)]
                for twist in ['trivial', 'quadratic']:
                    value = sum(entry_kernel[t] * (1 if twist == 'trivial' else chi[t])
                                for t in range(1, p))
                    assert value*value <= bound_squared
                    sums.append(dict(anchors=support, r=r, s=s, twist=twist, value=value, d=d))
                # The full real normal convolution operator checks every Mellin
                # character at once, including complex characters and equality.
                if p <= 17 or r+s in {2, R+1, 2*R}:
                    op = np.array([[entry_kernel[t] for t in row] for row in ratios], dtype=object)
                    full = np.array([[w[ratios[x-1][y-1]] * kernels[r][x, y] * kernels[s][x, y]
                                      for y in range(1, p)] for x in range(1, p)], dtype=object)
                    assert np.array_equal(op, full)
                    assert np.array_equal(op @ op.T, op.T @ op)
                    pivots = psd_certificate(bound_squared * np.eye(p-1, dtype=object) - op @ op.T)
                    certificates.append(dict(anchors=support, r=r, s=s, d=d,
                                             scaled_schur_pivots=pivots))
        for depth in sorted({1, 2, 3, R} & set(range(1, R+1))):
            U = depth*(depth+1)*(2*depth+1)//6
            c = (3*depth*depth-depth-2)//2
            for coeff in coeff_vectors(depth):
                norm = sum((x*x+y*y)*p**r for r, (x, y) in enumerate(coeff))
                abs_square = []
                for t in range(p):
                    real = sum(x*conv[r+1][t] for r, (x, y) in enumerate(coeff))
                    imag = sum(y*conv[r+1][t] for r, (x, y) in enumerate(coeff))
                    abs_square.append(real*real+imag*imag)
                C = (m+1-int(1 in support))*U
                for twist in ['trivial', 'quadratic']:
                    value = sum(w[t]*abs_square[t]*(1 if twist == 'trivial' else chi[t])
                                for t in range(1, p))
                    assert value*value <= C*C*p*norm*norm
                    aggregate_checks.append(dict(anchors=support, depth=depth, coefficients=coeff,
                                                 twist=twist, value=value, norm=norm))
                for signs in product([-1, 1], repeat=m):
                    cell = [t for t in range(1, p) if t not in support and
                            all(chi[(t-b) % p] == eps for b, eps in zip(support, signs))]
                    actual = sum(abs_square[t] for t in cell)
                    numerator = []
                    for t in range(p):
                        v = 1
                        for b, eps in zip(support, signs):
                            v *= 1+eps*chi[(t-b) % p]
                        numerator.append(v)
                    soft = sum(numerator[t]*abs_square[t] for t in range(1, p))
                    boundary = sum(numerator[t]*abs_square[t] for t in support)
                    assert (2**m)*actual == soft-boundary
                    assert all(numerator[t] in [0, 2**(m-1)] for t in support)
                    expansion = 0
                    for mask in range(2**m):
                        subtotal = 0
                        for t in range(1, p):
                            v = abs_square[t]
                            for j, (b, eps) in enumerate(zip(support, signs)):
                                if mask >> j & 1:
                                    v *= eps*chi[(t-b) % p]
                            subtotal += v
                        expansion += subtotal
                    assert expansion == soft
                    # Multiply (4) by p*2^m. The remaining sqrt(p) is cleared
                    # only after subtracting the nonnegative rational term.
                    alpha = c + ((2**(m-1))*(m+2-int(1 in support))-1)*U
                    beta = 2 + m*(2**(m-1))*U
                    deviation = abs((2**m)*actual-p*norm)
                    residual = max(deviation-beta*norm, 0)
                    assert residual*residual <= alpha*alpha*p*norm*norm
                    cell_checks.append(dict(anchors=support, signs=signs, depth=depth,
                                            coefficients=coeff, cardinality=len(cell), value=actual,
                                            boundary_numerator=boundary, norm=norm,
                                            cleared_error=[alpha, beta]))
    return dict(p=p, maximum_rank=R, kernel_entries=kernel_checks,
                pointwise_checks=point_checks, pair_sums=sums,
                mellin_operator_certificates=certificates,
                aggregate_checks=aggregate_checks, adjacency_cells=cell_checks)


def main():
    cases = []
    for p, R in [(5, 6), (13, 4), (17, 4), (29, 3), (41, 3)]:
        cases.append(one_field(p, R))
        print(f'passed p={p}, maximum rank={R}', flush=True)
    inputs = ['research/parallel5-anchor-aggregate-2026-09-04.md',
              'experiments/parallel5_anchor_aggregate_2026_09_04.py',
              'research/parallel4-kernel-aggregate-2026-09-04.md',
              'research/parallel3-necklace-2026-09-04.md',
              'experiments/parallel_necklace_2026_09_04.py',
              'experiments/parallel3_necklace_2026_09_04.py',
              'sources/katz-g2-hypergeometric.pdf',
              'sources/katz-gauss-kloosterman-monodromy.pdf']
    result = dict(status='Passed finite exact checks; uniform anchor-span theorem is in the proof note.',
                  arithmetic='Integer kernels, fraction-free positive-semidefinite certificates, integer-cleared quadratic and cell bounds.',
                  cases=cases,
                  input_sha256={name:sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs})
    dest = ROOT/'results/parallel5_anchor_aggregate_2026_09_04.json'
    dest.write_text(json.dumps(result, indent=2)+'\n')
    print(dest)


if __name__ == '__main__':
    main()
