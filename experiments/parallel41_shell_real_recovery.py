#!/usr/bin/env python3
"""Exact real-subfield recovery certificate, stdlib only.

The real generator was discovered by one bounded 16-dimensional lattice
reduction. This verifier uses only its coefficients and exact arithmetic;
it makes no shortest-vector or minimum-cofactor claim.
"""

from fractions import Fraction
from pathlib import Path
import json


def multiply(a, b):
    n = len(a)
    assert len(b) == n
    out = [0] * n
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            k = i + j
            out[k % n] += (-1 if k >= n else 1) * x * y
    return out


def conjugate(a):
    return [a[0]] + [-a[-j] for j in range(1, len(a))]


def norm_adjugate(a):
    n = len(a)
    if n == 1:
        return a[0], [1]
    e, o = a[::2], a[1::2]
    down, oo = multiply(e, e), multiply(o, o)
    for j, value in enumerate(oo):
        down[(j + 1) % (n // 2)] -= value * (
            -1 if j + 1 == n // 2 else 1
        )
    norm, below = norm_adjugate(down)
    lift = [0] * n
    lift[::2] = below
    opposite = [(-v if j % 2 else v) for j, v in enumerate(a)]
    inverse = multiply(opposite, lift)
    assert multiply(a, inverse) == [norm] + [0] * (n - 1)
    return norm, inverse


def determinant(matrix):
    """Fraction-free elimination, distinct from quadratic norm descent."""
    a = [row[:] for row in matrix]
    n, denominator, sign = len(a), 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            row = next(j for j in range(k + 1, n) if a[j][k])
            a[k], a[row] = a[row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % denominator == 0
                a[i][j] = numerator // denominator
        for i in range(k + 1, n):
            a[i][k] = 0
        denominator = pivot
    return sign * a[-1][-1]


def real_embed(coefficients):
    """Basis 1, X^j+X^-j, 1<=j<D, in degree N=2D."""
    d = len(coefficients)
    out = [0] * (2 * d)
    out[0] = coefficients[0]
    for j in range(1, d):
        out[j] = coefficients[j]
        out[2 * d - j] = -coefficients[j]
    return out


def real_norm(coefficients):
    d = len(coefficients)
    h = real_embed(coefficients)
    matrix = [[0] * d for _ in range(d)]
    for j in range(d):
        e = [0] * d
        e[j] = 1
        product = multiply(h, real_embed(e))
        assert product == real_embed(product[:d])
        for i in range(d):
            matrix[i][j] = product[i]
    return determinant(matrix)


def rational(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def main():
    p, n, N, D, g = 6700417, 64, 32, 16, 2
    assert pow(g, N, p) == p - 1
    F = [(pow(g, j, p) + p // 2) % p - p // 2 for j in range(N)]
    norm_F, _ = norm_adjugate(F)
    assert norm_F == p**31 * 1217
    tau = 1217
    V = Fraction(sum(v * v for v in F), p * p)
    assert V == Fraction(6120237, p) < 1

    real_coefficients = [1, 0, 1, 0, -1, -1, 1, 0, -1, 0, -1, 0, -1, 1, 1, 0]
    h = real_embed(real_coefficients)
    assert h == conjugate(h)
    assert real_norm(real_coefficients) == p
    norm_h, adj_h = norm_adjugate(h)
    assert norm_h == p * p
    assert all(c % p == 0 for c in adj_h)
    c = [v // p for v in adj_h]
    assert multiply(h, c) == [p] + [0] * (N - 1)
    assert c == conjugate(c)

    evaluate = lambda a, t: sum(v * pow(t, j, p) for j, v in enumerate(a)) % p
    roots = [pow(g, j, p) for j in range(1, n, 2)]
    h_values = [evaluate(h, t) for t in roots]
    c_values = [evaluate(c, t) for t in roots]
    assert h_values[0] == h_values[-1] == 0
    assert all(v != 0 for v in h_values[1:-1])
    assert c_values[0] != 0 and c_values[-1] != 0
    assert c_values[1:-1] == [0] * (N - 2)

    numerator = multiply(h, conjugate(F))
    assert all(v % p == 0 for v in numerator)
    f = [v // p for v in numerator]
    assert multiply(c, conjugate(f)) == F
    assert norm_adjugate(f)[0] == p * tau
    assert evaluate(f, pow(g, -1, p)) == 0

    # A second exact check of the general lambda^2*tau rule.
    h_large_real = [5, -2] + [0] * (D - 2)
    h_large = real_embed(h_large_real)
    assert real_norm(h_large_real) == 641 * p
    numerator_large = multiply(h_large, conjugate(F))
    assert all(v % p == 0 for v in numerator_large)
    f_large = [v // p for v in numerator_large]
    assert norm_adjugate(f_large)[0] == p * 641**2 * tau

    # Reverse recovery of the already known norm-(641*p) element.
    original = conjugate([2, -1] + [0] * (N - 2))
    assert norm_adjugate(original)[0] == p * 641
    raw = multiply(c, conjugate(original))
    assert norm_adjugate(raw)[0] == p**31 * 641
    centered = [(v + p // 2) % p - p // 2 for v in raw]
    raw_V = Fraction(sum(v * v for v in raw), p * p)
    centered_V = Fraction(sum(v * v for v in centered), p * p)
    assert raw_V == Fraction(56196475, p)
    assert centered_V == Fraction(20588347, p) < raw_V
    centered_norm = norm_adjugate(centered)[0]
    assert centered_norm % p**31 == 0
    centered_tau = centered_norm // p**31
    assert centered_tau == 178771841 > 641
    new_a = centered[0] % p
    assert new_a == 2922708
    assert centered == [
        (new_a * pow(g, j, p) + p // 2) % p - p // 2 for j in range(N)
    ]

    output = {
        "scope": "Exact real-subfield arithmetic recovery and nonmonotonicity of norm under coefficient centering; no minimum-cofactor or uniform shell claim.",
        "p": p, "n": n, "N": N, "real_degree": D, "g": g,
        "F": F, "V": rational(V), "tau": tau,
        "h_real_basis_coefficients": real_coefficients,
        "h": h, "real_norm_h": p, "complex_norm_h": p * p,
        "c_p_over_h": c,
        "h_at_odd_roots": h_values, "c_at_odd_roots": c_values,
        "recovered_f": f, "norm_recovered_f": p * tau,
        "recovered_cofactor": tau,
        "larger_h_real_basis_coefficients": h_large_real,
        "larger_real_cofactor": 641,
        "larger_recovered_cofactor": 641**2 * tau,
        "known_original_f": original,
        "known_original_cofactor": 641,
        "reverse_raw_F": raw,
        "reverse_raw_norm_defect": 641,
        "reverse_raw_squared_radius": rational(raw_V),
        "reverse_centered_F": centered,
        "reverse_centered_a": new_a,
        "reverse_centered_norm_defect": centered_tau,
        "reverse_centered_squared_radius": rational(centered_V),
        "discovery": "One bounded 16-dimensional reduction of the real evaluation-kernel lattice supplied h; no claim that the reduced vector or recovered cofactor is minimal.",
        "verification": "Exact negacyclic multiplication and recursive norms; separate 16x16 fraction-free determinant for each real norm; exact modular evaluations and rational comparisons.",
    }
    target = Path(__file__).resolve().parents[1] / "results/parallel41_shell_real_recovery_2026_09_06.json"
    target.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({
        "real_norm_h": p,
        "recovered_cofactor": tau,
        "scalar_reciprocal_cofactor_replaced": "1217^31",
        "centering_norm_defect_change": [641, centered_tau],
        "centering_squared_radius_change": [str(raw_V), str(centered_V)],
        "certificate": str(target),
    }, indent=2))


if __name__ == "__main__":
    main()
