#!/usr/bin/env python3
"""Exact finite principal-kernel shell certificate; stdlib only, no field scan.

This verifies one finite example. It proves no uniform period upper bound.
The analytic cosine, pi, and logarithm implications are in the companion note.
"""

from fractions import Fraction
from math import isqrt
from pathlib import Path
import json


def multiply(left, right):
    """Coefficient multiplication modulo X^d+1 over the integers."""
    d = len(left)
    assert len(right) == d
    product = [0] * d
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            quotient, remainder = divmod(i + j, d)
            product[remainder] += (-1 if quotient else 1) * x * y
    return product


def conjugate(value):
    """X -> X^(-1), with X^d=-1."""
    return [value[0]] + [-value[-j] for j in range(1, len(value))]


def norm_adjugate(value):
    """Quadratic norm descent, checking each lifted inverse identity."""
    d = len(value)
    if d == 1:
        return value[0], [1]
    even, odd = value[::2], value[1::2]
    down = multiply(even, even)
    odd_square = multiply(odd, odd)
    for j, coefficient in enumerate(odd_square):
        down[(j + 1) % (d // 2)] -= (
            coefficient * (-1 if j + 1 == d // 2 else 1)
        )
    norm, inverse_down = norm_adjugate(down)
    lift = [0] * d
    lift[::2] = inverse_down
    opposite = [(-x if j % 2 else x) for j, x in enumerate(value)]
    adjugate = multiply(opposite, lift)
    assert multiply(value, adjugate) == [norm] + [0] * (d - 1)
    return norm, adjugate


def rational(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def main():
    p, n, g = 215535361, 128, 25525303
    dimension = n // 2
    assert n & (n - 1) == 0
    assert p % n == 1
    assert all(p % divisor for divisor in range(2, isqrt(p) + 1))
    assert n**4 // 4 <= p <= n**4
    assert pow(g, n, p) == 1 and pow(g, n // 2, p) == p - 1

    f = [int(j in (0, 1, 19)) for j in range(dimension)]
    norm, adjugate = norm_adjugate(f)
    assert norm == p
    assert sum(c * pow(g, j, p) for j, c in enumerate(f)) % p == 0
    raw = multiply(conjugate(f), adjugate)
    assert multiply(f, raw) == [p * c for c in conjugate(f)]
    assert multiply(raw, conjugate(raw)) == [p * p] + [0] * (dimension - 1)
    assert sum(c * c for c in raw) == p * p

    centered = [(c + p // 2) % p - p // 2 for c in raw]
    assert centered == raw
    a = centered[0] % p
    assert a != 0
    inverse_generator = pow(g, -1, p)
    geometric = [
        (a * pow(inverse_generator, j, p) + p // 2) % p - p // 2
        for j in range(dimension)
    ]
    assert centered == geometric
    evaluate = lambda value, t: sum(
        c * pow(t, j, p) for j, c in enumerate(value)
    ) % p
    root_values = [evaluate(raw, pow(g, k, p)) for k in range(1, n, 2)]
    assert root_values[0] == a * dimension % p != 0
    assert root_values[1:] == [0] * (dimension - 1)

    V = Fraction(sum(c * c for c in centered), p * p)
    mean = Fraction(n * (p + 1), 24 * p)
    assert V == 1
    eta_lower = Fraction(n) - Fraction(1936, 49) * V
    assert eta_lower == Fraction(4336, 49) > 0
    log_upper = 15
    # e > 1+1+1/2+1/6 = 8/3, so this certifies log(p/n)<15.
    log_margin = n * 8**log_upper - p * 3**log_upper
    assert log_margin > 0
    amplitude_margin = eta_lower**2 - 4 * n * log_upper
    assert amplitude_margin == Fraction(361216, 2401) > 0

    certificate = {
        "scope": "One exact finite shell and literal C=2 amplitude obstruction; no uniform upper estimate.",
        "provenance": "Uses the previously identified norm-p trinomial from research/parallel30-triangle-remainder-2026-09-05.md; adjugate independently reconstructed here.",
        "p": p,
        "n": n,
        "N": dimension,
        "g": g,
        "inverse_generator": inverse_generator,
        "a": a,
        "a_coset_label": pow(a, n, p),
        "f": f,
        "norm_f": norm,
        "B": adjugate,
        "f_times_B": [p] + [0] * (dimension - 1),
        "A_raw": raw,
        "A_centered": centered,
        "A_times_conjugate_A": [p * p] + [0] * (dimension - 1),
        "sum_centered_squares": sum(c * c for c in centered),
        "geometric_form": "A_centered[j] = centered_residue(a*g^(-j) mod p)",
        "A_at_odd_roots_in_g_order": root_values,
        "V": rational(V),
        "vbar": rational(mean),
        "max_shell_deviation_lower": rational(mean - V),
        "eta_strict_lower_using_pi_lt_22_over_7": rational(eta_lower),
        "log_p_over_n_strict_upper": log_upper,
        "log_comparison_integer_margin": log_margin,
        "C_2_squared_comparison_margin": rational(amplitude_margin),
        "analytic_inputs": [
            "cos(x) >= 1-x^2/2 for every real x",
            "pi < 22/7",
            "e > 8/3",
        ],
        "finite_checks": [
            "primality by trial division through floor(sqrt(p))",
            "quartic window and exact dyadic generator order",
            "quadratic norm recursion and every lifted integer adjugate identity",
            "f*B=p and f*A=p*conjugate(f)",
            "A*conjugate(A)=p^2 as a full negacyclic coefficient identity",
            "all 64 coefficients already centered and have squared sum p^2",
            "all 64 coefficients equal the inverse-generator geometric form",
            "all 64 split-root values and exact scalar A(g)=a*N mod p",
            "exact rational inequalities implying eta>2*sqrt(n*log(p/n))",
        ],
    }
    output = Path(__file__).resolve().parents[1] / "results/parallel39_principal_shell_2026_09_06.json"
    output.write_text(json.dumps(certificate, indent=2) + "\n")
    print(json.dumps({
        "p": p, "n": n, "a": a, "V": str(V),
        "D_lower": str(mean - V), "eta_strict_lower": str(eta_lower),
        "literal_C_2_bound": "refuted at this finite quartic-window example",
        "certificate": str(output),
    }, indent=2))


if __name__ == "__main__":
    main()
