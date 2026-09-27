#!/usr/bin/env python3
"""Stdlib-only exact checks for the finite obstruction to a shell converse.

No scan over field cosets, numerical cosines, project build, or proof service.
The nonprincipality argument is ordinary mathematics in the companion note.
"""

from fractions import Fraction
from math import isqrt
from pathlib import Path
import json


def multiply(left, right):
    d = len(left)
    assert len(right) == d
    out = [0] * d
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            k = i + j
            out[k % d] += (-1 if k >= d else 1) * a * b
    return out


def conjugate(value):
    return [value[0]] + [-value[-j] for j in range(1, len(value))]


def norm_adjugate(value):
    d = len(value)
    if d == 1:
        return value[0], [1]
    even, odd = value[::2], value[1::2]
    down, odd_square = multiply(even, even), multiply(odd, odd)
    for j, coefficient in enumerate(odd_square):
        down[(j + 1) % (d // 2)] -= coefficient * (
            -1 if j + 1 == d // 2 else 1
        )
    norm, below = norm_adjugate(down)
    lift = [0] * d
    lift[::2] = below
    opposite = [(-x if j % 2 else x) for j, x in enumerate(value)]
    adjugate = multiply(opposite, lift)
    assert multiply(value, adjugate) == [norm] + [0] * (d - 1)
    return norm, adjugate


def rational(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def main():
    p, auxiliary_prime, n, g, a = 6700417, 641, 64, 2, 1
    N = n // 2
    for prime in (p, auxiliary_prime):
        assert all(prime % d for d in range(2, isqrt(prime) + 1))
        assert prime % n == 1
        assert pow(g, N, prime) == prime - 1
        assert pow(g, n, prime) == 1
    assert n**4 // 4 <= p <= n**4
    assert 2**N + 1 == p * auxiliary_prime

    F = [(a * pow(g, j, p) + p // 2) % p - p // 2 for j in range(N)]
    S = sum(c * c for c in F)
    V = Fraction(S, p * p)
    assert S == 41008140038829
    assert V == Fraction(6120237, 6700417) < 1
    evaluate = lambda value, t, prime: sum(
        c * pow(t, j, prime) for j, c in enumerate(value)
    ) % prime
    root_values = [evaluate(F, pow(g, j, p), p) for j in range(1, n, 2)]
    assert root_values == [0] * (N - 1) + [N]
    support_root = pow(g, -1, p)

    norm_F, adjugate_F = norm_adjugate(F)
    assert norm_F % p**(N - 1) == 0
    tau = norm_F // p**(N - 1)
    assert tau == 1217
    assert all(tau % d for d in range(2, isqrt(tau) + 1))
    assert tau < p and p % tau != 0
    divisor = p**(N - 2)
    assert all(c % divisor == 0 for c in adjugate_F)
    reciprocal = [c // divisor for c in adjugate_F]
    assert multiply(F, reciprocal) == [tau * p] + [0] * (N - 1)
    assert evaluate(reciprocal, support_root, p) == 0
    reciprocal_norm = norm_adjugate(reciprocal)[0]
    assert reciprocal_norm == p * tau**(N - 1)

    generator_relation = [2, -1] + [0] * (N - 2)
    assert norm_adjugate(generator_relation)[0] == p * auxiliary_prime
    assert evaluate(generator_relation, 2, p) == 0
    assert evaluate(generator_relation, 2, auxiliary_prime) == 0
    carry_numerator = multiply(conjugate(generator_relation), F)
    assert all(c % p == 0 for c in carry_numerator)
    carry = [c // p for c in carry_numerator]
    assert norm_adjugate(carry)[0] == auxiliary_prime * tau == 780097

    gram = multiply(F, conjugate(F))
    assert all(c % p == 0 for c in gram)
    assert gram != [p * p] + [0] * (N - 1)
    index = (auxiliary_prime - 1) // n
    assert index == 10 and auxiliary_prime < 26**2
    gauss_strict_upper = Fraction(1 + (index - 1) * 26, index)
    principal_strict_lower = Fraction(n) - Fraction(1936, 49)
    assert gauss_strict_upper == Fraction(47, 2)
    assert principal_strict_lower == Fraction(1200, 49)
    assert principal_strict_lower - gauss_strict_upper == Fraction(97, 98) > 0

    result = {
        "scope": "Exact finite V<1 counterexample to the implication that an evaluation kernel is principal; no uniform shell theorem.",
        "p": p, "n": n, "N": N, "g": g, "a": a,
        "support_root_of_F_mod_p": support_root,
        "F_centered": F,
        "F_at_odd_roots_in_g_order": root_values,
        "sum_centered_squares": S,
        "V": rational(V),
        "norm_F": norm_F,
        "tau": tau,
        "tau_is_prime": True,
        "F_times_conjugate_F_div_p": [c // p for c in gram],
        "scalar_reciprocal_T": reciprocal,
        "F_times_T": [tau * p] + [0] * (N - 1),
        "norm_T": reciprocal_norm,
        "least_positive_integer_m_for_integral_m_p_over_F": tau,
        "least_scalar_reciprocal_norm_cofactor": tau**(N - 1),
        "auxiliary_prime": auxiliary_prime,
        "generator_relation": generator_relation,
        "norm_generator_relation": p * auxiliary_prime,
        "carry_C_with_conjugate_2_minus_X_times_F_equal_p_C": carry,
        "norm_C": auxiliary_prime * tau,
        "auxiliary_subgroup_index": index,
        "auxiliary_Gauss_strict_upper": rational(gauss_strict_upper),
        "auxiliary_principal_kernel_strict_lower": rational(principal_strict_lower),
        "auxiliary_contradiction_margin": rational(Fraction(97, 98)),
        "analytic_inputs": [
            "Norm(f)=q in the order-64 evaluation kernel gives eta>=64-4*pi^2, by pass39",
            "pi<22/7",
            "multiplicative-character expansion gives M<=(1+(index-1)*sqrt(q))/index",
        ],
        "algebraic_inputs": [
            "If the p evaluation kernel is fR of norm p, then (2-X)/f is integral of norm 641",
            "A norm-641 quotient must vanish at 2 modulo 641 because f is invertible modulo 641",
            "Norm(m*p/F)=p*m^32/1217 is integral only if the prime 1217 divides m",
        ],
    }
    target = Path(__file__).resolve().parents[1] / "results/parallel40_shell_converse_2026_09_06.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "p": p, "n": n, "V": str(V), "tau": tau,
        "principal_kernel_converse": "refuted by the ordinary algebraic argument and exact inequalities",
        "least_scalar_reciprocal_cofactor": str(tau**(N - 1)),
        "certificate": str(target),
    }, indent=2))


if __name__ == "__main__":
    main()
