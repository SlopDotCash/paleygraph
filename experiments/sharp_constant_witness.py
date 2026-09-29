#!/usr/bin/env python3
"""Exact rational lower witness; floating outputs are labeled diagnostics only.

The proof that this polynomial bounds cosine, and the logarithm comparison,
are kernel checked in research/SharpConstantCounterexample.lean.
"""
from fractions import Fraction
from math import cos, fsum, log, pi, sqrt
import json
from pathlib import Path


def cosine_lower(u):
    c = 1 - u * u / 32
    d = 2 * c * c - 1
    return 2 * d * d - 1 if d >= 0 else Fraction(-1)


def main():
    p, n = 6700417, 64
    H = sorted({pow(2, j, p) for j in range(n)})
    assert len(H) == n and pow(2, n, p) == 1
    assert n ** 4 <= 4 * p and p <= n ** 4
    radii = [Fraction(63 * min(h, p - h), 10 * p) for h in H]
    assert all(0 <= u <= 4 for u in radii)
    lower = sum(map(cosine_lower, radii), Fraction(0))
    assert lower > 43
    # Combined with ln(p/n) < 12, this disproves M^2 <= 2n ln(p/n).
    assert 43 ** 2 > 2 * n * 12
    # Conversely, the independently computed M^2 <= 1970 implies M^2 < 2n ln p:
    # p^2 > 2^45 and ln 2 > 0.69 give 2n ln p > 1987.2.
    assert p ** 2 > 2 ** 45
    assert 2 * n * Fraction(45, 2) * Fraction(69, 100) > 1970
    # These approximate values do not participate in certificate checks.
    approx = fsum(cos(2 * pi * min(h, p - h) / p) for h in H)
    output = {
        "p": p, "n": n, "generator": 2, "frequency": 1,
        "subgroup": H,
        "exact_rational_cosine_lower": str(lower),
        "proved_period_real_lower_strict": 43,
        "proved_log_p_over_n_upper_strict": 12,
        "literal_natural_log_constant_sqrt2_refuted": True,
        "unspecified_constant_conjecture_refuted": False,
        "diagnostics_only": {
            "rational_lower_float": float(lower),
            "period_at_frequency_1_float": approx,
            "normalized_value_float": approx / sqrt(n * log(p / n)),
            "sqrt2_target_float": sqrt(2 * n * log(p / n)),
        },
    }
    path = Path(__file__).resolve().parents[1] / "results/sharp_constant_witness.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(f"Exact rational lower > 43; wrote {path}")
    print(json.dumps(output["diagnostics_only"], indent=2))


if __name__ == "__main__":
    main()
