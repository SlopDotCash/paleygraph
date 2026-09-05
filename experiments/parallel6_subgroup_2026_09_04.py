#!/usr/bin/env python3
"""Bounded exact certificates for the pass-six centered bilinear theorem.

Python standard library only. All inequality decisions use Fraction intervals;
decimal displays in the JSON are diagnostic, never used to decide a check.
This is a verifier of fixed cases and algebra, not a uniform-proof substitute.
"""

from fractions import Fraction as Q
from hashlib import sha256
from math import factorial, isqrt
from pathlib import Path
import json
import time


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/parallel6_subgroup_2026_09_04.json"
SCALE = 10**25
PI_SCALE = 10**30
CASES = [(5, 2), (5, 4), (13, 4), (17, 4), (17, 8), (17, 16),
         (41, 8), (73, 8)]
ORDERS = [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3),
          (3, 4), (4, 4)]


def encode(q):
    q = Q(q)
    return str(q.numerator) if q.denominator == 1 else str(q)


def floor_scaled(q, scale):
    return (q.numerator * scale) // q.denominator


def ceil_scaled(q, scale):
    return -((-q.numerator * scale) // q.denominator)


def outward(lo, hi, scale=SCALE):
    assert lo <= hi
    return Q(floor_scaled(lo, scale), scale), Q(ceil_scaled(hi, scale), scale)


def atan_interval(inv, terms):
    """Alternating-series bracket for arctan(1/inv), 0 < 1/inv < 1."""
    value = sum((Q((-1)**j, (2*j + 1) * inv**(2*j + 1))
                 for j in range(terms)), Q(0))
    next_term = Q((-1)**terms, (2*terms + 1) * inv**(2*terms + 1))
    return min(value, value + next_term), max(value, value + next_term)


def certified_pi():
    # tan(4 atan(1/5)) = 120/119, and tan(atan(1/239)) = 1/239.
    # The subtraction formula gives tan(4 atan(1/5)-atan(1/239))=1.
    # The angle is positive and < 4/5 < pi/2, identifying it with pi/4.
    tangent = (Q(120, 119) - Q(1, 239)) / (1 + Q(120, 119*239))
    assert tangent == 1
    a_lo, a_hi = atan_interval(5, 40)
    b_lo, b_hi = atan_interval(239, 12)
    lo, hi = outward(16*a_lo - 4*b_hi, 16*a_hi - 4*b_lo, PI_SCALE)
    assert Q(3) < lo < hi < Q(22, 7)
    return lo, hi


PI_LO, PI_HI = certified_pi()


def cosine_interval(k, p):
    """Rigorous interval for cos(2*pi*k/p), via Taylor and Lipschitz."""
    k %= p
    k = min(k, p-k)
    if not k:
        return Q(1), Q(1)
    angle_lo, angle_hi = Q(2*k, p)*PI_LO, Q(2*k, p)*PI_HI
    midpoint = (angle_lo + angle_hi) / 2
    term = Q(1)
    polynomial = term
    degree_half = 24
    for j in range(1, degree_half + 1):
        term *= -midpoint**2 / ((2*j-1)*(2*j))
        polynomial += term
    # Treat the Taylor polynomial as degree 2m+1 with zero odd coefficient.
    taylor_error = abs(midpoint)**(2*degree_half+2) / factorial(2*degree_half+2)
    angle_error = (angle_hi-angle_lo)/2  # cos is 1-Lipschitz.
    error = taylor_error + angle_error
    lo, hi = outward(polynomial-error, polynomial+error)
    return max(Q(-1), lo), min(Q(1), hi)


def sqrt_interval(q):
    """Exact outward bracket; integer comparisons certify both endpoints."""
    q = Q(q)
    assert q >= 0
    if q == 0:
        return Q(0), Q(0)
    root = isqrt((q.numerator*SCALE*SCALE)//q.denominator)
    lo = Q(root, SCALE)
    hi = lo if lo*lo == q else Q(root+1, SCALE)
    assert lo*lo <= q <= hi*hi
    return lo, hi


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def distinct_prime_factors(n):
    output = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            output.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        output.append(n)
    return output


def subgroup(p, n):
    assert is_prime(p) and (p-1) % n == 0
    factors = distinct_prime_factors(p-1)
    generator = next(g for g in range(2, p)
                     if all(pow(g, (p-1)//q, p) != 1 for q in factors))
    hgen = pow(generator, (p-1)//n, p)
    elements = sorted({pow(hgen, j, p) for j in range(n)})
    assert len(elements) == n and p-1 in elements
    return hgen, elements


def sum_walks(p, elements, highest=4):
    walk = [0]*p
    walk[0] = 1
    walks = {}
    for j in range(1, highest+1):
        nxt = [0]*p
        for x, count in enumerate(walk):
            for h in elements:
                nxt[(x+h) % p] += count
        walk = nxt
        assert sum(walk) == len(elements)**j
        walks[j] = walk
    return walks


def walk_statistics(p, n, walks):
    stats = {}
    for j, walk in walks.items():
        measure = [Q(c, n**j) for c in walk]
        origin = measure[0]
        mass = 1-origin
        centered = [measure[x]-mass/(p-1) for x in range(1, p)]
        variance = sum((c*c for c in centered), Q(0))
        energy = sum(c*c for c in walk)
        principal_removed = Q(energy) - Q(n**(2*j), p)
        assert Q(0) <= origin <= Q(1, n)
        assert sum(centered) == 0
        assert variance >= 0
        assert variance == (principal_removed/n**(2*j)
                            - Q(p, p-1)*(origin-Q(1, p))**2)
        stats[j] = {
            "u": origin, "m": mass, "V": variance, "E": energy,
            "E_centered": principal_removed, "measure": measure,
            "nu": centered,
        }
    return stats


def representative_frequencies(p, elements):
    unused = set(range(1, p))
    representatives = []
    while unused:
        a = min(unused)
        representatives.append(a)
        coset = {(a*h) % p for h in elements}
        assert coset <= unused
        unused -= coset
    assert len(representatives)*len(elements) == p-1
    return representatives


def absolute_period_interval(a, p, elements, cosines):
    n = len(elements)
    if n == p-1:
        # The sum of all nontrivial pth roots is exactly -1.
        return Q(1, n), Q(1, n)
    low = sum((cosines[(a*h) % p][0] for h in elements), Q(0))/n
    high = sum((cosines[(a*h) % p][1] for h in elements), Q(0))/n
    if low >= 0:
        return low, high
    if high <= 0:
        return -high, -low
    return Q(0), max(-low, high)


def check_centering_bilinear_coefficients(p, n, stats):
    """Check (9) modulo 1+X+...+X^(p-1), for every pair of orders.

    Additive root sums are encoded by their rational coefficient vector.
    A constant vector evaluates to zero. No approximate complex arithmetic.
    """
    checks = 0
    for r, s in ORDERS:
        left = [Q(0)]*p
        centered = [Q(0)]*p
        for x in range(1, p):
            for y in range(1, p):
                coefficient = (x*y) % p
                left[coefficient] += stats[r]["measure"][x]*stats[s]["measure"][y]
                centered[coefficient] += stats[r]["nu"][x-1]*stats[s]["nu"][y-1]
        # Difference B - centered + m_r m_s/(p-1) must evaluate to zero.
        difference = [left[j]-centered[j] for j in range(p)]
        difference[0] += stats[r]["m"]*stats[s]["m"]/(p-1)
        assert len(set(difference)) == 1
        checks += 1
    return checks


def check_field(p, n):
    generator, elements = subgroup(p, n)
    walks = sum_walks(p, elements)
    stats = walk_statistics(p, n, walks)
    representatives = representative_frequencies(p, elements)
    # All walks are invariant under H; this extends coset tests to all a != 0.
    for j, walk in walks.items():
        for h in elements:
            assert all(walk[(x*h) % p] == walk[x] for x in range(p))
    cosines = {k: cosine_interval(k, p) for k in range(p)}
    intervals = {a: absolute_period_interval(a, p, elements, cosines)
                 for a in representatives}
    minimum_margin = None
    core_checks = 0
    exact_equalities = 0
    pentagon_equalities = 0
    for a, (_, delta_hi) in intervals.items():
        assert 0 <= delta_hi <= 1
        for r, s in ORDERS:
            R, S = stats[r], stats[s]
            atom_removed = R["u"] + S["u"] - R["u"]*S["u"]
            atom_direct = sum((R["measure"][x]*S["measure"][y]
                               for x in range(p) for y in range(p)
                               if x == 0 or y == 0), Q(0))
            assert atom_removed == atom_direct
            lhs_hi = max(delta_hi**(r*s)-atom_removed, Q(0))**2
            root_lo, _ = sqrt_interval(p*R["V"]*S["V"])
            rhs_lo = (Q(p-n, n*(p-1))*R["m"]**2*S["m"]**2
                      + Q(n-1, n)*R["m"]*S["m"]*root_lo)
            if lhs_hi > rhs_lo:
                # The sharp F_5 case is an equality. Outward intervals cannot
                # prove equality; use the exact pentagon period instead.
                assert (p, n, a, r, s) == (5, 2, 2, 1, 1), (
                    p, n, a, r, s, float(lhs_hi-rhs_lo))
                assert R["u"] == S["u"] == 0
                assert R["V"] == S["V"] == Q(1, 4)
                delta_exact = (Q(1, 4), Q(1, 4))
                lhs_exact = q5_mul(delta_exact, delta_exact)
                rhs_exact = (Q(p-n, n*(p-1)), Q(n-1, n)*Q(1, 4))
                assert lhs_exact == rhs_exact
                margin = Q(0)
                pentagon_equalities += 1
            else:
                margin = rhs_lo-lhs_hi
            minimum_margin = margin if minimum_margin is None else min(margin, minimum_margin)
            exact_equalities += margin == 0 and n == p-1
            core_checks += 1
    centering_checks = check_centering_bilinear_coefficients(p, n, stats)
    return {
        "prime": p, "subgroup_size": n, "subgroup_generator": generator,
        "elements": elements, "frequency_coset_representatives": representatives,
        "orders": ORDERS,
        "core_certificates": core_checks,
        "strict_interval_certificates": core_checks-exact_equalities-pentagon_equalities,
        "covered_nonzero_frequency_order_pairs": (p-1)*len(ORDERS),
        "exact_full_group_equalities": exact_equalities,
        "exact_pentagon_equalities": pentagon_equalities,
        "centering_polynomial_certificates": centering_checks,
        "minimum_rational_margin": encode(minimum_margin),
        "minimum_margin_decimal_diagnostic": float(minimum_margin),
        "largest_period_interval_width": encode(max(hi-lo for lo, hi in intervals.values())),
        "moments": [
            {"order": j, "origin_count": walks[j][0],
             "origin_probability": encode(st["u"]), "energy": st["E"],
             "centered_energy": encode(st["E_centered"]),
             "nonzero_variance": encode(st["V"])}
            for j, st in stats.items()
        ],
    }


def exponent_ledger():
    def spread(j):
        return Q(j) if j <= 2 else Q(j+24, 9)

    def saving(r, s):
        return min(Q(1), max(Q(0), (spread(r)+spread(s)-4)/2))/(2*r*s)

    ledger = {(r, s): saving(r, s) for r in range(1, 65) for s in range(1, 65)}
    best = max(ledger.values())
    optimizers = [list(pair) for pair, value in ledger.items() if value == best]
    assert best == Q(1, 18) and optimizers == [[3, 3]]
    row1 = [saving(1, s) for s in range(1, 65)]
    assert max(row1) == Q(1, 42) and row1.index(max(row1))+1 == 21
    row2 = [saving(2, s) for s in range(1, 65)]
    assert max(row2) == Q(1, 24) and row2.index(max(row2))+1 == 3
    # Derivative signs of the two piecewise rational formulas establish tails.
    # r=1: (s-3)/(36s) has derivative 1/(12s^2) before s=21.
    # r=2: (s+6)/(72s) has derivative -1/(12s^2) before s=12.
    assert saving(1, 21) == Q(1, 42) and saving(2, 12) == Q(1, 48)
    assert Q(23, 24)-Q(17, 18) == Q(1, 72)
    assert Q(59, 12)-Q(44, 9) == Q(1, 36)
    higher = []
    for j in range(3, 33):
        n_exponent = (2*j-6)*Q(17, 18)+3
        b_exponent = (2*j-6)*Q(1, 18)+1
        two_exponent = Q(2*j-6, 6)
        assert n_exponent == Q(17*j-24, 9)
        assert b_exponent == Q(j+6, 9)
        assert two_exponent == Q(j-3, 3)
        higher.append({"order": j, "n_exponent": encode(n_exponent),
                       "B_exponent": encode(b_exponent),
                       "two_exponent": encode(two_exponent)})
    # Sixth-energy explicit constant:
    # 2/n <= sqrt((B+1)/n) iff 4 <= n(B+1); n>=2 and B>=1 suffice.
    # Four times (B+1) <= eight B, and 8^(1/18)=2^(1/6).
    assert 2*(1+1) == 4 and 4*(1+1) == 8*1
    assert Q(3, 18) == Q(1, 6)
    # Eighth-energy explicit constant: 1+sqrt(2)<3, square positive sides.
    assert 2 < (3-1)**2
    assert Q(1)-Q(1, 16) == Q(15, 16)
    return {
        "new_M_exponent": "17/18", "old_class_M_exponent": "23/24",
        "improvement_in_exponent": "1/72", "new_B_exponent": "1/18",
        "E4_feedback_exponent": "44/9", "feedback_exponents": higher,
        "finite_pairs_checked": len(ledger), "max_saving": encode(best),
        "maximizers": optimizers,
        "analytic_tail_reduction": "r,s>=3: 1/(2rs); r=1 peaks at s=21; r=2 peaks at s=3",
        "conditional_centered_E4_M_exponent": "15/16",
        "conditional_centered_E4_D_exponent": "1/16",
        "explicit_constants": ["sixth: 2^(1/6)", "eighth: 3^(1/4)"],
    }


def q5_mul(a, b):
    return (a[0]*b[0]+5*a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def q5_positive(pair):
    rational, radical = pair
    if radical >= 0 and rational >= 0:
        return rational > 0 or radical > 0
    if radical <= 0 and rational <= 0:
        return False
    return (5*radical*radical > rational*rational if radical > 0
            else rational*rational > 5*radical*radical)


def atom_counterexample():
    # Exact arithmetic in Q(sqrt(5)), represented by (rational, coefficient).
    delta = (Q(1, 4), Q(1, 4))
    delta_squared = q5_mul(delta, delta)
    uncorrected = q5_mul(delta_squared, delta_squared)
    corrected_base = (delta_squared[0]-Q(1, 2), delta_squared[1])
    assert q5_positive(corrected_base)
    corrected = q5_mul(corrected_base, corrected_base)
    rhs = (Q(3, 32), Q(1, 32))
    assert uncorrected == (Q(7, 32), Q(3, 32))
    assert corrected == (Q(3, 32), -Q(1, 32))
    assert q5_positive(tuple(x-y for x, y in zip(uncorrected, rhs)))
    assert q5_positive(tuple(x-y for x, y in zip(rhs, corrected)))
    return {
        "prime": 5, "subgroup": [1, 4], "frequency": 2, "orders": [1, 2],
        "delta": "(1+sqrt(5))/4", "correct_rhs": "(3+sqrt(5))/32",
        "uncorrected_lhs": "(7+3*sqrt(5))/32",
        "corrected_lhs": "(3-sqrt(5))/32",
        "claim_refuted": "Deleting zero-atom terms from the exact centered gate",
        "method": "Exact arithmetic and positivity in Q(sqrt(5))",
        "quartic_endpoint_scope": False,
    }


def main():
    started = time.perf_counter()
    fields = [check_field(p, n) for p, n in CASES]
    exponent_checks = exponent_ledger()
    counterexample = atom_counterexample()
    sources = [
        "research/parallel6-subgroup-2026-09-04.md",
        "experiments/parallel6_subgroup_2026_09_04.py",
        "research/parallel5-subgroup-2026-09-04.md",
        "research/parallel4-subgroup-2026-09-04.md",
    ]
    result = {
        "status": "all checks passed",
        "scope": "Fixed exact certificates; uniform theorem is proved in the note, not by enumeration",
        "fixed_prime_subgroup_case_count": len(fields),
        "distinct_prime_field_count": len({p for p, _ in CASES}),
        "walk_moment_count": 4*len(fields),
        "core_certificate_count": sum(f["core_certificates"] for f in fields),
        "strict_interval_certificate_count": sum(f["strict_interval_certificates"] for f in fields),
        "covered_nonzero_frequency_order_pairs": sum(f["covered_nonzero_frequency_order_pairs"] for f in fields),
        "centering_polynomial_certificate_count": sum(f["centering_polynomial_certificates"] for f in fields),
        "exact_full_group_equality_count": sum(f["exact_full_group_equalities"] for f in fields),
        "exact_pentagon_equality_count": sum(f["exact_pentagon_equalities"] for f in fields),
        "pi_certificate": {"lower": encode(PI_LO), "upper": encode(PI_HI),
                           "width": encode(PI_HI-PI_LO),
                           "method": "Machin identity and alternating arctangent intervals"},
        "cosine_certificate": {"Taylor_degree": 48, "outward_rounding_denominator": str(SCALE),
                               "remainder": "|midpoint|^50/50! plus half-angle-interval width"},
        "fields": fields,
        "exponent_ledger": exponent_checks,
        "atom_counterexample": counterexample,
        "source_sha256": {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},
        "elapsed_seconds": round(time.perf_counter()-started, 6),
    }
    RESULT.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key: result[key] for key in [
        "status", "fixed_prime_subgroup_case_count", "distinct_prime_field_count",
        "walk_moment_count", "core_certificate_count",
        "strict_interval_certificate_count",
        "covered_nonzero_frequency_order_pairs", "centering_polynomial_certificate_count",
        "exact_full_group_equality_count", "exact_pentagon_equality_count",
        "elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
