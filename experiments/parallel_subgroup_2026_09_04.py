"""Exact checks for the index-two projection bound; no floating arithmetic."""

from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def prime(p):
    return p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))


def subgroup(p, n):
    assert prime(p) and (p - 1) % n == 0 and n >= 4 and n & (n - 1) == 0
    for a in range(2, p):
        g = pow(a, (p - 1) // n, p)
        if pow(g, n // 2, p) == p - 1:
            h = [pow(g, i, p) for i in range(n)]
            assert len(set(h)) == n and pow(g, n, p) == 1
            return g, h[::2], h[1::2]
    raise AssertionError("No generator")


def conv(a, b, p):
    out = Counter()
    for x, ax in a.items():
        for y, by in b.items():
            out[(x + y) % p] += ax * by
    return out


def dot(a, b):
    return sum(v * b.get(x, 0) for x, v in a.items())


def fourth_check(p, n):
    g, kk, ll = subgroup(p, n)
    k = n // 2
    a, b = Counter(kk), Counter(ll)
    aa, ab, bb = conv(a, a, p), conv(a, b, p), conv(b, b, p)
    e, mixed, t = dot(aa, aa), dot(ab, ab), dot(aa, ab)
    assert mixed == dot(aa, bb)
    assert mixed <= e - k and 2 * mixed <= k**3
    assert all((aa.get(x, 0) - bb.get(x, 0)) % 2 == 1 for x in {(2 * h) % p for h in kk + ll})
    assert all(ab.get(g * x % p, 0) == v for x, v in ab.items())
    assert dot(aa, ab) == dot(bb, ab)
    parent = conv(a + b, a + b, p)
    ep = dot(parent, parent)
    assert ep == 2 * e + 6 * mixed + 8 * t
    old_rhs = (e - k * k) * mixed
    new_rhs2 = (e + mixed - 2 * k * k) * mixed
    assert 2 * t * t <= new_rhs2 <= 2 * old_rhs
    d = p - 1
    va2 = d * (e + mixed - 2 * k * k) - 2 * (k * k - k) ** 2
    vb = d * mixed - k**4
    centered = d * t - k * k * (k * k - k)
    assert va2 >= 0 and vb >= 0 and 2 * centered**2 <= va2 * vb
    return {
        "p": p, "parent_order": n, "child_order": k,
        "quartic_window_for_parent": n**4 <= 4 * p <= 4 * n**4,
        "E": e, "B": mixed, "T": t, "parent_energy": ep,
        "old_T_squared_upper": old_rhs,
        "new_T_squared_upper": str(Fraction(new_rhs2, 2)),
        "projection_saving": str(Fraction(2 * old_rhs - new_rhs2, 2)),
        "centered_gap_scaled": va2 * vb - 2 * centered**2,
    }


def depth_checks(p, n):
    g, kk, ll = subgroup(p, n)
    k, d = n // 2, p - 1
    walks = [Counter({0: 1})]
    for _ in range(12):
        walks.append(conv(walks[-1], Counter(kk), p))
    shifted = [Counter({g * x % p: v for x, v in r.items()}) for r in walks]

    def z(a, b):
        return dot(walks[a], shifted[b])

    checks = 0
    for q in (1, 2, 4, 8):
        odd_support = {x for x, v in walks[q].items() if v % 2}
        assert odd_support == {q * x % p for x in kk}
        assert dot(walks[q], walks[q]) - z(q, q) >= k
    for j in range(1, 9):
        for s in range(1, 5):
            f0, w0 = walks[j].get(0, 0), z(s, s)
            ej = dot(walks[j], walks[j])
            u_mass, w_mass = k**j - f0, k**(2 * s) - w0
            va2 = d * (ej + z(j, j) - 2 * f0**2) - 2 * u_mass**2
            vb = d * (z(2 * s, 2 * s) - w0**2) - w_mass**2
            cc = d * (z(j + s, s) - f0 * w0) - u_mass * w_mass
            assert va2 >= 0 and vb >= 0 and 2 * cc**2 <= va2 * vb
            # Independently construct the invariant convolution in this check.
            w = conv(walks[s], shifted[s], p)
            assert w.get(0, 0) == w0 and dot(walks[j], w) == z(j + s, s)
            assert dot(w, w) == z(2 * s, 2 * s)
            assert all(w.get(g * x % p, 0) == v for x, v in w.items())
            checks += 1
    return checks


def main():
    cases = []
    for p in range(5, 258, 4):
        if prime(p):
            n = 4
            while n <= 64 and (p - 1) % n == 0:
                cases.append((p, n))
                n *= 2
    cases += [(1153, 8), (18433, 16), (262657, 32), (6700417, 64), (67403009, 128)]
    fourth = [fourth_check(p, n) for p, n in cases]
    higher_cases = [(p, n) for p, n in cases if p <= 97 and n <= 32]
    higher_count = sum(depth_checks(p, n) for p, n in higher_cases)
    script = ROOT / "experiments/parallel_subgroup_2026_09_04.py"
    result = {
        "status": "Exact finite checks of proved projection inequalities; no uniform mixed-energy or moment bound proved.",
        "fourth_check_count": len(fourth),
        "fourth_cases_with_nonminimal_child_energy": sum(r["E"] > 3 * r["child_order"]**2 - 3 * r["child_order"] for r in fourth),
        "strict_projection_improvements": sum(Fraction(r["projection_saving"]) > 0 for r in fourth),
        "higher_depth_case_count": len(higher_cases),
        "higher_depth_identity_and_inequality_checks": higher_count,
        "dyadic_depth_parity_checks": 4 * len(higher_cases),
        "fourth_checks": fourth,
        "source_sha256": {str(script.relative_to(ROOT)): hashlib.sha256(script.read_bytes()).hexdigest()},
    }
    path = ROOT / "results/parallel_subgroup_2026_09_04.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("fourth_checks", "source_sha256")}, indent=2))
    print(json.dumps(fourth[-2:], indent=2))


if __name__ == "__main__":
    main()
