#!/usr/bin/env python3
"""Exact all-order weighted-model checks. These rows are NOT Paley kernels."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, isqrt, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def prime(p):
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    return all(p % d for d in range(3, isqrt(p)+1, 2))


def next_prime(p):
    while not prime(p):
        p += 1
    return p


def floor_root(p, q):
    low, high = 0, 1
    while high**q <= p:
        high *= 2
    while low+1 < high:
        mid = (low+high)//2
        if mid**q <= p:
            low = mid
        else:
            high = mid
    return low


def double_factorial(s):
    return prod(range(1, 2*s, 2))


@cache
def U(s, m):
    assert s >= 0 and m >= 0
    if s == 0:
        return Q(1)
    return Q(sum(comb(m, j)*(2*j-m)**(2*s) for j in range(m+1)), 2**m)


def V(s, m):
    return (U(s+1, m)-m*U(s, m))/2


def model(k, p):
    assert k >= 2 and p >= k**4
    L = k*k
    W = p-k-L
    alpha = Q(L+1, W)
    assert W > 0
    assert W-(k*k+1)*k*(k-1) >= k*k*(k-2) >= 0
    assert W-(k*k+1)*k >= k*(k-2)*(k*k+k+1) >= 0
    assert alpha*comb(k, 2) <= Q(1, 2) and alpha*k <= 1
    return k, p, L, W, alpha


def expected_corr(k, p, d):
    if d == 0:
        return Q(p)
    if d % 2:
        return Q(0)
    if d == 2:
        return Q(-1)
    return Q(k*k)


def subset_moment(k, p, m, s):
    assert 1 <= m <= k
    L = k*k
    return L*m**(2*s)+(p-L-m)*U(s, m)-(L+1)*V(s, m)+m*U(s, m-1)


def explicit_rows(k, p):
    k, p, L, W, alpha = model(k, p)
    rows = defaultdict(Q)
    for eps in product((-1, 1), repeat=k):
        e2 = sum(eps[i]*eps[j] for i in range(k) for j in range(i+1, k))
        rho = 1-alpha*e2
        assert Q(1, 2) <= rho <= Q(3, 2)
        rows[eps] += Q(W, 2**k)*rho
        COUNTS["explicit_cube_density_bounds"] += 1
    for sign in (-1, 1):
        rows[(sign,)*k] += Q(L, 2)
    for zero in range(k):
        for eps in product((-1, 1), repeat=k-1):
            row = eps[:zero]+(0,)+eps[zero:]
            rows[row] += Q(1, 2**(k-1))
    assert all(w > 0 for w in rows.values())
    return rows


def check_small_rows():
    reports = []
    for k in range(2, 9):
        p = next_prime(k**4)
        rows = explicit_rows(k, p)
        assert sum(rows.values()) == p
        for eps, w in rows.items():
            assert rows[tuple(-e for e in eps)] == w
        for i in range(k):
            assert sum(w for eps, w in rows.items() if not eps[i]) == 1
            assert sum(w*eps[i] for eps, w in rows.items()) == 0
            COUNTS["explicit_column_mass_and_mean"] += 2
            for j in range(k):
                assert sum(w*eps[i]*eps[j] for eps, w in rows.items()) == p*int(i == j)-1
                COUNTS["explicit_gram_entries"] += 1
        for d in range(k+1):
            for D in combinations(range(k), d):
                corr = sum(w*prod(eps[i] for i in D) for eps, w in rows.items())
                assert corr == expected_corr(k, p, d)
                COUNTS["explicit_all_subset_correlations"] += 1
        subsets = [D for m in range(1, k+1) for D in combinations(range(k), m)] if k <= 6 else [tuple(range(m)) for m in range(1, k+1)]
        for B in subsets:
            m = len(B)
            values = [(sum(eps[i] for i in B), w) for eps, w in rows.items()]
            for s in range(1, 5):
                moment = sum(w*f**(2*s) for f, w in values)
                assert moment == subset_moment(k, p, m, s)
                COUNTS["explicit_unsigned_moment_formulas"] += 1
                if s <= 2:
                    assert moment <= (double_factorial(s)+1)*p*m**s
                    COUNTS["explicit_unsigned_lower_moment_bounds"] += 1
            assert subset_moment(k, p, m, 1) == p*m-m*m
        vectors = set(product((-1, 0, 1), repeat=k)) if k <= 4 else {
            tuple(1 for _ in range(k)),
            tuple((-1)**i for i in range(k)),
            tuple(i-k//2 for i in range(k)),
            tuple(i % 3-1 for i in range(k)),
            *[tuple(int(i == j) for i in range(k)) for j in range(k)],
        }
        vectors.add(tuple(Q(1, i+1) for i in range(k)))
        for b in vectors:
            norm2 = sum(v*v for v in b)
            values = [(sum(v*e for v, e in zip(b, eps)), w) for eps, w in rows.items()]
            for s in (1, 2):
                moment = sum(w*f**(2*s) for f, w in values)
                assert moment <= (Q(3, 2)*double_factorial(s)+1)*p*norm2**s
                if s == 1:
                    assert moment == p*norm2-sum(b)**2
                COUNTS["explicit_real_coefficient_lower_bounds"] += 1
        canonical = [[list(eps), str(w)] for eps, w in sorted(rows.items())]
        reports.append({
            "k": k, "p": p, "rows": len(rows),
            "row_measure_sha256": sha256(json.dumps(canonical).encode()).hexdigest(),
            "coefficient_vectors": len(vectors), "unsigned_subsets": len(subsets),
            "sixth_ratio": str(subset_moment(k, p, k, 3)/(p*k**3)),
            "scope": "Actual rational weighted vectors, not prime-field character rows. Algebra fixtures include k below asymptotic threshold."
        })
        COUNTS["explicit_weighted_rows"] += len(rows)
    return reports


def choose(n, j):
    return comb(n, j) if 0 <= j <= n else 0


def elementary(F, N, d):
    a, b = (N+F)//2, (N-F)//2
    return sum((-1)**j*choose(b, j)*choose(a, d-j) for j in range(d+1))


def compressed_corr(k, p, d):
    k, p, L, W, alpha = model(k, p)
    cube_sum = Q(0)
    for j in range(k+1):
        F = 2*j-k
        rho = 1-alpha*Q(F*F-k, 2)
        cube_sum += Q(comb(k, j), 2**k)*rho*Q(elementary(F, k, d), comb(k, d))
    # Direct sum classes for zero rows, rather than assuming their cancellation.
    zero_sum = Q(0)
    if d < k:
        for j in range(k):
            F = 2*j-(k-1)
            zero_sum += Q(comb(k-1, j), 2**(k-1))*Q(elementary(F, k-1, d), comb(k-1, d))
        zero_sum *= k-d
    extreme = L if d % 2 == 0 else 0
    return extreme+W*cube_sum+zero_sum


def check_cube_moments():
    for m in range(65):
        for s in range(1, 13):
            assert V(s, m) >= 0
            assert U(s, m) <= double_factorial(s)*m**s
            COUNTS["binomial_pairing_and_covariance_checks"] += 2
    for k in range(2, 33):
        p = next_prime(k**4)
        for d in range(k+1):
            assert compressed_corr(k, p, d) == expected_corr(k, p, d)
            COUNTS["compressed_all_degree_correlations"] += 1


def check_large_slices():
    cases = [(3, k) for k in (8, 16, 32, 64, 128, 256)]
    cases += [(4, k) for k in (10, 16, 32)]
    cases += [(5, k) for k in (12, 16, 24)]
    cases += [(6, 14), (6, 16), (7, 16), (8, 18), (9, 20)]
    reports = []
    for r, target in cases:
        p = next_prime(target**(r+1))
        k = floor_root(p, r+1)
        assert prime(p) and k == target
        assert k >= 2*(r+1)
        model(k, p)
        assert k**(r+1) <= p < (k+1)**(r+1) < 2*k**(r+1)
        degrees = sorted(set(range(min(k, 16)+1)) | {k-1, k})
        for d in degrees:
            c = compressed_corr(k, p, d)
            assert c == expected_corr(k, p, d)
            if d:
                assert c*c <= (d-1)**2*p
                COUNTS["large_Weil_shaped_degree_bounds"] += 1
        for m in range(1, k+1):
            for s in range(1, r):
                moment = subset_moment(k, p, m, s)
                assert 0 <= moment <= (double_factorial(s)+1)*p*m**s
                COUNTS["large_unsigned_lower_moment_bounds"] += 1
        top = subset_moment(k, p, k, r)
        assert top >= k**(2*r+2)
        assert 2*top >= p*k**(r+1)
        lower_moments = {str(s): str(subset_moment(k, p, k, s)/(p*k**s)) for s in range(1, r)}
        hs = []
        for h in range(2, (r+1)//2+1):
            forbidden = (k-1)+sum((k-1)**(2*h-j) for j in range(1, h+1))
            assert p > h and forbidden < p
            assert forbidden < (h+1)*k**(2*h-1) < k**(2*h) <= p
            hs.append(h)
            COUNTS["B_h_label_existence_thresholds"] += 1
        reports.append({
            "r": r, "p": p, "k": k, "checked_degrees": degrees,
            "top_moment_ratio": str(top/(p*k**r)), "ratio_lower_bound": str(Q(k, 2)),
            "lower_moment_ratios": lower_moments, "allowed_relation_orders_checked": hs,
            "scope": "Exact compressed weighted measure on a genuine prime-size slice; no actual character kernel is evaluated."
        })
        COUNTS["large_prime_size_slices"] += 1
    return reports


def check_label_fixtures():
    reports = []
    for h in (2, 3, 4):
        r = 2*h-1
        for k in (4, 5):
            p = next_prime(k**(r+1))
            labels = [(h+1)**j for j in range(k)]
            assert h*max(labels) < p
            sums = [sum(t) % p for t in combinations_with_replacement(labels, h)]
            assert len(set(sums)) == len(sums) == comb(k+h-1, h)
            COUNTS["explicit_B_h_multiset_sums"] += len(sums)
            reports.append({"r": r, "h": h, "k": k, "p": p, "labels": labels, "distinct_h_sums": len(sums),
                            "scope": "Label-only finite fixture, below the asymptotic top-moment threshold; no link to sign rows is imposed."})
    return reports


def main():
    small = check_small_rows()
    print("Explicit weighted rows and coefficient checks passed.", flush=True)
    check_cube_moments()
    large = check_large_slices()
    labels = check_label_fixtures()
    inputs = [
        "research/parallel22-all-orders-obstruction-2026-09-05.md",
        "experiments/parallel22_all_orders_obstruction_2026_09_05.py",
        "research/parallel21-squarefree-moments-2026-09-05.md",
        "research/parallel19-relation-free-reduction-2026-09-05.md",
    ]
    data = {
        "status": "All exact checks passed. Weighted relaxation only; no actual Paley upper bound.",
        "arithmetic": "Python integers and Fraction; exact integer-root and trial-primality checks; no floating-point acceptance.",
        "counts": dict(COUNTS), "explicit_small_measures": small, "compressed_prime_slices": large,
        "explicit_relation_free_labels": labels,
        "input_sha256": {f:sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs},
        "limitations": [
            "The positive rational measure is not uniform sampling over p prime-field translates.",
            "Finite checks do not prove the asymptotic claims; those have a separate prose proof.",
            "Arbitrary real coefficient bounds are proved by cube density and pairings, not by enumerating all real vectors.",
            "B_h labels do not supply field-difference or multiplicative-character identities.",
            "All full classical, subgroup, spectral and official prize targets remain unproved.",
        ],
    }
    path = ROOT/"results/parallel22_all_orders_obstruction_2026_09_05.json"
    path.write_text(json.dumps(data, indent=2)+"\n")
    print(json.dumps({"status":data["status"],"counts":data["counts"],"output":str(path)},indent=2),flush=True)


if __name__ == "__main__":
    main()

