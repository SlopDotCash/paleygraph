#!/usr/bin/env python3
"""Exact checks for pass21. The large weighted fixtures are NOT Paley kernels."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import comb, factorial, isqrt, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def choose(n, j):
    return comb(n, j) if 0 <= j <= n else 0


def elementary(F, N, d):
    """Independent binomial expansion, not the recurrence."""
    assert abs(F) <= N and (N + F) % 2 == 0
    a, b = (N + F) // 2, (N - F) // 2
    return sum((-1) ** j * choose(b, j) * choose(a, d - j)
               for j in range(d + 1))


def hermite_discrete(F, N, d):
    if d == 0:
        return 1
    prev, cur = 1, F
    for j in range(1, d):
        prev, cur = cur, F * cur - j * (N - j + 1) * prev
    return cur


def kappa(F, N, d):
    assert 0 <= d <= N
    return Q(elementary(F, N, d), comb(N, d))


def check_pointwise():
    for N in range(129):
        for F in range(-N, N + 1, 2):
            for r in range(1, 9):
                h = hermite_discrete(F, N, 2 * r)
                assert h == factorial(2 * r) * elementary(F, N, 2 * r)
                COUNTS["recurrence_binomial_comparisons"] += 1
                z = F ** (2 * r)
                assert h >= -(8 * r * N) ** r
                assert abs(h) <= z + (8 * r * N) ** r
                assert z <= 2 ** r * h + 2 * (16 * r * N) ** r
                COUNTS["pointwise_inequalities"] += 3
            h4 = F**4 - (6*N - 8)*F**2 + 3*N*(N-2)
            h6 = (F**6 - (15*N-40)*F**4
                  + (45*N*N-210*N+184)*F**2 - 15*N*(N-2)*(N-4))
            assert h4 == hermite_discrete(F, N, 4)
            assert h6 == hermite_discrete(F, N, 6)
            COUNTS["explicit_polynomial_identities"] += 2


def chi(a, p):
    a %= p
    if not a:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def check_actual_fields():
    reports = []
    for p in (3, 5, 7, 11, 13):
        chars = [chi(a, p) for a in range(p)]
        cache = {}
        # Directly multiply field differences and then apply chi:
        # independent of signs, elementary coefficients, and recurrence.
        for d in (0, 2, 4, 6):
            for subset in combinations(range(p), d):
                cache[subset] = sum(chars[prod((x-c) % p for c in subset) % p]
                                    for x in range(p))
                COUNTS["direct_character_product_sums"] += 1
        cases = 0
        for n in range(min(6, p) + 1):
            for C in combinations(range(p), n):
                F = [sum(chars[(x-c) % p] for c in C) for x in range(p)]
                T = {d: sum(cache[q] for q in combinations(C, d))
                     for d in (2, 4, 6)}
                T[8] = 0
                assert T[2] == -comb(n, 2)
                for r in range(1, 5):
                    mom = sum(v ** (2*r) for v in F)
                    hsum = factorial(2*r) * T[2*r]
                    # Direct product aggregation versus row recurrence.
                    recsum = sum(hermite_discrete(F[x], n - int(x in C), 2*r)
                                 for x in range(p))
                    assert hsum == recsum
                    COUNTS["actual_aggregate_identities"] += 1
                    scale = p * n**r
                    assert hsum >= -(8*r)**r * scale
                    assert abs(hsum) <= mom + (8*r)**r * scale
                    assert mom <= 2**r * hsum + 2*(16*r)**r * scale
                    COUNTS["actual_aggregate_inequalities"] += 3
                a0 = 15*n**3 - 30*n**2 + 16*n
                a2 = 90*n*n - 300*n + 272
                a4 = 360*n - 960
                ledger = (p*a0 - a2*comb(n, 2) + a4*T[4] + 720*T[6]
                          - 15*sum(F[c]**4 for c in C)
                          - 15*sum(F[c]**2 for c in C) - n)
                assert ledger == sum(v**6 for v in F)
                COUNTS["actual_sixth_moment_ledgers"] += 1
                cases += 1
        reports.append({"p": p, "sets": cases, "max_set_size": min(6, p)})
    return reports


def prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def next_prime(n):
    while not prime(n):
        n += 1
    return n


def model(p):
    k = isqrt(isqrt(p))
    P = p - k
    L = isqrt(p)
    W = P - L
    e = (k - 1) % 2
    v = Q(P*k - k*e*e - L*k*k, W)
    assert W > 0 and k-2 <= v < k
    t = k % 2
    while (t+2)**2 <= v:
        t += 2
    b = t + 2
    lam = (v - t*t) / (b*b - t*t)
    assert 0 <= lam <= 1 and b <= k
    assert t*t <= v < b*b
    return dict(p=p, k=k, P=P, L=L, W=W, e=e, v=v, t=t, b=b, lam=lam)


def bulk(m, power):
    return (1-m["lam"])*m["t"]**power + m["lam"]*m["b"]**power


def model_moment(m, power):
    return m["L"]*m["k"]**power + m["W"]*bulk(m, power) + m["k"]*m["e"]**power


def correlation(m, d):
    k, L, W, e, t, b, lam = (m[j] for j in ("k", "L", "W", "e", "t", "b", "lam"))
    assert 0 <= d <= k
    if d % 2:
        return Q(0)
    zero = (k-d)*kappa(e, k-1, d) if d < k else Q(0)
    return L + W*((1-lam)*kappa(t, k, d) + lam*kappa(b, k, d)) + zero


def explicit_rows(m):
    """Enumerate actual weighted vectors for small fixtures; merge duplicates."""
    k = m["k"]
    rows = defaultdict(Q)

    def add_class(zero, F, mass):
        free = [i for i in range(k) if i != zero]
        N = len(free)
        a = (N+F)//2
        assert 0 <= a <= N and (N+F) % 2 == 0
        weight = Q(mass, comb(N, a))
        for plus in combinations(free, a):
            eps = [-1]*k
            if zero is not None:
                eps[zero] = 0
            for i in plus:
                eps[i] = 1
            rows[tuple(eps)] += weight

    for sign in (-1, 1):
        add_class(None, sign*k, Q(m["L"], 2))
        add_class(None, sign*m["t"], m["W"]*(1-m["lam"])/2)
        add_class(None, sign*m["b"], m["W"]*m["lam"]/2)
        for zero in range(k):
            add_class(zero, sign*m["e"], Q(1, 2))
    return {eps: w for eps, w in rows.items() if w}


def stringify_model(m):
    return {key: str(val) if isinstance(val, Q) else val for key, val in m.items()}


def check_explicit_models():
    reports = []
    for target in range(4, 9):
        p = next_prime(target**4)
        m = model(p)
        k = m["k"]
        assert k == target
        rows = explicit_rows(m)
        assert all(w > 0 for w in rows.values())
        assert sum(rows.values()) == p
        for eps, w in rows.items():
            assert rows[tuple(-v for v in eps)] == w
        for i in range(k):
            assert sum(w for eps, w in rows.items() if eps[i] == 0) == 1
            assert sum(w*eps[i] for eps, w in rows.items()) == 0
            COUNTS["explicit_column_zero_and_mean_checks"] += 2
            for j in range(k):
                assert sum(w*eps[i]*eps[j] for eps, w in rows.items()) == p*int(i == j)-1
                COUNTS["explicit_gram_entries"] += 1
        for d in range(min(k, 6) + 1):
            expected = correlation(m, d)
            for subset in combinations(range(k), d):
                total = sum(w*prod(eps[i] for i in subset) for eps, w in rows.items())
                assert total == expected
                COUNTS["explicit_distinct_subset_correlations"] += 1
        for power in (2, 4, 6):
            moment = sum(w*sum(eps)**power for eps, w in rows.items())
            assert moment == model_moment(m, power)
            COUNTS["explicit_weighted_moments"] += 1
        serialized = [[list(eps), str(w)] for eps, w in sorted(rows.items())]
        reports.append({
            "parameters": stringify_model(m), "distinct_rows": len(rows),
            "row_measure_sha256": sha256(json.dumps(serialized).encode()).hexdigest(),
            "correlations": {str(d): str(correlation(m, d)) for d in range(min(k, 6)+1)},
            "moments": {str(d): str(model_moment(m, d)) for d in (2, 4, 6)},
            "scope": "Explicit rational row enumeration; asymptotic k>=256 bounds are not asserted here."
        })
        COUNTS["explicit_weighted_rows"] += len(rows)
    return reports


def check_large_models():
    reports = []
    for target in (256, 257, 384, 512, 1024):
        p = next_prime(target**4)
        assert prime(p)
        m = model(p)
        k, e = m["k"], m["e"]
        assert k == target and k**4 <= p < (k+1)**4
        assert k*k <= m["L"] <= k*k+2*k
        assert model_moment(m, 2) == p*k-k*k
        assert correlation(m, 0) == p
        assert correlation(m, 2) == -1
        for d in (1, 3, 5):
            assert correlation(m, d) == 0
        for d in range(1, 7):
            assert correlation(m, d)**2 <= (d-1)**2*p
            COUNTS["large_individual_Weil_shaped_bounds"] += 1
        var = bulk(m, 4)-m["v"]**2
        assert var == (m["v"]-m["t"]**2)*(m["b"]**2-m["v"])
        assert 0 <= var <= 8*k
        assert (k-2)**2 <= bulk(m, 4) <= k*k+8*k
        assert bulk(m, 6) <= 4*k**3
        h4 = (1-m["lam"])*hermite_discrete(m["t"], k, 4) + m["lam"]*hermite_discrete(m["b"], k, 4)
        h6 = (1-m["lam"])*hermite_discrete(m["t"], k, 6) + m["lam"]*hermite_discrete(m["b"], k, 6)
        assert -2*k*k+2*k <= h4 <= -2*k*k+18*k-12 < 0
        assert abs(h6) <= 100*k**3
        assert abs(hermite_discrete(e, k-1, 6)) <= 100*k**3
        assert correlation(m, 4)**2 <= 4*p
        assert correlation(m, 6)**2 <= 25*p
        M4, M6 = model_moment(m, 4), model_moment(m, 6)
        assert M4 <= 3*p*k*k
        assert M6 >= m["L"]*k**6 >= k**8
        assert 2*M6 >= p*k**4
        # An explicit Sidon label fixture, checked independently by pair sums.
        q = next_prime(k)
        assert prime(q) and q % 2 == 1
        labels = [t+2*q*((t*t) % q) for t in range(k)]
        assert len(set(labels)) == k and 2*max(labels) < p
        pair_sums = set()
        for i in range(k):
            for j in range(i, k):
                v = (labels[i]+labels[j]) % p
                assert v not in pair_sums
                pair_sums.add(v)
                COUNTS["large_Sidon_unordered_pair_sums"] += 1
        assert len(pair_sums) == comb(k+1, 2)
        E2 = k + 4*comb(k, 2)
        assert E2 == 2*k*k-k
        reports.append({
            "parameters": stringify_model(m),
            "M4_over_p_k2": str(M4/(p*k*k)),
            "M6_over_p_k3": str(M6/(p*k**3)),
            "M6_ratio_lower_bound": str(Q(k, 2)),
            "correlations": {str(d): str(correlation(m, d)) for d in range(7)},
            "bulk_variance": str(var), "bulk_H4": str(h4), "bulk_H6": str(h6),
            "Sidon_labels": {"q": q, "labels": labels, "unordered_pair_sums": len(pair_sums), "E2": E2},
            "all_assertions_passed": True,
            "scope": "Compressed weighted model only. No actual large-field character kernel is evaluated."
        })
        COUNTS["large_weighted_fixtures"] += 1
    return reports


def check_constant_ledger():
    # Integer/rational checks of every numerical bound used in the proof
    # on a finite range. The prose proof, not this ledger, covers all k.
    for k in range(256, 4097):
        assert k**4-k**3-7*k > 0
        assert k**4-k*k-3*k > 0
        assert k*k-6*k+1 >= 0      # (sqrt(k)+1)^2 <= 2k
        assert k*k-8*k+4 >= 0      # (sqrt(k)+2)^2 <= 3k
        assert k*k-38*k-4 >= 0
        assert (k+1)**4 <= 2*k**4
        assert 79*k**3+120*k*k <= 100*k**3
        assert 1+15*k+45*k*k+15*k**3 <= 100*k**3
        assert 400*k <= 2*k*k and 200 < k*k
        for start, d in ((0, 4), (1, 4), (0, 6), (1, 6)):
            f = prod(k-j for j in range(start, start+d))
            assert f >= k**d * (1-Q(sum(range(start, start+d)), k))
            assert 2*f >= k**d
            COUNTS["factorial_product_bounds"] += 2
        COUNTS["constant_parameter_values"] += 1


def main():
    check_pointwise()
    print("Pointwise recurrence and inequalities passed.", flush=True)
    actual = check_actual_fields()
    print("Actual small-field identities passed.", flush=True)
    small = check_explicit_models()
    print("Explicit weighted sign measures passed.", flush=True)
    large = check_large_models()
    check_constant_ledger()
    inputs = [
        "research/parallel21-squarefree-moments-2026-09-05.md",
        "experiments/parallel21_squarefree_moments_2026_09_05.py",
        "research/parallel20-inversion-moments-2026-09-05.md",
        "research/parallel19-relation-free-reduction-2026-09-05.md",
        "research/parallel5-classical-2026-09-04.md",
        "sources/mcdonald-sahay-wyman-2210.03789v2.html",
    ]
    result = {
        "status": "All exact checks passed. Uniform character-moment and Paley targets remain unproved.",
        "arithmetic": "Python integers and Fraction; exact integer-square-root primality cutoffs; no floating acceptance.",
        "input_sha256": {s: sha256((ROOT/s).read_bytes()).hexdigest() for s in inputs},
        "counts": dict(COUNTS),
        "pointwise_scope": {"N": [0, 128], "r": [1, 8], "F": "All parity-compatible actual sign sums"},
        "actual_character_fields": actual,
        "explicit_weighted_models": small,
        "large_weighted_models": large,
        "constant_ledger_scope": "Finite k=256,...,4096; not a substitute for the all-k prose argument.",
        "limitations": [
            "Weighted models are not uniformly sampled prime-field character translates.",
            "Sidon labels impose no link between sign rows and field differences.",
            "No actual large-field character moment is calculated.",
            "Finite checks do not prove the missing one-sided aggregate upper bound.",
            "Independent mathematical review and formal verification remain outstanding.",
        ],
    }
    output = ROOT/"results/parallel21_squarefree_moments_2026_09_05.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "counts": result["counts"],
                      "output": str(output)}, indent=2), flush=True)


if __name__ == "__main__":
    main()

