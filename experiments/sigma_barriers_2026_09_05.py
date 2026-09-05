#!/usr/bin/env python3
"""Exact checks for research/sigma-barriers-2026-09-05.md.

Standard library only. Writes results/sigma_barriers_2026_09_05.json.
"""
import itertools
import json
import math
import os
import random
import sys
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def legendre_table(p):
    chi = [0] * p
    for x in range(1, p):
        chi[x * x % p] = 1
    for x in range(1, p):
        if chi[x] == 0:
            chi[x] = -1
    return chi


def check_b1(primes_small, primes_random, rng):
    """Shift orthogonality, second moment, Chung bound."""
    counts = {"shift": 0, "second_moment": 0, "chung": 0}
    for p in primes_small:
        chi = legendre_table(p)
        for c in range(1, p):
            s = sum(chi[x] * chi[(x + c) % p] for x in range(p))
            assert s == -1, (p, c, s)
            counts["shift"] += 1
        for mask in range(1, 1 << p):
            B = [i for i in range(p) if mask >> i & 1]
            tot = 0
            for x in range(p):
                f = sum(chi[(x + b) % p] for b in B)
                tot += f * f
            assert tot == len(B) * (p - len(B)), (p, B, tot)
            counts["second_moment"] += 1
    for p in primes_random:
        chi = legendre_table(p)
        for _ in range(20):
            nb = rng.randint(1, p - 1)
            na = rng.randint(1, p - 1)
            B = rng.sample(range(p), nb)
            A = rng.sample(range(p), na)
            F = [sum(chi[(x + b) % p] for b in B) for x in range(p)]
            assert sum(f * f for f in F) == nb * (p - nb)
            counts["second_moment"] += 1
            S = sum(F[a] for a in A)
            assert S * S <= na * nb * (p - nb), (p, A, B, S)
            counts["chung"] += 1
    return counts


def check_b2(primes, k_values, rng):
    """Exact tuple counts and inner sums for the Hoelder-Weil step."""
    out = []
    for p in primes:
        chi = legendre_table(p)
        for k in k_values:
            nb = min(p - 1, 6)
            B = rng.sample(range(p), nb)
            even_mult = 0
            weil_ok = 0
            max_nonsquare = 0
            for tup in itertools.product(B, repeat=2 * k):
                inner = 0
                for x in range(p):
                    v = 1
                    for b in tup:
                        v *= chi[(x + b) % p]
                        if v == 0:
                            break
                    inner += v
                cnt = Counter(tup)
                if all(m % 2 == 0 for m in cnt.values()):
                    even_mult += 1
                    # exactly p minus the number of distinct roots
                    assert inner == p - len(cnt), (p, tup, inner)
                else:
                    # Weil: |inner| <= (deg-1) sqrt(p) where deg = number of
                    # distinct roots with odd multiplicity (after removing
                    # square factors), plus zeros contribute at most len(cnt)
                    odd = sum(1 for m in cnt.values() if m % 2 == 1)
                    bound = (odd - 1) * math.sqrt(p) + (len(cnt) - odd)
                    assert abs(inner) <= bound + 1e-9, (p, tup, inner, bound)
                    weil_ok += 1
                    max_nonsquare = max(max_nonsquare, abs(inner))
            # the number of even-multiplicity ordered tuples is at most
            # (2k-1)!! nb^k; for k=2 it equals 3nb^2-2nb, exceeding k! nb^k
            dfact = math.factorial(2 * k) // (2 ** k * math.factorial(k))
            assert even_mult <= dfact * nb ** k, (p, k, nb, even_mult)
            if k == 2:
                assert even_mult == 3 * nb * nb - 2 * nb
            out.append({"p": p, "k": k, "|B|": nb, "even_multiplicity_tuples": even_mult,
                        "weil_checked": weil_ok, "max_abs_nonsquare_inner_sum": max_nonsquare,
                        "double_factorial_bound": dfact * nb ** k,
                        "k_factorial_Bk_is_violated": even_mult > math.factorial(k) * nb ** k})
    return out


def fp2_elements(p, nres):
    """F_{p^2} = F_p[t]/(t^2 - nres). Elements as (a, b) = a + b t."""
    return [(a, b) for a in range(p) for b in range(p)]


def fp2_mul(x, y, p, nres):
    a, b = x
    c, d = y
    return ((a * c + b * d * nres) % p, (a * d + b * c) % p)


def fp2_pow(x, e, p, nres):
    r = (1, 0)
    while e:
        if e & 1:
            r = fp2_mul(r, x, p, nres)
        x = fp2_mul(x, x, p, nres)
        e >>= 1
    return r


def check_b4(primes):
    """Quadratic character of F_{p^2} is +1 on F_p^*."""
    out = []
    for p in primes:
        chi = legendre_table(p)
        nres = next(x for x in range(2, p) if chi[x] == -1)
        q = p * p
        e = (q - 1) // 2
        allone = True
        for u in range(1, p):
            v = fp2_pow((u, 0), e, p, nres)
            if v != (1, 0):
                allone = False
        # also check that the character is nontrivial on F_{p^2}^*: some
        # element must map to -1 (t itself need not be a non-square, e.g. p=3)
        witness = None
        for z in fp2_elements(p, nres):
            if z == (0, 0):
                continue
            if fp2_pow(z, e, p, nres) == (p - 1, 0):
                witness = z
                break
        nontrivial = witness is not None
        assert allone and nontrivial, p
        out.append({"p": p, "q": q, "chi_q_on_Fp_star": "all +1",
                    "nonsquare_witness_a_plus_b_t": list(witness),
                    "two_set_sum_A=B=F_p": p * (p - 1)})
    return out


def check_b5(primes):
    """Interval case: with N = floor((n_p - 1)/2), S([1,N],[1,N]) = N^2."""
    out = []
    for p in primes:
        chi = legendre_table(p)
        n_p = next(x for x in range(1, p) if chi[x] == -1)
        N = (n_p - 1) // 2
        if N == 0:
            out.append({"p": p, "n_p": n_p, "N": 0, "S": 0, "N2": 0})
            continue
        S = sum(chi[(a + b) % p] for a in range(1, N + 1) for b in range(1, N + 1))
        assert S == N * N, (p, n_p, N, S)
        out.append({"p": p, "n_p": n_p, "N": N, "S": S, "N2": N * N})
    return out


def main():
    rng = random.Random(20260905)
    primes = primes_upto(500)
    small = [3, 5, 7, 11, 13]
    random_primes = [q for q in primes if 17 <= q <= 199]
    result = {
        "b1": check_b1(small, random_primes, rng),
        "b2": check_b2([7, 11, 13], [2, 3], rng),
        "b4": check_b4([q for q in primes if 3 <= q <= 60]),
        "b5": check_b5([q for q in primes if q >= 3]),
        "scope": "finite exact checks of the elementary facts in the barrier note; no asymptotic claim",
    }
    out = os.path.join(ROOT, "results", "sigma_barriers_2026_09_05.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=1)
    print("b1", result["b1"])
    print("b2 cases", len(result["b2"]), "b4 primes", len(result["b4"]), "b5 primes", len(result["b5"]))
    print("wrote", out)


if __name__ == "__main__":
    sys.exit(main())
