#!/usr/bin/env python3
"""Exact opposite-free sixth-energy decomposition and symmetric-set obstruction."""
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from math import isqrt, factorial
from pathlib import Path
import hashlib
import json

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'results/parallel21_subgroup_next_input_2026_09_05.json'
COUNTS = Counter()


def check(ok, name):
    assert ok, name
    COUNTS[name] += 1


def prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def factor(n):
    ans = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            ans.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        ans.append(n)
    return ans


def subgroup(p, n):
    check(prime(p) and (p-1) % n == 0, 'prime_subgroup_parameters')
    fs = factor(p-1)
    g = next(g for g in range(2, p) if all(pow(g, (p-1)//d, p) != 1 for d in fs))
    z = pow(g, (p-1)//n, p)
    S = sorted(pow(z, j, p) for j in range(n))
    check(len(set(S)) == n and pow(z, n, p) == 1, 'subgroup_order')
    return S


def tuple_weight(t):
    val = factorial(len(t))
    for m in Counter(t).values():
        val //= factorial(m)
    return val


def energies(S, p):
    r2 = Counter((a+b) % p for a in S for b in S)
    r3 = Counter()
    for x, w in r2.items():
        for a in S:
            r3[(x+a) % p] += w
    return r2, sum(v*v for v in r2.values()), sum(v*v for v in r3.values())


def primitive_six(S, p):
    """Independent enumeration: unordered triples, exact ordered weights, no opposite pair."""
    ix = {a: i for i, a in enumerate(S)}
    by_sum = defaultdict(list)
    for t in combinations_with_replacement(S, 3):
        mask = sum(1 << ix[a] for a in set(t))
        negmask = sum(1 << ix[-a % p] for a in set(t))
        if mask & negmask:
            continue
        by_sum[sum(t) % p].append((mask, negmask, tuple_weight(t)))
    ans = 0
    for s, left in by_sum.items():
        ns = -s % p
        if s > ns:
            continue
        fac = 1 if s == ns else 2
        for mask, negmask, w in left:
            for mask2, negmask2, w2 in by_sum.get(ns, []):
                if not (mask & negmask2):
                    ans += fac*w*w2
                    COUNTS['admissible_triple_pair_terms'] += 1
    return ans


def classify_four(S, p):
    counts = Counter()
    ss = set(S)
    for t in combinations_with_replacement(S, 4):
        if sum(t) % p:
            continue
        if any((-a % p) in t for a in t):
            continue
        key = ''.join(str(x) for x in sorted(Counter(t).values(), reverse=True))
        counts[key] += tuple_weight(t)
    check(set(counts) <= {'1111', '211', '31'}, 'four_multiplicity_classification')
    return counts


def inspect(S, p, name, is_subgroup=False, classify=False):
    n = len(S)
    check(len(set(S)) == n and 0 not in S and {-a % p for a in S} == set(S), 'symmetric_nonzero_set')
    r2, E2, E3 = energies(S, p)
    T4 = 3*n*n - 3*n
    T6 = 15*n**3 - 45*n*n + 40*n
    R6 = primitive_six(S, p)
    Sset = set(S)
    J = sum(3*a % p in Sset for a in S)
    U = sum(r2[2*a % p]-1 for a in S)
    rhs = T6 + (15*n-60)*(E2-T4) + 60*U - 30*J + R6
    check(E3 == rhs, 'general_exact_sixth_decomposition')
    check(T6 + R6 <= E3 <= T6 + 15*n*(E2-T4) + R6, 'opposite_pair_upper_lower_bounds')
    if is_subgroup:
        K = r2[2 % p]
        eps = int(3 in Sset)
        check(U == n*(K-1) and J == n*eps, 'subgroup_compression')
        check(E3 == T6+(15*n-60)*(E2-T4)+60*n*(K-1)-30*n*eps+R6, 'subgroup_exact_sixth_decomposition')
    else:
        K = eps = None
    if classify:
        cc = classify_four(S, p)
        check(sum(cc.values()) == E2-T4, 'primitive_four_total')
        check(cc['31'] == 4*J and cc['211'] == 6*U-12*J, 'primitive_four_repeated_values')
        check(2*(E3-T6-R6) == (30*n-120)*cc['1111']+(30*n-100)*cc['211']+(30*n-75)*cc['31'], 'weighted_extension_identity')
    return {'name': name, 'p': p, 'n': n, 'E2': E2, 'T4': T4, 'E3': E3, 'T6': T6,
            'opposite_free_six': R6, 'opposite_containing_excess': E3-T6-R6,
            'pair_representations_at_2': K, 'three_in_subgroup': eps,
            'sum_pair_representations_correction': U, 'three_overlap': J}


def symmetric_sidon(n):
    k = n//2
    q = next(q for q in range(k, 2*k) if prime(q))
    low = n**4//4
    p = low + ((1-low) % n)
    while not prime(p):
        p += n
    check(low <= p <= n**4 and (p-1) % n == 0, 'quartic_dyadic_prime_parameters')
    B = [8*q*q + t + 2*q*((t*t) % q) for t in range(k)]
    A = sorted(B + [p-b for b in B])
    pair_sums = [a+b for i, a in enumerate(B) for b in B[i:]]
    check(len(set(pair_sums)) == k*(k+1)//2, 'integer_sidon_pair_sums')
    check(6*max(B) < p and min(B)*2 > max(B), 'integer_no_wrap_and_sign_separation')
    r2B = Counter(a+b for a in B for b in B)
    r3B = Counter()
    for s, w in r2B.items():
        for a in B:
            r3B[s+a] += w
    E3B = sum(w*w for w in r3B.values())
    _, E2A, E3A = energies(A, p)
    check(E2A == 3*n*n-3*n, 'symmetric_sidon_intrinsic_fourth_energy')
    check(E3A == 20*E3B, 'symmetric_sidon_sixth_energy_identity')
    check(20*E3A >= n**4, 'symmetric_sidon_quartic_sixth_lower_bound')
    check(1 not in A, 'constructed_set_is_not_a_subgroup')
    out = {'n': n, 'k': k, 'q': q, 'p': p, 'positive_set': B, 'E2': E2A, 'E3': E3A,
           'E3_positive_set': E3B, 'T6': 15*n**3-45*n*n+40*n,
           'opposite_free_six_by_decomposition': E3A-(15*n**3-45*n*n+40*n),
           'E3_over_n_cubed': [E3A, n**3]}
    if n <= 32:
        direct = primitive_six(A, p)
        check(direct == out['opposite_free_six_by_decomposition'], 'symmetric_sidon_primitive_six_independent')
    return out


def main():
    sub = []
    for p in (5, 7, 13, 17, 29, 41, 73, 97, 193, 257):
        for n in (2, 4, 8, 16):
            if (p-1) % n == 0:
                sub.append(inspect(subgroup(p,n), p, f'H{n}@{p}', True, n <= 8))
    for p in (6700417, 6878593, 7041409, 7177601, 7204033, 7884353, 7987009,
              8019073, 9190913, 9877633, 10219457, 11127041, 12942337, 13640513, 14721281):
        sub.append(inspect(subgroup(p,64),p,f'H64@{p}',True,False))
    generic = []
    for p in (13, 17, 29):
        for pos in ((1,2,3), (1,2,4), (1,3,5)):
            S = sorted(set(pos) | {-x % p for x in pos})
            generic.append(inspect(S,p,f'symmetric{pos}@{p}',False,True))
    obstruction = [symmetric_sidon(n) for n in (16,32,64,128,256)]
    hashes = {}
    for rel in ('research/subgroup-target.md','research/parallel7-subgroup-2026-09-04.md',
                'research/sigma-subgroup-2026-09-05.md','research/coset-coherence.md',
                'research/cyclotomic-prime-average.md', 'experiments/parallel21_subgroup_next_input_2026_09_05.py'):
        hashes[rel] = hashlib.sha256((BASE/rel).read_bytes()).hexdigest()
    result = {'status': 'PASS', 'scope': 'Exact finite validation; asymptotic bounds proved separately in the note. No subgroup cancellation exponent is improved.',
              'checks': dict(COUNTS), 'subgroup_fixtures': sub, 'general_symmetric_fixtures': generic,
              'symmetric_sidon_obstruction': obstruction, 'input_sha256': hashes}
    OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':dict(COUNTS),'subgroup_fixtures':len(sub),'general_fixtures':len(generic),
                      'symmetric_sidon_obstruction':[{k:v for k,v in x.items() if k not in ('positive_set',)} for x in obstruction]},indent=2))

if __name__ == '__main__':
    main()
