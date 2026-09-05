#!/usr/bin/env python3
"""Exact sparse count of the distinct, opposite-free, fully unbalanced remainder."""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
import json
from math import comb, isqrt
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
CHECKS = Counter()

def check(ok, name):
    assert ok, name
    CHECKS[name] += 1

def prime(p):
    if p < 2 or p % 2 == 0: return p == 2
    return all(p % d for d in range(3, isqrt(p) + 1, 2))

def subgroup(p, n):
    check(prime(p) and n >= 4 and n & (n-1) == 0 and (p-1) % n == 0,
          'exact_prime_and_dyadic_divisor')
    for a in range(2, p):
        g = pow(a, (p-1)//n, p)
        if pow(g, n//2, p) != 1:
            H = sorted(pow(g, j, p) for j in range(n))
            check(len(set(H)) == n and pow(g, n, p) == 1, 'exact_subgroup_order')
            return H
    raise AssertionError('generator missing')

def free(word, p):
    return not any(-x % p in word for x in word)

def balance_count(word, p):
    count = 0
    for pair in combinations(range(1, 6), 2):
        side = {0, *pair}
        left = right = 1
        for i, a in enumerate(word):
            if i in side: left = left*a % p
            else: right = right*a % p
        count += (left + right) % p == 0
    return count

def sparse(p, n):
    started = time.monotonic()
    H = subgroup(p, n)
    Hset = set(H)
    pairs = list(combinations_with_replacement(H, 2))
    sums = defaultdict(list)
    shifted = defaultdict(list)
    r2 = Counter()
    for a, b in pairs:
        weight = 1 if a == b else 2
        sums[(a+b) % p].append((a,b,weight))
        shifted[(a*b-a-b) % p].append((a,b,weight))
        r2[(a+b) % p] += weight
    E2 = sum(c*c for c in r2.values())
    r4 = sum(c*r2.get((-2-s) % p, 0) for s,c in r2.items())
    w = {z: sum(row[2] for row in bucket) for z,bucket in shifted.items()}
    orbit_mass = Counter()
    # The kernel of z -> z^n is H. Raising to the subgroup index would
    # instead quotient by the complementary subgroup and gives wrong energies.
    exponent = n
    check(all(pow(h, exponent, p) == 1 for h in H), 'coset_label_kernel_contains_H')
    for z,c in w.items():
        if z: orbit_mass[pow(z, exponent, p)] += c
    E3 = n*n*w.get(0,0)**2 + n*sum(c*c for c in orbit_mass.values())
    check(sum(w.values()) == n*n, 'shifted_pair_mass')
    check(sum(orbit_mass.values()) + w.get(0,0) == n*n, 'orbit_partition_mass')
    check(w.get(0,0) == r2.get(1,0), 'zero_orbit_count')
    T4 = 3*n*n-3*n
    T6 = 15*n**3-45*n*n+40*n
    J6 = (15*n-60)*(E2-T4) + 60*n*(r2.get(2,0)-1) - 30*n*(3 in Hset)
    check(J6 >= 0, 'nonnegative_opposite_pair_remainder')

    # Mark a repeated positional pair, normalize its value to one, and divide
    # each word by its number of marked choices. Pair multiplicities are exact.
    repeated = Fraction(0)
    marked_all = 0
    repeated_patterns = Counter()
    for s, left in sums.items():
        for a,b,ab in left:
            for c,d,cd in sums.get((-2-s) % p, ()):
                word = (1,1,a,b,c,d)
                mass = ab*cd
                marked_all += mass
                if not free(word,p): continue
                multiplicities = Counter(word)
                nu = sum(comb(v,2) for v in multiplicities.values())
                weighted = Fraction(15*n*mass, nu)
                repeated += weighted
                repeated_patterns[tuple(sorted(multiplicities.values(),reverse=True))] += weighted
    check(marked_all == r4, 'independent_four_sum_convolution')
    check(repeated.denominator == 1, 'unmarked_repeated_count_integral')
    repeated = repeated.numerator
    check(repeated <= 15*n*(r4-6*n+8), 'previous_repeated_upper')
    check(all(v.denominator == 1 for v in repeated_patterns.values()), 'pattern_counts_integral')

    # One fixed balanced split: (ad,b,c,-bc,-a,-d). Summing all split
    # indicators counts every distinct opposite-free balanced word exactly once.
    fixed_balanced_normalized = 0
    balanced_samples = []
    for bucket in shifted.values():
        for a,d,ad_weight in bucket:
            for b,c,bc_weight in bucket:
                word = (a*d % p,b,c,-b*c % p,-a % p,-d % p)
                if len(set(word)) != 6 or not free(word,p): continue
                check(sum(word) % p == 0 and balance_count(word,p) == 1,
                      'actual_distinct_balanced_unique_split')
                fixed_balanced_normalized += ad_weight*bc_weight
                if len(balanced_samples) < 3: balanced_samples.append(word)
    balanced = 10*n*fixed_balanced_normalized
    shifted_energy = sum(c*c for c in w.values())
    X = shifted_energy - (6*n*n-9*n+4)
    check(0 <= balanced <= 10*n*X, 'previous_balanced_upper')
    D6 = E3-T6-J6-repeated-balanced
    check(D6 >= 0, 'nonnegative_distinct_unbalanced_remainder')
    check(D6 % (720*n) == 0, 'free_scaling_orbits_of_distinct_six_sets')
    check(balanced % (720*n) == 0, 'free_scaling_orbits_of_balanced_six_sets')
    return {'p':p,'n':n,'E2':E2,'E3':E3,'r4_at_2':r4,'X':X,'T6':T6,'J6':J6,
        'repeated_R6':repeated,'distinct_balanced_R6':balanced,
        'distinct_unbalanced_R6':D6, 'D6_scaling_orbits':D6//(720*n),
        'D6_over_n_cubed':str(Fraction(D6,n**3)),
        'repeated_patterns':{str(k):int(v) for k,v in sorted(repeated_patterns.items())},
        'fixed_balanced_normalized':fixed_balanced_normalized,
        'balanced_samples':balanced_samples,'seconds':time.monotonic()-started}

def direct(p,n):
    # Independent enumeration by unordered six-multisets, with factorial weights.
    H = subgroup(p,n)
    counts = Counter()
    facts = [1,1,2,6,24,120,720]
    for word in combinations_with_replacement(H,6):
        if sum(word) % p: continue
        multiplicities = Counter(word)
        denominator = 1
        for v in multiplicities.values(): denominator *= facts[v]
        weight = 720//denominator
        counts['E3'] += weight
        if all(multiplicities[x] == multiplicities[-x % p] for x in multiplicities):
            counts['T6'] += weight
        elif not free(word,p): counts['J6'] += weight
        elif len(multiplicities) < 6: counts['repeated_R6'] += weight
        elif balance_count(word,p): counts['distinct_balanced_R6'] += weight
        else: counts['distinct_unbalanced_R6'] += weight
    return counts

def main():
    # The unique-split statement does not need zero sum or subgroup membership.
    for p in [13,17,19,23]:
        for word in combinations(range(1,p),6):
            if free(word,p):
                check(balance_count(word,p) <= 1, 'all_distinct_opposite_free_unique_split')
    rows=[]
    for p,n in [(5,4),(17,8),(97,16),(73,4),(1153,8),(33713,16),(37201,16)]:
        row=sparse(p,n)
        expected=direct(p,n)
        for k in ['E3','T6','J6','repeated_R6','distinct_balanced_R6','distinct_unbalanced_R6']:
            check(row[k] == expected[k], 'independent_six_multiset_enumeration')
        rows.append(row)
    prior=json.loads((ROOT/'results/parallel25_subgroup_growing_orders_2026_09_05.json').read_text())
    for old in prior['cases']:
        row=sparse(old['p'],old['n'])
        for k in ['E3','repeated_R6','distinct_balanced_R6','distinct_unbalanced_R6']:
            check(row[k] == old['counts'].get(k,0), 'pass25_ordered_word_enumeration_agrees')
        check(row['E2']==old['E2'] and row['r4_at_2']==old['r4_at_2'], 'pass25_lower_counts_agree')
        rows.append(row)
        print(json.dumps({k:row[k] for k in ['p','n','E3','D6_scaling_orbits','seconds']}),flush=True)
    for n in [512,1024]:
        for anchor in [n**4//4,n**4//2,3*n**4//4]:
            p=anchor+(1-anchor)%n
            while not prime(p): p+=n
            row=sparse(p,n)
            rows.append(row)
            print(json.dumps({k:row[k] for k in ['p','n','E3','D6_scaling_orbits','seconds']}),flush=True)
    inputs=['experiments/parallel27_six_distinct_2026_09_05.py',
        'research/parallel23-subgroup-upper-2026-09-05.md',
        'research/parallel24-subgroup-unbalanced-2026-09-05.md',
        'research/parallel25-subgroup-growing-orders-2026-09-05.md',
        'results/parallel25_subgroup_growing_orders_2026_09_05.json']
    output={'status':'passed','checked_at_utc':datetime.now(timezone.utc).isoformat(),
        'checks':dict(CHECKS),'cases':rows,
        'input_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'scope':'Exact finite D6 counts via a proved sparse formula; no uniform asymptotic bound.'}
    (ROOT/'results/parallel27_six_distinct_2026_09_05.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(rows),'checks':dict(CHECKS)}),flush=True)

if __name__=='__main__': main()
