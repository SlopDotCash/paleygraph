#!/usr/bin/env python3
"""Exact checks for the triple-mass energy recovery theorem.

Literal counts validate finite instances. They do not prove the remaining
uniform triple-mass estimate or the Paley conjecture.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def pair_counts(a, b, p):
    return Counter((x + y) % p for x in a for y in b)


def energy(counts):
    return sum(c*c for c in counts.values())


def radical_upper_verified(actual_energy, base_square, N, masses):
    """Prove the bound using a rational LOWER bound for its radical RHS."""
    if not masses:
        assert actual_energy <= base_square
        return
    Q = 10**12
    lower_base = isqrt(base_square*Q*Q)
    lower_factor = isqrt(2*N*Q*Q)
    lower_sum = sum(isqrt(y*Q*Q) for y in masses)
    numerator = lower_base*Q + 2*lower_factor*lower_sum
    assert actual_energy*Q**4 <= numerator*numerator


def tower_check(case):
    p, N, g = case['p'], case['n'], case['generator']
    assert (p-1) % N == 0
    assert pow(g, N, p) == 1 and pow(g, N//2, p) == p-1
    levels, energies, pair_checks = [], {2: 6}, 0
    s = 4
    while s <= N:
        h = pow(g, N//s, p)
        H = [pow(h, j, p) for j in range(s)]
        K, L = H[::2], H[1::2]
        assert len(set(H)) == s and set(K).isdisjoint(L)
        mixed = pair_counts(K, L, p)
        child = pair_counts(K, K, p)
        parent = pair_counts(H, H, p)
        pair_checks += len(K)*len(L) + len(K)**2 + len(H)**2
        E, Ec, B = energy(parent), energy(child), energy(mixed)
        assert Ec == energies[s//2]
        roots = [pow(1+pow(h, j, p), s, p) for j in range(1, s//2, 2)]
        c = Counter(roots)
        Y = sum(v*(v-2) for v in c.values() if v >= 3)
        Z = sum((v-2)**2 for v in c.values() if v >= 3)
        singles = sum(v == 1 for v in c.values())
        expected = Counter()
        for v in c.values():
            expected[v] += s
        assert Counter(mixed.values()) == expected
        assert sum(max(v-2, 0)**2 for v in mixed.values()) == s*Z
        assert 0 <= Z <= Y
        T = sum(v*mixed.get(x, 0) for x, v in child.items())
        k = s//2
        assert E == 2*Ec + 6*B + 8*T
        assert B == 2*k*k + s*(Y-singles) and B >= s*Y
        # Subtract the two-representation baseline before Cauchy.
        excess_T = max(T-2*k*k, 0)
        assert excess_T**2 <= Ec*s*Z <= Ec*s*Y
        recurrence_excess = max(E-2*Ec-7*s*s-6*s*Y, 0)
        assert recurrence_excess**2 <= 64*Ec*s*Y
        energies[s] = E
        levels.append(dict(s=s,energy=E,child_energy=Ec,balanced_energy=B,
                           T=T,Y=Y,Z=Z,singletons=singles,all_passed=True))
        s *= 2
    assert energies[N] == case['parent_energy']
    assert levels[-1]['Y'] == case['triple_mass_Y']
    assert levels[-1]['balanced_energy'] == case['balanced_energy']
    R = sum(x['Y'] for x in levels)
    ell = len(levels)
    assert 3*N + 6*N*R <= energies[N]
    assert energies[N] <= 28*N*N + 16*N*ell*R
    radical_upper_verified(energies[N],14*N*N-25*N,N,[x['Y'] for x in levels])
    cutoffs = []
    for M in energies:
        upper_levels = [x for x in levels if x['s'] > M]
        remainder_mass = sum(x['Y'] for x in upper_levels)
        propagated = (N//M)*energies[M]
        assert propagated + 6*N*remainder_mass <= energies[N]
        base_square = propagated + 14*N*(N-M)
        radical_upper_verified(energies[N],base_square,N,[x['Y'] for x in upper_levels])
        cutoffs.append(M)
    # Compare the new square-root mass to the old weighted criterion.
    weighted = sum((Fraction(3,4)**(N.bit_length()-x['s'].bit_length()))
                   *Fraction(x['Y'],x['s']) for x in levels)
    # Cauchy gives (sum sqrt(Y))^2 <= 3*N*weighted. Its integer
    # expansion has cross radicals; a rational UPPER bound checks it.
    if R:
        Q=10**12
        upper_sum=sum(isqrt(x['Y']*Q*Q)+bool(x['Y']) for x in levels)
        assert Fraction(upper_sum**2,Q*Q) <= 3*N*weighted
    return dict(p=p,N=N,generator=g,levels=levels,total_triple_mass=R,
                energy=energies[N],cutoffs_checked=cutoffs,
                old_weighted_criterion=str(weighted),literal_pair_checks=pair_checks,
                quartic_window=case['quartic_window'],all_passed=True)


def abstract_checks():
    # Equal masses can saturate the logarithmic Cauchy loss using
    # individually admissible multiplicity lists. No field realization.
    rows=[]
    for m in [2,4,8,16]:
        N, c = 2**(3*m), 2**m
        levels=[]
        s=4*c
        while s<=N:
            assert s//4>=c
            levels.append(s)
            s*=2
        Y=c*(c-2)
        R=len(levels)*Y
        square_root_sum_squared=len(levels)**2*Y
        assert square_root_sum_squared==len(levels)*R
        rows.append(dict(m=m,N=N,cluster_size=c,active_levels=len(levels),
                         R=R,square_root_sum_squared=square_root_sum_squared,
                         realized_in_a_field=False))
    alpha,beta=Fraction(7,3),Fraction(49,20)
    cutoff=(alpha-1)/(beta-1)
    log_cutoff=Fraction(1,5)/(beta-1)
    assert cutoff==Fraction(80,87) and log_cutoff==Fraction(4,29)
    assert 1+cutoff*(beta-1)==alpha
    assert -log_cutoff*(beta-1)+Fraction(1,5)==0
    assert (beta-1)/2-Fraction(2,3)==Fraction(7,120)
    assert 1+2*((beta-1)/2)==beta
    return dict(equal_mass_profiles=rows,known_energy_cutoff_power=str(cutoff),
                known_energy_cutoff_log_divisor=str(log_cutoff),
                general_cutoff_formula='(alpha-1)/(beta-1)',all_passed=True)


def main():
    path=ROOT/'results/parallel48_collision_eliminants_2026_09_06.json'
    prior=json.loads(path.read_text())
    towers=[tower_check(case) for order in prior['orders'] for case in order['fields']]
    out=dict(scope='Two-sided energy recovery from first-precision tower triple mass; the uniform target estimate remains unproved.',
             towers=towers,abstract_checks=abstract_checks(),tower_cases=len(towers),
             level_checks=sum(len(t['levels']) for t in towers),
             cutoff_checks=sum(len(t['cutoffs_checked']) for t in towers),
             literal_pair_checks=sum(t['literal_pair_checks'] for t in towers),
             prior_certificate_sha256=sha256(path.read_bytes()).hexdigest(),all_passed=True)
    target=ROOT/'results/parallel50_tower_mass_2026_09_06.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['tower_cases','level_checks','cutoff_checks','literal_pair_checks','all_passed']}))


if __name__=='__main__':
    main()
