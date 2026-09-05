#!/usr/bin/env python3
"""Exact finite checks of alphabet-subspace restrictions on MCA challenges."""
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import random

ROOT = Path(__file__).resolve().parents[1]


class Field:
    def __init__(self, p, modulus):
        self.p = p
        self.modulus = modulus
        self.d = len(modulus) - 1
        assert modulus[-1] == 1
        self.q = p ** self.d
        self.digits = [self.unpack(i) for i in range(self.q)]
        self.add = [[self.pack([(a+b) % p for a,b in zip(x,y)])
                     for y in self.digits] for x in self.digits]
        self.neg = [self.pack([(-a) % p for a in x]) for x in self.digits]
        self.mul = [[self.multiply(x,y) for y in self.digits] for x in self.digits]
        # Every nonzero element is a unit: a direct field check of this
        # finite commutative polynomial quotient, without an assumed modulus.
        self.inv = [None] + [next((b for b in range(1,self.q)
                                  if self.mul[a][b] == 1), None)
                             for a in range(1,self.q)]
        assert all(x is not None for x in self.inv[1:])

    def unpack(self, x):
        out = []
        for _ in range(self.d):
            out.append(x % self.p)
            x //= self.p
        return out

    def pack(self, coeffs):
        return sum(x*self.p**i for i,x in enumerate(coeffs))

    def multiply(self, x, y):
        c = [0]*(2*self.d-1)
        for i,a in enumerate(x):
            for j,b in enumerate(y):
                c[i+j] = (c[i+j]+a*b) % self.p
        for i in range(len(c)-1,self.d-1,-1):
            a = c[i]
            for j,b in enumerate(self.modulus):
                c[i-self.d+j] = (c[i-self.d+j]-a*b) % self.p
        return self.pack(c[:self.d])


def bad_by_constant_interpolation(F, f, g, required):
    """Literal MCA witnesses for dimension-one codes: enumerate all S.

    A constant codeword explains a restriction exactly iff its values
    are equal. No ratio-set theorem is used in deciding this event.
    """
    n = len(f)
    out = set()
    # A witness larger than required has a required-sized witness with
    # the same failure only when tested directly; hence retain all sizes.
    for size in range(required,n+1):
        for S in combinations(range(n),size):
            if len({f[i] for i in S}) == 1 and len({g[i] for i in S}) == 1:
                continue
            for gamma in range(F.q):
                folded = {F.add[f[i]][F.mul[gamma][g[i]]] for i in S}
                if len(folded) == 1:
                    out.add(gamma)
    return out


def exhaustive_descent():
    prime = Field(3,[0,1])
    extension = Field(3,[1,0,1])
    words = list(product(range(3),repeat=3))
    comparisons = 0
    nonempty = 0
    for f,g in product(words,repeat=2):
        for required in range(1,4):
            bad_p = bad_by_constant_interpolation(prime,f,g,required)
            bad_q = bad_by_constant_interpolation(extension,f,g,required)
            assert bad_p == bad_q
            comparisons += 1
            nonempty += bool(bad_p)
    return {'p':3,'q':9,'modulus':[1,0,1],'n':3,'k':1,
            'received_pairs_exhausted':len(words)**2,
            'agreement_thresholds':[1,2,3],
            'exact_bad_set_comparisons':comparisons,
            'nonempty_bad_set_comparisons':nonempty}


def span(F, basis):
    out = set()
    for c in product(range(F.p),repeat=len(basis)):
        x = 0
        for a,b in zip(c,basis):
            x = F.add[x][F.mul[a][b]]
        out.add(x)
    return out


def subspace_checks():
    F = Field(3,[2,1,0,0,1])
    U = span(F,[1,3])
    V = span(F,[1,9])
    assert len(U) == len(V) == 9
    ratios = {F.mul[u][F.inv[v]] for u in U-{0} for v in V-{0}}
    allowed = ratios | {0}
    assert len(ratios) <= (len(U)-1)*(len(V)-1)//(F.p-1)
    for gamma in range(1,F.q):
        intersection = U & {F.mul[gamma][v] for v in V}
        assert (gamma in ratios) == (len(intersection)>1)
    rng = random.Random(20260904)
    records = []
    count = 0
    maximum = 0
    # Sample received words, but exhaust all challenges, agreement subsets,
    # and thresholds for each pair. This is finite evidence, not a proof.
    for _ in range(96):
        f = tuple(rng.choice(sorted(U)) for _ in range(3))
        g = tuple(rng.choice(sorted(V)) for _ in range(3))
        for required in range(1,4):
            bad = bad_by_constant_interpolation(F,f,g,required)
            assert bad <= allowed
            count += 1
            if len(bad)>maximum:
                maximum = len(bad)
                records = [{'f':f,'g':g,'required':required,'bad':sorted(bad)}]
    return {'p':F.p,'q':F.q,'modulus':F.modulus,'n':3,'k':1,
            'domain':[0,1,2],
            'U_basis':[1,3],'V_basis':[1,9],
            'U':sorted(U),'V':sorted(V),'ratio_set':sorted(ratios),
            'ratio_set_size':len(ratios),'counting_upper_bound':32,
            'received_pairs_sampled':96,'threshold_checks':count,
            'maximum_bad_count_observed':maximum,'maximum_example':records}


def degree_one_descent():
    """Enumerate actual degree-at-most-one codewords in both fields."""
    P = Field(5,[0,1])
    E = Field(5,[3,0,1])  # t^2=2
    domain = (0,1,2,3)
    subsets = [S for size in range(1,5) for S in combinations(range(4),size)]

    def restrictions(F):
        codes = [tuple(F.add[a][F.mul[b][x]] for x in domain)
                 for a,b in product(range(F.q),repeat=2)]
        return {S:{tuple(c[i] for i in S) for c in codes} for S in subsets}

    small = restrictions(P)
    large = restrictions(E)
    rng = random.Random(9202604)
    comparisons = 0
    for _ in range(64):
        f = tuple(rng.randrange(5) for _ in domain)
        g = tuple(rng.randrange(5) for _ in domain)
        for required in range(1,5):
            sets = []
            for F,table in [(P,small),(E,large)]:
                bad = set()
                for S in subsets:
                    if len(S)<required:
                        continue
                    if (tuple(f[i] for i in S) in table[S]
                        and tuple(g[i] for i in S) in table[S]):
                        continue
                    for gamma in range(F.q):
                        folded = tuple(F.add[f[i]][F.mul[gamma][g[i]]] for i in S)
                        if folded in table[S]:
                            bad.add(gamma)
                sets.append(bad)
            assert sets[0] == sets[1]
            comparisons += 1
    return {'p':5,'q':25,'modulus':[3,0,1],'domain':domain,'k':2,
            'base_codewords_enumerated':25,'extension_codewords_enumerated':625,
            'received_pairs_sampled':64,'exact_bad_set_comparisons':comparisons}


def run():
    result = {'status':'Finite checks only; uniform proofs are in the note.',
              'exhaustive_base_field_descent':exhaustive_descent(),
              'degree_one_polynomial_descent':degree_one_descent(),
              'alphabet_subspaces':subspace_checks()}
    p = 2130706433
    assert p**5 > 2**128
    result['pinned_profile'] = {'p':p,'q':p**6,
        'base_field_pair_mca_upper_numerator':1,
        'base_field_pair_mca_upper_denominator':p**5,
        'denominator_exceeds_2_to_128':True}
    paths = ['research/parallel-prize-subfield-2026-09-04.md',
             'experiments/parallel_prize_subfield_2026_09_04.py',
             'research/prize-reduction-audit.md',
             'research/official-profile-and-trace.md']
    result['source_sha256'] = {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    target = ROOT/'results/parallel_prize_subfield_2026_09_04.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    run()
