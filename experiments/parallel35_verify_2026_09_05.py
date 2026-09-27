#!/usr/bin/env python3
"""Exact checks for the single-degree comparison; no full conjecture claim."""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
from time import perf_counter
import json

ROOT = Path(__file__).resolve().parents[1]
CHECKS = Counter()


def check(ok, kind, data=None):
    if not ok:
        raise AssertionError((kind, data))
    CHECKS[kind] += 1


def rat(q):
    q = F(q)
    return dict(numerator=q.numerator, denominator=q.denominator)


def elementary(y, limit):
    e = [F(1)]+[F(0)]*limit
    for i, a in enumerate(y):
        for j in range(min(i+1, limit), 0, -1):
            e[j] += a*e[j-1]
    return e


def trim(p):
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p


def derivative(p):
    return trim([i*p[i] for i in range(1, len(p))] or [F(0)])


def remainder(a, b):
    a = a.copy()
    while len(a) >= len(b) and any(a):
        k, c = len(a)-len(b), a[-1]/b[-1]
        for j, v in enumerate(b):
            a[j+k] -= c*v
        trim(a)
    return a


def sturm(p):
    chain = [trim(p.copy()), derivative(p)]
    while any(chain[-1]):
        rem = [-v for v in remainder(chain[-2], chain[-1])]
        if not any(rem):
            break
        # Positive rescaling keeps signs and limits coefficient growth.
        scale = abs(rem[-1])
        chain.append([v/scale for v in rem])
    return chain


def sign(x):
    return (x > 0)-(x < 0)


def sign_at_radical(p, D, side):
    even = sum(v*D**(k//2) for k, v in enumerate(p) if k % 2 == 0)
    odd = side*sum(v*D**(k//2) for k, v in enumerate(p) if k % 2)
    if not even or not odd:
        return sign(even) if not odd else sign(odd)
    if sign(even) == sign(odd):
        return sign(even)
    return sign(even)*sign(even*even-odd*odd*D)


def changes(signs):
    nonzero = [x for x in signs if x]
    return sum(a != b for a, b in zip(nonzero, nonzero[1:]))


def eval_poly(p, x):
    answer = F(0)
    for v in reversed(p):
        answer = answer*x+v
    return answer


def row_case(y, s, label):
    N, r = len(y), 2*s
    n = 2*N
    total = sum(y)
    e = elementary(y, r)
    H = factorial(r)*e[r]
    main = total**r
    K = (256*s*n)**s
    check(H >= -K, 'pointwise lower')
    check(abs(H) <= 4**s*(main+(64*s*n)**s), 'pointwise absolute')
    check(main <= 16**s*H+2*(4096*s*n)**s, 'pointwise reverse')
    item = dict(N=N, s=s, label=label, sum=rat(total), H=rat(H), moment=rat(main),
                values=[dict(value=rat(v), count=m) for v, m in sorted(Counter(y).items())])
    if r <= N//2:
        mu = total/N
        z = [a-mu for a in y]
        ez = elementary(z, r)
        # h(X)=r! e_r(z_i+X/N), derived from the shifted product.
        h = [F(factorial(r)*comb(N-j, r-j), N**(r-j))*ez[j]
             for j in reversed(range(r+1))]
        a = F(prod(N-j for j in range(r)), N**r)
        check(h[-1] == a, 'normalization')
        check(F(1, 2**r) <= a <= 1, 'normalization bounds')
        check(eval_poly(h, total) == H, 'shifted coefficient identity')
        chain = sturm(h)
        D = 64*r*N
        variation_left = changes([sign_at_radical(q, D, -1) for q in chain])
        variation_right = changes([sign_at_radical(q, D, 1) for q in chain])
        distinct_degree = len(h)-len(chain[-1])
        check(variation_left-variation_right == distinct_degree,
              'exact Sturm root enclosure', (N, s, label))
        large = main >= (4*D)**s
        if large:
            check(H >= a*main/2**r, 'large-sum product lower')
        item.update(root_radius_squared=D, root_count=distinct_degree,
                    variation_left=variation_left, variation_right=variation_right,
                    normalization=rat(a), high_sum_branch=large,
                    shifted_polynomial=[rat(v) for v in h])
    else:
        check(n < 8*s and abs(H) <= n**r and main <= n**r, 'small-size branch')
    return item


def subgroup_counts():
    # Recheck the conclusion on prior exact counts, including dense sets
    # where pass33's epsilon condition is not small.
    records = []
    old = json.loads((ROOT/'results/parallel34_verification_2026_09_05.json').read_text())
    for case in old['cases']:
        p, n = case['p'], case['n']
        for row in case['rows'][1:]:
            s = row['s']
            B = F(row['centered_opposite_free']['numerator'], row['centered_opposite_free']['denominator'])
            T = F(row['energy'])-F(n**(2*s), p)
            check(B >= -(256*s*n)**s, 'finite-field lower')
            check(abs(B) <= 4**s*(T+(64*s*n)**s), 'finite-field absolute')
            check(T <= 16**s*B+2*(4096*s*n)**s, 'finite-field reverse')
            records.append(dict(p=p, n=n, s=s, B=rat(B), T=rat(T)))
    return records


def averaged_root_counterexample():
    p, H = 1153, [1, 75, 123, 140, 1013, 1030, 1078, 1152]
    check(all(p % d for d in range(2, 34)), 'witness primality')
    check(all(a*b % p in H for a in H for b in H), 'witness subgroup')
    check(8**4 <= 4*p <= 4*8**4, 'witness quartic window')
    reps = [1, 75, 123, 140]
    counts = [0]*5
    for signs in product([-1, 0, 1], repeat=4):
        if sum(a*b for a, b in zip(reps, signs)) % p == 0:
            counts[sum(v != 0 for v in signs)] += 1
    check(counts == [1, 0, 0, 0, 0], 'witness full signed enumeration')
    P = [p*counts[k]-comb(4, k)*2**k for k in range(5)]
    check(P == [1152, -8, -24, -32, -16], 'witness centered polynomial')
    # Translation t=1+2w gives 1153-t^4. Its two imaginary roots are
    # certified by this exact polynomial form, without numerical roots.
    translated = [sum(F(P[j], 2**j)*comb(j, k)*(-1)**(j-k)
                      for j in range(k, 5)) for k in range(5)]
    check(translated == [1153, 0, 0, 0, -1], 'witness translation')
    return dict(p=p, n=8, elements=H, signed_selections=81,
                unordered_zero_counts=counts, p_times_centered_coefficients=P,
                translated_coefficients=[int(v) for v in translated],
                real_roots=2, nonreal_roots=2,
                scope='The averaged generating polynomial need not be real-rooted; not a Paley counterexample.')


def newton_group_checks():
    def multiply(a, b, p):
        result = Counter()
        for u, v in a.items():
            for x, y in b.items():
                result[(u+x) % p] += v*y
        return {k:v for k,v in result.items() if v}
    previous = json.loads((ROOT/'results/parallel34_verification_2026_09_05.json').read_text())
    chosen = [(1153, [1,75,123,140,1013,1030,1078,1152])]
    chosen += [(c['p'],c['elements']) for c in previous['cases']
               if (c['p'],c['n']) in [(97,16),(33713,16)] or not c['subgroup']]
    output = []
    for p, A in chosen:
        N, limit = len(A)//2, 8
        reps = [a for a in A if a < (-a) % p]
        E = [{0:1}]+[{} for _ in range(limit)]
        powers = [Counter() for _ in range(limit+1)]
        for a in reps:
            pair = {a:1, (-a) % p:1}
            power = {0:1}
            for j in range(1,limit+1):
                power = multiply(power,pair,p)
                powers[j].update(power)
            for j in range(limit,0,-1):
                E[j] = dict(Counter(E[j])+Counter(multiply(E[j-1],pair,p)))
        for j in range(1,limit+1):
            formula = Counter({0:N*comb(j,j//2)}) if j % 2 == 0 else Counter()
            for l in range((j+1)//2):
                for a in A:
                    formula[((j-2*l)*a) % p] += comb(j,l)
            check(dict(formula) == dict(powers[j]), 'paired power sum identity')
            total = Counter()
            for k in range(1,j+1):
                for a,v in multiply(powers[k],E[j-k],p).items():
                    total[a] += (-1)**(k-1)*v
            check({a:v for a,v in total.items() if v} == {a:j*v for a,v in E[j].items()},
                  'group algebra Newton identity')
        output.append(dict(p=p,n=len(A),maximum_degree=limit,
                           unordered_zero_coefficients=[e.get(0,0) for e in E]))
    return output


def main():
    started = perf_counter()
    rows = []
    for N in [2, 4, 8, 16, 32, 64, 128]:
        patterns = [([F(2)]*N, 'constant positive'),
                    ([F(2 if j % 2 else -2) for j in range(N)], 'balanced endpoints'),
                    ([F(2 if j % 4 else -2) for j in range(N)], 'skew endpoints'),
                    ([F((7*j % 17)-8, 4) for j in range(N)], 'rational asymmetric')]
        for y, label in patterns:
            for s in [1, 2, 3, 4]:
                rows.append(row_case(y, s, label))
    # Reach the large-sum branch with nonconstant rows at higher degrees.
    for N, s in [(512, 1), (2048, 4)]:
        rows.append(row_case([F(2) if j % 5 else F(3, 2) for j in range(N)], s,
                             'large biased rational row'))
    fields = subgroup_counts()
    witness = averaged_root_counterexample()
    newton = newton_group_checks()
    for s in range(1, 33):
        check(1+16**s <= 2*16**s, 'reverse error absorption')
        for K in [1, 2, 5, 16]:
            check((16*K)**s+2*4096**s <= (16*K+8192)**s, 'single-degree constant')
    ten_constant = 90+1280**5
    check(ten_constant == 3435973836800090, 'ten-term lower constant')
    result = dict(status='single-degree finite checks passed; full targets unproved',
                  created_at_utc=datetime.now(timezone.utc).isoformat(),
                  duration_seconds=perf_counter()-started, check_counts=dict(CHECKS),
                  rational_rows=rows, finite_field_rows=fields, averaged_polynomial_witness=witness,
                  newton_group_cases=newton,
                  ten_term_lower=dict(constant=ten_constant,
                       conclusion='O_10 >= n^6 - (90+1280^5)n^5 for n>=90 and p<=n^4',
                       scope='A lower estimate only; no matching upper estimate or prime/subgroup existence claim.'),
                  input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in [
                      'experiments/parallel35_verify_2026_09_05.py',
                      'results/parallel34_verification_2026_09_05.json',
                      'sources/parallel35-single-degree/ravichandran-1609.04187v2.pdf']},
                  scope='Exact rational pointwise comparisons, Sturm enclosures, finite-field counts, and counterexample. The uniform proof uses the cited derivative-root theorem and is not Lean-verified or independently reviewed.')
    (ROOT/'results/parallel35_verification_2026_09_05.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],checks=sum(CHECKS.values()),
          rational_rows=len(rows),field_rows=len(fields),
          root_enclosures=CHECKS['exact Sturm root enclosure'],
          high_sum_cases=sum(x.get('high_sum_branch',False) for x in rows),
          duration_seconds=result['duration_seconds'],check_counts=dict(CHECKS)),indent=2))


if __name__ == '__main__':
    main()
