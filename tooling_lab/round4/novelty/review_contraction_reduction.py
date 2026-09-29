#!/usr/bin/env python3
"""Direct marked Walsh contractions and complete small subset checks."""
import sys
sys.dont_write_bytecode = True
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'marked_moments'/'contraction_reduction.py'
sys.path.insert(0, str(SOURCE.parent))


def load(path):
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def direct_moments(S, M, n, degrees):
    population = sorted(set(range(len(S)))-set(M))
    sums = {d: [0, 0] for d in degrees}
    tables = {d: {(a,b): sum((-1)**j*choose(b,j)*choose(a,d-j) for j in range(d+1))
                   for a in range(n+1) for b in range(n+1-a)} for d in degrees}
    count = 0
    for rest in combinations(population, n-len(M)):
        C = tuple(M)+rest
        rows = [([row[x] for x in C].count(1), [row[x] for x in C].count(-1)) for row in S]
        for d in degrees:
            value = sum(tables[d][a,b] for a,b in rows)
            sums[d][0] += value
            sums[d][1] += value*value
        count += 1
    return count, {d: (Fraction(t,count), Fraction(s,count)) for d,(t,s) in sums.items()}


def check_low_contractions(S, M):
    q = len(S)
    subsets = [I for d in range(4) for I in combinations(range(3), d)]
    f = {I: [0 if x in M else prod(S[x][M[i]] for i in I) for x in range(q)] for I in subsets}
    h = {I: sum(v) for I,v in f.items()}
    checked = 0
    for I in subsets:
        for J in subsets:
            if min(len(I), len(J)) > 1:
                continue
            L,R = (I,J) if len(I) <= len(J) else (J,I)
            toggle = lambda A,i: tuple(sorted(set(A) ^ {i}))
            if not L:
                expected = -sum(h[toggle(R,i)] for i in range(3))
            else:
                expected = -h[R]-sum(S[M[j]][M[L[0]]]*h[toggle(R,j)] for j in range(3))
            actual = sum(f[I][x]*S[x][y]*f[J][y] for x in range(q) for y in range(q))
            assert actual == expected, (M,I,J,actual,expected)
            checked += 1
    h2 = [sum(f[I][x] for I in subsets if len(I) == 2) for x in range(q)]
    h3 = f[(0,1,2)]
    Q = sum(h2[x]*S[x][y]*h3[y] for x in range(q) for y in range(q))
    return Q, checked


def main():
    source_hash = sha256(SOURCE.read_bytes()).hexdigest()
    reduction = load(SOURCE)
    base = load(HERE/'review_conditional_compiler.py')
    matrices = base.matrices()
    p = 29
    chi = [0 if x == 0 else (1 if pow(x,14,29) == 1 else -1) for x in range(29)]
    matrices['Paley29'] = [[chi[(y-x)%p] for x in range(p)] for y in range(p)]
    cases = [('Paley17',[0,1,2],3,[0,1,2,3]),
             ('Paley17',[0,1,3],4,[1,2,3,4]),
             ('Paley17',[0,1,2],6,[3,6]),
             ('Paley17',[0,1,3],7,[3,6]),
             ('Paley29',[1,4,0],7,[6]),
             ('Paley29',[1,9,0],7,[6]),
             ('Paley49',[0,1,7],4,[3,4]),
             ('Peisert49',[0,1,7],4,[3,4])]
    records = []
    contractions = completions = 0
    for name,marks,n,degrees in cases:
        S = matrices[name]
        counts = base.load(HERE.parent/'marked_moments'/'conditional_moments.py').matrix_counts(S,marks)
        Q, checks = check_low_contractions(S,marks)
        contractions += checks
        count, direct = direct_moments(S,marks,n,degrees)
        completions += count
        for degree in degrees:
            actual = reduction.reduce_three_marks(counts,n,degree)
            assert actual['Q'] == Q
            assert Fraction(*actual['second_moment']) == direct[degree][1]
            assert Fraction(*actual['mean']) == direct[degree][0]
            assert Fraction(*actual['variance']) == direct[degree][1]-direct[degree][0]**2
            cells = {tuple(cell['pattern']):cell['size'] for cell in counts['cells']}
            cell_input = reduction.reduce_from_cells(len(S),marks,cells,Q,n,degree)
            assert all(cell_input[key] == actual[key] for key in ('mean','second_moment','variance','cell_only_term','Q_coefficient'))
            walsh,c = reduction.walsh_decomposition(len(S),n,degree)
            for u in reduction.PATTERNS:
                for v in reduction.PATTERNS:
                    for permutation in ((1,0,2),(2,0,1)):
                        uu = tuple(u[i] for i in permutation)
                        vv = tuple(v[i] for i in permutation)
                        gap = reduction.coefficient(len(S),n,degree,u,v,1)-reduction.coefficient(len(S),n,degree,u,v,-1)
                        other = reduction.coefficient(len(S),n,degree,uu,vv,1)-reduction.coefficient(len(S),n,degree,uu,vv,-1)
                        assert gap == other
            if n == degree:
                assert actual['Q_coefficient'] == [0,1]
            records.append({'graph':name,'marks':marks,'n':n,'degree':degree,'completions':count,
                            'Q':Q,'cell_only_term':actual['cell_only_term'],
                            'Q_coefficient':actual['Q_coefficient'],'second_moment':actual['second_moment']})
        print(name,marks,n,'checked',len(degrees),flush=True)
    witness = [r for r in records if r['graph']=='Paley29']
    assert witness[0]['cell_only_term'] == witness[1]['cell_only_term'] == [941207,7475]
    assert witness[0]['Q_coefficient'] == witness[1]['Q_coefficient'] == [-2,7475]
    assert [r['Q'] for r in witness] == [-42,86]
    assert Fraction(*witness[0]['second_moment'])-Fraction(*witness[1]['second_moment']) == Fraction(256,7475)
    assert sha256(SOURCE.read_bytes()).hexdigest() == source_hash
    result = {'status':'passed','date':'2026-09-05','complete_subset_moment_checks':len(records),
              'complete_subsets_enumerated':completions,'direct_low_degree_contractions':contractions,
              'records':records,'candidate_sha256':source_hash,
              'reviewer_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'contraction_review.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
