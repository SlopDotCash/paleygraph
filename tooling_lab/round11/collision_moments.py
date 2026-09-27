#!/usr/bin/env python3
"""Exact moments of the pinning census, with a separate overlap count."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    path=HERE/'rank_three_review.json';source=json.loads(path.read_text());n,s=source['n'],source['s'];rows=[]
    def containment(size):return F(comb(n-size,s-size),comb(n,s)) if size<=s else F(0)
    for row in source['rows']:
        histogram=row['agreement_histogram'];total=sum(count for value,count in histogram)
        moments=[F(sum(value**power*count for value,count in histogram),total) for power in range(1,7)]
        edges=list(map(set,row['dependent_triples']));unions=Counter(len(a|b) for a,b in combinations(edges,2))
        first_dep=len(edges)*containment(3)
        second_dep=first_dep+2*sum(count*containment(size) for size,count in unions.items())
        target=comb(s,3)
        assert moments[0]==target-first_dep and moments[1]==target**2-2*target*first_dep+second_dep
        rows.append({'third_polynomial_coefficients':row['third_polynomial_coefficients'],
                     'raw_moments_1_through_6':[[x.numerator,x.denominator] for x in moments],
                     'variance':[(moments[1]-moments[0]**2).numerator,(moments[1]-moments[0]**2).denominator],
                     'unordered_dependent_triple_pair_union_sizes':[[a,b] for a,b in sorted(unions.items())]})
    equality=[a==b for a,b in zip(*(r['raw_moments_1_through_6'] for r in rows))]
    assert equality==[True,True,False,False,False,False]
    out={'status':'passed','scope':'Complete eleven-set distributions for the explicit two-space collision; moments do not imply a uniform bound.',
         'rows':rows,'moment_equality_1_through_6':equality,
         'source_sha256':{'collision_moments.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'rank_three_review.json':sha256(path.read_bytes()).hexdigest()}}
    (HERE/'collision_moments.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','moment_equality_1_through_6':equality}),flush=True)


if __name__=='__main__':main()
