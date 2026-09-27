#!/usr/bin/env python3
"""Cross-check older general conditioning through inclusion-exclusion."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import importlib.util
import json
from math import comb
from pathlib import Path
from preflight import literal
from two_insertions import from_matrix

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round4/marked_moments/conditional_moments.py'
spec=importlib.util.spec_from_file_location('old_conditioning',OLD)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def main():
    cases=[]
    for q,c,a in [(13,list(range(6)),[0,1]),(17,[0,1,2,3,4,6,10],[0,1])]:
        S=literal.literal_signs(q);base=set(c)-set(a);n=len(c)
        for degree in (2,6):
            first=second=0;count=0
            for size in range(3):
                for extra in combinations(a,size):
                    marks=sorted(base|set(extra));weight=(-1)**size*comb(q-len(marks),n-len(marks))
                    r=old.conditional_moments(old.matrix_counts(S,marks),n,degree)
                    first+=weight*F(*r['mean']);second+=weight*F(*r['second_moment']);count+=weight
            new=from_matrix(S,c,a,degree)
            assert count==new['distinct_final_sets'] and 2*first==new['target_sum'] and 2*second==new['target_square_sum']
            cases.append({'q':q,'n':n,'degree':degree,'older_marked_queries':4,'final_sets':count})
    controls=0
    q=17;c=[0,1,2,3,4,6,10];a=[0,1];S=literal.literal_signs(q)
    for d in range(7):
        initial=from_matrix(S,c,a,d)
        for multiplier,shift in [(1,5),(2,3),(3,7)]:
            cc=[(multiplier*x+shift)%q for x in c];aa=[(multiplier*x+shift)%q for x in a]
            r=from_matrix(S,cc,aa,d);flip=pow(multiplier,(q-1)//2,q)
            sign=1 if flip==1 else (-1)**d
            assert F(*r['mean'])==sign*F(*initial['mean']) and r['variance']==initial['variance']
            controls+=1
    out={'status':'passed','scope':'Older marked moment compiler with explicit exclusion of deleted points, plus affine controls in all degrees.',
         'cases':cases,'affine_controls':controls,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('interoperability_review.py','two_insertions.py','preflight.py')},
         'input_sha256':{'../round4/marked_moments/conditional_moments.py':sha256(OLD.read_bytes()).hexdigest(),
                         '../round7/local_edits/review_results.py':sha256((HERE.parent/'round7/local_edits/review_results.py').read_bytes()).hexdigest()}}
    (HERE/'interoperability_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('sha256')}),flush=True)


if __name__=='__main__':main()
