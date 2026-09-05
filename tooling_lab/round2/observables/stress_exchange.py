#!/usr/bin/env python3
"""Apply the frozen exchange diagnostic to specified arithmetic stress inputs."""
from hashlib import sha256
import json
from pathlib import Path
from exchange_spectrum import components,frac

HERE=Path(__file__).resolve().parent


def generator(p):
    factors=[];n=p-1;d=2
    while d*d<=n:
        if n%d==0:
            factors.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:factors.append(n)
    return next(g for g in range(2,p) if all(pow(g,(p-1)//q,p)!=1 for q in factors))


def main():
    cases=[]
    for p,n in [(61,6),(1297,6),(2437,7),(4129,8)]:
        g=generator(p);assert (p-1)%n==0
        subgroup=sorted(pow(g,i*((p-1)//n),p) for i in range(n))
        data=json.loads((HERE/f'data_{p}.json').read_text())
        choices=[('multiplicative_subgroup',subgroup),('arithmetic_progression',list(range(n))),
                 ('geometric_progression',sorted(pow(g,i,p) for i in range(n))),
                 ('random_cohort_max_T6',max(data['random'],key=lambda r:r['t6'])['C']),
                 ('random_cohort_min_T6',min(data['random'],key=lambda r:r['t6'])['C']),
                 ('three_anchor_max_T6',max(data['planted'],key=lambda r:r['t6'])['C']),
                 ('random_max_abs_triangle',max(data['random'],key=lambda r:abs(r['triangle']))['C'])]
        rows=[]
        for name,c in choices:
            values,av=components(c,p)
            rows.append({'name':name,'C':c,'T6':int(av[0]),'components':[frac(v) for v in values],
                         'highest_component':float(values[6]),
                         'nonconstant_lower_sum':float(sum(values[1:6])),
                         'one_swap_expected_T6':frac(av[1])})
        cases.append({'p':p,'n':n,'selected_stress_inputs':rows})
        print(p,[(r['name'],r['T6'],round(r['nonconstant_lower_sum'],6)) for r in rows],flush=True)
    pair_inputs=[[[0,1,2,3,7,9],[0,1,2,9,10,17]],[[0,1,2,4,38,52],[0,1,3,8,10,21]]]
    pairs=[]
    for left,right in pair_inputs:
        a,_=components(left,61);b,_=components(right,61);differences=[y-x for x,y in zip(a,b)]
        pairs.append({'p':61,'C':left,'D':right,'component_differences':[frac(v) for v in differences],
                      'target_difference':frac(sum(differences)),
                      'all_lower_components_equal':all(v==0 for v in differences[:6])})
    print('equal-feature pairs',pairs,flush=True)
    out={'status':'selected pointwise arithmetic tests; not exhaustive exceptional-set search',
         'cases':cases,'old_equal_feature_pairs':pairs,'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'stress_exchange_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
