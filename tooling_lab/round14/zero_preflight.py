#!/usr/bin/env python3
"""Detect a profile that cannot host a zero column, then reserve it explicitly."""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random
from arc_profile import profile
from dimension_preflight import horner,rank,review

HERE=Path(__file__).resolve().parent


def main():
    sourcepath=HERE.parent/'round12/hard_rank3.certificate.json';data=json.loads(sourcepath.read_text())['input'];data={**data,'dimension':3}
    p,n,d=data['p'],data['n'],data['dimension'];columns=[tuple(horner(b,x,p) for b in data['basis']) for x in data['domain']]
    zeros=[i for i,v in enumerate(columns) if not any(v)];nonzero=[i for i in range(n) if i not in zeros]
    naive=profile(n,data['s'],d)
    assert len(zeros)==1 and naive['queries']==8 and all(b>=d for b in naive['sizes'])
    # Every coordinate in this profile belongs to a queried d-subset. Any
    # occurrence of the actual zero column makes that subset dependent, so
    # no permutation of THIS profile can be an arc partition.
    fixed=profile(n-len(zeros),data['s']-len(zeros),d);assert fixed['queries']==10
    rng=random.Random(1410);ids=nonzero[:];found=None
    for attempt in range(1,1001):
        rng.shuffle(ids);blocks=[];offset=0
        for size in fixed['sizes']:
            blocks.append(sorted(ids[offset:offset+size]));offset+=size
        blocks.extend([[z] for z in zeros])
        if all(rank([columns[i] for i in Q],p)==d for B in blocks for Q in combinations(B,d)):
            found=blocks;break
    assert found is not None
    checked=review(data,found);assert checked['queries_rank_checked']==10 and checked['maximum_query_avoiding_set']==10
    out={'status':'passed','zero_coordinates':zeros,'naive_profile':naive,
         'naive_profile_infeasibility':'All blocks in the8-query profile have size at least3. Every coordinate then lies in a required independent triple, impossible for the actual zero column.',
         'reduced_profile':fixed,'search_attempts':attempt,'input':data,'blocks':found,'review':checked,
         'scope':'A structural obstruction to the chosen unconstrained arc profile, repaired by reserving zero coordinates. No global minimum over arbitrary query hypergraphs is claimed.',
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['zero_preflight.py','arc_profile.py','dimension_preflight.py']},
         'input_sha256':{'../round12/hard_rank3.certificate.json':sha256(sourcepath.read_bytes()).hexdigest()}}
    (HERE/'zero_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','naive_queries':8,'reserved_zero_queries':10,'attempts':attempt,'review':checked}),flush=True)


if __name__=='__main__':main()
