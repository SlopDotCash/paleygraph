#!/usr/bin/env python3
"""Avoid the n=d shortcut; test whether already available quartic data suffices."""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter
import local_edits as tool
from run_experiments import direct_signs, direct_target, direct_neighbours, check

HERE=Path(__file__).resolve().parent


def main():
    start=perf_counter();q=17;matrix=direct_signs(q)
    quartics={c:direct_target(matrix,c,4) for c in combinations(range(q),4)}
    rows=[]
    for n in (6,7,8):
        cache={c:direct_target(matrix,c,6) for c in combinations(range(q),n)}
        fibres={key:defaultdict(list) for key in ('histogram','histogram_and_quartic_deck','histogram_and_deletion_summaries')}
        count=0
        for tail in combinations(range(2,q),n-2):
            c=(0,1)+tail;r=tool.from_prime(q,list(c));oracle=direct_neighbours(matrix,c,6,cache);check(r,oracle)
            hist=tuple(map(tuple,r['signed_row_count_histogram']))
            deck=tuple(sorted(quartics[a] for a in combinations(c,4)))
            summaries=tuple(sorted((x['insertion_sum'],x['insertion_norm_squared'],r['internal_contraction'][a][a],sum(r['internal_contraction'][a])) for a,x in enumerate(r['deletion_records'])))
            for model,key in (('histogram',hist),('histogram_and_quartic_deck',(hist,deck)),('histogram_and_deletion_summaries',(hist,summaries))):
                fibres[model][key].append(r)
            count+=1
        models=[]
        for model,groups in fibres.items():
            gaps=[]
            for key,records in groups.items():
                low=min(records,key=lambda r:F(*r['neighbour_variance']))
                high=max(records,key=lambda r:F(*r['neighbour_variance']))
                gap=F(*high['neighbour_variance'])-F(*low['neighbour_variance'])
                if gap:gaps.append((gap,low,high))
            gaps.sort(key=lambda row:row[0],reverse=True)
            witness=None
            if gaps:
                gap,low,high=gaps[0]
                witness={'lower':low,'upper':high,'variance_gap':tool.frac(gap)}
            models.append({'model':model,'fibres':len(groups),'nontrivial_fibres':sum(len(v)>1 for v in groups.values()),
                           'ambiguous_variance_fibres':len(gaps),'witness':witness})
        row={'q':q,'n':n,'normalized_inputs':count,'all_set_targets_cached':len(cache),'direct_neighbour_targets':count*n*(q-n),'models':models}
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='models'}),flush=True)
        print(json.dumps([{k:v for k,v in m.items() if k!='witness'} for m in models]),flush=True)
    out={'status':'passed','scope':'Complete normalized actual-set scan at q17; quartic-deck and deletion-summary ablations. No sufficiency beyond this finite domain.',
         'cases':rows,'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('local_edits.py','run_experiments.py','ablate_information.py')},'elapsed_seconds':perf_counter()-start}
    (HERE/'ablation_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
