#!/usr/bin/env python3
"""Determine exactly which two-swap information the saved twins share."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'preflight_results.json').read_text());out=[]
    for pair in data['pairs']:
        members=[]
        for member in pair['members']:
            hist=member['histogram'];N=sum(n for x,n in hist);mean=F(*member['mean'])
            central=[sum(n*(F(x)-mean)**d for x,n in hist)/N for d in range(1,7)]
            members.append({'selected':member['selected'],'histogram':hist,
                            'central_moments_1_through_6':[[x.numerator,x.denominator] for x in central],
                            'fixed_deletion_mean_deck':sorted(F(*x['result']['mean']) for x in member['fibres']),
                            'fixed_deletion_variance_deck':sorted(F(*x['result']['variance']) for x in member['fibres'])})
        equality={k:members[0][k]==members[1][k] for k in ('histogram','fixed_deletion_mean_deck','fixed_deletion_variance_deck')}
        for r in members:
            for k in ('fixed_deletion_mean_deck','fixed_deletion_variance_deck'):r[k]=[[x.numerator,x.denominator] for x in r[k]]
        out.append({'members':members,'equal_features':equality})
    result={'status':'passed','pairs':out,'source_sha256':{'analyze_twins.py':sha256(Path(__file__).read_bytes()).hexdigest()},
            'input_sha256':{'preflight_results.json':sha256((HERE/'preflight_results.json').read_bytes()).hexdigest()}}
    (HERE/'twin_analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    for pair in out:print(json.dumps({'equal_features':pair['equal_features'],'members':[{k:v for k,v in m.items() if not k.endswith('_deck')} for m in pair['members']]}),flush=True)


if __name__=='__main__':main()
