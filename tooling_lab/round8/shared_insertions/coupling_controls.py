#!/usr/bin/env python3
"""Affine invariance and literal double-centering checks on the actual twins."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from covariance import query

HERE=Path(__file__).resolve().parent


def main():
    toy=json.loads((HERE/'results.json').read_text());contrasts=json.loads((HERE/'contrast_results.json').read_text())['cases']
    by_set={(r['q'],tuple(r['selected'])):r for r in contrasts};affine=entries=center_entries=0
    for pair in toy['pairs']:
        left,right=pair['members']
        assert left['flat_histogram']==right['flat_histogram']
        assert left['unlabelled_row_histograms']==right['unlabelled_row_histograms']
        assert left['diagonal_covariance_deck']==right['diagonal_covariance_deck']
        assert left['unlabelled_column_histograms']!=right['unlabelled_column_histograms']
        assert left['covariance_trace_squared']!=right['covariance_trace_squared']
        c1=by_set[(left['q'],tuple(left['selected']))];c2=by_set[(right['q'],tuple(right['selected']))]
        assert c1['total_covariance']==c2['total_covariance'] and c1['trace_covariance']==c2['trace_covariance']
        assert c1['contrast_covariance_trace_squared']!=c2['contrast_covariance_trace_squared']
        for member in pair['members']:
            q,n=member['q'],member['n'];m=q-n;c=member['selected'];E=member['delta_matrix']
            # Subtract the deletion mean for each insertion, then each row's
            # insertion mean. Form the covariance directly from these edits.
            centered=[[F(E[a][b])-F(sum(E[j][b] for j in range(n)),n) for b in range(m)] for a in range(n)]
            means=[sum(row)/m for row in centered]
            contrast=by_set[(q,tuple(c))]
            for a in range(n):
                for b in range(n):
                    value=sum((centered[a][j]-means[a])*(centered[b][j]-means[b]) for j in range(m))/m
                    assert value==F(contrast['contrast_covariance_numerator'][a][b],contrast['contrast_covariance_denominator'])
                    center_entries+=1
            base,_=query(q,c)
            for multiplier,shift in ((1,4),(3,4)):
                transformed=sorted((multiplier*x+shift)%q for x in c);actual,_=query(q,transformed)
                permutation=[transformed.index((multiplier*x+shift)%q) for x in c]
                assert actual['covariance_denominator']==base['covariance_denominator']
                for a in range(n):
                    for b in range(n):
                        assert base['covariance_numerator'][a][b]==actual['covariance_numerator'][permutation[a]][permutation[b]]
                        entries+=1
                affine+=1
    out={'status':'passed','affine_cases':affine,'affine_covariance_entries_checked':entries,'literal_centered_covariance_entries':center_entries,
         'scope':'Both translation and nonsquare-affine controls; literal projection replay and simpler-feature ablations on four saved twin inputs.',
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('coupling_controls.py','covariance.py','covariance_backend.cpp','covariance_backend')},
         'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('results.json','contrast_results.json')}}
    (HERE/'coupling_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if 'sha256' not in k}),flush=True)


if __name__=='__main__':main()
