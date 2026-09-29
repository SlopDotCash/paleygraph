#!/usr/bin/env python3
"""Independent standard-library readback, literal tiny witnesses and Gram comparisons.

No production module import. Completed larger calculations are bound by hashes
and checked entrywise against the separate synthetic-division implementation.
"""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def digest(path):return sha256(path.read_bytes()).hexdigest()


def literal_signs(q):
    def symbol(a):
        if a%q==0:return 0
        return 1 if pow(a%q,(q-1)//2,q)==1 else -1
    return [[symbol(x-y) for y in range(q)] for x in range(q)]


def elementary(values,d):
    result=0
    for chosen in combinations(values,d):
        term=1
        for value in chosen:term*=value
        result+=term
    return result


def target(signs,c,d):
    return sum(elementary([row[a] for a in c],d) for row in signs)


def check_algebra(record):
    q,n=record['q'],record['n'];c=record['selected'];d=record['degree']
    assert len(c)==n and len(set(c))==n and sorted(c)==c and 0<=min(c)<=max(c)<q
    assert record['neighbour_count']==n*(q-n)
    assert len(record['deletion_records'])==n and len(record['internal_contraction'])==n
    first=second=complete=inside=0
    for a,(row,deriv) in enumerate(zip(record['internal_contraction'],record['deletion_records'])):
        assert len(row)==n and deriv['point']==c[a]
        norm=q*deriv['insertion_norm_squared']-deriv['insertion_sum']**2
        assert norm==deriv['complete_transform_norm_squared']
        u=row[a];t=sum(row);w=sum(v*v for v in row)
        fd=-t-(q-n)*u;sd=norm-w+2*u*t+(q-n)*u*u
        assert deriv['delta_sum']==fd and deriv['delta_square_sum']==sd
        assert norm>=w and sd>=0
        first+=fd;second+=sd;complete+=norm;inside+=w
    assert complete==record['complete_insertion_energy'] and inside==record['excluded_internal_energy']
    mean=F(first,n*(q-n));sq=F(second,n*(q-n));T=record['target']
    assert F(*record['delta_mean'])==mean and F(*record['delta_second_moment'])==sq
    assert F(*record['neighbour_mean'])==T+mean
    assert F(*record['neighbour_second_moment'])==T*T+2*T*mean+sq
    assert F(*record['neighbour_variance'])==sq-mean*mean>=0


def check_gram(record,gram):
    for key in ('q','n','degree','selected','target','internal_contraction'):
        assert record[key]==gram[key],key
    assert [x['insertion_sum'] for x in record['deletion_records']]==gram['insertion_sums']
    assert [x['insertion_norm_squared'] for x in record['deletion_records']]==gram['insertion_norms_squared']


def direct_review(record):
    q,n,d=record['q'],record['n'],record['degree'];c=record['selected'];S=literal_signs(q)
    assert record['target']==target(S,c,d)
    for a in range(n):
        r=[elementary([S[x][c[b]] for b in range(n) if b!=a],d-1) for x in range(q)]
        assert sum(r)==record['deletion_records'][a]['insertion_sum']
        assert sum(x*x for x in r)==record['deletion_records'][a]['insertion_norm_squared']
        assert [sum(S[x][b]*r[x] for x in range(q)) for b in c]==record['internal_contraction'][a]
    hist={}
    for row in S:
        values=[row[a] for a in c];key=(values.count(0),values.count(1));hist[key]=hist.get(key,0)+1
    assert [[*key,count] for key,count in sorted(hist.items())]==record['signed_row_count_histogram']
    values=[target(S,sorted((set(c)-{a})|{b}),d) for a in c for b in set(range(q))-set(c)]
    mean=F(sum(values),len(values));second=F(sum(x*x for x in values),len(values))
    assert F(*record['neighbour_mean'])==mean
    assert F(*record['neighbour_second_moment'])==second
    assert F(*record['neighbour_variance'])==second-mean*mean
    return len(values),[min(values),max(values)]


def main():
    data=json.loads((HERE/'results.json').read_text())
    gram=json.loads((HERE/'gram_results.json').read_text())
    ablation=json.loads((HERE/'ablation_results.json').read_text())
    binding_count=0
    for saved in (data,gram,ablation):
        for group in ('source_sha256','input_sha256'):
            for name,value in saved.get(group,{}).items():assert digest(HERE/name)==value;binding_count+=1
    assert gram['status']=='passed' and ablation['status']=='passed'
    assert len(data['scale_cases'])==len(gram['scale_cases'])==8
    gram_map={(r['q'],tuple(r['selected'])):r for r in gram['scale_cases']}
    entries=deletions=0
    for record in data['scale_cases']:
        check_algebra(record);check_gram(record,gram_map[(record['q'],tuple(record['selected']))])
        entries+=record['n']**2;deletions+=record['n']
    witnesses=[data['collision_search']['witness']]
    quartic_pairs=0
    for case in ablation['cases']:
        witnesses.extend(m['witness'] for m in case['models'] if m['witness'])
        for model in case['models']:
            if model['model']=='histogram_and_quartic_deck' and model['witness']:
                pair=model['witness'];S=literal_signs(case['q'])
                decks=[sorted(target(S,list(c),4) for c in combinations(pair[side]['selected'],4)) for side in ('lower','upper')]
                assert decks[0]==decks[1]
                quartic_pairs+=1
    unique={}
    for pair in witnesses:
        low,high=pair['lower'],pair['upper']
        assert low['signed_row_count_histogram']==high['signed_row_count_histogram']
        assert low['target']==high['target'] and low['neighbour_mean']==high['neighbour_mean']
        assert F(*pair['variance_gap'])==F(*high['neighbour_variance'])-F(*low['neighbour_variance'])>0
        for row in (low,high):unique[(row['q'],tuple(row['selected']))]=row
    direct_values=0
    for row in unique.values():
        check_algebra(row);count,_=direct_review(row);direct_values+=count
    corruptions=0
    for key in ('target','neighbour_count','internal_contraction','delta_mean','neighbour_variance','complete_insertion_energy','deletion_records'):
        bad=deepcopy(data['scale_cases'][0])
        if key=='internal_contraction':bad[key][0][0]+=1
        elif key=='deletion_records':bad[key][0]['insertion_norm_squared']+=1
        elif isinstance(bad[key],list):bad[key][0]+=1
        else:bad[key]+=1
        try:check_algebra(bad);check_gram(bad,gram['scale_cases'][0])
        except AssertionError:corruptions+=1
        else:raise AssertionError(('corruption accepted',key))
    out={'status':'passed','scope':'Independent literal witness enumeration, exact rational algebra and entrywise separate-backend scale comparisons; no production import, no direct large-neighbour census.',
         'scale_cases':8,'scale_internal_entries_compared':entries,'scale_deletion_sums_and_norms_compared':deletions,
         'unique_toy_witness_sets':len(unique),'literal_toy_neighbour_targets':direct_values,'equal_quartic_deck_witness_pairs':quartic_pairs,'corrupted_records_rejected':corruptions,
         'prior_source_bindings_checked':binding_count,
         'source_sha256':{'review_results.py':digest(Path(__file__))},
         'input_sha256':{name:digest(HERE/name) for name in ('results.json','gram_results.json','ablation_results.json')}}
    (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if 'sha256' not in k}),flush=True)


if __name__=='__main__':main()
