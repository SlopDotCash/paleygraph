#!/usr/bin/env python3
"""Usable fast local-moment query plus exact one-sided finite-neighbour bounds."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from time import perf_counter

HERE=Path(__file__).resolve().parent


def enc(x):return [x.numerator,x.denominator]


def summarize(gram):
    q,n=gram['q'],gram['n'];V=gram['internal_contraction']
    assert 0<n<q and len(V)==n and all(len(row)==n for row in V)
    assert len(gram['insertion_sums'])==len(gram['insertion_norms_squared'])==n
    first=second=0
    for a,row in enumerate(V):
        total=q*gram['insertion_norms_squared'][a]-gram['insertion_sums'][a]**2
        inside=sum(x*x for x in row);u=row[a];t=sum(row)
        assert total>=inside
        first-=t+(q-n)*u
        second+=total-inside+2*u*t+(q-n)*u*u
    drift=F(first,n*(q-n));sq=F(second,n*(q-n));var=sq-drift*drift
    assert var>=0
    return {'q':q,'n':n,'degree':gram['degree'],'selected':gram['selected'],'target':gram['target'],
            'neighbour_count':n*(q-n),'delta_mean':enc(drift),'delta_second_moment':enc(sq),
            'neighbour_mean':enc(gram['target']+drift),'neighbour_variance':enc(var)}


def threshold_bound(summary, threshold, direction='upper'):
    """Cantelli, exactly rational, then integer floor for a finite uniform set.

    upper means T(C') >= threshold; lower means T(C') <= threshold.
    A threshold on the wrong side of the mean returns only the trivial bound.
    """
    if direction not in ('upper','lower'):raise ValueError('unknown direction')
    mean=F(*summary['neighbour_mean']);var=F(*summary['neighbour_variance'])
    if var<0 or summary['neighbour_count']<1:raise ValueError('invalid moments')
    threshold=F(threshold);gap=threshold-mean if direction=='upper' else mean-threshold
    probability=var/(var+gap*gap) if gap>0 else F(1)
    count=(summary['neighbour_count']*probability.numerator)//probability.denominator
    return {'event':f'T_after_swap {">=" if direction=="upper" else "<="} {threshold}',
            'threshold':enc(threshold),'direction':direction,'probability_upper':enc(probability),
            'neighbour_count_upper':count,'complement_count_lower':summary['neighbour_count']-count,
            'method':'Exact Cantelli bound and finite-count floor' if gap>0 else 'Trivial bound: threshold is on the other side of the mean'}


def query(q,selected,degree=6):
    selected=list(selected)
    if type(q)is not int or type(degree)is not int or any(type(x)is not int for x in selected):
        raise ValueError('integer q, degree and selected labels required')
    if not (5<=q<=10000000 and 1<=len(selected)<=64 and len(selected)<q and 0<=degree<=min(6,len(selected))):
        raise ValueError('outside prototype parameter bounds')
    if len(set(selected))!=len(selected) or any(x<0 or x>=q for x in selected):raise ValueError('invalid selected set')
    text=f'{q} {len(selected)} {degree}\n'+' '.join(map(str,selected))+'\n'
    start=perf_counter()
    run=subprocess.run([str(HERE/'gram_oracle')],input=text,text=True,capture_output=True)
    if run.returncode:raise ValueError(run.stderr.strip())
    gram=json.loads(run.stdout);gram['selected']=sorted(selected)
    assert (gram['q'],gram['n'],gram['degree'])==(q,len(selected),degree)
    out=summarize(gram)
    out['query_seconds']=perf_counter()-start
    out['source_sha256']={name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('query.py','gram_oracle.cpp','gram_oracle')}
    out['scope']='Exact local moments at the specified actual Paley set; no uniform bound or direct neighbour census.'
    return out,gram


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q',type=int,required=True)
    parser.add_argument('--selected',required=True,help='Comma-separated field labels')
    parser.add_argument('--degree',type=int,default=6)
    parser.add_argument('--upper-threshold',type=int)
    parser.add_argument('--lower-threshold',type=int)
    parser.add_argument('--full-certificate',action='store_true')
    args=parser.parse_args()
    out,gram=query(args.q,[int(x) for x in args.selected.split(',')],args.degree)
    out['threshold_bounds']=[threshold_bound(out,t,d) for t,d in ((args.upper_threshold,'upper'),(args.lower_threshold,'lower')) if t is not None]
    if args.full_certificate:out['contraction_certificate']=gram
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
