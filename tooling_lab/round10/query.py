#!/usr/bin/env python3
"""Exact prime-field query for every deletion-pair mean."""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent


def query(q,selected,degree=6):
    c=sorted(selected);n=len(c)
    if (type(q) is not int or type(degree) is not int or not (5<=q<=10_000_000 and 2<=n<=min(64,q-2)
        and 0<=degree<=min(6,n)) or len(set(c))!=n or any(type(x) is not int or not 0<=x<q for x in c)):
        raise ValueError('invalid field, selected set or degree')
    begin=time.monotonic()
    payload=' '.join(map(str,[q,n,degree,*c]))+'\n'
    p=subprocess.run([str(HERE/'pair_means_backend')],input=payload,text=True,capture_output=True)
    if p.returncode:raise ValueError(p.stderr.strip())
    row=json.loads(p.stdout);N=row['ordered_insertion_pairs_per_deletion']
    for pair in row['pairs']:
        f=F(pair['target_sum'],N);pair['mean']=[f.numerator,f.denominator]
    row['seconds']=time.monotonic()-begin
    return row


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q',type=int);parser.add_argument('--selected',type=int,nargs='+',required=True)
    parser.add_argument('--degree',type=int,default=6)
    a=parser.parse_args();print(json.dumps(query(a.q,a.selected,a.degree),indent=2))
