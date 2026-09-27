#!/usr/bin/env python3
"""Query exact shared-insertion covariance and its deletion-contrast projection."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from covariance import query
from contrast_analysis import analyze

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q',type=int,required=True)
    parser.add_argument('--selected',required=True,help='Comma-separated distinct field labels')
    parser.add_argument('--degree',type=int,default=6)
    parser.add_argument('--full-certificate',action='store_true')
    args=parser.parse_args()
    try:row,_=query(args.q,[int(x) for x in args.selected.split(',')],args.degree)
    except ValueError as exc:parser.error(str(exc))
    contrast=analyze(row)
    out={k:row[k] for k in ('q','n','degree','selected','target','neighbour_count','shared_insertions','marginal_neighbour_mean','marginal_neighbour_variance','covariance_trace_squared','seconds','scope')}
    for key in ('contrast_covariance_trace_squared','contrast_to_full_squared_norm_fraction'):out[key]=contrast[key]
    if args.full_certificate:out['covariance_certificate']=row;out['contrast_certificate']=contrast
    out['source_sha256']={name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('query_covariance.py','covariance.py','contrast_analysis.py','covariance_backend.cpp','covariance_backend')}
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
