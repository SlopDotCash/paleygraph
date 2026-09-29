#!/usr/bin/env python3
"""Independently verify a saved cover and emit a bound review report."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
from cover_review import verify

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('cover',type=Path)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    data=json.loads(args.cover.read_text());source=args.cover.parent/data['source']
    c=json.loads(source.read_text())['certificate'];result=verify(c,data['cover'])
    out={'status':'passed','name':data['name'],**result,
         'input_sha256':{args.cover.name:sha256(args.cover.read_bytes()).hexdigest(),data['source']:sha256(source.read_bytes()).hexdigest()},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['review_one.py','cover_review.py']}}
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
