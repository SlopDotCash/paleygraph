#!/usr/bin/env python3
"""Use a checked difference cover alongside a complete erasure query."""
import argparse
from pathlib import Path
import json
import sys
from cover_review import verify

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from erasure_query import query


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,required=True)
    parser.add_argument('--cover',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--candidate-budget',type=int,default=50000)
    parser.add_argument('--cyclic-start',type=int)
    args=parser.parse_args()
    try:
        data=json.loads(args.certificate.read_text());c=data.get('certificate',data)
        data=json.loads(args.cover.read_text());proof=data.get('cover',data)
        checked=verify(c,proof)
        out=query(c,json.loads(args.input.read_text()),args.candidate_budget,args.cyclic_start)
        out['difference_cover']={'status':'complete' if checked['universal_unique_completion'] else 'incomplete',
                                 **checked}
        out['universal_unique_completion']=out['universal_unique_completion'] or checked['universal_unique_completion']
        if out['status']=='complete' and out['universal_unique_completion']:assert out['count']<=1
    except (ValueError,AssertionError,TypeError,KeyError,IndexError) as error:
        parser.exit(2,'invalid input or certificate: '+str(error)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
