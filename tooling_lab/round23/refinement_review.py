#!/usr/bin/env python3
"""Check both the refined proof and preservation of its previously saved work."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
from cover_review import verify

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('files',nargs='+',type=Path)
    parser.add_argument('--output',type=Path,default=HERE/'refinement_review.json');args=parser.parse_args()
    inputs={};cases=[]
    for path in args.files:
        new=json.loads(path.read_text());name=new['name'].removesuffix('_refined');oldpath=path.parent/(name+'.json')
        old=json.loads(oldpath.read_text());source=path.parent/new['source'];c=json.loads(source.read_text())['certificate']
        for p in [path,oldpath]:inputs[p.name]=sha256(p.read_bytes()).hexdigest()
        inputs[new['source']]=sha256(source.read_bytes()).hexdigest()
        before,after=old['cover'],new['cover'];checked=verify(c,after)
        assert after['roots']==before['roots'] and after['cuts'][:len(before['cuts'])]==before['cuts']
        assert len(after['nodes'])>=len(before['nodes'])
        unresolved=set(before['unresolved_leaves']);preserved=0
        for i,node in enumerate(before['nodes']):
            current=after['nodes'][i]
            if i not in unresolved:assert node==current;preserved+=1
            else:
                assert node['input_box']==current['input_box']
                assert current['trace'][:len(node['trace'])]==node['trace']
        changes=after['refinement']
        assert changes=={'old_nodes':len(before['nodes']),'old_cuts':len(before['cuts']),
                         'old_unresolved_leaves':before['unresolved_leaves'],
                         'added_cuts':len(after['cuts'])-len(before['cuts']),
                         'added_nodes':len(after['nodes'])-len(before['nodes'])}
        result={'name':new['name'],**checked,'earlier_finished_nodes_preserved':preserved,**changes}
        cases.append(result);print(json.dumps(result),flush=True)
    out={'status':'passed','cases':cases,'input_sha256':inputs,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['refinement_review.py','cover_review.py']}}
    args.output.write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
