#!/usr/bin/env python3
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify


def main():
    c=json.loads((HERE/'consecutive36.json').read_text())['certificate'];prepared=certify(c)
    a,start=1234567,57;word=prepared[2]([a])[0];erased={(start+j)%c['N'] for j in range(36)}
    digits=[None if i in erased else v for i,v in enumerate(word)]
    path=HERE/'erasure36_input.json';path.write_text(json.dumps(digits,indent=2)+'\n')
    command=[sys.executable,str(HERE/'cover_query.py'),'--certificate',str(HERE/'consecutive36.json'),
             '--cover',str(HERE/'cover36.json'),'--input',str(path),'--cyclic-start',str(start)]
    result=subprocess.run(command,capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stderr;output=json.loads(result.stdout)
    assert output['difference_cover']['status']=='complete' and output['universal_unique_completion']
    assert output['status']=='complete' and output['count']==1 and output['completions']==[{'scalar':a,'digits':word}]
    (HERE/'erasure36_output.json').write_text(json.dumps(output,indent=2)+'\n')
    out={'status':'passed','cyclic_start':start,'erased_coordinates':36,'scalar':a,
         'candidate_box_size':output['candidate_box_size'],'scalar_candidates_checked':output['scalar_candidates_checked'],
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['consecutive36.json','cover36.json','erasure36_input.json']},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['query36_review.py','cover_query.py','../round22/conditioned_review.py']}}
    (HERE/'query36_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
