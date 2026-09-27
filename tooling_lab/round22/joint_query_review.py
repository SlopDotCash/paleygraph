#!/usr/bin/env python3
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys
from tempfile import TemporaryDirectory
from joint_query import verify_uniqueness

HERE=Path(__file__).resolve().parent


def main():
    c=json.loads((HERE/'metric_consecutive32.json').read_text())['certificate']
    proof=json.loads((HERE/'coupled_summary.json').read_text())
    command=[sys.executable,str(HERE/'joint_query.py'),'--certificate',str(HERE/'metric_consecutive32.json'),
             '--joint-certificate',str(HERE/'coupled_summary.json'),'--input',str(HERE/'erasure_input.json'),'--cyclic-start','57']
    run=subprocess.run(command,capture_output=True,text=True,timeout=120)
    assert run.returncode==0,run.stderr
    result=json.loads(run.stdout); original=json.loads((HERE/'erasure_output.json').read_text())
    assert result['universal_unique_completion'] and result['completions']==original['completions']
    assert result['count']==1 and result['completions'][0]['scalar']==1234567
    (HERE/'joint_output.json').write_text(json.dumps(result,indent=2)+'\n')
    bad=deepcopy(proof); bad['target_separators'].pop()
    try: verify_uniqueness(c,bad)
    except AssertionError: pass
    else: raise AssertionError('incomplete coverage accepted')
    out={'status':'passed','controls':['cyclic32_joint_cli','missing_target_certificate_rejected'],
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['metric_consecutive32.json','coupled_summary.json','erasure_input.json','erasure_output.json']},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['joint_query_review.py','joint_query.py','coupled_review.py','erasure_query.py']}}
    (HERE/'joint_query_review.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out),flush=True)


if __name__=='__main__': main()
