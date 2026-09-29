#!/usr/bin/env python3
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import subprocess
import sys
from cover_review import verify

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round22'))
from conditioned_review import certify,replay


def main():
    c=json.loads((HERE/'consecutive33.json').read_text())['certificate'];prepared=certify(c)
    proof=json.loads((HERE/'cover33_refined.json').read_text())['cover'];checked=verify(c,proof)
    data=json.loads((HERE/'query_experiments.json').read_text());count=0
    for record in data['cyclic_queries']:
        start=record['start'];digits=record['input'];N=c['N'];d=len(c['erased'])
        known={j:digits[(j+start)%N]*(1 if j+start<N else -1) for j in range(d,N)}
        assert list(map(list,known.items()))==record['rotated_known']
        expected=replay(c,known,50000,prepared)
        scalars=sorted(x['scalar']*pow(c['g'],-start,c['p'])%c['p'] for x in expected['completions'])
        expected['completions']=[{'scalar':a,'digits':w} for a,w in zip(scalars,prepared[2](scalars))]
        assert expected==record['output'] and scalars==[1234567];count+=expected['scalar_candidates_checked']
    controls=[]
    with TemporaryDirectory(prefix='round23-query-') as directory:
        temp=Path(directory)
        def run(cert,cover,digits,start=None,budget=50000):
            for n,v in [('cert',cert),('cover',cover),('input',digits)]:
                (temp/(n+'.json')).write_text(json.dumps(v))
            command=[sys.executable,str(HERE/'cover_query.py'),'--certificate',str(temp/'cert.json'),
                     '--cover',str(temp/'cover.json'),'--input',str(temp/'input.json'),'--candidate-budget',str(budget)]
            if start is not None:command+=['--cyclic-start',str(start)]
            return subprocess.run(command,capture_output=True,text=True,timeout=120)
        digits=json.loads((HERE/'erasure_input.json').read_text());result=run(c,proof,digits,57)
        assert result.returncode==0,result.stderr;output=json.loads(result.stdout)
        assert output['difference_cover']['status']=='complete' and output['universal_unique_completion']
        assert output['completions']==data['cyclic_queries'][6]['output']['completions']
        (HERE/'erasure_output.json').write_text(json.dumps(output,indent=2)+'\n');controls.append('cyclic33_verified_cli')
        # Uniqueness certification does not promise enumeration within budget.
        source=json.loads((HERE/'consecutive33.json').read_text())
        r=next(r for r in source['records'] if r['kind']=='actual' and r['output']['candidate_box_size']>1)
        values=dict(r['known']);result=run(c,proof,[values.get(i) for i in range(N)],budget=1)
        assert result.returncode==0;out=json.loads(result.stdout)
        assert out['status']=='budget_exceeded' and out['universal_unique_completion'] and 'count' not in out and 'completions' not in out
        controls.append('uniqueness_with_explicit_enumeration_budget')
        tiny=json.loads((HERE/'../round22/small_p17_basic_erase4.json').read_text())['certificate']
        incomplete=json.loads((HERE/'cover_p17_basic_erase4.json').read_text())['cover']
        result=run(tiny,incomplete,[None]*4);assert result.returncode==0
        out=json.loads(result.stdout);assert out['count']==17 and out['difference_cover']['status']=='incomplete' and not out['universal_unique_completion']
        controls.append('incomplete_cover_preserves17_completions')
        bad=deepcopy(proof);bad['roots'].pop();result=run(c,bad,digits,57)
        assert result.returncode==2;controls.append('missing_root_rejected')
    out={'status':'passed','cyclic_queries_checked':len(data['cyclic_queries']),'scalar_candidates_checked':count,
         'controls':controls,'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest()
             for n in ['consecutive33.json','cover33_refined.json','query_experiments.json','erasure_input.json','../round22/small_p17_basic_erase4.json','cover_p17_basic_erase4.json']},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['query_review.py','cover_query.py','cover_review.py','../round22/conditioned_review.py']}}
    (HERE/'query_review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
