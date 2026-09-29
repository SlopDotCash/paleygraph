#!/usr/bin/env python3
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from erasure_query import query
from conditioned_review import certify, replay

HERE = Path(__file__).resolve().parent


def main():
    source = HERE/'metric_consecutive32.json'; data = json.loads(source.read_text()); c = data['certificate']
    prepared = certify(c); records = []; anchor = 1234567; word = prepared[2]([anchor])[0]
    for start in range(c['N']):
        erased = {(start+j) % c['N'] for j in range(len(c['erased']))}
        digits = [None if j in erased else v for j,v in enumerate(word)]
        out = query(c, digits, start=start)
        rotated = {j:digits[(j+start) % c['N']]*(1 if j+start < c['N'] else -1) for j in range(len(c['erased']), c['N'])}
        expected = replay(c, rotated, 50000, prepared)
        scalars = sorted(x['scalar']*pow(c['g'], -start, c['p']) % c['p'] for x in expected['completions'])
        expected['completions'] = [{'scalar':a,'digits':w} for a,w in zip(scalars,prepared[2](scalars))]
        assert {k:out[k] for k in expected} == expected
        assert [v['scalar'] for v in out['completions']] == [anchor]
        records.append({'cyclic_start':start, 'input':digits, 'output':out})
    sample = records[57]
    (HERE/'erasure_input.json').write_text(json.dumps(sample['input'], indent=2)+'\n')
    (HERE/'erasure_output.json').write_text(json.dumps(sample['output'], indent=2)+'\n')
    controls = []
    with TemporaryDirectory(prefix='round22-query-') as temp:
        temp = Path(temp)
        def run(cert, digits, budget=50000, start=None):
            (temp/'cert.json').write_text(json.dumps(cert)); (temp/'input.json').write_text(json.dumps(digits))
            cmd=[sys.executable,str(HERE/'erasure_query.py'),'--certificate',str(temp/'cert.json'),'--input',str(temp/'input.json'),'--candidate-budget',str(budget)]
            if start is not None: cmd += ['--cyclic-start', str(start)]
            return subprocess.run(cmd,capture_output=True,text=True,timeout=45)
        p=run(c,sample['input'],start=57); assert p.returncode == 0 and json.loads(p.stdout) == sample['output']
        controls.append({'name':'cyclic32_cli', 'status':'passed'})
        large=json.loads((HERE/'conditioned_random32.json').read_text()); known=dict(large['records'][0]['known'])
        digits=[known.get(j) for j in range(64)]; p=run(large['certificate'],digits)
        out=json.loads(p.stdout); assert p.returncode == 0 and out['status']=='budget_exceeded' and 'count' not in out and 'completions' not in out
        controls.append({'name':'incomplete_has_no_false_list', 'status':'passed'})
        tiny=json.loads((HERE/'small_p17_basic_erase4.json').read_text()); p=run(tiny['certificate'],[None]*4)
        assert p.returncode == 0 and json.loads(p.stdout)['count']==17
        controls.append({'name':'all17_completions', 'status':'passed'})
        for label in ['wrong_length','boolean_digit','wrong_mask','zero_budget','tampered_lambda']:
            cert=deepcopy(c); digits=sample['input'][:]; budget=50000
            if label=='wrong_length': digits.pop()
            if label=='boolean_digit': digits[next(i for i,v in enumerate(digits) if v is not None)]=True
            if label=='wrong_mask': digits[next(i for i,v in enumerate(digits) if v is None)]=0
            if label=='zero_budget': budget=0
            if label=='tampered_lambda': cert['conditioning']['cuts'][0]['lambda'][0][0]+=1
            p=run(cert,digits,budget,57); assert p.returncode==2
            controls.append({'name':label,'status':'rejected'})
    # Exhibit why shrinking around the old center would exclude a true word.
    regression=None
    regression_source=HERE/'small_p41_general_erase1.json'
    small=json.loads(regression_source.read_text()); c=small['certificate']; prepared=certify(c)
    for anchor in range(c['p']):
        actual=prepared[2]([anchor])[0]
        known={i:actual[i] for i in c['conditioning']['known_coordinates']}
        rhs=-sum(c['projection_weights'][i]*y for i,y in known.items()) % c['k']
        x0=[rhs//c['kernel_gcd']*v for v in c['bezout_lift']]
        for j,cut in enumerate(c['conditioning']['cuts']):
            oldcenter=-sum(F(*c['inverse'][i][j])*x0[i] for i in range(len(c['erased'])))
            z=sum(F(*c['inverse'][i][j])*(actual[u]-x0[i]) for i,u in enumerate(c['erased']))
            assert z.denominator==1
            if abs(z-oldcenter)>F(*cut['radius']):
                shift=sum(F(*v)*known[i] for i,v in zip(c['conditioning']['known_coordinates'],cut['lambda']))
                assert abs(z-oldcenter-shift)<=F(*cut['radius'])
                regression={'source':regression_source.name,'anchor_scalar':anchor,'coordinate':j,'known':list(map(list,known.items())),
                            'integer_coordinate':int(z),'old_center':[oldcenter.numerator,oldcenter.denominator],
                            'center_shift':[shift.numerator,shift.denominator],'radius':cut['radius']}
                break
        if regression: break
    assert regression
    out={'status':'passed','cyclic_controls':records,'controls':controls,'omitted_center_regression':regression,
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['metric_consecutive32.json','conditioned_random32.json','small_p17_basic_erase4.json','small_p41_general_erase1.json']},
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['query_review.py','erasure_query.py','conditioned_erasure.py','conditioned_review.py','metric_review.py']}}
    (HERE/'query_review.json').write_text(json.dumps(out,separators=(',', ':'))+'\n')
    print(json.dumps({'cyclic_controls':len(records),'cli_controls':len(controls),'center_regression_scalar':regression['anchor_scalar']}),flush=True)


if __name__=='__main__': main()
