#!/usr/bin/env python3
"""Independent query/word review and coefficient-based lifting of all branches."""
from hashlib import sha256
import json
from pathlib import Path
from conditional_verifier import validate_conditional,evaluate
from arc_review import review_run

HERE=Path(__file__).resolve().parent


def review_decoding(c,result,replay=True):
    assert result['schema']=='conditional_arc_pruning_v1' and result['certificate_review']==validate_conditional(c)
    data=c['plan']['input'];p,n,k=(data[x] for x in ('p','n','k'));received=c['plan']['received'];outside=c['plan']['outside_coordinates'];w=c['plan']['projective_direction'];pivot=c['plan']['pivot']
    assert len(result['branches'])==len(c['parts'])
    reviewed=[];lifts=[];outputs=[];whole=0;seen=set()
    for index,(part,branch) in enumerate(zip(c['parts'],result['branches'])):
        assert branch['transcript']['received']==[received[i] for i in outside]
        if replay:reviewed.append(review_run(part['certificate'],branch))
        for j,row in enumerate(branch['filter']['outputs']):
            if part['kind']=='light':
                params=row['parameters'][:];value=sum(params[t]*w[t] for t in range(3))%p
                if value in part['excluded_values']:
                    lifts.append({'branch':index,'branch_output':j,'parameters':params,'status':'excluded_heavy_value','equation_value':value});continue
            else:
                # Solve the single affine relation using canonical free-index
                # order derived here, rather than trusting producer lift data.
                free=[t for t in range(3) if t!=pivot];params=[0]*3
                for t,a in zip(free,row['parameters']):params[t]=a
                value=part['value'];params[pivot]=(value-sum(params[t]*w[t] for t in free))*pow(w[pivot],-1,p)%p
            assert tuple(params) not in seen;seen.add(tuple(params))
            coefficients=[0]*k
            for t,a in enumerate(data['origin']):coefficients[t]=a
            for a,b in zip(params,data['basis']):
                for t,v in enumerate(b):coefficients[t]=(coefficients[t]+a*v)%p
            word=[evaluate(coefficients,x,p) for x in data['domain']];whole+=n
            assert [word[i] for i in outside]==row['codeword']
            support=[i for i in range(n) if word[i]==received[i]];accepted=len(support)>=data['s']
            lifts.append({'branch':index,'branch_output':j,'parameters':params,'status':'accepted' if accepted else 'insufficient_total_agreement','equation_value':value,'agreement_count':len(support)})
            if accepted:outputs.append({'parameters':params,'coefficients':coefficients,'codeword':word,'agreement_support':support})
    outputs.sort(key=lambda r:r['parameters']);assert outputs==result['outputs'] and lifts==result['lifts']
    b=result['branches'];expected={'query_systems_solved':sum(r['transcript']['queries_executed'] for r in b),
        'branch_candidates_generated':sum(len(r['transcript']['candidates']) for r in b),
        'branch_candidate_coordinate_evaluations':sum(r['filter']['candidate_coordinate_evaluations'] for r in b),
        'branch_output_word_coordinate_evaluations':sum(r['filter']['output_word_coordinate_evaluations'] for r in b),
        'final_whole_word_coordinate_evaluations':whole,'branch_full_scan_coordinate_evaluations':sum(r['filter']['full_scan_coordinate_evaluations'] for r in b)}
    expected['total_explicit_coordinate_evaluations']=expected['branch_candidate_coordinate_evaluations']+expected['branch_output_word_coordinate_evaluations']+whole
    assert expected==result['metrics']
    return {'status':'passed','branches':reviewed,'whole_word_lifts_checked':len(seen),'whole_word_coordinates_checked':whole,'outputs_checked':len(outputs),'query_replay_performed':replay}


def main():
    source=HERE/'conditional_decoding.json';produced=json.loads(source.read_text());records=[];bindings={source.name:sha256(source.read_bytes()).hexdigest()}
    for row in produced['cases']:
        name=row['case'];cpath=HERE/f'{name}.conditional.json';rpath=HERE/row['file'];oldpath=HERE.parent/'round15'/f'{name}.run.json'
        c=json.loads(cpath.read_text());result=json.loads(rpath.read_text());old=json.loads(oldpath.read_text())
        checked=review_decoding(c,result)
        assert c['plan']['received']==old['transcript']['received']
        baseline=[{k:o[k] for k in ('parameters','coefficients','codeword','agreement_support')} for o in old['filter']['outputs']]
        assert result['outputs']==baseline
        record={'case':name,'review':checked,'metrics':result['metrics'],
                'baseline_queries':old['transcript']['queries_executed'],'baseline_coordinate_evaluations':old['filter']['total_explicit_coordinate_evaluations'],
                'baseline_complete_supplied_space_list_matched':True};records.append(record)
        for path in [cpath,rpath,oldpath]:bindings[str(path.relative_to(HERE)) if path.parent==HERE else f'../round15/{path.name}']=sha256(path.read_bytes()).hexdigest()
        print(json.dumps(record),flush=True)
    out={'status':'passed','scope':'Every branch query independently solved by integer Cramer determinants and every branch candidate word fully evaluated. All affine lifts and final lists checked, including comparison with frozen complete supplied-space baselines.',
         'cases':records,'queries_independently_replayed':sum(b['queries_independently_replayed'] for r in records for b in r['review']['branches']),
         'candidate_words_fully_evaluated':sum(b['candidate_words_fully_evaluated'] for r in records for b in r['review']['branches']),
         'candidate_coordinate_values_checked':sum(b['full_word_coordinate_values_checked'] for r in records for b in r['review']['branches']),
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['decoding_review.py','conditional_verifier.py','../round14/arc_review.py','../round14/arc_verifier.py']},'input_sha256':bindings}
    (HERE/'decoding_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
