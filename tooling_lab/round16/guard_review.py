#!/usr/bin/env python3
"""Complete independent replay of singleton-reserving conditional decodings."""
from hashlib import sha256
import json
from pathlib import Path
from decoding_review import review_decoding

HERE=Path(__file__).resolve().parent


def main():
    source=HERE/'guard_results.json';results=json.loads(source.read_text());rows=[];bindings={source.name:sha256(source.read_bytes()).hexdigest()}
    for row in results['cases']:
        name=row['case'];paths=[HERE/row['certificate'],HERE/row['decoded'],HERE/f'{name}.decoded.json',HERE.parent/'round15'/f'{name}.run.json']
        c,result,conditional,baseline=[json.loads(p.read_text()) for p in paths]
        checked=review_decoding(c,result)
        assert c['plan']['received']==baseline['transcript']['received']
        assert result['outputs']==conditional['outputs']==[{k:r[k] for k in ('parameters','coefficients','codeword','agreement_support')} for r in baseline['filter']['outputs']]
        assert all(any(len(B)==1 for B in part['certificate']['blocks']) for part in c['parts'])
        record={'case':name,'review':checked,'metrics':result['metrics'],'first_conditional_metrics':conditional['metrics'],
                'baseline_queries':baseline['transcript']['queries_executed'],'baseline_coordinate_evaluations':baseline['filter']['total_explicit_coordinate_evaluations'],
                'both_complete_supplied_space_lists_matched':True};rows.append(record);print(json.dumps(record),flush=True)
        for path in paths:bindings[path.name if path.parent==HERE else f'../round15/{path.name}']=sha256(path.read_bytes()).hexdigest()
    out={'status':'passed','scope':'Every guarded branch query independently replayed and every candidate fully evaluated; lifts and final lists match both checked supplied-space baselines. Runtime and globally optimal test ordering are not claimed.',
         'cases':rows,'queries_independently_replayed':sum(b['queries_independently_replayed'] for r in rows for b in r['review']['branches']),
         'candidate_words_fully_evaluated':sum(b['candidate_words_fully_evaluated'] for r in rows for b in r['review']['branches']),
         'candidate_coordinate_values_checked':sum(b['full_word_coordinate_values_checked'] for r in rows for b in r['review']['branches']),
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['guard_review.py','decoding_review.py','conditional_verifier.py','../round14/arc_review.py','../round14/arc_verifier.py']},'input_sha256':bindings}
    (HERE/'guard_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
