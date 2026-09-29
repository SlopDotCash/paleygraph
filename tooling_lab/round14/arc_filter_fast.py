#!/usr/bin/env python3
"""Same transcript-bound decisions, with a shared coordinate-testing order.

The reference rebuilds and sorts all absent-block indices for each candidate.
This version stores the order once and stops its iteration as soon as the
agreement bound rejects a candidate. It does not change the mathematics.
"""
from math import comb


def filter_candidates(prepared,transcript):
    data=prepared['data'];p,n,k,d,s=(data[x] for x in ('p','n','k','dimension','s'));blocks=prepared['blocks'];received=transcript['received'];columns=prepared['columns'];origin=prepared['origin']
    base_caps=[min(len(B),d-1) for B in blocks];base=sum(base_caps)
    order=sorted(range(len(blocks)),key=lambda j:(len(blocks[j]),j))
    outputs=[];decisions=[];tests=materialized=0
    for index,row in enumerate(transcript['candidates']):
        params=row['parameters'];positive={v['block']:v for v in row['positive_blocks']};support=[];upper=base
        for j,record in positive.items():
            mask=record['support_mask'];size=mask.bit_count();assert record['count']==comb(size,d)
            upper+=size-base_caps[j];support.extend(i for t,i in enumerate(blocks[j]) if mask>>t&1)
        initial=upper;tested=[]
        for j in order:
            if upper<s:break
            if j in positive:continue
            B=blocks[j];cap=base_caps[j];matches=0
            for done,i in enumerate(B,1):
                agrees=(origin[i]+sum(a*b for a,b in zip(params,columns[i])))%p==received[i]
                tested.append([i,int(agrees)]);tests+=1
                if agrees:matches+=1;support.append(i)
                next_cap=min(base_caps[j],matches+len(B)-done);upper+=next_cap-cap;cap=next_cap
                if upper<s:break
            assert matches<=base_caps[j]
        accepted=upper>=s
        if accepted:
            assert upper==len(support)>=s
            word=[(a+sum(c*v for c,v in zip(params,column)))%p for a,column in zip(origin,columns)];materialized+=n
            support.sort();assert [i for i in range(n) if word[i]==received[i]]==support
            coefficients=[((data['origin'][j] if j<len(data['origin']) else 0)+sum(params[i]*(data['basis'][i][j] if j<len(data['basis'][i]) else 0) for i in range(d)))%p for j in range(k)]
            outputs.append({'parameters':params,'coefficients':coefficients,'codeword':word,'agreement_support':support,'first_query_index':row['first_query_index']})
        decisions.append({'candidate_index':index,'initial_agreement_upper':initial,'final_agreement_upper':upper,'tested_coordinates':tested,'accepted':accepted})
    return {'schema':'arc_transcript_pruning_v1','outputs':outputs,'decisions':decisions,
            'candidate_coordinate_evaluations':tests,'output_word_coordinate_evaluations':materialized,
            'total_explicit_coordinate_evaluations':tests+materialized,
            'full_scan_coordinate_evaluations':len(transcript['candidates'])*n,
            'scope':'All query labels, including absence from each block, give deterministic agreement bounds. Filtering is exact inside the supplied space, conditional on a checked cover and complete query transcript.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    from arc_pruning import prepare,collect
    if len(sys.argv)!=3:raise SystemExit('Usage: arc_filter_fast.py certificate.json received.json')
    prepared=prepare(json.loads(Path(sys.argv[1]).read_text()));transcript=collect(prepared,json.loads(Path(sys.argv[2]).read_text()))
    print(json.dumps({'transcript':transcript,'filter':filter_candidates(prepared,transcript)},separators=(',',':')))
