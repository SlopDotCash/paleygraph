#!/usr/bin/env python3
"""Scientific summary of checked query counts, not an end-to-end benchmark."""
from hashlib import sha256
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    family=json.loads((HERE/'family_results.json').read_text())
    actual=json.loads((HERE/'actual_results.json').read_text())
    previous=json.loads((HERE.parent/'round12/pruning_results.json').read_text())
    assert family['spaces_covered']==2894 and [r['coverage_review']['queries_rank_checked'] for r in actual['cases']]==[10,2070]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    xs,ys=zip(*family['query_count_histogram']);axes[0].bar(xs,ys,color='#285E8E')
    axes[0].set(xticks=xs,xlabel='Queries in checked cover',ylabel='Actual rank-three spaces',title='All 2,894 small spaces')
    small_random=previous['small_runs'][0]['run']['sampling_budget']['trials'];assert small_random==432
    axes[1].bar(['Sampled\n2⁻⁴⁶ failure bound','Deterministic'],[small_random,10],color=['#AC7041','#285E8E'])
    axes[1].set(ylabel='Queries per received word',title='Hard small space')
    for i,v in enumerate([small_random,10]):axes[1].text(i,v+8,str(v),ha='center')
    vals=[previous['large_universal_baseline_trial_budget']['trials'],previous['large_certificate_trial_budget']['trials'],2070]
    axes[2].bar(['Degree\nbaseline','Instance\nbound','Deterministic'],vals,color=['#C5A484','#AC7041','#285E8E'])
    axes[2].set(ylabel='Queries per received word',title='Length 1,024; sampled bounds 2⁻⁴²')
    for i,v in enumerate(vals):axes[2].text(i,v+30,str(v),ha='center')
    axes[2].set_ylim(0,2350);axes[1].set_ylim(0,480)
    fig.suptitle('Query counts exclude construction and verification costs',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(HERE/f'overview.{ext}',dpi=170)
    out={'status':'passed','scope':'Query counts only, with different guarantee types shown explicitly. These bars are not end-to-end runtime comparisons.',
         'source_sha256':{'plot_results.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['family_results.json','actual_results.json','../round12/pruning_results.json']},
         'artifact_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['overview.png','overview.pdf']}}
    (HERE/'plot_metadata.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
