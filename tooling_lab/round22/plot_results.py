#!/usr/bin/env python3
from hashlib import sha256
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    initial=json.loads((HERE/'conditioned_summary.json').read_text())['cases'][0]
    metric=json.loads((HERE/'metric_summary.json').read_text())['cases'][0]
    joint=json.loads((HERE/'coupled_review.json').read_text())
    caps=[initial['baseline_cap'],initial['universal_cap'],metric['universal_cap']]
    counts=[joint['sign_representatives'],joint['remaining_after_pair_cuts'],joint['unresolved_targets']+joint['continuous_feasible_targets']]
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    colors=['#687789','#3479a8','#16857a']
    axes[0].bar(range(3),caps,color=colors,width=.58)
    axes[0].set_yscale('log'); axes[0].set_ylim(100,120000)
    axes[0].set_xticks(range(3),['Centered\nbounds','Known-digit\ncuts','Changed metric\nand cuts'])
    axes[0].set_ylabel('Uniform candidate combinations (log scale)')
    axes[0].set_title('Complete decoder: 32 consecutive erasures')
    for x,y in enumerate(caps): axes[0].text(x,y*1.14,f'{y:,}',ha='center',fontsize=12)
    axes[1].plot(range(3),counts,'o-',color='#16857a',linewidth=2.5,markersize=8)
    axes[1].set_yscale('symlog',linthresh=1); axes[1].set_ylim(-.1,40000)
    axes[1].set_yticks([0,1,10,100,1000,10000],['0','1','10','100','1,000','10,000'])
    axes[1].set_xticks(range(3),['Individual\nbounds','156 pair\ncuts','86 targeted\nseparators'])
    axes[1].set_ylabel('Nonzero difference directions, up to sign')
    axes[1].set_title('Universal uniqueness: all directions excluded')
    for x,y in enumerate(counts): axes[1].annotate(f'{y:,}',(x,y),xytext=(0,12),textcoords='offset points',ha='center',fontsize=12)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.15)
    fig.suptitle('Exact finite certificates at p = 2,013,265,921, N = 64',fontsize=15)
    fig.savefig(HERE/'results.png',dpi=180); fig.savefig(HERE/'results.pdf')
    out={'status':'produced','candidate_caps':caps,'difference_sign_representatives':counts,
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['conditioned_summary.json','metric_summary.json','coupled_review.json']},
         'source_sha256':{'plot_results.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'plot_metadata.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__': main()
