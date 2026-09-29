#!/usr/bin/env python3
from hashlib import sha256
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    sizes=[9841,1328602,3690562];nodes=[131,1337,1867]
    orbit=json.loads((HERE/'orbit_large.json').read_text())['output'];counts=orbit['survivors_after_coordinates'][:20]
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained');x=np.arange(3)
    axes[0].bar(x-.18,sizes,.35,color='#697989',label='Difference directions, up to sign')
    axes[0].bar(x+.18,nodes,.35,color='#16857a',label='Verified tree nodes')
    axes[0].set_yscale('log');axes[0].set_ylim(50,2e7);axes[0].set_xticks(x,['32 erased','33 erased','36 erased'])
    axes[0].set_title('Complete unique-completion certificates');axes[0].set_ylabel('Count (log scale)')
    axes[0].legend(loc='upper left',fontsize=8,frameon=False)
    for i,(s,n) in enumerate(zip(sizes,nodes)):
        axes[0].text(i-.18,s*1.18,f'{s:,}',ha='center',fontsize=8)
        axes[0].text(i+.18,n*1.18,f'{n:,}',ha='center',fontsize=8)
    axes[1].plot(range(1,21),counts,'o-',color='#3479a8',markersize=4)
    axes[1].set_yscale('symlog',linthresh=1);axes[1].set_ylim(-.1,2e8)
    axes[1].set_xticks([1,5,10,15,20]);axes[1].set_xlabel('Fixed-coordinate constraints applied')
    axes[1].set_ylabel('Surviving scalar candidates');axes[1].set_title('A 40-erasure lattice candidate has no actual pair')
    axes[1].annotate('27,904,523',(1,counts[0]),xytext=(6,8),textcoords='offset points',fontsize=9)
    axes[1].annotate('0',(20,0),xytext=(-8,12),textcoords='offset points',fontsize=12)
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15)
    fig.suptitle('Exact finite tooling at p = 2,013,265,921, N = 64',fontsize=15)
    fig.savefig(HERE/'results.png',dpi=180);fig.savefig(HERE/'results.pdf')
    out={'status':'produced','complete_difference_domain_sizes':sizes,'tree_nodes':nodes,'orbit_survivors':counts,
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['cover_review.json','refinement33_review.json','cover36_review.json','orbit_large.json','orbit_review.json']},
         'source_sha256':{'plot_results.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'plot_metadata.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
