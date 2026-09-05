#!/usr/bin/env python3
"""Publication-style summary of saved exact round5 evidence."""
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    inputs={name:json.loads((HERE/name).read_text()) for name in
            ('third_moment/results.json','trace_invariants/residue_results.json','piece_discovery/results.json')}
    third=inputs['third_moment/results.json'];pieces=inputs['piece_discovery/results.json']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'bold','axes.labelcolor':'#263343',
                         'text.color':'#263343','xtick.color':'#263343','ytick.color':'#263343'})
    fig,axs=plt.subplots(2,2,figsize=(13.8,8.5))
    fig.subplots_adjust(left=.08,right=.97,top=.87,bottom=.12,wspace=.29,hspace=.52)
    fig.suptitle('Round 5 · exact compression, hidden-piece discovery, and limits',x=.08,y=.975,ha='left',fontsize=19,weight='bold')
    fig.text(.08,.924,'Independent finite checks; global moments and ordinary codeword lists, without prize-proof claims.',fontsize=11)
    colors=['#74859c','#37879b','#b66939']
    ax=axs[0,0];x=np.arange(2);w=.24
    rows=[next(r for r in third['scale_cases'] if r['q']==q) for q in (65537,1000033)]
    refined=[next(r for r in inputs['trace_invariants/residue_results.json']['comparisons']
                  if (r['q'],r['n'])==(q,n)) for q,n in ((65537,16),(1000033,31))]
    for i,(name,values) in enumerate((('Direct inventory',[r['direct_inventory_coefficient_evaluations'] for r in rows]),
                                     ('Four families',[r['coefficient_evaluations'] for r in rows]),
                                     ('Seven statistics',[r['coefficient_evaluations'] for r in refined]))):
        bars=ax.bar(x+(i-1)*w,values,w,label=name,color=colors[i])
        ax.bar_label(bars,padding=3,fontsize=10)
    ax.set_yscale('log');ax.set_ylim(8,4200);ax.set_xticks(x,['q = 65,537\nn = 16','q = 1,000,033\nn = 31'])
    ax.set_ylabel('Coefficient evaluations (log scale)');ax.set_title('A  Constant algebra budget',loc='left',pad=12)
    ax.legend(frameon=False,fontsize=9,ncol=3,loc='upper left',bbox_to_anchor=(-.02,1.02))

    x=np.arange(3);w=.32
    for panel,key,title,ylabel in ((axs[0,1],'standardized_third_float','B  Extra moment distinguishes actual twins','Standardized third moment'),
                                    (axs[1,0],'maximum','C  That ordering does not order maxima','Maximum of T₆ over all n-column sets')):
        for i,graph in enumerate(('Paley49','Peisert49')):
            values=[next(r[key] for r in third['twins'] if r['n']==n and r['graph']==graph) for n in (6,7,8)]
            bars=panel.bar(x+(i-.5)*w,values,w,label=graph,color=colors[i+1])
            if key=='maximum':panel.bar_label(bars,padding=3,fontsize=10)
        panel.set_xticks(x,['n = 6','n = 7','n = 8']);panel.set_ylabel(ylabel)
        panel.set_title(title,loc='left',pad=12);panel.axhline(0,color='#9aabba',lw=.7)
        panel.legend(frameon=False,fontsize=10,loc='lower left' if key=='standardized_third_float' else 'upper left')
        if key=='maximum':panel.set_ylim(0,162)

    ax=axs[1,1]
    large=pieces['large'];three=next(r for r in large if r['result']['minimum_polynomial_cover_size_certified']==3)
    two=next(r for r in large if r['result']['minimum_polynomial_cover_size_certified']==2)
    caps=[three['chunk_baseline']['root_count_cap'],three['result']['static_certificate']['outside_piece_agreement_cap'],
          two['result']['static_certificate']['outside_piece_agreement_cap']]
    threshold=three['result']['s'];assert threshold==two['result']['s']
    bars=ax.bar(range(3),caps,color=colors,width=.62);ax.bar_label(bars,padding=3,fontsize=11)
    ax.axhline(threshold,color='#303a49',ls='--',lw=1.2)
    ax.text(2.45,threshold+22,f'Agreement threshold {threshold}',ha='right',fontsize=10)
    ax.set_xticks(range(3),['Chunk cover\n16 pieces','Blind discovery\n3 pieces','Blind discovery\n2 pieces'])
    ax.set_ylim(0,1150);ax.set_ylabel('Agreement cap for an unlisted polynomial')
    ax.set_title('D  Discovery makes a list certificate possible',loc='left',pad=12)
    fig.text(.08,.035,'D uses structured interleaved inputs with n = 1024, k = 64. All 16 saved hard coset stacks remain unresolved by discovery.\nA excludes inventory construction. B–C use complete order49 censuses through n = 8 (450,978,066 subsets per graph).',fontsize=10,color='#526173')
    for ax in axs.flat:
        ax.grid(axis='y',alpha=.13);ax.set_axisbelow(True)
    for ext in ('png','pdf'):fig.savefig(HERE/f'overview.{ext}',dpi=180,facecolor='white')
    meta={'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
          'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in inputs},
          'scope':'Saved exact evidence displayed; floating values only for standardized moments and plotting.'}
    (HERE/'plot_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')


if __name__=='__main__':main()
