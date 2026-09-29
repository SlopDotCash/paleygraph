#!/usr/bin/env python3
"""Scientific figure generated only from completed exact round7 outputs."""
from hashlib import sha256
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    names=('local_edits/ablation_results.json','local_edits/threshold_results.json')
    data={name:json.loads((HERE/name).read_text()) for name in names}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'bold','text.color':'#263343',
                         'axes.labelcolor':'#263343','xtick.color':'#263343','ytick.color':'#263343'})
    fig,(a,b)=plt.subplots(1,2,figsize=(13,5.7))
    fig.subplots_adjust(left=.065,right=.975,top=.73,bottom=.24,wspace=.27)
    fig.suptitle('Round 7 · exact local fluctuations retain the chosen input',x=.065,y=.97,ha='left',fontsize=19,weight='bold')
    fig.text(.065,.885,'New arithmetic implementation of established conditional-variance tools; independently checked through q = 6,700,417.',fontsize=10.5)
    colors=['#37879b','#b66939','#74859c'];rows=data[names[0]]['cases'];x=np.arange(3)
    for j,label in enumerate(('Signed row histogram','Also unordered quartic deck','Also deletion summaries')):
        bars=a.bar(x+(j-1)*.24,[r['models'][j]['ambiguous_variance_fibres'] for r in rows],.24,color=colors[j],label=label)
        a.bar_label(bars,padding=3,fontsize=9)
    a.set_ylim(0,75);a.set_xticks(x,[f'n = {r["n"]}' for r in rows]);a.set_ylabel('Fibres with different local variances')
    a.set_title('A  Available summaries can still lose information',loc='left',pad=13,fontsize=11.5)
    a.legend(frameon=False,fontsize=9,loc='upper left')
    bounds=data[names[1]]['scale_threshold_bounds']
    for label,name,color in (('Arithmetic progression','arithmetic_progression',colors[1]),('One seeded uniform set','seeded_uniform',colors[0])):
        chosen=[r for r in bounds if r['family']==name]
        b.plot([r['q'] for r in chosen],[100*r['half_magnitude_persistence_lower_float'] for r in chosen],'-o',color=color,label=label)
    b.set_xscale('log');b.set_ylim(-3,108);b.set_xlabel('Prime order q (log scale)')
    b.set_ylabel('Neighbours retaining sign and >½ magnitude (%)')
    b.set_title('B  One-step persistence on selected critical-size inputs',loc='left',pad=13,fontsize=11.5)
    b.legend(frameon=False,fontsize=9,loc='upper left')
    for ax in (a,b):ax.grid(axis='y',alpha=.13);ax.set_axisbelow(True)
    fig.text(.065,.105,'A: complete normalized Paley17 scans; the six-point shortcut fails at n = 7 and 8. Zero observed ambiguity is not a general sufficiency theorem.\nB: exact local moments plus classical Cantelli; n = 6, 16, 31, 50. These bounds concern individual one-swap neighbourhoods, not all sets or multiple steps.',fontsize=9.5,color='#526173')
    for ext in ('png','pdf'):fig.savefig(HERE/f'overview.{ext}',dpi=180,facecolor='white')
    (HERE/'plot_metadata.json').write_text(json.dumps({'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in names},
        'scope':'Exact saved counts and rational-bound display values; no interpolation or uncertainty bands claimed.'},indent=2)+'\n')


if __name__=='__main__':main()
