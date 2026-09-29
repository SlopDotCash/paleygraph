#!/usr/bin/env python3
"""Show the exact coupling separation and the explicitly different reference model."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    names=('shared_insertions/contrast_results.json','shared_insertions/reference_results.json')
    data={name:json.loads((HERE/name).read_text()) for name in names}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'bold','text.color':'#263343',
                         'axes.labelcolor':'#263343','xtick.color':'#263343','ytick.color':'#263343'})
    fig,(a,b)=plt.subplots(1,2,figsize=(13,5.7));fig.subplots_adjust(left=.075,right=.97,top=.73,bottom=.24,wspace=.30)
    fig.suptitle('Round 8 · keep the insertion shared across deletion choices',x=.075,y=.97,ha='left',fontsize=19,weight='bold')
    fig.text(.075,.885,'Exact covariance distinguishes actual inputs whose complete edit marginals agree; a separate control removes the common response.',fontsize=10)
    colors=['#37879b','#b66939'];twins=data[names[0]]['cases'][:2]
    values=[F(*x['contrast_covariance_trace_squared']) for x in twins]
    bars=a.bar(range(2),[float(x) for x in values],color=colors,width=.5)
    a.bar_label(bars,padding=4,labels=[f'{v.numerator:,}/{v.denominator:,}' for v in values],fontsize=10)
    a.set_ylim(0,21000);a.set_xticks(range(2),['10','15']);a.set_xlabel('Seventh point (common points: 0, 1, 2, 3, 4, 6)')
    a.set_ylabel('Squared Frobenius norm of projected covariance')
    a.set_title('A  Equal marginal laws, different edit coupling',loc='left',pad=13,fontsize=11.5)
    a.text(.03,.94,'Paley17 · n = 7 · common deletion mode removed',transform=a.transAxes,fontsize=9.5)
    rows=data[names[1]]['comparisons']
    for family,label,color in (('arithmetic_progression','Arithmetic progression',colors[1]),('seeded_uniform','One seeded uniform set',colors[0])):
        selected=[r for r in rows if r['family']==family]
        b.plot([r['q'] for r in selected],[r['actual_over_reference_float'] for r in selected],'-o',color=color,label=label)
    b.axhline(1,color='#74859c',ls='--',lw=1);b.set_xscale('log');b.set_ylim(.5,1.42)
    b.set_xlabel('Prime order q (log scale)');b.set_ylabel('Contrast fraction / independent-sign reference')
    b.set_title('B  Compare against an explicit common-mode baseline',loc='left',pad=13,fontsize=11.5)
    b.legend(frameon=False,loc='upper left',fontsize=9)
    for ax in (a,b):ax.grid(axis='y',alpha=.13);ax.set_axisbelow(True)
    fig.text(.075,.10,'A: all one-swap values and the unlabelled per-deletion distributions agree; the per-insertion distributions differ. This toy is outside n⁴ ≤ q.\nB: reference = 25 / [25 + (n−1)(n−5)²] for independent signs. It is not a Paley law; these are eight selected inputs at n = 6, 16, 31, 50.',fontsize=9.5,color='#526173')
    for ext in ('png','pdf'):fig.savefig(HERE/f'overview.{ext}',dpi=180,facecolor='white')
    (HERE/'plot_metadata.json').write_text(json.dumps({'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
       'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in names},
       'scope':'Saved exact fractions and finite selected-input comparisons; no fitted Paley asymptotic.'},indent=2)+'\n')


if __name__=='__main__':main()
