#!/usr/bin/env python3
"""Scientific figure generated only from saved exact results."""
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    collision=json.loads((HERE/'rank_three_review.json').read_text());scale=json.loads((HERE/'rank_three_scale.json').read_text())
    actual=json.loads((HERE/'actual_results.json').read_text());colors=['#9a3b42','#216a83']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8.4),layout='constrained')
    ax=axes[0,0]
    for row,color,label in zip(collision['rows'],colors,['span(1, X, X⁶+X²)','span(1, X, X⁴+X²)']):
        h=dict(row['agreement_histogram']);xs=range(146,164);ys=[h.get(x,0)/4368 for x in xs]
        ax.plot(xs,ys,marker='o',markersize=3,color=color,label=label)
        ax.axvline(row['minimum'],color=color,ls=':',alpha=.7)
    ax.set(title='Same Tutte polynomial, different worst pinning',xlabel='Informative triples inside an eleven-set',ylabel='Fraction of all 4,368 eleven-sets')
    ax.legend(fontsize=9);ax.text(.03,.94,'Equal means and variances; minima 147 and 148',transform=ax.transAxes,va='top',fontsize=9)
    ax=axes[0,1]
    h0=dict(collision['rows'][0]['agreement_histogram']);h1=dict(collision['rows'][1]['agreement_histogram']);xs=np.arange(147,164)
    delta=np.array([h0.get(int(x),0)-h1.get(int(x),0) for x in xs]);ax.bar(xs,delta,color=[colors[0] if v>=0 else colors[1] for v in delta],width=.75)
    ax.axhline(0,color='#666666',lw=.7);ax.set(title='What the global summaries discard',xlabel='Informative triples',ylabel='Count difference: sextic − quartic')
    ax=axes[1,0]
    for offset,label,color in [(-.10,'quartic',colors[1]),(.10,'sextic',colors[0])]:
        rows=[r for r in scale['rows'] if label in r['name']]
        for j,r in enumerate(rows):
            low,high=np.array(r['interval'])/comb(r['n'],3);x=j+offset
            ax.plot([x,x],[low,high],lw=3,color=color);ax.plot(x,low,'o',ms=4,color=color);ax.plot(x,high,'_',ms=8,color=color)
        ax.plot([],[],color=color,lw=3,label=label)
    rows=scale['rows'][::2];ax.plot(range(5),[r['universal_degree_bound']/comb(r['n'],3) for r in rows],color='#777777',marker='x',ls='--',label='Degree-bound baseline')
    ax.set_xticks(range(5),[str(r['n']) for r in rows]);ax.set(title='Rank-three query: exact points and intervals',xlabel='Number of coordinates n',ylabel='Worst-case success probability')
    ax.legend(fontsize=9);ax.set_ylim(0,.38)
    ax=axes[1,1];labels=['F17, n=16','F65537, n=1,024'];values=[];baseline=[]
    for key in ('small','large'):
        r=actual[key];values.append(r['result']['minimum_injective_pairs']/comb(r['n'],2));baseline.append(comb(r['s']-r['k']+2,2)/comb(r['n'],2))
    x=np.arange(2);ax.bar(x-.17,baseline,.34,color='#bbbbbb',label='Degree-bound baseline');bars=ax.bar(x+.17,values,.34,color='#216a83',label='Exact declared-space probability')
    for b,v in zip(bars,values):ax.text(b.get_x()+b.get_width()/2,v+.005,f'{v:.4f}',ha='center',fontsize=9)
    ax.set_xticks(x,labels);ax.set(title='Rank-two tool on saved actual clusters',ylabel='Worst-case pair success probability',ylim=(0,.20));ax.legend(fontsize=9,loc='upper left')
    fig.suptitle('Round 11 · joint incidence changes an adversarial query',fontsize=16)
    fig.supxlabel('Declared polynomial spaces only. No general decoder, prize proof or historical originality certification.',fontsize=10)
    for extension in ('png','pdf'):fig.savefig(HERE/f'overview.{extension}',dpi=180)
    plt.close(fig)
    out={'status':'passed','source_sha256':{'plot_results.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('rank_three_review.json','rank_three_scale.json','actual_results.json')},
         'figure_sha256':{f'overview.{extension}':sha256((HERE/f'overview.{extension}').read_bytes()).hexdigest() for extension in ('png','pdf')}}
    (HERE/'plot_metadata.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
