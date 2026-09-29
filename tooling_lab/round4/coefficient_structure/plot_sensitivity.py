#!/usr/bin/env python3
"""Scientific plot of exact coefficient signs and fixed-size scaling."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from closed_coefficient import contraction_coefficient_degree6
HERE=Path(__file__).resolve().parent

def main():
    data=json.loads((HERE/'all_n_validation.json').read_text());cases=[r for r in data['cases'] if r['q'] in [29,49,101,1297]]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(11,8),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.1]);top=fig.add_subplot(gs[0,:]);left=fig.add_subplot(gs[1,0]);right=fig.add_subplot(gs[1,1])
    colors={-1:'#c85a43',0:'#c6cbd0',1:'#347db2'}
    for row,case in enumerate(cases):
        q=case['q']
        for interval in case['sign_intervals']:
            a=(interval['from_n']-.5)/q;b=min(1,(interval['through_n']+.5)/q)
            top.barh(row,b-a,left=a,height=.68,color=colors[interval['sign']],linewidth=0)
    top.set(yticks=range(len(cases)),yticklabels=[f'q = {c["q"]}' for c in cases],xlim=(0,1),xlabel='Set size n/q',title='Exact coefficient sign: multiple reversals, not a monotone sensitivity')
    top.legend(handles=[Patch(color=colors[s],label=label) for s,label in [(-1,'negative'),(1,'positive'),(0,'zero')]],loc='upper center',bbox_to_anchor=(.5,1.01),ncols=3,frameon=False,fontsize=8)
    top.grid(axis='x',alpha=.15)
    ns=np.arange(6,30);cs=[float(contraction_coefficient_degree6(29,int(n))) for n in ns]
    left.plot(ns,cs,marker='o',ms=3,color='#344d66');left.axhline(0,color='#777',lw=.8)
    left.scatter([7,8],[cs[1],cs[2]],s=70,facecolors='none',edgecolors='#c85a43')
    left.set(xlabel='n, with q = 29',ylabel='c(29,n,6)',title='The verified ablation changes sign at n7 → n8')
    left.grid(alpha=.15)
    qs=[1297,2017,4001,10009,65537,1000033];ratios=[-q**3*float(contraction_coefficient_degree6(q,31))/(28*27*26*25) for q in qs]
    right.semilogx(qs,ratios,marker='o',color='#347db2');right.axhline(1,color='#777',ls='--',lw=.8)
    right.set(xlabel='q (log scale); n fixed at 31',ylabel='−q³ c(q,31,6) / (28)₄',title='Fixed n: the normalized coefficient tends to 1')
    right.grid(alpha=.15)
    fig.suptitle('Set-size sensitivity of the three-mark contraction coefficient',fontsize=15)
    fig.supxlabel('c is the coefficient of Q = h₂ᵀ S h₃ in a conditional second moment. Exact signs; numerical values only for display.\nThese are conditional-average diagnostics, not bounds on exceptional sets.',fontsize=9)
    fig.savefig(HERE/'coefficient_sensitivity.png',dpi=180);fig.savefig(HERE/'coefficient_sensitivity.svg')
    (HERE/'plot_metadata.json').write_text(json.dumps({'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'numpy':np.__version__,'matplotlib':matplotlib.__version__,'python':sys.executable,'fixed_n_ratio_samples':list(zip(qs,ratios))},indent=2)+'\n')
if __name__=='__main__':main()
