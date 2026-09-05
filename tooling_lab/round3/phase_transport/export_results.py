#!/usr/bin/env python3
"""Export concise comparisons and a scientific figure from saved outputs.
Does not rerun the numerical experiments. Plot runtime is recorded separately.
"""
import csv,hashlib,json,math,platform,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    pilot=json.loads((HERE/'results.json').read_text())
    revised=json.loads((HERE/'refinement.json').read_text())
    validation=json.loads((HERE/'validation.json').read_text())
    rows=[]
    for p,r in zip(pilot['cases'],revised['cases']):
        assert (p['p'],p['side'])==(r['p'],r['side'])
        b=p['baseline'];s=r['baseline'];nulls=[t['off_axes_quartic_sum'] for t in r['S3_equivariant_shuffles']]
        row={'p':p['p'],'side':p['side'],'rayleigh_H':p['input_rayleigh_H'],
             'ordinary_overlap_max':b['ordinary_translation_max'],
             'full_phase_overlap_max':b['all_nonidentity_max'],
             'full_peak_a':b['peak_location'][0],'full_peak_b':b['peak_location'][1],
             'vertical_line_energy':s['vertical_line'],
             'largest_line_is_vertical':b['largest_line_slope']==p['p'],
             'off_axes_quartic_sum':s['off_axes_quartic_sum'],
             'S3_null_quartic_min':min(nulls),'S3_null_quartic_median':float(np.median(nulls)),
             'S3_null_quartic_max':max(nulls),
             'largest_nonvertical_line_is_horizontal':abs(s['nonvertical_line_max']-s['horizontal_line'])<1e-10,
             'chirp_c1_rayleigh_H':r['chirp_orbit_fixed_operator_rayleighs'][1]['rayleigh_H_same_operator']}
        rows.append(row)
    out={'status':'generic tested phase summaries did not provide a new operator-specific separator',
         'normalization':'unit-norm zero-extended f; A[a,b]=sum conj(f[x])*f[x+a]*e_p(b*(x+a/2))',
         'cases':rows,'all_six_full_peaks_on_amplitude_axis':all(r['full_peak_a']==0 for r in rows),
         'all_six_largest_lines_are_amplitude_axis':all(r['largest_line_is_vertical'] for r in rows),
         'five_of_six_nonvertical_maxima_are_old_translation_axis':sum(r['largest_nonvertical_line_is_horizontal'] for r in rows)==5,
         'max_direct_phase_check_error':max(r['direct_complex_max_error'] for r in validation['cases']),
         'limitations':['six selected near-edge witnesses, only p401 negative exceeds sqrt(3)/2',
                        'eight nulls per case, exploratory rather than calibrated statistical comparison',
                        'b=0 checked with exact integer correlations; b!=0 numerical',
                        'phase-plane invariants cannot alone classify fixed-operator Rayleigh behavior; no rejection of all phase information'],
         'solver_runtime':pilot['runtime'],
         'plot_runtime':{'python':sys.executable,'version':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__}}
    (HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n')
    with (HERE/'summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,9.5),layout='constrained')
    planes=np.load(HERE/'phase_plane_401_negative.npz')
    for ax,key,title in zip(axs[0],['original','chirped'],['Original p=401 negative witness','Same vector times the chirp exp(2πix²/p)']):
        data=abs(planes[key]).copy();data[0,0]=np.nan
        im=ax.imshow(data.T,origin='lower',interpolation='nearest',extent=(-.5,400.5,-.5,400.5),vmin=0,vmax=.48,cmap='magma',aspect='equal')
        ax.axvline(0,color='#00d9d4',lw=1.2,label='a=0: amplitudes only')
        ax.axhline(0,color='#90ee90',lw=1.2,label='b=0: ordinary translations')
        ax.scatter([0],[121],s=65,facecolors='none',edgecolors='white',linewidths=1.1,clip_on=False)
        ax.set(xlabel='Translation a',ylabel='Modulation b',title=title)
        ax.set_xticks([0,100,200,300,400]);ax.set_yticks([0,100,200,300,400])
    axs[0,0].legend(loc='upper right',fontsize=8,framealpha=.85)
    cb=fig.colorbar(im,ax=axs[0],shrink=.85,pad=.02);cb.set_label('|A(a,b)|; identity entry omitted')

    ax=axs[1,0];positions=np.arange(len(rows));labels=[str(r['p'])+(' +' if r['side']=='positive' else ' −') for r in rows]
    for i,r in enumerate(revised['cases']):
        vals=[x['off_axes_quartic_sum'] for x in r['S3_equivariant_shuffles']]
        jitter=np.linspace(-.15,.15,len(vals))
        ax.scatter(i+jitter,vals,color='#8a99a6',s=24,alpha=.85,label='S3-preserving nulls' if i==0 else None)
    ax.scatter(positions,[r['off_axes_quartic_sum'] for r in rows],color='#bc3f25',marker='D',s=45,label='Original witnesses',zorder=3)
    ax.set(xticks=positions,xticklabels=labels,ylabel='Σ a≠0,b≠0 |A(a,b)|⁴',title='After removing axes and preserving known symmetry')
    ax.legend(frameon=False,fontsize=8);ax.grid(axis='y',alpha=.18)

    ax=axs[1,1];width=.34
    ax.bar(positions-width/2,[r['rayleigh_H'] for r in rows],width,color='#334b66',label='Original')
    ax.bar(positions+width/2,[r['chirp_c1_rayleigh_H'] for r in rows],width,color='#d39343',label='Chirp c=1; same H')
    ax.axhline(0,color='#888',lw=.7)
    ax.set(xticks=positions,xticklabels=labels,ylabel='⟨f,Hf⟩',title='Same global ambiguity summaries, different quotient')
    ax.legend(frameon=False,fontsize=8);ax.grid(axis='y',alpha=.18)
    fig.suptitle('Phase diagnostics: stronger-looking peaks were not new spectral evidence',fontsize=15)
    fig.supxlabel('Unit norm f; A(a,b)=Σₓ conjugate(f(x)) f(x+a) exp(2πib(x+a/2)/p). Half is computed modulo p.\nSix selected witnesses; eight nulls each. Nonzero-phase values are numerical; no spectral upper bound.',fontsize=9)
    fig.savefig(HERE/'phase_controls.png',dpi=180)
    fig.savefig(HERE/'phase_controls.svg')
    print(json.dumps({'cases':len(rows),'phase_plot':str(HERE/'phase_controls.png'),'plot_runtime':out['plot_runtime']}))


if __name__=='__main__':main()
