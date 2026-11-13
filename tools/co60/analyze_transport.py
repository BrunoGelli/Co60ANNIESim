#!/usr/bin/env python3
"""Validate/plot ACTUAL standalone P1 ROOT; truth selections, never detector triggers."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import sys
import ROOT
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load(path):
    f=ROOT.TFile.Open(str(path)); run=f.Get('runs')
    if not run or run.GetEntries()!=1: raise ValueError('incomplete ROOT: missing runs completion record')
    run.GetEntry(0)
    meta={k:str(getattr(run,k)) for k in ('run_uuid','mode','metadata_json','scope','version')}
    for k in ('requested','completed','seed'): meta[k]=int(getattr(run,k))
    for k in ('cut_mm','source_cut_mm','mass_Q_MeV','decay_threshold_s','nuclide_half_life_threshold_ns'):
        meta[k]=float(getattr(run,k))
    for k in ('synthetic','correlated','primary_prompt_clock'): meta[k]=bool(getattr(run,k))
    meta['region_names']=[str(x) for x in run.region_names]
    volumes={int(v.volume_id):dict(name=str(v.name),material=str(v.material),region=str(v.region_name)) for v in f.Get('volumes')}
    f.Close()
    cols=['decay_id','source_x_mm','source_y_mm','source_z_mm','source_valid','closure_MeV',
          'released_total_MeV','neutrino_emitted_MeV','escaped_neutrino_MeV','escaped_visible_MeV',
          'emitted_gamma_MeV','emitted_beta_MeV','total_edep_MeV','edep_regions','edep_bgo_MeV','edep_wbls_MeV',
          'cos_gamma_pair','n_gamma','n_beta','unknown_ancestry','bgo_contained_gammas','beta_daughter_excitation_MeV']
    events=ROOT.RDataFrame('events',str(path)).AsNumpy(cols)
    gamma=ROOT.RDataFrame('gamma_transport',str(path)).AsNumpy(['decay_id','root_gamma_id','initial_MeV',
          'edep_regions','escaped_MeV','dir_x','dir_y','dir_z','time_ns','first_volume_id','first_process','first_valid','bgo_interactions'])
    return dict(path=path,meta=meta,events=events,gamma=gamma,volumes=volumes)

def interval(k,n):
    # Wilson score 95% intervals include both zero and unit counts.
    z=stats.norm.ppf(.975); p=k/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den; width=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return dict(count=int(k),n=int(n),efficiency=float(p),ci95_low=float(max(0,center-width)),
                ci95_high=float(min(1,center+width)),rate_Hz=float(13061*p),
                rate_ci95_Hz=[float(13061*max(0,center-width)),float(13061*min(1,center+width))])

def validate(data,strict):
    m,e,g=data['meta'],data['events'],data['gamma']; n=m['completed']
    if m['synthetic'] or n!=m['requested'] or not np.array_equal(e['decay_id'],np.arange(n)):
        raise ValueError('not an all-generated-decay real run')
    for k,v in e.items():
        if k!='edep_regions' and not np.isfinite(v).all(): raise ValueError('nonfinite '+k)
    if not np.all(e['source_valid']==1) or not np.all(e['unknown_ancestry']==0):
        raise ValueError('source navigation or ancestry failure')
    if np.max(np.abs(e['closure_MeV']))>1e-5: raise ValueError('visible energy closure > 10 eV')
    if np.max(np.abs(e['neutrino_emitted_MeV']-e['escaped_neutrino_MeV']))>1e-8:
        raise ValueError('neutrino energy loss')
    inputs=json.loads(m['metadata_json'])
    for p,h in {**inputs['input_sha256'],**inputs['nuclear_data_sha256']}.items():
        if sha(p)!=h: raise ValueError('input/data changed '+p)
    center=inputs['source_active_center_mm']; r=inputs['active_radius_mm']; half=inputs['active_half_length_mm']
    r2=((e['source_x_mm']-center[0])**2+(e['source_y_mm']-center[1])**2)/r**2
    uz=(e['source_z_mm']-center[2]+half)/(2*half)
    if np.any(r2>1+1e-12) or np.any(uz<0) or np.any(uz>1): raise ValueError('source outside active cylinder')
    source_p=[float(stats.kstest(r2,'uniform').pvalue),float(stats.kstest(uz,'uniform').pvalue)]
    matrix=np.stack(e['edep_regions']).astype(float)
    gamma_matrix=np.stack(g['edep_regions']).astype(float)
    if not np.allclose(matrix.sum(axis=1),e['total_edep_MeV'],atol=1e-10): raise ValueError('region sum mismatch')
    gc=g['initial_MeV']-gamma_matrix.sum(axis=1)-g['escaped_MeV']
    if np.max(np.abs(gc))>1e-5: raise ValueError('gamma-lineage energy closure > 10 eV')
    # Check complete aggregate deposit joins and per-event energy, including zero records.
    deps=ROOT.RDataFrame('depositions',str(data['path'])).AsNumpy(['decay_id','edep_MeV'])
    joined=np.bincount(deps['decay_id'].astype(int),weights=deps['edep_MeV'],minlength=n)
    if not np.allclose(joined,e['total_edep_MeV'],atol=1e-10): raise ValueError('deposition join mismatch')
    ids=list(zip(g['decay_id'],g['root_gamma_id']))
    if len(ids)!=len(set(ids)) or not np.isin(g['decay_id'],e['decay_id']).all(): raise ValueError('gamma ID join')
    # Dominant pair opening angle, allowing evaluated internal conversion/minor branches.
    u=e['cos_gamma_pair']; u=u[np.abs(u)<=1]
    if len(u)<.98*n: raise ValueError('missing dominant two-gamma cascade')
    bins=np.linspace(-1,1,21); h,_=np.histogram(u,bins)
    def anti(u): return (u+u**3/24+u**5/120)/2.1
    probability=np.diff(anti(bins)); chi2=np.sum((h-len(u)*probability)**2/(len(u)*probability))
    p=float(stats.chi2.sf(chi2,19))
    leg2=(3*u*u-1)/2
    leg2_mean=float(leg2.mean()); leg2_se=float(leg2.std(ddof=1)/np.sqrt(len(leg2)))
    main=(np.abs(g['initial_MeV']-1.173239)<.001)|(np.abs(g['initial_MeV']-1.332514)<.001)
    contained=np.bincount(g['decay_id'].astype(int),weights=(main & (np.abs(gamma_matrix[:,0]-g['initial_MeV'])<.001)),minlength=n)
    isotropy_p={}
    for en in (1.173239,1.332514):
        select=np.abs(g['initial_MeV']-en)<.001
        isotropy_p[str(en)+'_z']=float(stats.kstest((g['dir_z'][select]+1)/2,'uniform').pvalue)
        phi=(np.arctan2(g['dir_y'][select],g['dir_x'][select])+np.pi)/(2*np.pi)
        isotropy_p[str(en)+'_phi']=float(stats.kstest(phi,'uniform').pvalue)
    # Energy ledger checked against independent installed Geant4 ground-state masses.
    mass_residual=e['released_total_MeV']-m['mass_Q_MeV']
    if m['mode']=='full':
        if not np.all(e['n_beta']==1): raise ValueError('missing beta/Co decay')
        if g['time_ns'].max()>1 or g['time_ns'].min()<0: raise ValueError('prompt clock failure')
        if np.max(np.abs(mass_residual))>1e-4: raise ValueError('mass-Q release discrepancy >100 eV')
        if np.count_nonzero(main)<1.98*n: raise ValueError('emitted nuclear gamma energies failed')
        excitation=e['beta_daughter_excitation_MeV']
        endpoint=np.where(np.abs(excitation-2.505753)<1e-5,.317057,1.490296)
        if np.any(e['emitted_beta_MeV']>endpoint+1e-6): raise ValueError('beta energy outside evaluated endpoint')
    if m['mode']=='pair':
        allowed=np.isclose(g['initial_MeV'],1.173,atol=1e-12,rtol=0)|np.isclose(g['initial_MeV'],1.332,atol=1e-12,rtol=0)
        if not np.all(e['n_gamma']==2) or not allowed.all() or not np.all(e['n_beta']==0):
            raise ValueError('benchmark particle/energy contract failed')
    dominant=int(np.count_nonzero(np.abs(e['beta_daughter_excitation_MeV']-2.505753)<1e-5))
    minor=int(np.count_nonzero(np.abs(e['beta_daughter_excitation_MeV']-1.332514)<1e-5))
    if strict:
        if min(source_p)<1e-5: raise ValueError('source uniformity test failed')
        if min(isotropy_p.values())<1e-5: raise ValueError('single-gamma marginal anisotropy')
        if p<1e-4 or leg2_mean<2*leg2_se: raise ValueError('angular correlation validation failed')
        if m['mode']=='full' and abs(minor/n-.0012)>5*np.sqrt(.0012*.9988/n):
            raise ValueError('evaluated minor beta branching inconsistent')
    b=e['edep_bgo_MeV']; w=e['edep_wbls_MeV']
    bgo_interaction=np.zeros(n,dtype=bool)
    np.logical_or.at(bgo_interaction,g['decay_id'].astype(int),g['bgo_interactions']>0)
    definitions={
      'BGO_gamma_interaction':bgo_interaction,
      'BGO_any_Edep':b>0,
      'BGO_Edep_gt_10keV':b>.010,
      'WbLS_Edep_gt_10keV':w>.010,
      'BGO_WbLS_both_gt_10keV':(b>.010)&(w>.010),
      'BGO_at_least_one_gamma_contained_1keV_tolerance':contained>0,
      'BGO_1173_window_10keV':np.abs(b-1.173)<.010,
      'BGO_1332_window_10keV':np.abs(b-1.332)<.010,
      'BGO_sum_window_10keV':np.abs(b-2.505)<.010}
    result=dict(status='PASS',kind='ACTUAL_STANDALONE_TRUTH_TRANSPORT',root=str(data['path']),sha256=sha(data['path']),
        mode=m['mode'],n=n,seed=m['seed'],source_center_mm=center,source_radius_mm=r,active_half_length_mm=half,
        max_visible_closure_MeV=float(np.abs(e['closure_MeV']).max()),
        max_gamma_lineage_closure_MeV=float(np.abs(gc).max()),
        installed_mass_Q_MeV=m['mass_Q_MeV'],max_release_vs_mass_Q_MeV=float(np.abs(mass_residual).max()) if m['mode']=='full' else None,
        zero_bgo_events=int(np.count_nonzero(b==0)),zero_wbls_events=int(np.count_nonzero(w==0)),
        zero_both_events=int(np.count_nonzero((b==0)&(w==0))),
        source_uniformity_KS_p=source_p,
        gamma_energies_MeV_quantiles=np.quantile(g['initial_MeV'],[0,.01,.5,.99,1]).tolist(),
        gamma_time_ns_quantiles=np.quantile(g['time_ns'],[0,.5,.99,1]).tolist(),
        dominant_pair_events=len(u),angular_chi2=float(chi2),angular_ndf=19,angular_p=p,
        angular_P2_mean=leg2_mean,angular_P2_SE=leg2_se,
        single_gamma_direction_KS_p=isotropy_p,
        stored_vs_mainline_containment_differences=int(np.count_nonzero(contained!=e['bgo_contained_gammas'])),
        beta_branch_counts=dict(excitation_2505_753keV=dominant,excitation_1332_514keV=minor),
        efficiencies={key:interval(int(mask.sum()),n) for key,mask in definitions.items()},
        mean_region_edep_MeV=dict(zip(m['region_names'],matrix.mean(axis=0).tolist())),
        strict_generator_tests=strict,
        interpretation='truth energy/interactions only; not BGO readout or DAQ efficiencies; provisional material fallbacks')
    return result


def figures(full,pair,out):
    if out.exists(): raise ValueError('choose a NEW analysis output directory')
    out.mkdir(parents=True)
    summary={d['meta']['mode']:validate(d,True) for d in (full,pair)}
    plt.rcParams.update({'font.size':11,'savefig.dpi':170,'pdf.fonttype':42})
    saved=[]
    label=f"full Co-60, N={full['meta']['completed']:,} • provisional P0 geometry • no optics"
    def save(fig,name,caption):
        fig.suptitle(label,fontsize=12); fig.savefig(out/(name+'.png'),bbox_inches='tight')
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight'); plt.close(fig)
        saved.append((name,caption))
    def hist(ax,data,bins,label_,n=None):
        counts,edges=np.histogram(data,bins); centers=(edges[1:]+edges[:-1])/2
        den=1 if n is None else n
        ax.stairs(counts/den,edges,label=label_)
        ax.errorbar(centers,counts/den,yerr=np.sqrt(counts)/den,fmt='none',alpha=.4,lw=.7)
    e=full['events']; n=len(e['decay_id']); b=e['edep_bgo_MeV']; w=e['edep_wbls_MeV']; bins=np.linspace(0,2.85,144)
    fig,ax=plt.subplots(figsize=(10,6)); hist(ax,b,bins,'All generated decays')
    ax.set(xlabel='BGO true deposited energy [MeV]',ylabel='Decays / bin',yscale='log',ylim=(.5,None))
    for en in (1.173,1.332,2.505): ax.axvline(en,ls='--',alpha=.4)
    ax.legend(); save(fig,'bgo_edep','Unconditional event spectrum, including zero deposits. Error bars are sqrt(N) per bin.')
    masks={'All decays':np.ones(n,dtype=bool),'BGO > 0.1 MeV':b>.1,
           'BGO 1.173 ± 0.010 MeV':np.abs(b-1.173)<.01,
           'BGO 1.332 ± 0.010 MeV':np.abs(b-1.332)<.01,
           'BGO sum 2.505 ± 0.010 MeV':np.abs(b-2.505)<.01}
    fig,axs=plt.subplots(1,2,figsize=(15,6))
    for title,mask in masks.items():
        count=int(mask.sum())
        hist(axs[0],w[mask],bins,title+f' (n={count})',n)
        if count: hist(axs[1],w[mask],bins,title,count)
    for ax in axs:
        ax.set(xlabel='WbLS true deposited energy [MeV]',yscale='log',ylim=(1e-5,None)); ax.legend(fontsize=8)
    axs[0].set_ylabel('Counts / all generated decays / bin'); axs[1].set_ylabel('Counts / selected decays / bin')
    save(fig,'wbls_conditioned','Left preserves all-decay normalization; right normalizes each truth-energy selection separately. '
         'No PE thresholds, BGO readout resolution or tagger trigger are modeled.')
    fig,ax=plt.subplots(figsize=(9,7)); h=ax.hist2d(b,w,bins=np.linspace(0,2.85,144),norm=matplotlib.colors.LogNorm(),cmin=1)
    ax.plot([0,2.505],[2.505,0],'--',color='white',label='rounded gamma-pair sum (reference only)')
    ax.set(xlabel='BGO Edep [MeV]',ylabel='WbLS Edep [MeV]'); fig.colorbar(h[3],ax=ax,label='Decays / cell'); ax.legend()
    save(fig,'bgo_wbls_joint','Actual full-decay events. The sum line is a benchmark reference, not a constraint for beta-containing full decays.')
    matrix=np.stack(e['edep_regions']).astype(float); escaped=e['escaped_visible_MeV']; nu=e['neutrino_emitted_MeV']
    names=full['meta']['region_names']+['escaped visible','neutrino']; budget=np.column_stack([matrix,escaped,nu])
    fig,axs=plt.subplots(1,2,figsize=(16,7))
    axs[0].bar(range(len(names)),budget.mean(axis=0),yerr=budget.std(axis=0,ddof=1)/np.sqrt(n))
    axs[0].set(xticks=range(len(names)),xticklabels=names,ylabel='Mean energy per decay [MeV]'); axs[0].tick_params(axis='x',rotation=70)
    closure=e['closure_MeV']*1e6
    axs[1].hist(closure,bins=60); axs[1].set(xlabel='Visible ledger residual [eV]',ylabel='Decays')
    axs[1].text(.04,.95,f"max |residual|={np.abs(closure).max():.3g} eV\n"
                f"max |release − G4 mass Q|={summary['full']['max_release_vs_mass_Q_MeV']*1e6:.3g} eV",
                transform=axs[1].transAxes,va='top')
    save(fig,'energy_budget','Region deposits plus outgoing kinetic/annihilation energy and neutrino accounting. '
         'Error bars on means are event-sample standard errors. Escape is defined at the outer tank surface; no external reentry.')
    fig,ax=plt.subplots(figsize=(11,6)); hist(ax,escaped,np.linspace(0,2.85,115),'Escaped visible energy',n)
    hist(ax,nu,np.linspace(0,2.85,115),'Neutrino energy',n)
    ax.set(xlabel='Energy per decay [MeV]',ylabel='Counts / generated decays / bin',yscale='log'); ax.legend()
    save(fig,'escape_energy','Separates invisible neutrino energy from particles leaving the modeled tank. Escape includes e+ annihilation reservoir if relevant.')
    fig,axs=plt.subplots(1,2,figsize=(16,6)); keys=['BGO_gamma_interaction','BGO_Edep_gt_10keV',
        'BGO_at_least_one_gamma_contained_1keV_tolerance','BGO_WbLS_both_gt_10keV']
    vals=[summary['full']['efficiencies'][k] for k in keys]
    labels=['BGO gamma\ninteraction','BGO deposit\n>10 keV','≥1 gamma fully in BGO\nwithin 1 keV','BGO & WbLS\neach >10 keV']
    for ax,scale,ylab in zip(axs,(1,13061),('Probability / original decay','Truth-selection rate at 13,061 Bq [s⁻¹]')):
        y=np.array([v['efficiency'] for v in vals]); lower=np.array([v['ci95_low'] for v in vals]); upper=np.array([v['ci95_high'] for v in vals])
        ax.bar(range(4),y*scale,yerr=np.array([y-lower,upper-y])*scale,capsize=5)
        ax.set(xticks=range(4),xticklabels=labels,ylabel=ylab)
    save(fig,'truth_efficiencies_rates','Wilson 95% statistical intervals; all generated decays form the denominator. '
         'Containment uses gamma ancestry, not an event-level photopeak window. These are not readout/tagging/DAQ efficiencies or operating rates.')
    fig,axs=plt.subplots(2,3,figsize=(16,10))
    for d,color in ((full,'#087f8c'),(pair,'#b06a27')):
        ev=d['events']; gam=d['gamma']; count=d['meta']['completed']; mode=d['meta']['mode']
        hist(axs[0,0],ev['edep_bgo_MeV'],bins,mode,count); hist(axs[0,1],ev['edep_wbls_MeV'],bins,mode,count)
        u=ev['cos_gamma_pair']; u=u[np.abs(u)<=1]; hist(axs[0,2],u,np.linspace(-1,1,21),mode,len(u))
        axs[1,0].hist(gam['initial_MeV'],bins=300,range=(0,2.6),histtype='step',label=mode)
        axs[1,1].hist(ev['emitted_beta_MeV'],bins=80,range=(0,1.5),histtype='step',label=mode)
        axs[1,2].hist(gam['time_ns'],bins=100,range=(0,.1),histtype='step',label=mode)
    u=np.linspace(-1,1,200); axs[0,2].plot(u,(1+u*u/8+u**4/24)/2.1*.1,'k--',label='dominant-cascade model × bin width')
    for ax,xlab in zip(axs.flat,('BGO Edep [MeV]','WbLS Edep [MeV]','Gamma opening-angle cosine',
                      'Emitted gamma energy [MeV]','Emitted beta energy [MeV]','Prompt-relative gamma time [ns]')):
        ax.set_xlabel(xlab); ax.legend(fontsize=8)
    axs[0,0].set_yscale('log'); axs[0,1].set_yscale('log'); axs[1,0].set_yscale('log'); axs[1,1].set_yscale('log')
    axs[0,0].set_ylabel('Counts / generated decay / bin'); axs[0,1].set_ylabel('Counts / generated decay / bin')
    axs[0,2].set_ylabel('Counts / dominant pair / bin')
    save(fig,'generator_comparison','Full mode uses evaluated beta/atomic/cascade machinery with native correlations. '
         'Pair mode emits correlated 1.173 and 1.332 MeV gammas at t=0, no beta/atomic emissions. Lower panels show raw counts; '
         'their sample sizes are recorded in the summaries.')
    fig,axs=plt.subplots(1,3,figsize=(15,5)); center=np.array(summary['full']['source_center_mm'])
    axs[0].hist2d(e['source_x_mm']-center[0],e['source_y_mm']-center[1],bins=40)
    axs[0].set(xlabel='Source x − active center [mm]',ylabel='Source y − active center [mm]',aspect='equal')
    axs[1].hist(e['source_z_mm']-center[2],bins=40); axs[1].set(xlabel='Source z − active center [mm]',ylabel='Decays')
    axs[2].bar(['all events','zero BGO','zero WbLS','zero both'],[n,int((b==0).sum()),int((w==0).sum()),int(((b==0)&(w==0)).sum())])
    axs[2].set_ylabel('Original decay rows retained'); axs[2].tick_params(axis='x',rotation=20)
    save(fig,'source_and_zero_retention','Actual source samples in the provisional active region; zero-energy responses remain in events. '
         'Cup center is (0,750,0) mm; active center is (0,750,3.9) mm, not the cup origin.')
    (out/'physics_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (out/'truth_efficiencies.csv').open('w') as f:
        rows=[dict(mode=mode,selection=k,**v) for mode,s in summary.items() for k,v in s['efficiencies'].items()]
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    body='<h1>P1 standalone Co-60 truth transport</h1><p><b>Actual simulated events; provisional P0 geometry/materials. '
    body+='No optical tracking, BGO readout or RAT integration. Rates are truth selections, not DAQ predictions.</b></p>'
    body+='<p>All counts/intervals: <a href="physics_summary.json">physics_summary.json</a>; '
    body+='<a href="truth_efficiencies.csv">truth_efficiencies.csv</a>.</p>'
    for name,caption in saved:
        body+=f'<h2>{name.replace("_"," ")}</h2><a href="{name}.pdf">PDF</a><figure><img src="{name}.png">'
        body+='<figcaption>'+html.escape(caption)+'</figcaption></figure>'
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>P1 truth transport</title>'
        '<style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:20px}img{width:100%}figure{margin:20px 0}</style>'+body)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    v=sub.add_parser('validate'); v.add_argument('root',type=Path); v.add_argument('--out',type=Path,required=True)
    v.add_argument('--strict',action='store_true')
    plot=sub.add_parser('plot'); plot.add_argument('--full',type=Path,required=True); plot.add_argument('--pair',type=Path,required=True)
    plot.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.command=='validate':
        result=validate(load(a.root),a.strict); a.out.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
    else:
        result=figures(load(a.full),load(a.pair),a.out); print('PASS: actual physics plots in',a.out)
