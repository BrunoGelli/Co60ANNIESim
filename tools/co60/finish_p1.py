#!/usr/bin/env python3
"""Write P1 validation figures, run manifest and report from actual private outputs."""
import argparse
import hashlib
import html
import json
import re
from pathlib import Path
import subprocess
import sys
import ROOT
import numpy as np
from scipy import stats,integrate
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from analyze_transport import interval

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main(a):
    root=a.results.resolve(); out=root/'analysis'
    if not (out/'index.html').exists(): raise ValueError('run main physics plots first')
    summaries={name:json.loads((root/(name+'_validation.json')).read_text()) for name in
               ('full_10000','pair_10000','full_100000','pair_100000','full_cut005_10000')}
    if not all(s['status']=='PASS' and s['strict_generator_tests'] for s in summaries.values()):
        raise ValueError('all strict validations must pass before reporting')
    full=summaries['full_100000']; pair=summaries['pair_100000']; n=full['n']
    f=ROOT.TFile.Open(str(root/'full_100000.root')); runs=f.Get('runs'); runs.GetEntry(0)
    metadata=json.loads(str(runs.metadata_json)); mode=str(runs.mode); f.Close()
    input_dir=root/('inputs_v2' if (root/'inputs_v2').exists() else 'inputs')
    gam=ROOT.RDataFrame('gamma_transport',str(root/'full_100000.root')).AsNumpy(
        ['initial_MeV','decay_id','time_ns'])
    time_arrays=[]; energy_means={}
    for energy in (1.173239,1.332514):
        mask=np.abs(gam['initial_MeV']-energy)<.001
        times=np.full(n,np.nan); times[gam['decay_id'][mask].astype(int)]=gam['time_ns'][mask]
        time_arrays.append(times)
        en=gam['initial_MeV'][mask]
        energy_means[str(energy)]=dict(n=len(en),mean_MeV=float(en.mean()),min_MeV=float(en.min()),max_MeV=float(en.max()))
    delta=time_arrays[1]-time_arrays[0]; delta=delta[np.isfinite(delta)]
    tau=.000735/np.log(2) # Evaluated 1332.514-keV Ni-60 level: half-life 7.35e-13 s.
    if np.any(delta<0): raise ValueError('cascade clock reverses')
    timing_p=float(stats.kstest(delta,'expon',args=(0,tau)).pvalue)
    expected_p2=integrate.quad(lambda u:((3*u*u-1)/2)*(1+u*u/8+u**4/24)/2.1,-1,1)[0]
    control=ROOT.RDataFrame('gamma_transport',str(root/'ni60_timing_5000.root')).AsNumpy(['initial_MeV','time_ns'])
    stationary=control['time_ns'][np.abs(control['initial_MeV']-1.332514)<.001]
    stationary_p=float(stats.kstest(stationary,'expon',args=(0,tau)).pvalue)
    if stationary_p<1e-5: raise ValueError('stationary Ni-60 evaluated lifetime check failed')
    plt.rcParams.update({'font.size':10,'savefig.dpi':170,'pdf.fonttype':42})
    fig,axs=plt.subplots(2,3,figsize=(16,10),layout='constrained')
    for ax,en,reference in zip(axs[0,:2],(1.173239,1.332514),(1.173,1.332)):
        mask=np.abs(gam['initial_MeV']-en)<.001
        lo=min(reference,en)-.00005; hi=max(reference,en)+.00008
        ax.hist(gam['initial_MeV'][mask]*1000,bins=np.linspace(lo*1000,hi*1000,120),histtype='step')
        ax.axvline(en*1000,ls='--',color='#777',label='nuclear level difference')
        ax.axvline(reference*1000,ls=':',color='#a56622',label='rounded benchmark')
        ax.set(xlabel='Emitted gamma energy [keV]',ylabel='Full-decay gamma count',yscale='log'); ax.legend(fontsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.ticklabel_format(axis='x',useOffset=False,style='plain')
    branch=full['beta_branch_counts']; values=[branch['excitation_2505_753keV'],branch['excitation_1332_514keV']]
    eff=[interval(k,n) for k in values]; y=np.array([v['efficiency'] for v in eff])
    err=np.array([y-np.array([v['ci95_low'] for v in eff]),np.array([v['ci95_high'] for v in eff])-y])
    axs[0,2].errorbar(range(2),y,yerr=err,fmt='o',capsize=4,label='sample, Wilson 95%')
    axs[0,2].plot(range(2),[.9988,.0012],'x',label='installed RDM data')
    axs[0,2].set(xticks=[0,1],xticklabels=['Ni* 2505.753 keV','Ni* 1332.514 keV'],ylabel='Beta branch fraction',yscale='log')
    axs[0,2].legend(fontsize=8)
    axs[1,0].hist(delta,bins=80,range=(0,.01),density=True,histtype='step',label='actual cascade delay')
    axs[1,0].hist(stationary,bins=80,range=(0,.01),density=True,histtype='step',label='stationary Ni-60 control')
    t=np.linspace(0,.01,300); axs[1,0].plot(t,stats.expon.pdf(t,scale=tau),'--',label='evaluated level lifetime')
    axs[1,0].set(xlabel='t(1332 gamma) − t(1173 gamma) [ns]',ylabel='Probability density [1/ns]')
    axs[1,0].legend(fontsize=8)
    for i,s in enumerate((full,pair)):
        axs[1,1].errorbar(i,s['angular_P2_mean'],yerr=1.96*s['angular_P2_SE'],fmt='o',capsize=5)
    axs[1,1].axhline(expected_p2,ls='--',label='dominant-cascade expectation')
    axs[1,1].axhline(0,ls=':',color='#777',label='independent isotropic expectation')
    axs[1,1].set(xticks=[0,1],xticklabels=['evaluated full','correlated pair'],ylabel='Mean Legendre P₂, 95% statistical interval')
    axs[1,1].legend(fontsize=8)
    axs[1,2].axis('off')
    text=(f"Actual full decays: {n:,}\n"
          f"Pair benchmark: {pair['n']:,}\n\n"
          f"Angular chi² / ndf: {full['angular_chi2']:.2f} / {full['angular_ndf']}\n"
          f"Angular p: {full['angular_p']:.3f}\n"
          f"Stationary Ni lifetime KS p: {stationary_p:.3f}\n"
          'Recoil cascade delay: timing limitation\n\n'
          'Full: native evaluated branching, beta,\ninternal conversion / atomic relaxation.\n'
          'Benchmark: two correlated photons only.\nNo activity clock or optical response.')
    axs[1,2].text(.05,.95,text,va='top',fontsize=12)
    fig.suptitle('Actual generator validation • prompt clock, evaluated source and angular correlations',fontsize=15)
    for ext in ('png','pdf'): fig.savefig(out/('generator_validation.'+ext),bbox_inches='tight')
    plt.close(fig)
    # Cut comparison: conservative independent-run statistical scale; same seed is not common event history.
    base=summaries['full_10000']; cut=summaries['full_cut005_10000']
    keys=['BGO_Edep_gt_10keV','WbLS_Edep_gt_10keV','BGO_WbLS_both_gt_10keV',
          'BGO_at_least_one_gamma_contained_1keV_tolerance']
    comparisons={}
    fig,ax=plt.subplots(figsize=(12,6),layout='constrained')
    for i,(s,label) in enumerate(((base,'bulk 0.010 mm'),(cut,'bulk 0.005 mm'))):
        vals=[s['efficiencies'][k] for k in keys]; y=np.array([v['efficiency'] for v in vals])
        err=np.array([y-np.array([v['ci95_low'] for v in vals]),np.array([v['ci95_high'] for v in vals])-y])
        ax.errorbar(np.arange(4)+(i-.5)*.12,y,yerr=err,fmt='o',capsize=4,label=label)
    for k in keys:
        v1=base['efficiencies'][k]; v2=cut['efficiencies'][k]
        p1,p2=v1['efficiency'],v2['efficiency']; se=np.sqrt(p1*(1-p1)/base['n']+p2*(1-p2)/cut['n'])
        comparisons[k]=dict(delta=p2-p1,combined_binomial_SE=float(se),delta_over_SE=float((p2-p1)/se))
    ax.set(xticks=range(4),xticklabels=['BGO >10 keV','WbLS >10 keV','BGO & WbLS >10 keV','≥1 main gamma contained'],
           ylabel='Truth probability / decay, Wilson 95% intervals',
           title='Production-cut spot check: 10,000 full decays each; source-region cut remains 0.001 mm')
    ax.legend(); fig.text(.06,-.02,'No systematic convergence certificate: this is one factor-two bulk-cut comparison at finite statistics.')
    for ext in ('png','pdf'): fig.savefig(out/('cut_comparison.'+ext),bbox_inches='tight')
    plt.close(fig)
    natural=ROOT.RDataFrame('gamma_transport',str(root/'natural_clock_10.root')).AsNumpy(['time_ns'])['time_ns']
    if len(natural)==0 or min(natural)<1e12: raise ValueError('native Co-60 clock/threshold diagnostic failed')
    validation=dict(status='PASS_WITH_TIMING_LIMITATION' if timing_p<1e-5 else 'PASS',energy_lines=energy_means,cascade_delay_mean_ns=float(delta.mean()),
        expected_cascade_delay_mean_ns=float(tau),cascade_delay_KS_p=timing_p,
        stationary_Ni60_mean_ns=float(stationary.mean()),stationary_Ni60_KS_p=stationary_p,
        timing_limitation='native recoil transport delay differs from a stationary exponential; needs investigation before timing/optical use',
        expected_P2=expected_p2,cut_comparison=comparisons,
        natural_clock_gamma_time_range_s=[float(min(natural)/1e9),float(max(natural)/1e9)],
        note='natural-clock sample is a threshold diagnostic only; precision-losing absolute times are not used in physics predictions')
    (root/'additional_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    existing=(out/'index.html').read_text().split('<!--P1_EXTRA-->')[0].split('<h2>Additional generator and cut checks</h2>')[0]
    extra='<!--P1_EXTRA--><h2>Additional generator and cut checks</h2>'
    extra+='<p><b>Timing limitation:</b> native recoil-cascade delays differ from the stationary exponential. '
    extra+='The stationary Ni-60 lifetime control passes. No timing selections enter the energy/rate results.</p>'
    for name in ('generator_validation','cut_comparison'):
        extra+=f'<a href="{name}.pdf">PDF</a><img src="{name}.png">'
    extra+='<p><a href="../P1_REPORT.md">P1 report</a> · <a href="../additional_validation.json">Additional checks</a></p>'
    (out/'index.html').write_text(existing+extra)
    # Concise report: quantities are model truth, never physical trigger/DAQ rates.
    lines=['# P1 standalone Co-60 transport — provisional validation\n',
      'Actual Geant4 11.4.0 transport; no RAT integration, optical photons, BGO readout, pileup or DAQ model. '
      'Truth rates below use 13,061 decays/s. Activity does not affect generation.\n',
      '## Completed runs\n',
      'Full evaluated RDM and the correlated gamma-pair benchmark each passed 10,000-event validation, then 100,000 events. '
      'Seeds: full 601001/601002; pair 602001/602002. A 10,000-event half-bulk-cut run (seed 601001) and '
      '10-event natural-clock diagnostic (600103) also completed, plus 5,000 stationary Ni-60 timing-control events '
      '(600104). ROOT files/logs are in this private directory.\n',
      '## Geometry and materials\n',
      'The accepted P0 manifest/CSG is retained: complete holder assembly, 10 mm steel separator, 50×50 mm BGO '
      'with provisional 1 mm bore radius, octagonal SuperSANDI WbLS and surrounding Gd-water/tank shells. '
      'All 18 transport material/navigation probes and sampled overlap checks pass. Source samples are uniform '
      'inside r=1 mm, thickness=0.01 mm active region at **(0,750,3.9) mm**; (0,750,0) mm is the cup center. '
      'The disk is provisionally 2 mm polystyrene with 25.4 mm diameter; radial/axial gaps are 0.25/0.1 mm. '
      'These dimensions and placement are not measured hardware.\n',
      'Uploaded material densities/compositions, including WbLS, Gd-water and BGO, are exported unchanged. '
      'Missing-core definitions require explicit fallbacks: NIST stainless steel, air, water, Pyrex PMT glass, '
      'acrylic/Delrin/vacuum; foam is assumed 0.03 g/cm³ polystyrene, PPO uses C15H11NO. '
      'The uploaded nylon fractions are interpreted as mass fractions and preserved even though their '
      'hydrogen content deserves review. Full tank PMT/support/black-sheet solids and the hall are omitted. '
      'Thus these are conditional predictions for this provisional standalone model, not final ANNIE efficiencies.\n',
      '## Source and energy validation\n',
      'Native G4VRadioactiveDecay uses the installed RadioactiveDecay6.1.2, PhotonEvaporation6.1.2, '
      'ENSDFSTATE3.0 and EMLOW8.8 data. The default long-decay threshold is one year; configured threshold '
      'is 1e15 s. A primary-only prompt adapter bypasses Co-60 waiting before creating daughters, without '
      'changing its physical lifetime (~2.399839e8 s mean). Ni-60 picosecond levels are retained using '
      'a 1e-6 ns nuclide half-life threshold. Native correlations/internal conversion/atomic relaxation are enabled. '
      'The benchmark emits correlated 1.173 and 1.332 MeV photons in one event (sum **2.505 MeV**); it has no beta/atomic branch.\n',
      f"Full sample: {full['beta_branch_counts']['excitation_2505_753keV']:,} dominant and "
      f"{full['beta_branch_counts']['excitation_1332_514keV']} minor beta branches (expected fractions 0.9988/0.0012). "
      f"Angular P₂={full['angular_P2_mean']:.6f} ± {full['angular_P2_SE']:.6f} (one standard error); "
      f"dominant-cascade expectation {expected_p2:.6f}; angular goodness-of-fit p={full['angular_p']:.3f}. "
      f"Stationary Ni-60 lifetime KS p={stationary_p:.3f}. Direction/source uniformity checks pass.\n",
      f"**Timing limitation:** the stationary Ni-60 control has mean delay {stationary.mean():.7f} ns "
      f"(expected {tau:.7f} ns), but the full recoiling cascade gives mean separation {delta.mean():.7f} ns "
      f"and rejects that simple stationary exponential (KS p={timing_p:.3g}). "
      'The native in-flight/stopping/at-rest time handling needs investigation before any fine timing or optical inference. '
      'All emissions remain prompt on the nanosecond scale; no timing cuts are used in the energy/efficiency results. '
      'No lifetime values were fitted or retuned to hide this difference.\n',
      f"Visible transport residual max={full['max_visible_closure_MeV']*1e6:.4g} eV; "
      f"gamma-lineage residual max={full['max_gamma_lineage_closure_MeV']*1e6:.4g} eV. "
      f"Independent release vs installed bare-ion mass-Q residual max={full['max_release_vs_mass_Q_MeV']*1e6:.4g} eV "
      '(acceptance 100 eV). The installed bare-ion mass-Q is '
      f"{full['installed_mass_Q_MeV']:.9f} MeV; this is not silently equated to a neutral-atom mass difference. "
      'Native beta/atomic/mass-model approximations remain; no evaluated data or masses were retuned. '
      'Escape is scored at the outer tank boundary, with neutrinos separate and an e+ annihilation-energy reservoir. '
      'External hall backscatter/reentry is omitted.\n',
      f"All {n:,} original decay rows are present, including {full['zero_bgo_events']:,} zero-BGO, "
      f"{full['zero_wbls_events']:,} zero-WbLS and {full['zero_both_events']:,} zero-both records. "
      'Stable UUID/decay joins, per-volume/root deposits, nuclear emissions, first interactions and gamma containment are persisted in flat native ROOT.\n',
      '## Full-decay truth results\n',
      '| Selection | Count / 100,000 | Probability, 95% Wilson interval | Rate at 13,061 Bq, 95% interval [s⁻¹] |',
      '|---|---:|---:|---:|']
    names={'BGO_gamma_interaction':'BGO nuclear-gamma interaction (including Rayleigh)',
      'BGO_Edep_gt_10keV':'BGO Edep >10 keV',
      'BGO_at_least_one_gamma_contained_1keV_tolerance':'≥1 main cascade gamma contained in BGO (1 keV tolerance)',
      'BGO_WbLS_both_gt_10keV':'BGO and WbLS each >10 keV'}
    for key,label in names.items():
        v=full['efficiencies'][key]
        lines.append(f"| {label} | {v['count']:,} | {v['efficiency']:.5f} [{v['ci95_low']:.5f}, {v['ci95_high']:.5f}] "
          f"| {v['rate_Hz']:.1f} [{v['rate_ci95_Hz'][0]:.1f}, {v['rate_ci95_Hz'][1]:.1f}] |")
    maxshift=max(abs(v['delta_over_SE']) for v in comparisons.values())
    lines.extend(['\nStatistical intervals exclude geometry/material/activity/physics-model uncertainty. '
      'Containment follows the original gamma and its descendants, not a photopeak selected on summed event energy. '
      'No BGO trigger/tagger/DAQ operating-rate claim is implied.\n',
      f"The factor-two bulk-cut check (0.010→0.005 mm; source cut 0.001 mm unchanged) changes the four listed "
      f"probabilities by at most {maxshift:.2f} combined binomial standard errors. This is a spot check, "
      'not a precision convergence certificate. Full/pair comparison and conditional WbLS spectra are in the HTML/PNG/PDF package.\n',
      '## Evidence and remaining work\n',
      '[Physics plots](analysis/index.html), `*_validation.json`, `additional_validation.json`, '
      f'`{input_dir.name}/transport_provenance.json`, `run_manifest.json`, and `logs/` preserve values, '
      'units, seeds, material fallbacks, code/data hashes and commands. Early failed smoke outputs are debugging artifacts '
      'without completion metadata and are excluded from analysis. The adapter now retains both its copied track '
      'and the required non-owning step pointer through Geant4 step updates.\n',
      'Before final ANNIE predictions: resolve missing-core material definitions; verify placement, bore/source thickness/active distribution; '
      'include remaining detector internals and validate full geometry; investigate source-model/mass/atomic details and broader cuts. '
      'Matching RAT optical/PMT validation and a BGO response model remain later work. Original archive and tuned optics are untouched; no public push occurred.\n',
      'Source references: [Geant4 11.4 photon evaporation](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/hadronic/Evaporation/SimulationPhotonEvaporation.html), '
      '[RDM](https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/hadronic/Radioactive/radDecay.html), '
      '[dominant Co-60 angular-correlation benchmark](https://legend-exp.github.io/remage/validation/v0.24.0/gammacorr.html). '
      'Installed source/data hashes, not web documentation alone, pin actual behavior.\n'])
    (root/'P1_REPORT.md').write_text('\n'.join(lines))
    repo=Path(__file__).resolve().parents[2]
    source_files=[repo/'transport/co60/transport.cc',repo/'transport/co60/CMakeLists.txt',
      repo/'integration/annie_ratpac/ManifestGeometry.hh',*list((repo/'tools/co60').glob('*.py')),repo/'tools/co60/run_p1.sh']
    completed=['full_10000','pair_10000','full_100000','pair_100000','full_cut005_10000','natural_clock_10','ni60_timing_5000']
    commands=[]
    for name in completed:
        if name=='natural_clock_10': mode,N,seed,cut='natural_clock_test',10,600103,.01
        elif name=='ni60_timing_5000': mode,N,seed,cut='ni60_timing_test',5000,600104,.01
        else:
            mode='pair' if name.startswith('pair') else 'full'; N=100000 if name.endswith('100000') else 10000
            seed=(602002 if N==100000 else 602001) if mode=='pair' else (601002 if N==100000 else 601001)
            cut=.005 if 'cut005' in name else .01
        commands.append([str(root/'build/co60_transport'),str(input_dir),mode,str(N),str(seed),str(cut),str(root/(name+'.root')),'1'])
    artifacts=[root/(name+'.root') for name in completed]+list(root.glob('*_validation.json'))+[
               root/'additional_validation.json',root/'build/co60_transport',root/'co60_transport_reference',root/'transport_reference.cc']
    artifacts=[p for p in artifacts if p.exists()]
    input_files=[p for p in input_dir.iterdir() if p.is_file()]
    archive=a.app.parent/'ANNIE_RATPAC2.tar.gz'
    optics=a.app/'ratdb/ANNIE/OPTICS_SANDI_wbls_Nov22.ratdb'
    model=json.loads((input_dir/'geometry_manifest.json').read_text())
    if sha(optics)!=model['static_optics']['sha256']: raise ValueError('tuned optics changed')
    preserve={str(Path(p)):sha(Path(p))==h for p,h in metadata['input_sha256'].items()}
    if not all(preserve.values()): raise ValueError('authoritative input changed')
    manifest=dict(status=validation['status'],standalone_transport_commands=commands,source_sha256={str(p.relative_to(repo)):sha(p) for p in source_files},
        actual_input_sha256={str(p):sha(p) for p in input_files},output_sha256={str(p):sha(p) for p in artifacts},
        authoritative_input_preservation=preserve,archive_sha256=sha(archive),
        unchanged_optics_sha256=sha(optics),
        current_branch=subprocess.check_output(['git','-C',str(repo),'branch','--show-current'],text=True).strip(),
        source_head=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
        native_RAT_runtime_used=False,full_optical_tracking=False,activity_in_generation=False,
        host_kernel=subprocess.check_output(['uname','-srm'],text=True).strip())
    prefix=Path(subprocess.check_output(['geant4-config','--prefix'],text=True).strip())
    cache=prefix.parent/'build/CMakeCache.txt'
    if cache.exists():
        match=re.search(r'^CMAKE_HOME_DIRECTORY:INTERNAL=(.*)$',cache.read_text(),re.M)
        if match:
            src=Path(match.group(1)); manifest['geant4_source_directory']=str(src)
            wanted=[src/'source/processes/hadronic/models/radioactive_decay/src'/f for f in
                    ('G4VRadioactiveDecay.cc','G4BetaMinusDecay.cc','G4ITDecay.cc')]
            wanted.append(src/'source/processes/hadronic/models/de_excitation/photon_evaporation/src/G4PhotonEvaporation.cc')
            manifest['inspected_geant4_source_sha256']={str(p):sha(p) for p in wanted if p.exists()}
    (root/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(validation['status']+': report and additional checks:',root/'P1_REPORT.md')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--results',type=Path,required=True)
    p.add_argument('--app',type=Path,required=True); main(p.parse_args())
