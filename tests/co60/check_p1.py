#!/usr/bin/env python3
"""Final P1 artifact and guard checks. Never starts additional transport events."""
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import subprocess
import sys
from PIL import Image
import ROOT

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main(root):
    manifest=json.loads((root/'run_manifest.json').read_text())
    extra=json.loads((root/'additional_validation.json').read_text())
    assert manifest['status']=='PASS_WITH_TIMING_LIMITATION'
    assert extra['stationary_Ni60_KS_p']>1e-5 and extra['cascade_delay_KS_p']<1e-5
    for section in ('source_sha256','actual_input_sha256','output_sha256'):
        for name,h in manifest[section].items():
            path=Path(name) if section!='source_sha256' else Path(__file__).resolve().parents[2]/name
            assert sha(path)==h,path
    assert manifest['archive_sha256']=='0045f8f8a116190d82ad4eb89bc869f4c7803c19877a68fd97ac1bedbeeb03dd'
    assert all(manifest['authoritative_input_preservation'].values())
    assert not manifest['activity_in_generation'] and not manifest['full_optical_tracking']
    for stem,n in [('full_10000',10000),('pair_10000',10000),('full_100000',100000),('pair_100000',100000),
                   ('full_cut005_10000',10000),('natural_clock_10',10),('ni60_timing_5000',5000)]:
        file=ROOT.TFile.Open(str(root/(stem+'.root')))
        run=file.Get('runs'); assert run.GetEntries()==1; run.GetEntry(0)
        assert int(run.completed)==int(run.requested)==file.Get('events').GetEntries()==n
        assert not run.synthetic and run.correlated
        for name in ('events','volumes','emissions','gamma_transport','depositions','decay_vertices','escapes'):
            assert file.Get(name).InheritsFrom('TTree')
        file.Close()
        log=(root/'logs'/(stem+'.log')).read_text()
        assert 'navigation probes passed=18' in log and 'SUCCESS ROOT complete' in log
        assert 'segmentation violation' not in log
    class Links(HTMLParser):
        def __init__(self): super().__init__(); self.paths=[]
        def handle_starttag(self,tag,attrs):
            self.paths.extend(v for k,v in attrs if k in ('href','src'))
    links=Links(); links.feed((root/'analysis/index.html').read_text())
    for path in links.paths: assert (root/'analysis'/path).is_file(),path
    pdfs=list((root/'analysis').glob('*.pdf')); assert len(pdfs)==10
    for pdf in pdfs:
        assert pdf.read_bytes().startswith(b'%PDF-') and pdf.stat().st_size>1000
        Image.open(pdf.with_suffix('.png')).verify()
    full=json.loads((root/'full_100000_validation.json').read_text())
    assert full['zero_bgo_events']==64723 and full['zero_wbls_events']==8214
    assert full['stored_vs_mainline_containment_differences']==0
    # Existing outputs and independently isotropic production are rejected BEFORE simulation.
    binary=root/'build/co60_transport'; original=root/'full_100000.root'; oldsha=sha(original)
    cmd=[str(binary),str(root/'inputs_v2'),'full','1','42','.01',str(original),'1']
    result=subprocess.run(cmd,capture_output=True,text=True)
    assert result.returncode!=0 and 'refusing existing ROOT output' in result.stderr
    assert sha(original)==oldsha
    rejected=root/'SHOULD_NOT_EXIST.root'; assert not rejected.exists()
    cmd[-2]=str(rejected); cmd[-1]='0'
    result=subprocess.run(cmd,capture_output=True,text=True)
    assert result.returncode!=0 and 'correlations must be enabled' in result.stderr and not rejected.exists()
    report=dict(status='PASS_ARTIFACTS_WITH_DECLARED_TIMING_LIMITATION',native_completed_runs=7,
        png_pdf_pairs=10,main_full_event_count=100000,zero_BGO_events_retained=64723,
        zero_WbLS_events_retained=8214,root_overwrite_guard='PASS',independent_pair_rejection='PASS',
        hashes_verified=True,validator_sha256=sha(Path(__file__)),command=[sys.executable,*sys.argv])
    (root/'final_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main(Path(sys.argv[1]).resolve())
