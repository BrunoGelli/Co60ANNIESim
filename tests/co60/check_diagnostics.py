#!/usr/bin/env python3
"""Verify diagnostic artifacts against their private inputs, no detector simulation."""
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/co60'))
from ratdb_subset import blocks


def main(out):
    summary=json.loads((out/'diagnostic_summary.json').read_text())
    provenance=json.loads((out/'provenance.json').read_text())
    for name,h in provenance['inputs_sha256'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h,name
    assert summary['input_preservation']['status']=='PASS'
    assert summary['geometry']['status']=='PASS'
    assert summary['geometry']['navigator_sections']==7
    assert summary['geometry']['point_checks']==6
    assert summary['scan']['cases']==9 and summary['scan']['recorded_overlap_rejection']=='PASS'
    assert summary['root']['status']=='PASS_SYNTHETIC_ONLY'
    assert summary['root']['events']==3 and summary['root']['zero_pe_events']==2
    assert summary['root']['all_event_ids']==[0,1,2]
    assert summary['root']['pmt_pe_join_counts']==[0,1,0]
    assert summary['pmts']['channels']==132 and not summary['pmts']['matching_factory_verified']
    assert not summary['optics']['runtime_properties_verified']

    class Links(HTMLParser):
        def __init__(self): super().__init__(); self.paths=[]
        def handle_starttag(self,tag,attrs):
            for key,value in attrs:
                if key in ('href','src'): self.paths.append(value)
    links=Links(); links.feed((out/'index.html').read_text())
    for name in links.paths:
        assert (out/name).is_file(),name
    pdfs=list(out.glob('*.pdf'))
    assert len(pdfs)==10, 'expected ten paired diagnostic figures'
    for path in pdfs:
        assert path.read_bytes().startswith(b'%PDF-') and path.stat().st_size>1000
        png=path.with_suffix('.png'); Image.open(png).verify()
        assert png.name in links.paths and path.name in links.paths

    optics_input=next(Path(p) for p in provenance['inputs_sha256'] if p.endswith('/OPTICS_SANDI_wbls_Nov22.ratdb'))
    original=[b for b in blocks(optics_input) if b['index'].startswith('wbls')]
    assert original==json.loads((out/'static_wbls_tables.json').read_text()),'optical arrays changed'
    with (out/'pmt_channels.csv').open() as f: pmts=list(csv.DictReader(f))
    assert len(pmts)==len({r['channel'] for r in pmts})==132
    offset=summary['pmts']['local_offset_prediction_mm']
    for r in pmts:
        for k,d in zip('xyz',offset):
            assert abs(float(r['table_'+k+'_mm'])-float(r['predicted_'+k+'_mm'])-d)<1e-9
        assert abs(np.linalg.norm([float(r['dir_'+k]) for k in 'xyz'])-1)<1e-7
    with (out/'cup_scan.csv').open() as f: scans=list(csv.DictReader(f))
    assert len(scans)==9
    for r in scans:
        assert r['result']=='PASS' and r['probes_passed']=='18' and r['overlap_detected']=='False'
        depth=float(r['disk_thickness_mm'])+2*float(r['axial_clearance_each_mm'])
        assert abs(float(r['pocket_depth_mm'])-depth)<1e-10
        assert abs(float(r['steel_bottom_mm'])+depth-10)<1e-10
        assert abs(float(r['steel_wall_mm'])+12.7+float(r['radial_clearance_mm'])-25)<1e-10
    # Confirm no material-trace gaps or duplicated trace columns.
    requests={r[0]:r for r in map(str.split,(out/'section_requests.tsv').read_text().splitlines())}
    traces={k:{} for k in requests}
    with (out/'sections.csv').open() as f:
        for r in csv.DictReader(f):
            key=(float(r['u0_mm']),float(r['u1_mm']))
            traces[r['section']].setdefault(key,[]).append(r)
    count=0
    for name,columns in traces.items():
        request=requests[name]; assert len(columns)==int(request[9])
        for key,intervals in columns.items():
            expected=float(request[7])
            for r in intervals:
                assert abs(float(r['v0_mm'])-expected)<1e-8
                assert float(r['v1_mm'])>expected
                expected=float(r['v1_mm']); count+=1
            assert abs(expected-float(request[8]))<1e-8
    assert count==summary['geometry']['traced_intervals']
    active=[r for col in traces['cup'].values() for r in col if r['volume']=='Co60_active']
    assert active and all(abs(float(r['v1_mm'])-float(r['v0_mm'])-.01)<2e-7 for r in active)
    cap=[r for col in traces['nylon_cap'].values() for r in col if r['volume']=='superSANDIvessel']
    # The uploaded polyhedra has a 0.00025 mm radial-to-cap taper. Away from
    # the axis its end cap is slightly thicker than the nominal 0.025 mm.
    assert cap
    for r in cap:
        mid=(float(r['u0_mm'])+float(r['u1_mm']))/2
        face_distance=abs(mid)*np.cos(np.pi/8)
        expected=.025+.00025*face_distance/(900-.025)
        assert abs(float(r['v1_mm'])-float(r['v0_mm'])-expected)<1e-8
    result=dict(status='PASS',diagnostic_png_pdf_pairs=len(pdfs),
                validation_command=[sys.executable,*sys.argv],
                validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                input_hashes_verified=len(provenance['inputs_sha256']),
                contiguous_navigator_intervals=count,optical_arrays_unchanged=True,
                all_132_channel_predictions_checked=True,zero_response_fixture_events=2,
                active_layer_and_nylon_cap_thickness_checked=True)
    (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__': main(Path(sys.argv[1]).resolve())
