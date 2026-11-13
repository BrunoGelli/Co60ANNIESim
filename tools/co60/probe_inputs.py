#!/usr/bin/env python3
"""Create deterministic material probes from the PRIVATE generated geometry."""
import json
import sys
from pathlib import Path

out=Path(sys.argv[1])
m=json.loads((out/'geometry_manifest.json').read_text())
s=m['vessel_shape']; v=m['vessel']['nylon_vessel_position']
(out/'vessel.txt').write_text(' '.join(map(str,[s['height_mm'],s['radius_parameter_mm'],
    s['nylon_thickness_mm'],s['gradient_mm'],s['faces'],*v]))+'\n')
positions={'superSANDIDetMed':v}
vols={b['index']:b for b in m['volumes']}
for b in m['volumes']:
    positions[b['index']]=[x+y for x,y in zip(positions[b['mother']],b['position'])]
rows=[]
def probe(label, volume, delta=(0,0,0)):
    xyz=[x+y for x,y in zip(positions[volume],delta)]
    rows.append([label,*xyz,volume,vols[volume]['material']])
probe('source_active','Co60_active')
probe('plastic_disk','Co60_disk',(5,0,0))
probe('radial_air_clearance','Co60_pocket',(12.7+m['parameters']['radial_clearance_mm']/2,0,0))
probe('axial_air_clearance','Co60_pocket',(0,0,m['parameters']['disk_thickness_mm']/2+
                                                        m['parameters']['axial_clearance_mm']/2))
probe('steel_wall','AmBe_placer',(20,0,0))
probe('steel_bottom','AmBe_placer',(0,0,-5+m['cup_bottom_mm']/2))
probe('separator_10mm','AmBe_steel_disk')
probe('BGO_interior','AmBe_BGO',(10,0,0))
# Hole is air inherited from housing-air mother.
bgo=positions['AmBe_BGO']
rows.append(['BGO_bore',*bgo,'AmBe_housing_air','air'])
probe('PMT_window','AmBe_PMT_glass',(0,0,-48.5))
probe('PMT_vacuum','AmBe_PMT_vacuum')
probe('PVC_wall','AmBe_housing_body',(30,0,50))
mat=m['vessel']['detection_medium_material']
rows.extend([
 ['WbLS_origin',0,0,0,'superSANDIDetMed',mat],
 ['WbLS_by_holder',100,750,0,'superSANDIDetMed',mat],
 ['WbLS_below_cap',0,0,v[2]+s['height_mm']/2-.04,'superSANDIDetMed',mat],
 ['nylon_cap',0,0,v[2]+s['height_mm']/2-.0125,'superSANDIvessel','superSANDI_nylon'],
 ['water_above_cap',0,0,v[2]+s['height_mm']/2+.01,'detector','water_gdS_0p2'],
 ['water_outside',1200,0,0,'detector','water_gdS_0p2']])
(out/'probes.tsv').write_text('\n'.join(' '.join(map(str,row)) for row in rows)+'\n')
# A complete transform/units record for review; no PMTINFO offsets are reapplied.
(out/'global_centers.json').write_text(json.dumps(positions,indent=2)+'\n')
