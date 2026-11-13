#!/usr/bin/env python3
"""Generate a PRIVATE variant from the supplied dirty application, never upstream.

Outputs contain collaboration geometry: destination must be outside this repo.
No optical tables are written. Unknown dimensions are mandatory CLI inputs.
"""
import argparse
import hashlib
import json
import copy
import math
from pathlib import Path
from ratdb_subset import blocks, dump_blocks


def prepare(app, out, p):
    repo = Path(__file__).resolve().parents[2]
    if out.resolve().is_relative_to(repo):
        raise ValueError('private generated geometry must be outside the public repository')
    for value in vars(p).values():
        if isinstance(value,float) and not math.isfinite(value):
            raise ValueError('nonfinite geometry parameter')
    if not all(math.isfinite(x) for x in p.source_position_mm):
        raise ValueError('nonfinite position')
    if not (0 < p.disk_thickness_mm + 2*p.axial_clearance_mm < 10):
        raise ValueError('disk and clearances must leave a bottom in the 10 mm cup')
    if not (0 < p.active_radius_mm <= 12.7 and
            0 < p.active_thickness_mm < p.disk_thickness_mm and
            abs(p.active_z_mm) + p.active_thickness_mm/2 < p.disk_thickness_mm/2):
        raise ValueError('active region must be strictly inside disk')
    if not (0 < p.radial_clearance_mm < 12.3 and p.axial_clearance_mm > 0 and
            0 < p.bgo_bore_radius_mm < 25):
        raise ValueError('positive clearances and nonzero BGO bore required')
    base = app/'ratdb/ANNIE'
    detector = blocks(base/'ANNIE_run2_config14_scan.geo')
    inner = next(b for b in detector if b['index'] == 'InnerStructure_Holders_Blacksheets')
    inner['enable_superSANDI'] = 1
    vessel = next(b for b in detector if b['name'] == 'SuperSANDIGeo')
    vessel['co60_geometry_enabled'] = 1
    vessel['co60_geometry_manifest'] = str((out/'assembly.tsv').resolve())
    original = blocks(base/'ANNIE_AmBesource_v2.geo')
    assembly = [copy.deepcopy(b) for b in original if b.get('enable', 1) and
                b['index'] not in ('AmBe_case', 'AmBeSource')]
    byid = {b['index']: b for b in assembly}
    housing = byid['AmBe_housing_body']
    housing['mother'] = 'superSANDIDetMed'
    # Mother frame is vessel-local; user position is in detector/tank frame.
    housing['position'] = [x-y for x,y in zip(p.source_position_mm, vessel['nylon_vessel_position'])]
    byid['AmBe_BGO']['r_min'] = p.bgo_bore_radius_mm
    # A solid steel cup envelope with AIR pocket daughter displaces steel legally.
    # Cup spans z=[-5,5] in housing frame; recess opens at +5, towards shield/BGO.
    depth = p.disk_thickness_mm + 2*p.axial_clearance_mm
    def tube(index, mother, material, radius, half, z):
        return dict(name='GEO', index=index, enable=1, valid_begin=[0,0], valid_end=[0,0],
                    mother=mother, type='tube', material=material, r_min=0., r_max=radius,
                    size_z=half, position=[0.,0.,z])
    pocket = tube('Co60_pocket', 'AmBe_placer', 'air', 12.7+p.radial_clearance_mm,
                  depth/2, 5-depth/2)
    disk = tube('Co60_disk', 'Co60_pocket', p.plastic_material, 12.7,
                p.disk_thickness_mm/2, 0.)
    active = tube('Co60_active', 'Co60_disk', p.plastic_material, p.active_radius_mm,
                  p.active_thickness_mm/2, p.active_z_mm)
    assembly.extend([pocket, disk, active])
    # Keep references as supplied for comparison; origin remains original cup center.
    # Active center is explicitly exported; do not mistake housing origin for it.
    out.mkdir(parents=True, exist_ok=True)
    (out/'ANNIE_config14_SuperSANDI_Co60.geo').write_text(dump_blocks(detector))
    (out/'ANNIE_Co60_holder.geo').write_text(dump_blocks(assembly))
    ordered, pending = [], assembly[:]
    known = {'superSANDIDetMed'}
    while pending:
        ready = [b for b in pending if b['mother'] in known]
        if not ready:
            raise ValueError('invalid or cyclic mother hierarchy')
        for b in ready:
            ordered.append(b); known.add(b['index']); pending.remove(b)
    # Shared loader input: name mother material type x y z N z/rmin/rmax triples.
    lines = []
    for b in ordered:
        if b['type'] == 'tube':
            zs = [-b['size_z'], b['size_z']]
            ri, ro = [b.get('r_min',0)]*2, [b['r_max']]*2
        elif b['type'] == 'revolve':
            zs, ri, ro = b['z'], b['r_min'], b['r_max']
        else:
            raise ValueError('unsupported assembly solid')
        fields = [b['index'], b['mother'], b['material'], b['type'], b.get('surface','-'),
                  *b.get('position',[0,0,0]), len(zs)]
        fields.extend(v for triple in zip(zs,ri,ro) for v in triple)
        lines.append(' '.join(map(str,fields)))
    (out/'assembly.tsv').write_text('\n'.join(lines)+'\n')
    optics_path = base/'OPTICS_SANDI_wbls_Nov22.ratdb'
    optics = next(b for b in blocks(optics_path) if b['index']==vessel['detection_medium_material'])
    assert optics['LIGHT_YIELD'] == 110.0, 'unexpected supplied tuning; stop'
    cup_center_z = p.source_position_mm[2]
    active_center_z = cup_center_z + 5-depth/2+p.active_z_mm
    manifest = dict(status='PROVISIONAL P0 GEOMETRY ONLY; no optical validation',
                    parameters=vars(p), vessel=vessel,
                    vessel_shape=dict(height_mm=2800., radius_parameter_mm=900.,
                                      nylon_thickness_mm=.025, gradient_mm=.00025, faces=8),
                    static_optics=dict(material=optics['index'], LIGHT_YIELD=optics['LIGHT_YIELD'],
                                       sha256=hashlib.sha256(optics_path.read_bytes()).hexdigest(),
                                       resolved_in_RAT=False),
                    baseline_bgo_r_min_mm=by_baseline(original,'AmBe_BGO')['r_min'],
                    cup_bottom_mm=10-depth, cup_wall_mm=25-12.7-p.radial_clearance_mm,
                    active_center_detector_mm=[*p.source_position_mm[:2], active_center_z],
                    volumes=ordered)
    (out/'geometry_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def by_baseline(records, index):
    return next(b for b in records if b['index']==index)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--app', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    for name in ('disk-thickness-mm','radial-clearance-mm','axial-clearance-mm',
                 'active-radius-mm','active-thickness-mm','active-z-mm','bgo-bore-radius-mm'):
        parser.add_argument('--'+name, type=float, required=True)
    parser.add_argument('--plastic-material', required=True, help='provisional RAT material name')
    parser.add_argument('--source-position-mm', nargs=3, type=float, required=True,
                        help='cup center in detector frame; +z points to BGO/PMT')
    args = parser.parse_args()
    app, out = args.app, args.out
    del args.app, args.out
    prepare(app, out, args)
