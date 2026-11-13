#!/usr/bin/env python3
"""Export PRIVATE transport inputs from the accepted P0 model and uploaded materials.

Optics are never modified/loaded. Missing tuned-core material definitions get
explicit provisional fallbacks; they are recorded, not represented as exact RAT.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from ratdb_subset import blocks

def main(a):
    out=a.out.resolve(); repo=Path(__file__).resolve().parents[2]
    if out.exists() or out.is_relative_to(repo): raise ValueError('use NEW private output directory')
    out.mkdir(parents=True)
    db=a.app/'ratdb/ANNIE'; g=a.p0/'geometry'
    manifest=json.loads((g/'geometry_manifest.json').read_text())
    for name in ('assembly.tsv','vessel.txt','probes.tsv','geometry_manifest.json'):
        shutil.copy2(g/name,out/name)
    shutil.copy2(a.visualization/'envelopes.tsv',out/'envelopes.tsv')
    tables={}; provenance={}; input_files=list(g.glob('*'))+[a.visualization/'envelopes.tsv']
    for name in ('MATERIALS_ANNIE.ratdb','MATERIALS_SANDI_wbls_Nov22.ratdb','MATERIALS_SuperSANDI.ratdb'):
        path=db/name; input_files.append(path)
        for b in blocks(path):
            tables[b['index']]=b; provenance[b['index']]=str(path)
    # NIST-based material FALLBACKS for absent tuned-core definitions, no silent placeholders.
    fallback={
      'air':('G4_AIR',0.),'water':('G4_WATER',1.),'glass':('G4_Pyrex_Glass',0.),
      'stainless_steel':('G4_STAINLESS-STEEL',0.),'pmt_vacuum':('G4_Galactic',0.),
      'acrylic_black':('G4_PLEXIGLASS',0.),'delrin_black':('G4_POLYOXYMETHYLENE',0.),
      'polystyrene_foam':('G4_POLYSTYRENE',a.foam_density_g_cm3)}
    # PPO C15H11NO; explicit chemical-formula fallback, negligible mass fraction but recorded.
    atom={'C':12.0107,'H':1.00794,'N':14.0067,'O':15.9994}
    mass=sum(atom[k]*n for k,n in [('C',15),('H',11),('N',1),('O',1)])
    ppo={k:atom[k]*n/mass for k,n in [('C',15),('H',11),('N',1),('O',1)]}
    symbols=dict(Hydrogen='H',Carbon='C',Oxygen='O',Nitrogen='N',Gadolinium='Gd',
                 Sulphur='S',Bismuth='Bi',Germanium='Ge',Silicon='Si',Chlorine='Cl')
    used=[]; lines=[]; built=set()
    def build(name):
        if name in built: return
        if name in tables:
            b=tables[name]; ne=b['nelements']; nm=b['nmaterials']
            if ne and nm: raise ValueError('mixed elemental/material table needs explicit support')
            pairs=list(zip(b['elements'],b['elemprop'])) if ne else list(zip(b['materials'],b['matprop']))
            if not ne:
                for sub,_ in pairs: build(sub)
            else: pairs=[(symbols[k],v) for k,v in pairs]
            if abs(sum(v for _,v in pairs)-1)>.001: raise ValueError('bad fraction sum '+name)
            lines.append(' '.join(map(str,['E' if ne else 'M',name,b['density'],len(pairs),
                                          *[v for pair in pairs for v in pair]])))
            used.append(dict(name=name,status='uploaded density/composition',source=provenance[name],table=b))
        elif name=='ppo':
            lines.append(' '.join(map(str,['E',name,1.06,len(ppo),*[v for pair in ppo.items() for v in pair]])))
            used.append(dict(name=name,status='PROVISIONAL missing-core fallback',formula='C15H11NO',density_g_cm3=1.06))
        elif name in fallback:
            nist,density=fallback[name]
            lines.append(f'B {name} {density} {nist}')
            used.append(dict(name=name,status='PROVISIONAL missing-core fallback',nist=nist,
                             density_override_g_cm3=density or None))
        else: raise ValueError('unresolved transport material '+name)
        built.add(name)
    materials={b['material'] for b in manifest['volumes']}|{'water_gdS_0p2','superSANDI_nylon',
                    manifest['vessel']['detection_medium_material'],'air','acrylic_black','stainless_steel'}
    for name in sorted(materials): build(name)
    (out/'materials.tsv').write_text('\n'.join(lines)+'\n')
    positions={'superSANDIDetMed':manifest['vessel']['nylon_vessel_position']}
    for b in manifest['volumes']: positions[b['index']]=[x+y for x,y in zip(positions[b['mother']],b['position'])]
    active=next(b for b in manifest['volumes'] if b['index']=='Co60_active')
    (out/'source.tsv').write_text(' '.join(map(str,[*positions['Co60_active'],active['r_max'],active['size_z']]))+'\n')
    input_files=[p for p in input_files if p.is_file()]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    data=dict(status='PROVISIONAL standalone P1; no optics; missing core explicitly approximated',
        material_definitions=used,source_active_center_mm=positions['Co60_active'],
        active_radius_mm=active['r_max'],active_half_length_mm=active['size_z'],
        geometry_parameters=manifest['parameters'],input_sha256={str(p.resolve()):sha(p) for p in input_files},
        omissions=['tank PMT solids','internal supports/black sheets','hall; escape scored at tank outside boundary'],
        material_semantics='elemprop and matprop interpreted as mass fractions; provided values preserved')
    datasets={r.split()[0]:dict(environment=r.split()[1],path=r.split()[2])
              for r in subprocess.check_output(['geant4-config','--datasets'],text=True).splitlines()}
    for key in ('RadioactiveDecay','PhotonEvaporation','G4EMLOW','G4ENSDFSTATE'):
        if not Path(datasets[key]['path']).is_dir(): raise ValueError('missing dataset '+key)
    nuclear_files=[Path(datasets['RadioactiveDecay']['path'])/'z27.a60',
                   Path(datasets['PhotonEvaporation']['path'])/'z28.a60',
                   Path(datasets['G4ENSDFSTATE']['path'])/'ENSDFSTATE.dat']
    data['datasets']=datasets
    data['nuclear_data_sha256']={str(p):sha(p) for p in nuclear_files}
    data['software']={cmd:subprocess.check_output([cmd,'--version'],text=True).splitlines()[0]
                      for cmd in ('geant4-config','root-config','c++','cmake')}
    (out/'transport_provenance.json').write_text(json.dumps(data,indent=2)+'\n')
    print('Private transport inputs:',out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('app','p0','visualization','out'): p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--foam-density-g-cm3',type=float,required=True,help='explicit unmeasured fallback')
    main(p.parse_args())
