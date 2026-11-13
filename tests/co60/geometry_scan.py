#!/usr/bin/env python3
"""Nine explicitly provisional geometry cases, plus an overlap rejection check.

No transport or optical simulation. Private outputs must stay outside public repo.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--app',type=Path,required=True)
p.add_argument('--probe',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
repo=Path(__file__).resolve().parents[2]
assert not a.out.resolve().is_relative_to(repo)
results=[]
for thickness in (1.,2.,3.):
    for clearance in (.1,.25,.5):
        case=a.out/f'PROVISIONAL_t{thickness}_c{clearance}'
        cmd=[sys.executable,str(repo/'tools/co60/prepare_geometry.py'),'--app',str(a.app),
             '--out',str(case),'--disk-thickness-mm',str(thickness),
             '--radial-clearance-mm',str(clearance),'--axial-clearance-mm','.1',
             '--active-radius-mm','1','--active-thickness-mm','.01','--active-z-mm','0',
             '--bgo-bore-radius-mm','1','--plastic-material','polystyrene',
             '--source-position-mm','0','750','0']
        subprocess.run(cmd,check=True)
        subprocess.run([sys.executable,str(repo/'tools/co60/probe_inputs.py'),str(case)],check=True)
        run=[str(a.probe),str(case/'assembly.tsv'),str(case/'vessel.txt'),str(case/'probes.tsv')]
        r=subprocess.run(run,capture_output=True,text=True)
        (case/'navigation.log').write_text(r.stdout+r.stderr)
        if r.returncode: raise RuntimeError(str(case))
        results.append(dict(case=case.name,status='PASS',kind='GEOMETRY_ONLY',
                            command=run,prepare_command=cmd))
# Enlarge BGO in one copied manifest so it intersects PMT/steel. Must fail overlap check.
case=a.out/'PROVISIONAL_t2.0_c0.25'
lines=(case/'assembly.tsv').read_text().splitlines()
for i,line in enumerate(lines):
    fields=line.split()
    if fields[0]=='AmBe_BGO':
        fields[9]='-26.0'; fields[12]='26.0'
        lines[i]=' '.join(fields)
bad=case/'INTENTIONALLY_OVERLAPPING.tsv'
bad.write_text('\n'.join(lines)+'\n')
r=subprocess.run([str(a.probe),str(bad),str(case/'vessel.txt'),str(case/'probes.tsv')],
                 capture_output=True,text=True)
(case/'overlap_rejection.log').write_text(r.stdout+r.stderr)
assert r.returncode==1 and 'Overlap is detected' in r.stdout+r.stderr
summary=dict(status='PASS',kind='GEOMETRY_ONLY_PROVISIONAL_SCAN',cases=results,
             overlap_rejection='PASS')
(a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS: nine provisional geometry cases; intentional overlap rejected')
