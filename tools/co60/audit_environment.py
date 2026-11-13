#!/usr/bin/env python3
"""Read-only audit. Writes private provenance outside the public repo; no full env dump."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from datetime import datetime, timezone

def capture(argv):
    try:
        r=subprocess.run(argv,capture_output=True,text=True,timeout=30)
        return dict(command=argv,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr)
    except (OSError,subprocess.TimeoutExpired) as e:
        return dict(command=argv,error=str(e))

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for buf in iter(lambda:f.read(1024*1024),b''): h.update(buf)
    return h.hexdigest()

def repository(path):
    commands={name:capture(['git','-C',str(path),*args]) for name,args in {
        'head':['rev-parse','HEAD'],'status':['status','--porcelain=v1','-uall'],
        'branch':['branch','--show-current'],'diff':['diff','--binary']}.items()}
    names=capture(['git','-C',str(path),'ls-files','-co','--exclude-standard'])['stdout'].splitlines()
    hashes={name:sha(path/name) for name in sorted(set(names)) if (path/name).is_file()}
    return dict(path=str(path.resolve()),git=commands,sha256=hashes)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--app',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--public-repo',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.resolve().is_relative_to(a.public_repo.resolve()):
        raise ValueError('audit includes private source diffs; write outside public repo')
    commands=[['c++','--version'],['cmake','--version'],['root-config','--version'],
              ['geant4-config','--version'],['root-config','--features'],['geant4-config','--datasets']]
    paths=['/exp/annie/app/users/augusthy/RATPAC2/TunedRATPAC',
           '/pnfs/annie/scratch/users/augusthy/input','/cvmfs','/ratpac-setup']
    result=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
        archive=dict(path=str(a.archive.resolve()),sha256=sha(a.archive)),
        baseline=repository(a.baseline),app=repository(a.app),public=repository(a.public_repo),
        tools={t:shutil.which(t) for t in ('singularity','apptainer','annie','root-config','geant4-config')},
        runtime_paths={s:dict(exists=Path(s).exists(),readable=os.access(s,os.R_OK)) for s in paths},
        versions=[capture(c) for c in commands],
        geant4_data_env={k:v for k,v in os.environ.items() if k.startswith('G4') and k.endswith('DATA')},
        shipped_cache=(a.baseline/'build/CMakeCache.txt').read_text(),
        shipped_executable_sha256=sha(a.baseline/'install/bin/annie'),
        shipped_ldd=capture(['ldd',str(a.baseline/'install/bin/annie')]))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print('Private audit written:',a.out)
