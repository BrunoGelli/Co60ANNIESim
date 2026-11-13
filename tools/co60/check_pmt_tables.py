#!/usr/bin/env python3
"""Compare archived versus regenerated PMTINFO without applying extra offsets.

Raw-coordinate bounds are diagnostic only: the missing tuned factory must resolve
any global offsets/rotations before asserting PMT placements are valid.
"""
import json
import math
from pathlib import Path
import sys
from ratdb_subset import blocks

original, regenerated=map(lambda p:{b['name']:b for b in blocks(Path(p))},sys.argv[1:3])
active={b['pos_table'] for b in blocks(Path(sys.argv[3]))
        if b.get('type')=='pmtarray' and b.get('enable',1)}
report={'status':'STATIC_ONLY; full factory alignment UNVERIFIED','tables':[],'samples':[]}
all_ids=[]
for name,b in original.items():
    if name not in active: continue
    other=regenerated[name]
    ids=b['channel_number']; all_ids.extend(ids)
    assert ids==other['channel_number']
    for k in ('x','y','z','dir_x','dir_y','dir_z'):
        assert len(b[k])==len(ids)
        assert all(math.isfinite(x) for x in b[k])
    delta=max(abs(x-y) for k in ('x','y','z') for x,y in zip(b[k],other[k]))
    dir_error=max(abs(math.sqrt(sum(b[k][i]**2 for k in ('dir_x','dir_y','dir_z')))-1)
                  for i in range(len(ids)))
    report['tables'].append(dict(name=name,n_channels=len(ids),has_efficiency_corr='efficiency_corr' in b,
        max_position_delta_mm=delta,
        max_direction_norm_error=dir_error,raw_z_range_mm=[min(b['z']),max(b['z'])]))
    if any(s in name for s in ('bottomgrid','topgrid','side_ring_1')):
        report['samples'].append(dict(name=name,channel=ids[0],
            raw_position_mm=[b[k][0] for k in ('x','y','z')],
            direction=[b[k][0] for k in ('dir_x','dir_y','dir_z')]))
assert len(all_ids)==len(set(all_ids))
report['n_unique_channels']=len(all_ids)
print(json.dumps(report,indent=2))
