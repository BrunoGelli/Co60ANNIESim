#!/usr/bin/env python3
"""Synthetic schema I/O only. Three fabricated records, NO detector simulation."""
import copy
import json
from pathlib import Path
import sys
import uuid
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/co60'))
from root_schema import SCHEMA,write,read,validate

def blank(table, **values):
    row={k: '' if t=='S' else 0. for k,t in SCHEMA[table].items()}
    for k,t in SCHEMA[table].items():
        if t=='L': row[k]=0
    row.update(values)
    return row

def main(out):
    run_uuid=str(uuid.uuid4())
    rows={k:[] for k in SCHEMA}
    rows['runs']=[blank('runs',run_uuid=run_uuid,schema_version=1,generator_mode='SYNTHETIC_IO_TEST',
        is_synthetic=1,source_activity_Bq=13061.,seed=6012,n_generated=3,
        geometry_sha256='UNVALIDATED',optics_sha256='UNVALIDATED',
        provenance_json=json.dumps({'transport':'NONE','purpose':'schema smoke'}),
        units_json=json.dumps({'length':'mm','time':'ns prompt-relative','energy':'MeV','charge':'PE'}))]
    rows['events']=[blank('events',run_uuid=run_uuid,decay_id=i,rat_mc_id=100+i,
                         n_rat_evs=[0,2,0][i],weight=1.,n_detected_pe=int(i==1)) for i in range(3)]
    # Synthetic split-trigger case and an untriggered case with deposit; both retained.
    rows['events'][1]['edep_wbls_MeV']=.25
    rows['events'][2]['edep_steel_MeV']=.01
    rows['volumes']=[blank('volumes',run_uuid=run_uuid,volume_id=1,name='fixture_wbls',
                          material='fixture',region='wbls',mother='fixture_world'),
                     blank('volumes',run_uuid=run_uuid,volume_id=2,name='fixture_steel',
                          material='fixture',region='steel',mother='fixture_world')]
    rows['channels']=[blank('channels',run_uuid=run_uuid,pmt_id=4,channel_id=401,role='tank')]
    rows['gamma_transport']=[blank('gamma_transport',run_uuid=run_uuid,decay_id=1,
                                   root_gamma_id=7,track_id=9,initial_MeV=1.173,dir_z=1.,
                                   first_interaction_valid=1,first_volume_id=1,first_process='fixture')]
    rows['depositions']=[blank('depositions',run_uuid=run_uuid,decay_id=i,deposit_id=0,
        volume_id=1 if i==1 else 2,material='fixture',track_id=11,parent_id=9,
        root_gamma_id=7 if i==1 else -1,pdg=11,process='fixture',
        edep_MeV=.25 if i==1 else .01) for i in (1,2)]
    rows['pmt_pe']=[blank('pmt_pe',run_uuid=run_uuid,decay_id=1,pe_id=0,pmt_id=4,
        channel_id=401,role='tank',hit_time_ns=12.,front_end_time_ns=13.,charge_PE=1.,
        wavelength_nm=400.,process='fixture',root_gamma_id=-1,creation_volume_id=-1)]
    write(out,rows)
    loaded=read(out)
    assert loaded==rows
    for mutation in ('lost_zero','duplicate_event','orphan_pe','nan','duplicate_pe','orphan_gamma'):
        bad=copy.deepcopy(rows)
        if mutation=='lost_zero': bad['events'].pop(0)
        elif mutation=='duplicate_event': bad['events'][2]['decay_id']=1
        elif mutation=='orphan_pe': bad['pmt_pe'][0]['channel_id']=99
        elif mutation=='nan': bad['events'][0]['source_x_mm']=float('nan')
        elif mutation=='duplicate_pe': bad['pmt_pe'].append(copy.deepcopy(bad['pmt_pe'][0]))
        elif mutation=='orphan_gamma': bad['depositions'][0]['root_gamma_id']=999
        try: validate(bad)
        except ValueError: pass
        else: raise AssertionError('accepted '+mutation)
    print(json.dumps({'status':'PASS','kind':'SYNTHETIC_IO_ONLY',
                      'entries':{k:len(v) for k,v in loaded.items()},
                      'decay_ids':[r['decay_id'] for r in loaded['events']],
                      'n_rat_evs':[r['n_rat_evs'] for r in loaded['events']],
                      'pe_time_ns':loaded['pmt_pe'][0]['hit_time_ns'],
                      'negative_cases':6},indent=2))

if __name__=='__main__': main(Path(sys.argv[1]))
