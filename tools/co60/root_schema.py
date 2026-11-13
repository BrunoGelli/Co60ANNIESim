"""Native flat ROOT P0 schema. Transport adapter is deferred until tuned core access.

All numeric columns are scalar leaves; string columns are std::string. No RAT DS
or ToolAnalysis dictionaries. No ADC waveforms, QE filter, or physics generator.
"""
from array import array
import math
from pathlib import Path

# D=double, L=signed 64 bit, S=std::string; names carry physical units.
SCHEMA = {
 'runs': dict(run_uuid='S',schema_version='L',generator_mode='S',is_synthetic='L',
              source_activity_Bq='D',seed='L',n_generated='L',geometry_sha256='S',
              optics_sha256='S',provenance_json='S',units_json='S'),
 'events': dict(run_uuid='S',decay_id='L',rat_mc_id='L',n_rat_evs='L',weight='D',
                source_x_mm='D',source_y_mm='D',source_z_mm='D',
                emitted_visible_MeV='D',edep_bgo_MeV='D',edep_wbls_MeV='D',
                edep_other_water_MeV='D',edep_steel_MeV='D',edep_plastic_MeV='D',
                edep_other_MeV='D',escaped_MeV='D',residual_MeV='D',
                n_detected_pe='L',truth_valid='L',ancestry_valid='L'),
 'depositions': dict(run_uuid='S',decay_id='L',deposit_id='L',volume_id='L',
                     material='S',track_id='L',parent_id='L',root_gamma_id='L',
                     pdg='L',process='S',x_mm='D',y_mm='D',z_mm='D',time_ns='D',edep_MeV='D'),
 'gamma_transport': dict(run_uuid='S',decay_id='L',root_gamma_id='L',track_id='L',
                         parent_id='L',initial_MeV='D',dir_x='D',dir_y='D',dir_z='D',
                         first_interaction_valid='L',first_volume_id='L',first_process='S',
                         first_x_mm='D',first_y_mm='D',first_z_mm='D',first_time_ns='D',
                         edep_bgo_MeV='D',edep_wbls_MeV='D',edep_other_MeV='D',escaped_MeV='D'),
 'pmt_pe': dict(run_uuid='S',decay_id='L',pe_id='L',pmt_id='L',channel_id='L',role='S',
                hit_time_ns='D',front_end_time_ns='D',charge_PE='D',wavelength_nm='D',
                process='S',root_gamma_id='L',creation_volume_id='L',ancestry_valid='L'),
 'volumes': dict(run_uuid='S',volume_id='L',name='S',material='S',region='S',mother='S'),
 'channels': dict(run_uuid='S',pmt_id='L',channel_id='L',role='S',x_mm='D',y_mm='D',z_mm='D')
}


def validate(rows):
    if len(rows['runs']) != 1:
        raise ValueError('one runs row required per file/chunk')
    run=rows['runs'][0]; uuid=run['run_uuid']
    events=rows['events']
    ids=[r['decay_id'] for r in events]
    if len(ids)!=run['n_generated'] or sorted(ids)!=list(range(run['n_generated'])):
        raise ValueError('every generated decay must have exactly one events row')
    for table, records in rows.items():
        for row in records:
            if set(row)!=set(SCHEMA[table]):
                raise ValueError('missing or extra columns in '+table)
            if row['run_uuid']!=uuid:
                raise ValueError('run UUID mismatch')
            if 'decay_id' in row and row['decay_id'] not in ids:
                raise ValueError('orphan decay foreign key')
            for name,typ in SCHEMA[table].items():
                if typ=='D' and not math.isfinite(row[name]):
                    raise ValueError('nonfinite '+name)
    for table,keys in [('depositions',('decay_id','deposit_id')),
                       ('gamma_transport',('decay_id','root_gamma_id')),
                       ('pmt_pe',('decay_id','pe_id')),('volumes',('volume_id',)),
                       ('channels',('pmt_id',))]:
        values=[tuple(r[k] for k in keys) for r in rows[table]]
        if len(set(values))!=len(values):
            raise ValueError('duplicate primary key in '+table)
    volumes={r['volume_id'] for r in rows['volumes']}
    channels={(r['pmt_id'],r['channel_id'],r['role']) for r in rows['channels']}
    gammas={(r['decay_id'],r['root_gamma_id']) for r in rows['gamma_transport']}
    for table in ('depositions','pmt_pe'):
        for r in rows[table]:
            if r['root_gamma_id']>=0 and (r['decay_id'],r['root_gamma_id']) not in gammas:
                raise ValueError('orphan gamma lineage')
    for r in rows['depositions']:
        if r['volume_id'] not in volumes or r['edep_MeV']<0:
            raise ValueError('invalid deposition')
    for r in rows['pmt_pe']:
        if (r['pmt_id'],r['channel_id'],r['role']) not in channels:
            raise ValueError('orphan PMT channel')
    for event in events:
        if event['n_detected_pe']!=sum(r['decay_id']==event['decay_id'] for r in rows['pmt_pe']):
            raise ValueError('PE count mismatch')


def write(path, rows):
    import ROOT
    validate(rows)
    if Path(path).exists():
        raise FileExistsError(path)
    root=ROOT.TFile.Open(str(path),'CREATE')
    if not root or root.IsZombie():
        raise OSError(path)
    try:
        for table,columns in SCHEMA.items():
            tree=ROOT.TTree(table,table)
            buffers={}
            for name,typ in columns.items():
                buf=ROOT.std.string() if typ=='S' else array('d' if typ=='D' else 'q',[0])
                buffers[name]=buf
                if typ=='S': tree.Branch(name,buf)
                else: tree.Branch(name,buf,name+'/'+typ)
            for row in rows[table]:
                for name,value in row.items():
                    if columns[name]=='S': buffers[name].assign(value)
                    else: buffers[name][0]=value
                tree.Fill()
            tree.Write()
    finally:
        root.Close()


def read(path):
    import ROOT
    root=ROOT.TFile.Open(str(path),'READ')
    if not root or root.IsZombie():
        raise OSError(path)
    try:
        rows={}
        for table,columns in SCHEMA.items():
            tree=root.Get(table)
            if not tree or not tree.InheritsFrom('TTree'):
                raise ValueError('missing TTree '+table)
            if {b.GetName() for b in tree.GetListOfBranches()}!=set(columns):
                raise ValueError('branch mismatch '+table)
            rows[table]=[{name: (str(getattr(r,name)) if typ=='S' else
                         float(getattr(r,name)) if typ=='D' else int(getattr(r,name)))
                         for name,typ in columns.items()} for r in tree]
    finally:
        root.Close()
    validate(rows)
    return rows
