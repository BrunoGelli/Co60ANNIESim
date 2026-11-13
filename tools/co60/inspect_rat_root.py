#!/usr/bin/env python3
"""Inspect an ACTUAL ANNIE flat ROOT run; no simulated results are supplied here.

Report values and multiplicity. Legacy output rows count EVs, not decays; compare
meta.dsentries to the generated MC count. This cannot establish splitter semantics
without a full outroot cross-check from the matching core.
"""
import argparse
import json
import ROOT

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('file')
p.add_argument('--expected-mc',type=int,required=True)
p.add_argument('--require-pe',action='store_true')
a=p.parse_args()
f=ROOT.TFile.Open(a.file,'READ')
if not f or f.IsZombie(): raise OSError(a.file)
report={'file':a.file,'keys':[{ 'name':k.GetName(),'class':k.GetClassName()} for k in f.GetListOfKeys()]}
for name in ('meta','output'):
    tree=f.Get(name)
    if not tree or not tree.InheritsFrom('TTree'): raise ValueError('missing '+name)
    report[name]={'entries':tree.GetEntries(),'branches':[b.GetName() for b in tree.GetListOfBranches()]}
meta=f.Get('meta')
counts=[int(r.dsentries) for r in meta]
report['dsentries_values']=counts
if counts != [a.expected_mc]: raise ValueError('MC count mismatch: '+str(counts))
output=f.Get('output')
samples=[]; total_pe=0; untriggered=0
for row in output:
    npe=int(row.mcpecount); total_pe+=npe
    untriggered+=int(row.evid==-1)
    if len(samples)<5:
        samples.append(dict(evid=int(row.evid),subev=int(row.subev),mcpecount=npe,
                            mcnhits=int(row.mcnhits),mcke=float(row.mcke),
                            pmt_ids=list(row.mcPMTID) if hasattr(row,'mcPMTID') else [],
                            pe_times=list(row.mcPEHitTime)[:10] if hasattr(row,'mcPEHitTime') else []))
report.update(sample_rows=samples,sum_row_mcpecount=total_pe,untriggered_rows=untriggered,
              warning='MC hits repeat across EV rows; sum_row_mcpecount is NOT an all-decay total')
if a.require_pe and total_pe==0: raise ValueError('no detected PE in this run')
print(json.dumps(report,indent=2))
f.Close()
