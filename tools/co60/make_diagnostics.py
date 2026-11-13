#!/usr/bin/env python3
"""Private local P0 plots: Geant4 sections, static RATDB, synthetic ROOT I/O.

No transport, optical response, runtime factory alignment, or public upload.
All figures/data go outside this public repository to a NEW output directory.
"""
import argparse
import csv
import hashlib
import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch
import numpy as np

from ratdb_subset import blocks
from root_schema import read as read_root

REPO=Path(__file__).resolve().parents[2]
COLORS={'water_gdS_0p2':'#b8d8f1','wbls1pct_ly95_gdS0p2pct_SANDI':'#76d6ce',
        'superSANDI_nylon':'#993ac7','stainless_steel':'#76818c','pvc_white':'#e1d9c6',
        'air':'#fafafa','BGO_scint':'#326dc9','glass':'#a0d7b2','pmt_vacuum':'#e9f2eb',
        'polystyrene_foam':'#3d454e','delrin_black':'#282f37','polystyrene':'#f2ab52',
        'acrylic_black':'#32323c'}
PROVISIONAL='#c06408'
LOCKED='#006b65'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_save(path, data):
    path.write_text(json.dumps(data,indent=2)+'\n')


class Suite:
    def __init__(self,a):
        self.a=a
        self.out=a.out.resolve()
        if self.out.is_relative_to(REPO) or self.out.exists():
            raise ValueError('choose a NEW private output directory outside public repository')
        self.out.mkdir(parents=True)
        self.figures=[]
        self.sources={}
        self.results={}
        self.p0=a.p0.resolve()
        self.app=a.app.resolve()
        self.g=self.p0/'geometry'
        self.m=json.loads((self.g/'geometry_manifest.json').read_text())
        self.vols={b['index']:b for b in self.m['volumes']}
        self.centers={'superSANDIDetMed':np.array(self.m['vessel']['nylon_vessel_position'])}
        for b in self.m['volumes']:
            self.centers[b['index']]=self.centers[b['mother']]+b['position']
        self.db=self.app/'ratdb/ANNIE'
        self.detector=blocks(self.db/'ANNIE_run2_config14_scan.geo')
        self.inputs=[self.g/'geometry_manifest.json',self.g/'assembly.tsv',self.g/'vessel.txt',
                     self.g/'probes.tsv',self.db/'ANNIE_run2_config14_scan.geo',
                     self.db/'ANNIE_AmBesource_v2.geo',self.app/'src/GeoANNIEFactory.cc',
                     self.db/'OPTICS_SANDI_wbls_Nov22.ratdb',
                     self.db/'PMTINFO_run2_config14_scan.ratdb']
        self.inputs.extend(self.p0.glob('*.root'))
        self.original_hashes={str(p):sha(p) for p in self.inputs}
        plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                             'savefig.dpi':170,'pdf.fonttype':42})

    def run(self,cmd,log):
        with (self.out/log).open('w') as f:
            f.write('$ '+subprocess.list2cmdline(list(map(str,cmd)))+'\n'); f.flush()
            subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True)

    def save(self,fig,name,caption):
        fig.savefig(self.out/(name+'.png'),bbox_inches='tight')
        fig.savefig(self.out/(name+'.pdf'),bbox_inches='tight')
        plt.close(fig)
        self.figures.append(dict(name=name,caption=caption))

    def geometry(self):
        # Assert private source constants agree with the exported vessel geometry.
        source=(self.app/'src/GeoANNIEFactory.cc').read_text()
        for text in ('nv_height = 2.8 * CLHEP::m','vessel_radius = 0.9 * CLHEP::m',
                     'nylon_thickness = 25 * CLHEP::um','gradient = 0.25 * CLHEP::um',
                     'G4SubtractionSolid("detectionMedium",nv_outer,nylonVessel)'):
            if text not in source: raise ValueError('re-audit changed vessel construction: '+text)
        if self.m['vessel_shape']!=dict(height_mm=2800.,radius_parameter_mm=900.,
                nylon_thickness_mm=.025,gradient_mm=.00025,faces=8):
            raise ValueError('vessel manifest differs from source construction')
        baseline={b['index']:b for b in blocks(self.db/'ANNIE_AmBesource_v2.geo')}
        for name in ('AmBe_BGO','AmBe_steel_disk','AmBe_placer','AmBe_PMT_glass',
                     'AmBe_PMT_vacuum','AmBe_housing_body','AmBe_housing_air'):
            for key in ('type','r_max','size_z','z','material'):
                if key in baseline[name] and self.vols[name][key]!=baseline[name][key]:
                    raise ValueError('unexpected change in uploaded dimensions: '+name+' '+key)
        # Export actual nested tank/liner/detector tube dimensions, in detector frame.
        # Hall/world and runtime-built PMTs/structure are intentionally NOT substituted.
        lines=[]
        dims={b['index']:b for b in self.detector}
        for name in ('tank','liner','detector'):
            b=dims[name]; mother='detector_frame_world' if name=='tank' else b['mother']
            if name!='tank' and b['position']!=[0,0,0]:
                raise ValueError('nonzero nested envelope offset needs explicit handling')
            lines.append(f"{name} {mother} {b['material']} tube - 0 0 0 2 "
                         f"{-b['size_z']} 0 {b['r_max']} {b['size_z']} 0 {b['r_max']}")
        (self.out/'envelopes.tsv').write_text('\n'.join(lines)+'\n')
        self.dims=dims
        c=self.centers['AmBe_housing_body']; bx=self.centers['AmBe_BGO']
        # name u-axis v-axis fixed-axis fixed umin umax vmin vmax columns
        views=[('assembly',0,2,1,c[1],c[0]-45,c[0]+45,c[2]-43,c[2]+348,901),
               ('cup',0,2,1,c[1],c[0]-30,c[0]+30,c[2]-7,c[2]+18,1201),
               ('bore',0,1,2,bx[2],bx[0]-27,bx[0]+27,bx[1]-27,bx[1]+27,1001),
               ('detector_yz',1,2,0,0,-1600,1600,-2050,2050,1001),
               ('detector_xy',0,1,2,self.m['active_center_detector_mm'][2],
                -1600,1600,-1600,1600,1001),
               ('holder_boundary',1,2,0,c[0],c[1]-100,c[1]+260,c[2]-55,c[2]+360,901),
               ('nylon_cap',0,2,1,0,-1,1,1374.94,1375.03,101)]
        self.views={row[0]:row for row in views}
        (self.out/'section_requests.tsv').write_text('\n'.join(' '.join(map(str,r)) for r in views)+'\n')
        self.run(['cmake','-S',REPO/'tests/co60','-B',self.out/'build'],'cmake.log')
        self.run(['cmake','--build',self.out/'build','--target','geometry_sections','-j2'],'build.log')
        self.run([self.out/'build/geometry_sections',self.g/'assembly.tsv',self.g/'vessel.txt',
                  self.out/'envelopes.tsv',self.out/'section_requests.tsv',self.out/'sections.csv'],
                 'sections.log')
        bound_line=next(line for line in (self.out/'sections.log').read_text().splitlines()
                        if line.startswith('VESSEL_OUTER_BOUNDS_LOCAL_MM '))
        bounds=list(map(float,bound_line.split()[1:]))
        self.sections={v[0]:[] for v in views}
        with (self.out/'sections.csv').open() as f:
            for row in csv.DictReader(f):
                for k in ('u0_mm','u1_mm','v0_mm','v1_mm'): row[k]=float(row[k])
                if row['v1_mm']<=row['v0_mm']: raise ValueError('nonpositive section interval')
                self.sections[row['section']].append(row)
        # Validate representative points using actual traced boundary intervals.
        checks=[('cup',c[0],c[2]+3.9,'Co60_active'),
                ('cup',c[0]+5,c[2]+3.9,'Co60_disk'),
                ('cup',c[0]+20,c[2],'AmBe_placer'),
                ('assembly',bx[0],bx[2],'AmBe_housing_air'),
                ('assembly',bx[0]+10,bx[2],'AmBe_BGO'),
                ('nylon_cap',0,1374.9875,'superSANDIvessel')]
        # Active center is parameter-derived, not fixed to the default pocket location.
        checks[0]=(checks[0][0],self.m['active_center_detector_mm'][0],
                   self.m['active_center_detector_mm'][2],checks[0][3])
        checks[1]=(checks[1][0],c[0]+5,self.centers['Co60_disk'][2],checks[1][3])
        for view,u,v,name in checks:
            found={r['volume'] for r in self.sections[view]
                   if r['u0_mm']<=u<r['u1_mm'] and r['v0_mm']<=v<r['v1_mm']}
            if found!={name}: raise ValueError(f'section check failed: {view,u,v,name,found}')
        self.results['geometry']=dict(status='PASS',navigator_sections=len(views),
             traced_intervals=sum(map(len,self.sections.values())),point_checks=len(checks),
             vessel_outer_local_bounds_mm=bounds,
             scope='tank/liner/liquid envelopes, exact vessel CSG and manifest assembly; no tank PMTs/structure')
        json_save(self.out/'section_checks.json',self.results['geometry'])
        self.cutaway()
        self.geometry_details()
        self.detector_views()

    def draw_section(self,ax,name,local=(0,0),aspect='equal'):
        polys=[]; colors=[]
        for r in self.sections[name]:
            u0,u1=r['u0_mm']-local[0],r['u1_mm']-local[0]
            v0,v1=r['v0_mm']-local[1],r['v1_mm']-local[1]
            polys.append([(u0,v0),(u1,v0),(u1,v1),(u0,v1)])
            colors.append('#b83729' if r['volume']=='Co60_active' else COLORS[r['material']])
        # Boundary coordinates along the traced axis are Geant4 exact to stepping epsilon;
        # transverse columns are sampled. Rasterizing these shapes keeps PDFs compact.
        ax.add_collection(PolyCollection(polys,facecolors=colors,edgecolors='none',rasterized=True))
        view=self.views[name]
        ax.set(xlim=(view[5]-local[0],view[6]-local[0]),ylim=(view[7]-local[1],view[8]-local[1]))
        ax.set_aspect(aspect)

    @staticmethod
    def dim(ax,p1,p2,text,kind='baseline',offset=(0,0)):
        color={'baseline':'#202832','provisional':PROVISIONAL,'fixed':LOCKED}[kind]
        p1=np.array(p1); p2=np.array(p2); d=np.array(offset)
        ax.plot([p1[0],p1[0]+d[0]],[p1[1],p1[1]+d[1]],color=color,lw=.6)
        ax.plot([p2[0],p2[0]+d[0]],[p2[1],p2[1]+d[1]],color=color,lw=.6)
        ax.annotate('',p1+d,p2+d,arrowprops=dict(arrowstyle='<->',color=color,lw=.9,
                    linestyle='--' if kind=='provisional' else '-'))
        mid=(p1+p2)/2+d
        ax.text(mid[0],mid[1],text,ha='center',va='center',fontsize=9,color=color,
                rotation=90 if p1[0]==p2[0] else 0,
                bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1.2))

    def cutaway(self):
        m=self.m; p=m['parameters']; c=self.centers['AmBe_housing_body']
        fig=plt.figure(figsize=(15,11)); gs=fig.add_gridspec(2,3,width_ratios=[1.05,1.5,1.1],
                                                         height_ratios=[1.05,1],wspace=.42,hspace=.28)
        ax=fig.add_subplot(gs[:,0]); self.draw_section(ax,'assembly',(c[0],c[2]))
        ax.set(xlabel='Holder-local x [mm]',ylabel='Holder-local z [mm]',title='Axial cut, through center\n+z toward passive PMT')
        for name,label in [('AmBe_foam_bottom','foam pad'),('AmBe_steel_disk','steel separator'),
                ('AmBe_BGO','BGO'),('AmBe_PMT_glass','R1398 envelope'),('AmBe_PMT_base','PMT base')]:
            z=self.centers[name][2]-c[2]
            ax.annotate(label,xy=(12,z),xytext=(44,z),ha='left',fontsize=8,
                        arrowprops=dict(arrowstyle='-',color='#333',lw=.6))
        body=self.vols['AmBe_housing_body']; zmin,zmax=min(body['z']),max(body['z'])
        self.dim(ax,(-35,zmin),(-35,zmax),f'{zmax-zmin:g} mm',offset=(-5,0))
        self.dim(ax,(-35,290),(35,290),'OD 70 mm')
        self.dim(ax,(-25,270),(25,270),'bore 50 mm')
        self.dim(ax,(-25,15),(-25,65),'50 mm',offset=(-6,0))
        self.dim(ax,(-14,65),(-14,163),'98 mm',offset=(-5,0))
        ax.text(-42,344,'Tapered PVC cap\n(actual revolve)',fontsize=8,va='top')
        ax=fig.add_subplot(gs[0,1:]); self.draw_section(ax,'cup',(c[0],c[2]))
        ax.set(xlim=(-33,36),ylim=(-11,23),xlabel='Holder-local x [mm]',ylabel='Holder-local z [mm]',
               title='Recess / disk / separator: actual mother–daughter displacement')
        self.dim(ax,(-25,-5),(-25,5),'cup 10 mm',offset=(-4,0))
        self.dim(ax,(25,5),(25,15),'shield 10 mm',offset=(4,0))
        cz=self.centers['Co60_disk'][2]-c[2]
        pocket=self.vols['Co60_pocket']; pr=pocket['r_max']; depth=2*pocket['size_z']
        self.dim(ax,(-pr,5),(pr,5),f'pocket Ø{2*pr:g} mm','provisional',offset=(0,3))
        self.dim(ax,(-12.7,cz-p['disk_thickness_mm']/2),(12.7,cz-p['disk_thickness_mm']/2),
                 'disk Ø25.4 mm','fixed',offset=(0,-3.2))
        ax.annotate(f"disk t={p['disk_thickness_mm']:g} mm; radial gap={p['radial_clearance_mm']:g} mm\n"
                    f"axial gap={p['axial_clearance_mm']:g} mm each; recess depth={depth:g} mm",
                    xy=(11,cz),xytext=(0,-10),ha='center',fontsize=9,color=PROVISIONAL,
                    arrowprops=dict(arrowstyle='->',color=PROVISIONAL))
        ax.annotate(f"active region: r={p['active_radius_mm']:g} mm, t={p['active_thickness_mm']:g} mm",
                    xy=(0,self.m['active_center_detector_mm'][2]-c[2]),xytext=(0,20.7),
                    ha='center',color=PROVISIONAL,fontsize=9,arrowprops=dict(arrowstyle='->',color=PROVISIONAL))
        ax=fig.add_subplot(gs[1,1]); bx=self.centers['AmBe_BGO']
        self.draw_section(ax,'bore',(bx[0],bx[1]))
        ax.set(xlim=(-29,29),ylim=(-31,31),xlabel='BGO-local x [mm]',ylabel='BGO-local y [mm]',
               title='BGO transverse section\nAir-filled central bore')
        self.dim(ax,(-25,0),(25,0),'Ø50 mm',offset=(0,-28))
        rad=p['bgo_bore_radius_mm']
        ax.annotate(f'Provisional bore r={rad:g} mm\nØ{2*rad:g} mm',xy=(rad,0),xytext=(0,21),
                    ha='center',fontsize=10,color=PROVISIONAL,arrowprops=dict(arrowstyle='->',color=PROVISIONAL))
        ax=fig.add_subplot(gs[1,2]); ax.axis('off')
        rows=[('Uploaded model','black, solid dimensions'),('Locked disk diameter','teal, solid dimension'),
              ('New provisional choices','orange / dashed dimensions'),
              ('Steel bottom',f"{m['cup_bottom_mm']:g} mm (derived)"),
              ('Steel radial wall',f"{m['cup_wall_mm']:g} mm (derived)"),
              ('Source plastic',p['plastic_material']+' (assumed)'),
              ('Cup detector center',str(tuple(p['source_position_mm']))+' mm'),
              ('Active detector center',str(tuple(m['active_center_detector_mm']))+' mm')]
        y=.98
        for i,(label,value) in enumerate(rows):
            ax.text(0,y,label,weight='bold',fontsize=10,color=PROVISIONAL if i>=2 else '#202832')
            ax.text(0,y-.045,value,fontsize=9,wrap=True); y-=.115
        fig.suptitle('Co-60 holder cutaway • geometry only • no hardware certification',fontsize=16,y=.99)
        fig.text(.05,.035,'Solid sections: Geant4 navigator intervals from the private manifest. '
                 'Uploaded model dimensions may themselves be estimates.\n'
                 'R1398 is passive; BGO light yield remains zero. PTFE wrap is a skin surface (no finite thickness drawn). '
                 'Active layer thickness is resolved along the trace axis.',fontsize=9)
        self.save(fig,'holder_cutaway','Dimensioned manifest-derived cutaway. Black = uploaded ANNIE model; '
                  'teal = locked source diameter; orange = provisional redesign/placement. No optical response.')

    def geometry_details(self):
        c=self.centers['AmBe_housing_body']; records=[]
        for b in self.m['volumes']:
            tube=b['type']=='tube'
            zlo,zhi=(-b['size_z'],b['size_z']) if tube else (min(b['z']),max(b['z']))
            ro=[b['r_max']] if tube else b['r_max']; ri=[b.get('r_min',0)] if tube else b['r_min']
            center=self.centers[b['index']]-c
            status='provisional' if b['index'].startswith('Co60_') else 'uploaded model'
            if b['index']=='AmBe_BGO': status='uploaded; bore provisional'
            if b['index']=='AmBe_placer': status='uploaded envelope; recess new'
            records.append(dict(volume=b['index'],material=b['material'],shape=b['type'],
                envelope_OD_mm=' / '.join(f'{2*r:g}' for r in dict.fromkeys(ro)),
                intrinsic_ID_mm=' / '.join(f'{2*r:g}' for r in dict.fromkeys(ri)),
                axial_length_mm=zhi-zlo,z_min_holder_mm=center[2]+zlo,z_max_holder_mm=center[2]+zhi,
                dimension_provenance=status))
        with (self.out/'assembly_dimensions.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
        fig,ax=plt.subplots(figsize=(17,8)); ax.axis('off')
        cols=['volume','envelope_OD_mm','intrinsic_ID_mm','axial_length_mm','z_min_holder_mm',
              'z_max_holder_mm','dimension_provenance']
        labels=['Volume','Envelope OD','Intrinsic ID','Axial length','z min','z max','Provenance']
        vals=[[f'{r[k]:g}' if isinstance(r[k],(int,float)) else r[k] for k in cols] for r in records]
        table=ax.table(cellText=vals,colLabels=labels,loc='center',cellLoc='left',
                       colWidths=[.24,.10,.09,.10,.08,.08,.23])
        table.auto_set_font_size(False); table.set_fontsize(10); table.scale(1,2)
        for i,r in enumerate(records):
            for j in range(len(cols)):
                table[i+1,j].set_facecolor('#fff0d9' if 'provisional' in r['dimension_provenance']
                                         or 'new' in r['dimension_provenance'] else '#eef2f6')
        ax.set_title('Complete modeled assembly dimensions [mm] • intrinsic solids, before daughter displacement',
                     fontsize=15,pad=25)
        fig.text(.08,.06,'Envelope OD / intrinsic ID are solid parameters; e.g. the PVC body has intrinsic ID=0 '
                 'and its 50 mm air daughter makes the bore.\n'
                 'The steel cup also has intrinsic ID=0; the air-pocket daughter forms its recess. '
                 'Uploaded PMT/foam dimensions may be estimates.\n'
                 'PMT window = 1.5 mm; glass radial wall = 0.8 mm, derived from the nested uploaded glass/vacuum solids. '
                 'Surface wrapping has no geometric thickness.',fontsize=10)
        self.save(fig,'assembly_dimensions','Complete manifest dimensions and axial bounds in the holder frame. '
                  'Orange rows include new assumptions; original solids and daughter-created openings are distinguished. '
                  'The 25.4 mm source disk diameter is locked; its other dimensions are provisional.')
        fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        p=self.m['parameters']; center=self.centers['Co60_active']; disk=self.centers['Co60_disk']
        self.draw_section(axs[0],'cup',aspect='auto')
        axs[0].set(xlim=(center[0]-1.5,center[0]+1.5),ylim=(center[2]-.025,center[2]+.025),
                   title='Active-layer zoom (unequal scales)',xlabel='Detector x [mm]',ylabel='Detector z [mm]')
        axs[0].text(center[0],center[2],f"active t={p['active_thickness_mm']:g} mm",ha='center',va='center',fontsize=9)
        self.draw_section(axs[1],'cup',aspect='auto')
        edge=disk[0]+12.7
        axs[1].set(xlim=(edge-.2,edge+p['radial_clearance_mm']+.2),ylim=(disk[2]-.5,disk[2]+.5),
                   title=f"Radial air gap = {p['radial_clearance_mm']:g} mm",xlabel='Detector x [mm]',ylabel='Detector z [mm]')
        self.draw_section(axs[2],'cup',aspect='auto')
        top=disk[2]+p['disk_thickness_mm']/2
        axs[2].set(xlim=(disk[0]+2,disk[0]+3),ylim=(top-.1,top+p['axial_clearance_mm']+.1),
                   title=f"Axial air gap = {p['axial_clearance_mm']:g} mm",xlabel='Detector x [mm]',ylabel='Detector z [mm]')
        for ax in axs: ax.ticklabel_format(useOffset=False,style='plain')
        fig.suptitle('Actual navigator-boundary close-ups • provisional source geometry • unequal axis scales',fontsize=14)
        self.save(fig,'source_clearance_details','Close-ups of the existing Geant4 section data, with no enlarged solids. '
                  'Trace-axis boundaries resolve the active layer and axial gap; transverse sampling is about 0.05 mm. '
                  'Colors: red active region, orange plastic, off-white air, gray steel. No physical fit certification.')

    def detector_views(self):
        fig,axs=plt.subplots(1,2,figsize=(15,9),layout='constrained')
        for ax,name,title in zip(axs,('detector_yz','detector_xy'),
                                 ('Axial section x=0','Transverse section at active z')):
            self.draw_section(ax,name)
            ax.set(xlabel=('Detector y [mm]' if name.endswith('yz') else 'Detector x [mm]'),
                   ylabel=('Detector z [mm]' if name.endswith('yz') else 'Detector y [mm]'),title=title)
            xy=(750,140) if name.endswith('yz') else (0,750)
            ax.annotate('Co-60 holder\nprovisional placement',xy=xy,xytext=(1000,650),
                        color=PROVISIONAL,arrowprops=dict(arrowstyle='->',color=PROVISIONAL))
        for ax in axs:
            ax.legend(handles=[Patch(color=COLORS[n],label=label) for n,label in [
                ('wbls1pct_ly95_gdS0p2pct_SANDI','WbLS'),('water_gdS_0p2','Gd-water'),
                ('superSANDI_nylon','nylon (25 µm; subpixel here)'),
                ('stainless_steel','steel'),('acrylic_black','liner')]],loc='lower left',fontsize=8)
        fig.suptitle('ANNIE detector-envelope sections • actual Geant4 tube / octagonal CSG solids\n'
                     'Detector frame; tank PMTs, internal supports, black sheets and hall omitted',fontsize=14)
        self.save(fig,'detector_sections','Full detector-envelope sections, not a full runtime detector export. '
                  'Tank/liner/detector radii and lengths are from uploaded config14; nylon/WbLS use the uploaded factory CSG. '
                  'No invented circular substitute for the octagon. Thin boundaries are subpixel at this scale.')
        fig=plt.figure(figsize=(14,9)); gs=fig.add_gridspec(1,2,width_ratios=[1.5,1],wspace=.3)
        ax=fig.add_subplot(gs[0]); self.draw_section(ax,'holder_boundary')
        ax.set(xlabel='Detector y [mm]',ylabel='Detector z [mm]',
               title='Holder near octagonal WbLS boundary, x=0')
        ax.text(665,325,'WbLS',color='#065e59',fontsize=13)
        ax.text(990,255,'Gd-water',color='#245984',fontsize=11,ha='center',rotation=90)
        ax=fig.add_subplot(gs[1]); self.draw_section(ax,'nylon_cap',aspect='auto')
        ax.set(xlabel='Detector x [mm]',ylabel='Detector z [mm]',title='Nylon cap zoom (unequal axis scales)')
        ax.ticklabel_format(axis='y',style='plain',useOffset=False)
        ax.text(0,1374.958,'WbLS',ha='center'); ax.text(0,1374.988,'nylon 0.025 mm',ha='center')
        ax.text(0,1375.015,'Gd-water',ha='center')
        fig.suptitle('Boundary zooms from Geant4 navigation • stand-in material identities only',fontsize=14)
        self.save(fig,'boundary_zooms','Boundary/holder and 25 µm nylon cap sections. '
                  'The right panel intentionally uses unequal scales to resolve the actual thin layer; thickness was not enlarged in geometry.')

    def scans(self):
        summary_path=self.p0/'scans/summary.json'; self.inputs.append(summary_path)
        summary=json.loads(summary_path.read_text()); records=[]
        for record in summary['cases']:
            folder=self.p0/'scans'/record['case']
            path=folder/'geometry_manifest.json'; self.inputs.append(path)
            m=json.loads(path.read_text()); p=m['parameters']
            log=(folder/'navigation.log').read_text()
            self.inputs.append(folder/'navigation.log')
            probes=[line for line in log.splitlines() if line.startswith('PROBE ') and 'label' not in line]
            overlap='Overlap is detected' in log
            if len(probes)!=18 or any(line.split()[-1]!='1' for line in probes):
                raise ValueError('scan navigation failed '+str(folder))
            passed='GEOMETRY_ONLY PASS' in log and not overlap
            records.append(dict(case=folder.name,disk_thickness_mm=p['disk_thickness_mm'],
                radial_clearance_mm=p['radial_clearance_mm'],axial_clearance_each_mm=p['axial_clearance_mm'],
                pocket_diameter_mm=25.4+2*p['radial_clearance_mm'],
                pocket_depth_mm=p['disk_thickness_mm']+2*p['axial_clearance_mm'],
                steel_bottom_mm=m['cup_bottom_mm'],steel_wall_mm=m['cup_wall_mm'],
                disk_edge_to_air_bore_mm=25-12.7,probes_passed=len(probes),overlap_detected=overlap,
                result='PASS' if passed else 'FAIL',navigation_log_sha256=sha(folder/'navigation.log')))
        if len(records)!=9 or len({(r['disk_thickness_mm'],r['radial_clearance_mm']) for r in records})!=9:
            raise ValueError('expected nine distinct scan cases')
        with (self.out/'cup_scan.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
        self.scan_records=records
        fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        ts=sorted({r['disk_thickness_mm'] for r in records}); cs=sorted({r['radial_clearance_mm'] for r in records})
        for ax,key,title in zip(axs,('steel_bottom_mm','steel_wall_mm','pocket_depth_mm'),
                                ('Remaining steel bottom [mm]','Remaining radial wall [mm]','Recess depth [mm]')):
            values=np.array([[next(r[key] for r in records if r['disk_thickness_mm']==t and
                          r['radial_clearance_mm']==c) for c in cs] for t in ts])
            ax.imshow(values,cmap='YlGnBu',vmin=0,vmax=max(1,values.max()))
            for i,t in enumerate(ts):
                for j,c in enumerate(cs): ax.text(j,i,f'{values[i,j]:g}\nPASS 18/18',ha='center',va='center',
                    color='white' if values[i,j]>.6*values.max() else '#17212b',fontsize=10)
            ax.set(xticks=range(3),xticklabels=cs,yticks=range(3),yticklabels=ts,
                   xlabel='Provisional radial gap [mm]',ylabel='Provisional disk thickness [mm]',title=title)
        fig.suptitle('Nine provisional cup geometries • overlap / navigation evidence, no transport')
        self.save(fig,'cup_scan','Every cell: no overlaps detected in recorded 20,000-sample/placement checks; '
                  '18/18 material probes passed. Axial gap is 0.1 mm per face in every case. '
                  'Clearances are analytic dimensions, not measured hardware or minimum clearance to every ANNIE structure.')
        negative_logs=list((self.p0/'scans').glob('*/overlap_rejection.log'))
        if len(negative_logs)!=1: raise ValueError('missing intentional-overlap evidence')
        negative=negative_logs[0].read_text()
        if 'Overlap is detected' not in negative or 'GEOMETRY_ONLY FAIL' not in negative:
            raise ValueError('intentional overlap was not rejected')
        self.inputs.append(negative_logs[0])
        self.results['scan']=dict(status='PASS' if all(r['result']=='PASS' for r in records) else 'FAIL',
                                 cases=len(records),recorded_overlap_rejection=summary['overlap_rejection'],
                                 rejection_log_sha256=sha(negative_logs[0]))

    def optics(self):
        path=self.db/'OPTICS_SANDI_wbls_Nov22.ratdb'
        if sha(path)!=self.m['static_optics']['sha256']:
            raise ValueError('optics does not match geometry manifest hash')
        tables=[b for b in blocks(path) if b['index'].startswith('wbls')]
        if len(tables)!=2: raise ValueError('expected two uploaded WbLS tables')
        json_save(self.out/'static_wbls_tables.json',tables)
        fig,axs=plt.subplots(2,3,figsize=(16,10),layout='constrained')
        labels=['WbLS (no Gd)','WbLS + Gd (selected material)']
        for i,(b,label) in enumerate(zip(tables,labels)):
            style=dict(color=['#965ba8','#087e8b'][i],ls=['-','--'][i],label=label,lw=1.5)
            for ax,key,title,ylabel in [(axs[0,0],'ABSLENGTH','Absorption length','Native RATDB length (mm convention)'),
                (axs[0,1],'RINDEX','Refractive index','RINDEX (dimensionless)'),
                (axs[0,2],'SCINTILLATION','Scintillation emission','Raw dy_dwavelength table value')]:
                x=np.array(b[key+'_value1']); y=np.array(b[key+'_value2'])
                if len(x)!=len(y) or not np.isfinite([*x,*y]).all(): raise ValueError('invalid optical array')
                # Emission input is not necessarily sorted; plot samples in wavelength order, no clipping.
                order=np.argsort(x,kind='stable'); ax.plot(x[order],y[order],**style)
                ax.set(title=title,xlabel='Wavelength [nm]',ylabel=ylabel); ax.grid(alpha=.2)
            axs[0,0].set_yscale('log')
            axs[1,0].bar(i,b['LIGHT_YIELD'],color=style['color'],width=.6)
            axs[1,0].text(i,b['LIGHT_YIELD']+2,str(b['LIGHT_YIELD']),ha='center')
            axs[1,1].bar(np.arange(3)+(i-.5)*.32,b['SCINTWAVEFORM_value1'],width=.32,
                         color=style['color'],label=label)
            axs[1,2].plot(range(1,4),b['SCINTWAVEFORM_value2'],'o'+style['ls'],color=style['color'],label=label)
        for ax in axs[0]: ax.legend(fontsize=8)
        axs[0,2].axhline(0,color='#777',lw=.7)
        negative=sum(v<0 for v in tables[0]['SCINTILLATION_value2'])
        axs[0,2].text(.03,.9,f'Undoped input: {negative} negative samples\nPreserved, not clipped',
                       transform=axs[0,2].transAxes,fontsize=8)
        axs[1,0].set(title='Native LIGHT_YIELD parameter',ylabel='Raw RAT parameter; not PE/MeV',
                      xticks=[0,1],xticklabels=['no Gd','Gd selected'],ylim=(0,130))
        axs[1,1].set(title='Signed SCINTWAVEFORM_value1',ylabel='Raw time parameter (ns convention)',
                      xticks=range(3),xticklabels=['component 1','component 2','component 3'])
        for i,v in enumerate(tables[1]['SCINTWAVEFORM_value1']):
            axs[1,1].text(i,v-5,f'{v:g}',ha='center',va='top')
        axs[1,1].set_ylim(-140,5)
        axs[1,2].set(title='SCINTWAVEFORM_value2 weights',xlabel='Component',ylabel='Raw weight',
                     yscale='log',xticks=[1,2,3])
        for i,v in enumerate(tables[1]['SCINTWAVEFORM_value2']):
            axs[1,2].annotate(str(v),(i+1,v),xytext=(5,6),textcoords='offset points')
        axs[1,2].text(.03,.93,f"SCINT_RISE_TIME = {tables[1]['SCINT_RISE_TIME']} (both tables)\n"
                      'Signed values are not reinterpreted as a waveform',transform=axs[1,2].transAxes,fontsize=9,va='top')
        fig.suptitle('Authoritative uploaded WbLS optical tables • STATIC RATDB ONLY\n'
                     'Not resolved Geant4 material properties; no optical photon simulation',fontsize=15)
        self.save(fig,'wbls_static_optics','Values parsed from the untouched uploaded tuning file. '
                  'Array values and signs are preserved; emission points are ordered by wavelength for display only. '
                  'Units follow RAT length/time/wavelength conventions, pending matching-core verification. '
                  'LIGHT_YIELD 110 is a native parameter, not calibrated PE/MeV. '
                  'SCINT_RISE_TIME=0.25 and alpha time arrays match the general time arrays; no pulse shape was synthesized.')
        self.results['optics']=dict(status='STATIC_ONLY',sha256=sha(path),tables=[b['index'] for b in tables],
            light_yield=[b['LIGHT_YIELD'] for b in tables],undoped_negative_emission_samples=negative,
            runtime_properties_verified=False)

    def pmts(self):
        active={b['pos_table'] for b in self.detector if b.get('type')=='pmtarray' and b.get('enable',1)}
        tables=[b for b in blocks(self.db/'PMTINFO_run2_config14_scan.ratdb') if b['name'] in active]
        byid={b['index']:b for b in self.detector}; offset=np.zeros(3); current='detector'
        while current:
            b=byid[current]; offset+=b.get('position',[0,0,0]); current=b.get('mother','')
        channels=[]
        for b in tables:
            for i,channel in enumerate(b['channel_number']):
                pos=np.array([b[k][i] for k in ('x','y','z')]); direction=[b[k][i] for k in ('dir_x','dir_y','dir_z')]
                pred=pos-offset
                channels.append(dict(table=b['name'],table_row=i,channel=channel,
                    **dict(zip(('table_x_mm','table_y_mm','table_z_mm'),pos)),
                    **dict(zip(('dir_x','dir_y','dir_z'),direction)),
                    **dict(zip(('predicted_x_mm','predicted_y_mm','predicted_z_mm'),pred))))
        if len(channels)!=132 or len({r['channel'] for r in channels})!=132:
            raise ValueError('expected 132 unique stored channel numbers')
        with (self.out/'pmt_channels.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(channels[0])); w.writeheader(); w.writerows(channels)
        self.channels=channels
        fig,axs=plt.subplots(2,3,figsize=(16,11),layout='constrained')
        for row,prefix in enumerate(('table','predicted')):
            for col,(u,v) in enumerate((('x','y'),('x','z'),('y','z'))):
                ax=axs[row,col]
                for b in tables:
                    rows=[r for r in channels if r['table']==b['name']]
                    x=[r[f'{prefix}_{u}_mm'] for r in rows]; y=[r[f'{prefix}_{v}_mm'] for r in rows]
                    dx=[150*r['dir_'+u] for r in rows]; dy=[150*r['dir_'+v] for r in rows]
                    ax.quiver(x,y,dx,dy,angles='xy',scale_units='xy',scale=1,width=.004)
                    ax.scatter(x,y,s=10)
                ax.set(xlabel=f'{prefix} {u} [mm]',ylabel=f'{prefix} {v} [mm]',aspect='equal',
                       title=('Raw PMTINFO table' if row==0 else 'Presentation-described prediction (UNVERIFIED)'))
                ax.grid(alpha=.2)
        fig.suptitle(f'132 PMT coordinate/direction diagnostics • arrow length 150 mm for display\n'
                     f'Prediction: table − ancestor-translation sum {tuple(np.round(offset,3))} mm; '
                     'no measured runtime positions',fontsize=13)
        self.save(fig,'pmt_frames','Upper row: untouched PMTINFO coordinates. Lower row: diagnostic prediction using '
                  'the translation-only GetLocalOffset convention described on presentation pages 9–12. '
                  'Neither lower row nor directions have been verified against the matching factory. '
                  'Scan/glass-front corrections are already in the table and are not applied twice.')
        fig,axs=plt.subplots(5,2,figsize=(16,23),layout='constrained')
        for ax,b in zip(axs.flat,tables):
            rows=[r for r in channels if r['table']==b['name']]
            x=np.array([r['table_x_mm'] for r in rows]); y=np.array([r['table_y_mm'] for r in rows])
            ax.scatter(x,y,s=18,color='#086e80')
            ax.quiver(x,y,[120*r['dir_x'] for r in rows],[120*r['dir_y'] for r in rows],
                      angles='xy',scale_units='xy',scale=1,color='#086e80',width=.005)
            for r in rows: ax.annotate(str(r['channel']),(r['table_x_mm'],r['table_y_mm']),
                                      xytext=(4,5),textcoords='offset points',fontsize=9)
            if all(abs(r['dir_z'])>.99 for r in rows):
                ax.text(.02,.02,f"dir_z = {rows[0]['dir_z']:+g} (out/in of XY plane)",transform=ax.transAxes)
            ax.set(title=b['name'].replace('PMTINFO_config14_scan_','')+f' • {len(rows)} channels',
                   xlabel='Raw table x [mm]',ylabel='Raw table y [mm]',aspect='equal')
            ax.margins(.2); ax.grid(alpha=.15)
        fig.suptitle('All 132 stored channel identifiers • RAW TABLE FRAME • directions shown\n'
                     'Labels are channel_number, not an inferred internal RAT PMT index',fontsize=15)
        self.save(fig,'pmt_identifiers','Every active stored channel is labeled in its own table group; '
                  'XY arrows represent dir_x/dir_y (120 mm display length). Top/bottom dir_z is stated explicitly. '
                  'Complete original and predicted coordinates/directions are in pmt_channels.csv.')
        self.results['pmts']=dict(status='STATIC_ONLY',channels=132,local_offset_prediction_mm=offset.tolist(),
            maximum_direction_norm_error=max(abs(np.linalg.norm([r['dir_x'],r['dir_y'],r['dir_z']])-1) for r in channels),
            matching_factory_verified=False)

    def root(self):
        files=list(self.p0.glob('*schema_smoke.root'))
        if len(files)!=1: raise ValueError('expected one synthetic ROOT schema smoke file')
        rows=read_root(files[0]); run=rows['runs'][0]
        if run['is_synthetic']!=1 or run['generator_mode']!='SYNTHETIC_IO_TEST':
            raise ValueError('diagnostics require explicitly synthetic fixtures')
        events=rows['events']; ids=[r['decay_id'] for r in events]
        pe=[r['n_detected_pe'] for r in events]
        joins=[sum(p['decay_id']==i for p in rows['pmt_pe']) for i in ids]
        if pe!=joins or len(events)!=run['n_generated']: raise ValueError('ROOT join/count mismatch')
        fig,axs=plt.subplots(1,3,figsize=(15,5),layout='constrained')
        axs[0].bar(ids,[1]*len(ids),color=['#a7c7e8' if p==0 else '#0c8190' for p in pe])
        for i in ids: axs[0].text(i,1.03,'present',ha='center')
        axs[0].set(title='One events row for every fixture ID',xlabel='Synthetic decay_id',ylabel='Row count',
                   xticks=ids,ylim=(0,1.4),yticks=[0,1])
        axs[1].bar(np.array(ids)-.18,pe,width=.36,label='events.n_detected_pe')
        axs[1].bar(np.array(ids)+.18,joins,width=.36,label='joined pmt_pe rows')
        axs[1].set(title='Zero-PE fixtures stay in events',xlabel='Synthetic decay_id',ylabel='Fabricated PE count',xticks=ids)
        axs[1].legend(fontsize=8)
        axs[2].bar(ids,[r['n_rat_evs'] for r in events],color='#8a73a7')
        axs[2].set(title='Hypothetical EV multiplicity',xlabel='Synthetic decay_id',ylabel='Fabricated n_rat_evs',xticks=ids)
        fig.suptitle('SYNTHETIC ROOT I/O FIXTURES ONLY • NOT PHYSICS OR RAT RESPONSE',fontsize=15)
        self.save(fig,'synthetic_root','Reads the existing native ROOT file with PyROOT and validates its schema/joins. '
                  'All three fabricated events remain, including both zero-PE entries and the hypothetical split-EV entry. '
                  'This does not establish actual RAT all-decay retention, triggering, PE efficiency or energy response.')
        json_save(self.out/'synthetic_root_summary.json',dict(input=str(files[0]),sha256=sha(files[0]),
            entries={k:len(v) for k,v in rows.items()},events=events))
        self.results['root']=dict(status='PASS_SYNTHETIC_ONLY',events=len(events),
            zero_pe_events=sum(p==0 for p in pe),all_event_ids=ids,pmt_pe_join_counts=joins)

    def presentation(self):
        path=self.a.presentation.resolve()
        if not path.is_file(): raise FileNotFoundError(path)
        self.inputs.append(path)
        self.run(['pdftotext','-layout',path,self.out/'presentation_text.txt'],'presentation.log')
        text=(self.out/'presentation_text.txt').read_text(); pages=text.split('\f')
        if len(pages)<23: raise ValueError('presentation page count changed')
        # Embed source evidence only in the PRIVATE report; no public repo copy.
        for page in (9,11,21,22,23):
            self.run(['pdftoppm','-f',page,'-l',page,'-scale-to','1500','-png','-singlefile',
                      path,self.out/f'presentation_p{page}'],f'presentation_p{page}.log')
        self.results['presentation']=dict(path=str(path),sha256=sha(path),pages=23,
            comparison='p9/11 PMT offset; p21/22 original holder; p23 SuperSANDI validation requirement')

    def report(self):
        self.results['input_preservation']={'status':'PASS'}
        for path,h in self.original_hashes.items():
            if sha(path)!=h: raise ValueError('input was modified '+path)
        inputs={str(p.resolve()):sha(p) for p in self.inputs}
        json_save(self.out/'diagnostic_summary.json',self.results)
        code=list((REPO/'tools/co60').glob('*.py'))+[REPO/'tests/co60/geometry_sections.cc',
                   REPO/'tests/co60/CMakeLists.txt',
                   REPO/'integration/annie_ratpac/ManifestGeometry.hh']
        json_save(self.out/'provenance.json',dict(inputs_sha256=inputs,
            code_sha256={str(p.relative_to(REPO)):sha(p) for p in code},
            public_head=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),
            command=[sys.executable,*sys.argv],numpy=np.__version__,matplotlib=matplotlib.__version__,
            software={tool:subprocess.check_output([tool,'--version'],text=True).splitlines()[0]
                      for tool in ('root-config','geant4-config','cmake','c++')}))
        intro='''<h1>Local Co-60 P0 visualization and diagnostics</h1>
<p><b>Geometry/static-table inspection and synthetic I/O only. No Co-60 transport,
optical photon response, simulated efficiency, or final runtime alignment.</b></p>
<p>All figures are private artifacts. The original application/tuning files are unchanged.
Geometry sections follow actual Geant4 navigator boundaries through manifest-built
solids; materials are identity-only stand-ins without optical tables.</p>
<p>Full-detector views show the uploaded tank/liner/liquid envelopes, SuperSANDI
octagonal CSG, and holder in detector frame. They omit hall, tank PMT shapes,
internal supports and black sheets. They are not a complete RAT geometry export.</p>
<p>Geant4 section traces advance in one axis to actual boundaries (epsilon 1e-7 mm);
the perpendicular axis is sampled in 101–1201 columns. Therefore narrow features
can be missed between columns. The 0.01 mm active layer and 0.025 mm cap are resolved
along the trace axis. Full-scale thin boundaries are subpixel, not absent.</p>
<p>Black dimension labels refer to the authoritative uploaded ANNIE model, whose
own estimates are not hardware metrology. Teal marks the locked disk diameter.
Orange/dashed labels mark new provisional dimensions or placement.</p>
'''
        parts=[intro,'<h2>Results</h2><pre>'+html.escape(json.dumps(self.results,indent=2))+'</pre>']
        for f in self.figures:
            name=f['name']; parts.append(f'<h2>{html.escape(name.replace("_"," "))}</h2>'
                f'<a href="{name}.pdf">PDF</a> · <a href="{name}.png">PNG</a>'
                f'<figure><img src="{name}.png"><figcaption>{html.escape(f["caption"])}</figcaption></figure>')
        parts.append('<h2>Nine provisional cup configurations</h2><p>All dimensions below are mm; '
            'gaps are surface clearances within the modeled cup, not fit tolerances certified for hardware.</p>')
        cols=['disk_thickness_mm','radial_clearance_mm','axial_clearance_each_mm','pocket_diameter_mm',
              'pocket_depth_mm','steel_bottom_mm','steel_wall_mm','probes_passed','result']
        parts.append('<table><tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in cols)+'</tr>')
        parts.extend('<tr>'+''.join('<td>'+html.escape(f'{r[c]:g}' if isinstance(r[c],float) else str(r[c]))+
                     '</td>' for c in cols)+'</tr>' for r in self.scan_records)
        parts.append('</table><p><a href="cup_scan.csv">Complete scan CSV</a> · '
                     '<a href="pmt_channels.csv">All 132 PMT channels and directions</a> · '
                     '<a href="sections.csv">Geant4 section intervals</a> · '
                     '<a href="static_wbls_tables.json">Static optical tables</a> · '
                     '<a href="provenance.json">Input/code hashes and command</a></p>')
        parts.append('''<h2>Comparison to original ANNIE presentation</h2>
<p>The supplied 23-page PDF was inspected locally. Pages 21–22 describe the original
380 mm × 70 mm PVC holder, 50 × 50 mm BGO, 10 × 50 mm separator and steel placer.
Those dimensions agree with the retained modeled stack. The original 7.8 mm AmBe
capsule is removed; the disk/recess and bore are explicit P0 changes. The slide is
not a bore measurement. The uploaded v2 file is solid; 1 mm bore radius remains
user-approved provisional input. R1398 dimensions and some foam/plastic details
are estimated in the uploaded model, not newly measured here.</p>
<p>Pages 9 and 11 explain translation-only GetLocalOffset=(0,1198.6,-3339) mm and
the scan/glass-front corrections. The PMT diagnostic displays the predicted
subtraction once. The matching core is still required to verify actual placements,
rotations, channel mapping and shape alignment. Page 23 specifically calls for
SuperSANDI geometry checks; this suite supplies local shape evidence, not its final
runtime certificate. Its dimensions are taken from the factory, not inferred from
a slide rendering.</p>''')
        for page in (9,11,21,22,23):
            parts.append(f'<details><summary>Original presentation page {page}</summary>'
                         f'<img src="presentation_p{page}.png"></details>')
        parts.append('''<h2>Still requires the matching RAT stack</h2>
<p>Actual baseline meta/output values; no-hit and split-EV retention; resolved
WbLS material/MPT and time semantics; full PMT/structure/holder alignment and
overlaps; optical wrapping/coupling behavior. No transport or optical response
has been implemented. See diagnostic_summary.json for precisely scoped checks.</p>''')
        style='body{font:16px system-ui;max-width:1250px;margin:32px auto;padding:0 24px;color:#182630} '
        style+='img{width:100%;height:auto} figure{margin:24px 0} figcaption{font-size:14px} '
        style+='table{border-collapse:collapse;font-size:13px} td,th{border:1px solid #bbb;padding:6px} '
        style+='pre{white-space:pre-wrap;background:#eef3f5;padding:16px} details{margin:12px 0}'
        (self.out/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
            '<title>Co-60 local P0 diagnostics</title><style>'+style+'</style><body>'+''.join(parts)+'</body></html>')
        (self.out/'README.md').write_text('# Private local P0 diagnostics\n\nOpen [index.html](index.html) for all figures, '
            'comparisons and limitations. Every diagnostic figure has PNG and PDF.\n\n'
            'No physics/optical response was simulated. Full runtime alignment and MPT resolution remain blocked.\n\n'
            'Reproduce with the command in `provenance.json`, choosing a NEW private output directory. '
            'Build/section logs and hash manifests are retained. See the public-project '
            '`docs/co60/DIAGNOSTICS.md` for the CLI.\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--app',type=Path,required=True,help='preserved authoritative uploaded application')
    p.add_argument('--p0',type=Path,required=True,help='previous passed wrapper output')
    p.add_argument('--presentation',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    s=Suite(p.parse_args())
    for label,method in [('Geant4 sections',s.geometry),('scan summary',s.scans),('static optics',s.optics),
                         ('PMT coordinates',s.pmts),('synthetic ROOT',s.root),('presentation comparison',s.presentation),
                         ('HTML report',s.report)]:
        print(label,flush=True); method()
    print('PASS: private diagnostic suite:',s.out/'index.html',flush=True)


if __name__=='__main__': main()
