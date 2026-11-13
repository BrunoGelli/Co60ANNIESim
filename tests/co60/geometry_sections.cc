// Static geometry diagnostic only: trace navigator boundaries, no particle transport.
// Read PRIVATE shape parameters and assembly with the same adapter as P0 tests.
#include "ManifestGeometry.hh"
#include <G4Box.hh>
#include <G4Polyhedra.hh>
#include <G4SubtractionSolid.hh>
#include <G4NistManager.hh>
#include <G4Navigator.hh>
#include <G4GeometryManager.hh>
#include <iostream>
#include <iomanip>
#include <set>

int main(int argc,char** argv) {
  try {
    if (argc!=6) throw std::runtime_error(
      "usage: geometry_sections assembly.tsv vessel.txt envelopes.tsv requests.tsv sections.csv");
    // Material identities only. No MPT, optical tracking or density-dependent result.
    std::set<std::string> names;
    std::set<std::string> surfaceNames;
    for (int arg:{1,3}) {
      std::ifstream file(argv[arg]); std::string line;
      while (std::getline(file,line)) {
        std::istringstream row(line); std::string name,mother,material,type,surface;
        if (!(row>>name>>mother>>material>>type>>surface)) throw std::runtime_error("bad input manifest");
        names.insert(material); if (surface!="-") surfaceNames.insert(surface);
      }
    }
    names.insert("air"); names.insert("superSANDI_nylon");
    names.insert("wbls1pct_ly95_gdS0p2pct_SANDI");
    auto* base=G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    for (const auto& name:names) new G4Material(name,base->GetDensity(),base);
    std::map<std::string,G4OpticalSurface*> surfaces;
    for (const auto& name:surfaceNames) surfaces[name]=new G4OpticalSurface(name);
    auto get=[](const std::string& n){return G4Material::GetMaterial(n);};
    auto* worldLV=new G4LogicalVolume(new G4Box("frame",20000,20000,20000),get("air"),"frame");
    auto* world=new G4PVPlacement(nullptr,{},worldLV,"detector_frame_world",nullptr,false,0);
    auto envelopes=co60::BuildManifest(argv[3],worldLV,"detector_frame_world",surfaces);
    G4LogicalVolume* detector=nullptr;
    for (auto* pv:envelopes) if (pv->GetName()=="detector") detector=pv->GetLogicalVolume();
    if (!detector) throw std::runtime_error("missing detector envelope");
    double h,r,t,g,x,y,z; int faces;
    std::ifstream vessel(argv[2]);
    if (!(vessel>>h>>r>>t>>g>>faces>>x>>y>>z)) throw std::runtime_error("bad vessel parameters");
    // Exact construction path in supplied GeoANNIEFactory::ConstructSuperSANDI.
    double zs[]={-h/2,-h/2+t,-h/2+t+g,h/2-t-g,h/2-t,h/2};
    double ri[]={0,0,r-t,r-t,0,0},ro[]={r,r,r,r,r,r};
    auto* nylon=new G4Polyhedra("nylonVessel",0,360*deg,faces,6,zs,ri,ro);
    double zi[]={-h/2,h/2},ii[]={0,0},oo[]={r,r};
    auto* outer=new G4Polyhedra("nv_inner",0,360*deg,faces,2,zi,ii,oo);
    auto* liquid=new G4SubtractionSolid("detectionMedium",outer,nylon);
    G4ThreeVector low,high; outer->BoundingLimits(low,high);
    std::cout<<std::setprecision(15)<<"VESSEL_OUTER_BOUNDS_LOCAL_MM "
             <<low.x()<<' '<<low.y()<<' '<<low.z()<<' '
             <<high.x()<<' '<<high.y()<<' '<<high.z()<<'\n';
    auto* nylonLV=new G4LogicalVolume(nylon,get("superSANDI_nylon"),"nv_logic");
    auto* liquidLV=new G4LogicalVolume(liquid,get("wbls1pct_ly95_gdS0p2pct_SANDI"),"detmed_logic");
    new G4PVPlacement(nullptr,{x,y,z},nylonLV,"superSANDIvessel",detector,false,0);
    new G4PVPlacement(nullptr,{x,y,z},liquidLV,"superSANDIDetMed",detector,false,0);
    co60::BuildManifest(argv[1],liquidLV,"superSANDIDetMed",surfaces);
    G4GeometryManager::GetInstance()->CloseGeometry();
    G4Navigator nav; nav.SetWorldVolume(world);
    std::ifstream requests(argv[4]); std::ofstream out(argv[5]);
    if (!requests || !out) throw std::runtime_error("cannot open requests/output");
    out<<"section,u0_mm,u1_mm,v0_mm,v1_mm,volume,material\n"<<std::setprecision(15);
    std::string name; int u,v,w,n; double fixed,umin,umax,vmin,vmax;
    size_t intervals=0;
    while (requests>>name>>u>>v>>w>>fixed>>umin>>umax>>vmin>>vmax>>n) {
      if (std::set<int>{u,v,w}!=std::set<int>{0,1,2} || n<1 || umin>=umax || vmin>=vmax)
        throw std::runtime_error("invalid section request");
      const double du=(umax-umin)/n,eps=1e-7;
      for (int col=0;col<n;++col) {
        double boundary=vmin;
        G4ThreeVector direction; direction[v]=1;
        for (int step=0;boundary<vmax-eps;++step) {
          if (step>1000) throw std::runtime_error("navigator failed to advance");
          G4ThreeVector point; point[u]=umin+(col+.5)*du; point[v]=boundary+eps; point[w]=fixed;
          auto* pv=nav.LocateGlobalPointAndSetup(point,nullptr,false);
          if (!pv) throw std::runtime_error("section outside world");
          double safety=0;
          double length=nav.ComputeStep(point,direction,vmax-point[v],safety);
          if (length<eps/10) throw std::runtime_error("zero navigation step at "+name);
          const double next=std::min(vmax,point[v]+length);
          out<<name<<','<<umin+col*du<<','<<umin+(col+1)*du<<','<<boundary<<','<<next<<','
             <<pv->GetName()<<','<<pv->GetLogicalVolume()->GetMaterial()->GetName()<<'\n';
          boundary=next; ++intervals;
        }
      }
      std::cout<<"SECTION "<<name<<" columns="<<n<<" boundary_epsilon_mm="<<eps<<'\n';
    }
    G4GeometryManager::GetInstance()->OpenGeometry();
    std::cout<<"PASS: "<<intervals<<" navigator intervals; geometry only; optical-free material stand-ins\n";
    return 0;
  } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 2; }
}
