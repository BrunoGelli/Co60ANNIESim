// P0 navigation only. No particle transport, no optical or full ANNIE validation.
// Assembly geometry uses the SAME loader as the private RAT adapter.
#include "ManifestGeometry.hh"
#include <G4Box.hh>
#include <G4Polyhedra.hh>
#include <G4SubtractionSolid.hh>
#include <G4NistManager.hh>
#include <G4Navigator.hh>
#include <G4GeometryManager.hh>
#include <Randomize.hh>
#include <iostream>
#include <iomanip>

int main(int argc,char** argv) {
  try {
    if (argc!=4) throw std::runtime_error("usage: geometry_probe assembly.tsv vessel.txt probes.tsv");
    CLHEP::HepRandom::setTheSeed(6012);
    // Geometry-only stand-ins, explicitly named like private materials for navigation.
    // NO material property table is installed. No optics/transport claim is possible.
    auto* nist=G4NistManager::Instance();
    const std::map<std::string,std::string> materials{
      {"water_gdS_0p2","G4_WATER"},{"wbls1pct_ly95_gdS0p2pct_SANDI","G4_WATER"},
      {"superSANDI_nylon","G4_NYLON-6-6"},{"pvc_white","G4_POLYVINYL_CHLORIDE"},
      {"air","G4_AIR"},{"stainless_steel","G4_STAINLESS-STEEL"},
      {"BGO_scint","G4_BGO"},{"polystyrene_foam","G4_POLYSTYRENE"},
      {"glass","G4_Pyrex_Glass"},{"pmt_vacuum","G4_Galactic"},
      {"delrin_black","G4_POLYOXYMETHYLENE"},{"polystyrene","G4_POLYSTYRENE"}};
    for (const auto& [alias,nistName]: materials) {
      auto* mat=nist->FindOrBuildMaterial(nistName);
      if (!mat) throw std::runtime_error("missing NIST material "+nistName);
      new G4Material(alias,mat->GetDensity(),mat);
    }
    auto get=[](const std::string& name){return G4Material::GetMaterial(name);};
    std::map<std::string,G4OpticalSurface*> surfaces;
    for (const auto& name:{"ptfe","pvc_white"}) surfaces[name]=new G4OpticalSurface(name);
    auto* worldLV=new G4LogicalVolume(new G4Box("world",2500,2500,2500),get("air"),"world");
    auto* world=new G4PVPlacement(nullptr,{},worldLV,"world",nullptr,false,0);
    auto* tankLV=new G4LogicalVolume(new G4Tubs("tank",0,1516.6,1973.8,0,360*deg),
                                    get("water_gdS_0p2"),"detector");
    auto* tank=new G4PVPlacement(nullptr,{},tankLV,"detector",worldLV,false,0);
    double h,r,t,g,x,y,z; int faces;
    std::ifstream vessel(argv[2]);
    if (!(vessel>>h>>r>>t>>g>>faces>>x>>y>>z)) throw std::runtime_error("bad vessel parameters");
    // Same polyhedra/subtraction construction as uploaded factory; parameters from manifest.
    double zs[]={-h/2,-h/2+t,-h/2+t+g,h/2-t-g,h/2-t,h/2};
    double ri[]={0,0,r-t,r-t,0,0}, ro[]={r,r,r,r,r,r};
    auto* nylon=new G4Polyhedra("nylonVessel",0,360*deg,faces,6,zs,ri,ro);
    double zi[]={-h/2,h/2}, ii[]={0,0}, oo[]={r,r};
    auto* outer=new G4Polyhedra("nv_inner",0,360*deg,faces,2,zi,ii,oo);
    auto* liquid=new G4SubtractionSolid("detectionMedium",outer,nylon);
    auto* nylonLV=new G4LogicalVolume(nylon,get("superSANDI_nylon"),"nv_logic");
    auto* liquidLV=new G4LogicalVolume(liquid,get("wbls1pct_ly95_gdS0p2pct_SANDI"),"detmed_logic");
    auto* nylonPV=new G4PVPlacement(nullptr,{x,y,z},nylonLV,"superSANDIvessel",tankLV,false,0);
    auto* liquidPV=new G4PVPlacement(nullptr,{x,y,z},liquidLV,"superSANDIDetMed",tankLV,false,0);
    auto placements=co60::BuildManifest(argv[1],liquidLV,"superSANDIDetMed",surfaces);
    placements.insert(placements.begin(),{tank,nylonPV,liquidPV});
    bool failed=false;
    for (auto* pv:placements) failed=pv->CheckOverlaps(20000,0,true,10)||failed;
    G4GeometryManager::GetInstance()->CloseGeometry();
    G4Navigator navigator; navigator.SetWorldVolume(world);
    std::ifstream probes(argv[3]);
    std::string label,expectedPV,expectedMaterial;
    std::cout<<"PROBE label x_mm y_mm z_mm volume material pass\n"<<std::setprecision(12);
    while (probes>>label>>x>>y>>z>>expectedPV>>expectedMaterial) {
      auto* pv=navigator.LocateGlobalPointAndSetup({x,y,z},nullptr,false);
      std::string name=pv?std::string(pv->GetName()):"OUTSIDE";
      std::string mat=pv?std::string(pv->GetLogicalVolume()->GetMaterial()->GetName()):"NONE";
      bool ok=name==expectedPV && mat==expectedMaterial;
      std::cout<<"PROBE "<<label<<' '<<x<<' '<<y<<' '<<z<<' '<<name<<' '<<mat<<' '<<ok<<'\n';
      failed=failed||!ok;
    }
    std::cout<<"GEOMETRY_ONLY "<<(failed?"FAIL":"PASS")<<"; stand-in materials; no ANNIE PMTs or optics\n";
    G4GeometryManager::GetInstance()->OpenGeometry();
    return failed?1:0;
  } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 2; }
}
