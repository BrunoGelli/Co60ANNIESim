// Standalone P1 truth transport. No RAT code, optical physics or tuned optics.
// Geometry/material inputs are private files; native ROOT has no experiment dictionaries.
#include "ManifestGeometry.hh"
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VPhysicsConstructor.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4DecayPhysics.hh>
#include <G4EmParameters.hh>
#include <G4VRadioactiveDecay.hh>
#include <G4PhysicsListHelper.hh>
#include <G4NuclearLevelData.hh>
#include <G4DeexPrecoParameters.hh>
#include <G4NuclideTable.hh>
#include <G4GenericIon.hh>
#include <G4Ions.hh>
#include <G4IonTable.hh>
#include <G4Gamma.hh>
#include <G4ParticleGun.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4Event.hh>
#include <G4Track.hh>
#include <G4Step.hh>
#include <G4VProcess.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Polyhedra.hh>
#include <G4SubtractionSolid.hh>
#include <G4Navigator.hh>
#include <G4TransportationManager.hh>
#include <G4Region.hh>
#include <G4ProductionCuts.hh>
#include <G4Version.hh>
#include <Randomize.hh>
#include <G4RandomDirection.hh>
#include <TFile.h>
#include <TTree.h>
#include <TUUID.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <map>
#include <array>
#include <memory>
#include <cmath>

namespace fs=std::filesystem;
using namespace CLHEP;
const std::vector<std::string> regions={"bgo","wbls","gd_water","steel","source_plastic",
 "pvc","foam","glass","nylon","pmt_base","liner","vacuum_air","other"};
int region(const std::string& volume,const std::string& mat) {
  if (volume=="AmBe_BGO") return 0;
  if (volume=="superSANDIDetMed") return 1;
  if (volume=="detector") return 2;
  if (mat=="stainless_steel") return 3;
  if (volume=="Co60_disk" || volume=="Co60_active") return 4;
  if (mat=="pvc_white") return 5;
  if (mat=="polystyrene_foam") return 6;
  if (mat=="glass") return 7;
  if (mat=="superSANDI_nylon") return 8;
  if (mat=="delrin_black") return 9;
  if (mat=="acrylic_black") return 10;
  if (mat=="air" || mat=="pmt_vacuum") return 11;
  return 12;
}
std::string contents(const fs::path& p) {
  std::ifstream f(p); if (!f) throw std::runtime_error("cannot open "+p.string());
  return std::string(std::istreambuf_iterator<char>(f),{});
}
bool sourceIon(const G4Track& t) {
  auto* ion=dynamic_cast<const G4Ions*>(t.GetDefinition());
  return t.GetParentID()==0 && ion && ion->GetAtomicNumber()==27 && ion->GetAtomicMass()==60
         && ion->GetExcitationEnergy()==0;
}

// G4 11.4 DecayAnalog ignores preassigned time and resamples the at-rest lifetime.
// This PRIMARY-only adapter bypasses that wait before any daughters are made.
// No particle lifetime, daughter clock, evaluated channel, or angular kernel is changed.
class PromptCoDecay final:public G4VRadioactiveDecay {
 public:
  explicit PromptCoDecay(bool prompt):G4VRadioactiveDecay("RadioactiveDecay"),prompt_(prompt) {
    std::cout<<"RDM default threshold_s="<<GetThresholdForVeryLongDecayTime()/s<<'\n';
    SetThresholdForVeryLongDecayTime(1.e15*s); SetARM(true);
    std::cout<<"RDM configured threshold_s="<<GetThresholdForVeryLongDecayTime()/s<<'\n';
  }
  G4VParticleChange* DecayIt(const G4Track& track,const G4Step& step) override {
    if (prompt_ && sourceIon(track)) {
      // ParticleChange retains its current-track pointer through UpdateStep.
      // Keep the copy alive until the next primary decay, beyond this call.
      promptTrack_=std::make_unique<G4Track>(track); auto& prompt=*promptTrack_;
      // G4Track's copy constructor deliberately clears this non-owning pointer;
      // ParticleChange::Initialize requires it to initialize step bookkeeping.
      prompt.SetStep(track.GetStep());
      prompt.SetTrackID(track.GetTrackID()); prompt.SetParentID(track.GetParentID());
      prompt.SetGlobalTime(0); prompt.SetLocalTime(0); prompt.SetProperTime(0);
      // DecayAnalog samples the long parent wait only when status==fStopButAlive.
      prompt.SetTrackStatus(fAlive);
      return G4VRadioactiveDecay::DecayIt(prompt,step);
    }
    return G4VRadioactiveDecay::DecayIt(track,step);
  }
 protected:
  G4double GetMeanLifeTime(const G4Track& track,G4ForceCondition* force) override {
    return prompt_ && sourceIon(track)?0.:G4VRadioactiveDecay::GetMeanLifeTime(track,force);
  }
 private: bool prompt_; std::unique_ptr<G4Track> promptTrack_;
};
class SourceDecayPhysics final:public G4VPhysicsConstructor {
 public:
  explicit SourceDecayPhysics(bool prompt):G4VPhysicsConstructor("Co60RDM"),prompt_(prompt) {}
  void ConstructParticle() override { G4GenericIon::GenericIon(); }
  void ConstructProcess() override {
    G4EmParameters::Instance()->SetAuger(true);
    G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true);
    G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(new PromptCoDecay(prompt_),G4GenericIon::GenericIon());
  }
 private: bool prompt_;
};
class Physics final:public G4VModularPhysicsList {
 public:
  Physics(double cut,bool prompt,bool correlated):cut_(cut) {
    SetVerboseLevel(0);
    RegisterPhysics(new G4EmStandardPhysics_option4(0));
    RegisterPhysics(new G4DecayPhysics(0));
    // Resolve the evaluated picosecond Ni-60 levels instead of assigning zero
    // lifetimes to levels below the default 1 ns nuclide-table threshold.
    G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(1.e-6*ns);
    auto* p=G4NuclearLevelData::GetInstance()->GetParameters();
    p->SetStoreICLevelData(true); p->SetIsomerProduction(true); p->SetInternalConversionFlag(true);
    p->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()/std::log(2.));
    p->SetCorrelatedGamma(correlated); p->SetTwoJMAX(20);
    RegisterPhysics(new SourceDecayPhysics(prompt));
    G4EmParameters::Instance()->SetFluo(true);
  }
  void SetCuts() override { SetDefaultCutValue(cut_); SetCutsWithDefault(); }
 private: double cut_;
};

class Detector final:public G4VUserDetectorConstruction {
 public:
  explicit Detector(fs::path input):input_(std::move(input)) {}
  G4VPhysicalVolume* Construct() override {
    auto* nist=G4NistManager::Instance();
    std::istringstream input(contents(input_/"materials.tsv")); std::string line;
    while (std::getline(input,line)) {
      std::istringstream row(line); char kind; std::string name; double density; row>>kind>>name>>density;
      if (kind=='B') {
        std::string base; row>>base; auto* m=nist->FindOrBuildMaterial(base);
        if (!m) throw std::runtime_error("missing NIST material "+base);
        new G4Material(name,density>0?density*g/cm3:m->GetDensity(),m);
      } else {
        int n; row>>n; auto* m=new G4Material(name,density*g/cm3,n);
        for (int i=0;i<n;++i) {
          std::string component; double fraction; row>>component>>fraction;
          if (kind=='E') m->AddElement(nist->FindOrBuildElement(component),fraction);
          else m->AddMaterial(G4Material::GetMaterial(component),fraction);
        }
      }
      if (!row) throw std::runtime_error("invalid material row "+line);
    }
    auto get=[](const std::string& name){return G4Material::GetMaterial(name);};
    auto* worldLV=new G4LogicalVolume(new G4Box("frame",5000,5000,5000),get("pmt_vacuum"),"frame");
    auto* world=new G4PVPlacement(nullptr,{},worldLV,"detector_frame_world",nullptr,false,0);
    // Optical surfaces remain geometry labels only. No optical processes are registered.
    std::map<std::string,G4OpticalSurface*> surfaces;
    for (const auto& name:{"ptfe","pvc_white"}) surfaces[name]=new G4OpticalSurface(name);
    auto placed=co60::BuildManifest((input_/"envelopes.tsv").string(),worldLV,"detector_frame_world",surfaces);
    G4LogicalVolume* tankLiquid=nullptr;
    for (auto* pv:placed) if (pv->GetName()=="detector") tankLiquid=pv->GetLogicalVolume();
    if (!tankLiquid) throw std::runtime_error("missing detector mother");
    double h,r,t,gradient,x,y,z; int faces;
    std::istringstream vessel(contents(input_/"vessel.txt")); vessel>>h>>r>>t>>gradient>>faces>>x>>y>>z;
    if (!vessel) throw std::runtime_error("bad vessel parameters");
    double zs[]={-h/2,-h/2+t,-h/2+t+gradient,h/2-t-gradient,h/2-t,h/2};
    double ri[]={0,0,r-t,r-t,0,0},ro[]={r,r,r,r,r,r};
    auto* nylon=new G4Polyhedra("nylonVessel",0,360*deg,faces,6,zs,ri,ro);
    double zi[]={-h/2,h/2},ii[]={0,0},oo[]={r,r};
    auto* outer=new G4Polyhedra("nv_inner",0,360*deg,faces,2,zi,ii,oo);
    auto* medium=new G4SubtractionSolid("detectionMedium",outer,nylon);
    auto* nvLV=new G4LogicalVolume(nylon,get("superSANDI_nylon"),"nv_logic");
    auto* liquidLV=new G4LogicalVolume(medium,get("wbls1pct_ly95_gdS0p2pct_SANDI"),"detmed_logic");
    placed.push_back(new G4PVPlacement(nullptr,{x,y,z},nvLV,"superSANDIvessel",tankLiquid,false,0));
    placed.push_back(new G4PVPlacement(nullptr,{x,y,z},liquidLV,"superSANDIDetMed",tankLiquid,false,0));
    auto assembly=co60::BuildManifest((input_/"assembly.tsv").string(),liquidLV,"superSANDIDetMed",surfaces);
    placed.insert(placed.end(),assembly.begin(),assembly.end());
    for (auto* pv:placed) if (pv->CheckOverlaps(20000,0,false,5)) throw std::runtime_error("geometry overlap "+pv->GetName());
    // Source-region production cut resolves the provisional 10 um active layer.
    auto* sourceRegion=new G4Region("source_region"); auto* cuts=new G4ProductionCuts;
    cuts->SetProductionCut(sourceCut); sourceRegion->SetProductionCuts(cuts);
    for (auto* pv:assembly) if (pv->GetName()=="Co60_disk") sourceRegion->AddRootLogicalVolume(pv->GetLogicalVolume());
    all=placed; all.insert(all.begin(),world);
    // Dump actual transport material densities/element mass fractions for audit.
    for (const auto* mat:*G4Material::GetMaterialTable()) {
      std::cout<<"MATERIAL "<<mat->GetName()<<" density_g_cm3="<<mat->GetDensity()/(g/cm3);
      for (size_t i=0;i<mat->GetNumberOfElements();++i)
        std::cout<<' '<<(*mat->GetElementVector())[i]->GetSymbol()<<':'<<mat->GetFractionVector()[i];
      std::cout<<'\n';
    }
    return world;
  }
  void Probes() {
    G4Navigator nav; nav.SetWorldVolume(all[0]);
    std::istringstream probes(contents(input_/"probes.tsv"));
    std::string label,expected,material; double x,y,z; int count=0;
    while (probes>>label>>x>>y>>z>>expected>>material) {
      auto* pv=nav.LocateGlobalPointAndSetup({x,y,z},nullptr,false);
      if (!pv || pv->GetName()!=expected || pv->GetLogicalVolume()->GetMaterial()->GetName()!=material)
        throw std::runtime_error("material probe failed "+label);
      ++count;
    }
    std::cout<<"TRANSPORT material/navigation probes passed="<<count<<'\n';
    if (count!=18) throw std::runtime_error("expected 18 probes");
  }
  std::vector<G4VPhysicalVolume*> all;
  double sourceCut=.001*mm;
 private: fs::path input_;
};

struct Gamma {
  int root_gamma_id=0,parent_id=0;
  double initial_MeV=0,dir_x=0,dir_y=0,dir_z=0,time_ns=0,x_mm=0,y_mm=0,z_mm=0;
  std::vector<double> edep_regions=std::vector<double>(regions.size(),0);
  double escaped_MeV=0;
  bool first_valid=false; int first_volume_id=-1,bgo_interactions=0;
  std::string first_process;
  double first_x_mm=0,first_y_mm=0,first_z_mm=0,first_time_ns=0,first_path_mm=0;
};
struct Deposit { double energy=0,x=0,y=0,z=0,first=0; int steps=0; };
struct State {
  Long64_t decay_id=0;
  double source_x_mm=0,source_y_mm=0,source_z_mm=0;
  double released_total_MeV=0,neutrino_emitted_MeV=0,emitted_gamma_MeV=0,emitted_beta_MeV=0;
  double escaped_visible_MeV=0,escaped_neutrino_MeV=0,total_edep_MeV=0,closure_MeV=0;
  double edep_bgo_MeV=0,edep_wbls_MeV=0,edep_gd_water_MeV=0,cos_gamma_pair=2;
  double beta_daughter_excitation_MeV=-1;
  int n_gamma=0,n_beta=0,n_nuclear_vertices=0,unknown_ancestry=0,n_tracks=0,source_valid=0;
  int source_ion_track_id=0;
  int bgo_contained_gammas=0;
  std::vector<double> edep_regions=std::vector<double>(regions.size(),0);
  std::map<int,int> lineage;
  std::map<int,Gamma> gammas;
  std::map<std::pair<int,int>,Deposit> deposits;
};

class Writer {
 public:
  Writer(const fs::path& path,const std::string& mode_,int n,long seed_,double cut_,const fs::path& input,bool correlation_)
    :mode(mode_),requested(n),seed(seed_),cut_mm(cut_),correlated(correlation_) {
    if (fs::exists(path)) throw std::runtime_error("refusing existing ROOT output");
    file.reset(TFile::Open(path.c_str(),"CREATE"));
    if (!file || file->IsZombie()) throw std::runtime_error("cannot create ROOT output");
    run_uuid=TUUID().AsString(); metadata_json=contents(input/"transport_provenance.json");
    events=new TTree("events","one row per generated decay, including zero-response events");
    events->Branch("run_uuid",&run_uuid);
#define E(field) events->Branch(#field,&state.field)
    E(decay_id); E(source_x_mm); E(source_y_mm); E(source_z_mm); E(released_total_MeV);
    E(neutrino_emitted_MeV); E(emitted_gamma_MeV); E(emitted_beta_MeV); E(escaped_visible_MeV);
    E(escaped_neutrino_MeV); E(total_edep_MeV); E(closure_MeV); E(edep_bgo_MeV); E(edep_wbls_MeV);
    E(edep_gd_water_MeV); E(edep_regions); E(cos_gamma_pair); E(beta_daughter_excitation_MeV);
    E(n_gamma); E(n_beta); E(n_nuclear_vertices); E(unknown_ancestry); E(n_tracks); E(source_valid); E(bgo_contained_gammas);
#undef E
    emissions=new TTree("emissions","nuclear-decay daughters or benchmark primaries; no transport secondaries");
    bindCommon(emissions); emissions->Branch("track_id",&track_id); emissions->Branch("parent_id",&parent_id);
    emissions->Branch("root_gamma_id",&root_id); emissions->Branch("pdg",&pdg); emissions->Branch("energy_MeV",&energy);
    emissions->Branch("time_ns",&time); emissions->Branch("x_mm",&x); emissions->Branch("y_mm",&y); emissions->Branch("z_mm",&z);
    emissions->Branch("dir_x",&dx); emissions->Branch("dir_y",&dy); emissions->Branch("dir_z",&dz);
    emissions->Branch("excitation_MeV",&excitation); emissions->Branch("creator",&process);
    gammaTree=new TTree("gamma_transport","nuclear gamma roots with descendant-deposit lineage");
    bindCommon(gammaTree);
#define G(field) gammaTree->Branch(#field,&gamma.field)
    G(root_gamma_id); G(parent_id); G(initial_MeV); G(dir_x); G(dir_y); G(dir_z); G(time_ns); G(x_mm); G(y_mm); G(z_mm);
    G(edep_regions); G(escaped_MeV); G(first_valid); G(first_volume_id); G(first_process);
    G(first_x_mm); G(first_y_mm); G(first_z_mm); G(first_time_ns); G(first_path_mm); G(bgo_interactions);
#undef G
    deposits=new TTree("depositions","non-optical deposits aggregated by decay/root/physical volume");
    bindCommon(deposits); deposits->Branch("root_gamma_id",&root_id); deposits->Branch("volume_id",&volume_id);
    deposits->Branch("edep_MeV",&energy); deposits->Branch("n_steps",&steps);
    deposits->Branch("weighted_x_mm",&x); deposits->Branch("weighted_y_mm",&y); deposits->Branch("weighted_z_mm",&z);
    deposits->Branch("first_time_ns",&time);
    escapes=new TTree("escapes","first exit from outer tank; track killed to define open boundary");
    bindCommon(escapes); escapes->Branch("root_gamma_id",&root_id); escapes->Branch("track_id",&track_id);
    escapes->Branch("pdg",&pdg); escapes->Branch("accounted_MeV",&energy); escapes->Branch("time_ns",&time);
    vertices=new TTree("decay_vertices","evaluated RDM energy release; prompt primary clock only");
    bindCommon(vertices); vertices->Branch("track_id",&track_id); vertices->Branch("pdg",&pdg);
    vertices->Branch("excitation_MeV",&excitation); vertices->Branch("release_MeV",&energy);
    vertices->Branch("time_ns",&time);
  }
  void bindCommon(TTree* tree) { tree->Branch("run_uuid",&run_uuid); tree->Branch("decay_id",&state.decay_id); }
  void SetupVolumes(Detector& detector) {
    auto* tree=new TTree("volumes","transport physical-volume/material/region dictionary");
    std::string name,material,region_name; double density; int region_id;
    tree->Branch("run_uuid",&run_uuid); tree->Branch("volume_id",&volume_id); tree->Branch("name",&name);
    tree->Branch("material",&material); tree->Branch("density_g_cm3",&density);
    tree->Branch("region_id",&region_id); tree->Branch("region_name",&region_name);
    int index=0;
    for (auto* pv:detector.all) {
      name=pv->GetName(); auto* m=pv->GetLogicalVolume()->GetMaterial(); material=m->GetName();
      density=m->GetDensity()/(g/cm3); region_id=region(name,material); region_name=regions[region_id];
      volume_id=index++; volumes[pv]=volume_id; volumeRegions[volume_id]=region_id; tree->Fill();
    }
    tree->Write(); delete tree;
  }
  void EndEvent() {
    state.total_edep_MeV=std::accumulate(state.edep_regions.begin(),state.edep_regions.end(),0.);
    state.edep_bgo_MeV=state.edep_regions[0]; state.edep_wbls_MeV=state.edep_regions[1];
    state.edep_gd_water_MeV=state.edep_regions[2];
    state.closure_MeV=state.released_total_MeV-state.neutrino_emitted_MeV-state.total_edep_MeV-state.escaped_visible_MeV;
    Gamma *lo=nullptr,*hi=nullptr;
    for (auto& [id,g]:state.gammas) {
      if (std::abs(g.initial_MeV-1.173239)<.001) lo=&g;
      if (std::abs(g.initial_MeV-1.332514)<.001) hi=&g;
      const bool mainLine=std::abs(g.initial_MeV-1.173239)<.001 || std::abs(g.initial_MeV-1.332514)<.001;
      if (mainLine && std::abs(g.edep_regions[0]-g.initial_MeV)<.001) ++state.bgo_contained_gammas;
      gamma=g; gammaTree->Fill();
    }
    if (lo && hi) state.cos_gamma_pair=lo->dir_x*hi->dir_x+lo->dir_y*hi->dir_y+lo->dir_z*hi->dir_z;
    for (auto& [key,d]:state.deposits) {
      root_id=key.first; volume_id=key.second; energy=d.energy; steps=d.steps;
      x=d.x/d.energy; y=d.y/d.energy; z=d.z/d.energy; time=d.first; deposits->Fill();
    }
    events->Fill(); ++completed;
    if (completed%10000==0) std::cout<<"COMPLETED "<<completed<<'\n';
  }
  void Finish() {
    file->cd(); events->Write(); emissions->Write(); gammaTree->Write(); deposits->Write(); escapes->Write(); vertices->Write();
    TTree runs("runs","standalone P1 truth metadata");
    double activity_Bq=13061; std::string version=G4Version,em="G4EmStandardPhysics_option4",unit="mm,ns,MeV; prompt-relative";
    std::string scope="P0 envelopes/assembly; PMTs/supports omitted; missing-core material fallbacks; no optics";
    double source_cut_mm=std::min(.001,cut_mm),decay_threshold_s=1.e15;
    double nuclide_half_life_threshold_ns=1.e-6;
    double mass_Q_MeV=(G4IonTable::GetIonTable()->GetIon(27,60,0.)->GetPDGMass()
                      -G4IonTable::GetIonTable()->GetIon(28,60,0.)->GetPDGMass()-electron_mass_c2)/MeV;
    bool primary_prompt_clock=mode!="natural_clock_test";
    int schema_version=2;
    auto region_names=regions; bool synthetic=false;
#define R(field) runs.Branch(#field,&field)
    R(run_uuid); R(mode); R(requested); R(completed); R(seed); R(cut_mm); R(activity_Bq); R(metadata_json);
    R(version); R(em); R(unit); R(scope); R(region_names); R(synthetic);
    R(correlated); R(primary_prompt_clock); R(source_cut_mm); R(decay_threshold_s);
    R(nuclide_half_life_threshold_ns); R(mass_Q_MeV);
    R(schema_version);
#undef R
    runs.Fill(); runs.Write(); file->Close();
  }
  State state;
  std::map<const G4VPhysicalVolume*,int> volumes; std::map<int,int> volumeRegions;
  std::string run_uuid,mode,metadata_json,process;
  TTree *events,*emissions,*gammaTree,*deposits,*escapes,*vertices;
  int track_id=0,parent_id=0,root_id=-1,pdg=0,volume_id=0,steps=0;
  double energy=0,time=0,x=0,y=0,z=0,dx=0,dy=0,dz=0,excitation=0;
 private:
  std::unique_ptr<TFile> file;
  Gamma gamma;
  int requested,completed=0; long seed; double cut_mm; bool correlated;
};

class Generator final:public G4VUserPrimaryGeneratorAction {
 public:
  Generator(Writer& w,const fs::path& input):w_(w),gun_(1) {
    std::istringstream f(contents(input/"source.tsv")); double x,y,z;
    f>>x>>y>>z>>radius_>>half_; center_={x,y,z}; if (!f) throw std::runtime_error("bad source region");
  }
  void GeneratePrimaries(G4Event* event) override {
    w_.state=State{}; auto& s=w_.state; s.decay_id=event->GetEventID();
    double r=radius_*std::sqrt(G4UniformRand()),phi=twopi*G4UniformRand();
    auto point=center_+G4ThreeVector(r*std::cos(phi),r*std::sin(phi),(2*G4UniformRand()-1)*half_);
    s.source_x_mm=point.x(); s.source_y_mm=point.y(); s.source_z_mm=point.z();
    auto* nav=G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    auto* pv=nav->LocateGlobalPointAndSetup(point,nullptr,false);
    if (!pv || pv->GetName()!="Co60_active") throw std::runtime_error("sampled source outside active volume");
    s.source_valid=1; gun_.SetParticlePosition(point); gun_.SetParticleTime(0);
    if (w_.mode=="full" || w_.mode=="natural_clock_test" || w_.mode=="ni60_timing_test") {
      auto* ion=w_.mode=="ni60_timing_test"?G4IonTable::GetIonTable()->GetIon(28,60,1332.514*keV)
                                           :G4IonTable::GetIonTable()->GetIon(27,60,0.);
      gun_.SetParticleDefinition(ion); gun_.SetParticleCharge(0.); gun_.SetParticleEnergy(0.);
      gun_.SetParticleMomentumDirection({0,0,1}); gun_.GeneratePrimaryVertex(event);
    } else {
      auto n1=G4RandomDirection(); double u;
      do { u=2*G4UniformRand()-1; } while (G4UniformRand()*(1+1./8+1./24)>1+u*u/8+u*u*u*u/24);
      double az=twopi*G4UniformRand(); auto a=n1.orthogonal().unit(); auto b=n1.cross(a);
      auto n2=u*n1+std::sqrt(1-u*u)*(std::cos(az)*a+std::sin(az)*b);
      gun_.SetParticleDefinition(G4Gamma::Gamma()); gun_.SetParticleEnergy(1.173*MeV);
      gun_.SetParticleMomentumDirection(n1); gun_.GeneratePrimaryVertex(event);
      gun_.SetParticleEnergy(1.332*MeV); gun_.SetParticleMomentumDirection(n2); gun_.GeneratePrimaryVertex(event);
      s.released_total_MeV=2.505;
    }
  }
 private: Writer& w_; G4ParticleGun gun_; G4ThreeVector center_; double radius_=0,half_=0;
};
class EventAction final:public G4UserEventAction {
 public: explicit EventAction(Writer& w):w_(w) {} void EndOfEventAction(const G4Event*) override { w_.EndEvent(); }
 private: Writer& w_;
};
class Tracking final:public G4UserTrackingAction {
 public: explicit Tracking(Writer& w):w_(w) {}
  void PreUserTrackingAction(const G4Track* t) override {
    auto& s=w_.state; int id=t->GetTrackID(),parent=t->GetParentID(),pdg=t->GetDefinition()->GetPDGEncoding();
    if (sourceIon(*t)) s.source_ion_track_id=id;
    bool nuclear=t->GetCreatorProcess() && t->GetCreatorProcess()->GetProcessName()=="RadioactiveDecay";
    bool benchmark=parent==0 && w_.mode=="pair";
    int root=-1;
    if (pdg==22 && (nuclear || benchmark)) {
      root=id; Gamma gamma; gamma.root_gamma_id=id; gamma.parent_id=parent;
      gamma.initial_MeV=t->GetKineticEnergy()/MeV;
      auto dir=t->GetMomentumDirection(),pos=t->GetPosition();
      gamma.dir_x=dir.x(); gamma.dir_y=dir.y(); gamma.dir_z=dir.z(); gamma.time_ns=t->GetGlobalTime()/ns;
      gamma.x_mm=pos.x()/mm; gamma.y_mm=pos.y()/mm; gamma.z_mm=pos.z()/mm;
      s.gammas[id]=gamma; ++s.n_gamma; s.emitted_gamma_MeV+=gamma.initial_MeV;
    } else if (nuclear || parent==0) {
      root=std::abs(pdg)==12?-4:(t->GetDefinition()->GetParticleType()=="nucleus"?-3:-2);
    } else {
      auto it=s.lineage.find(parent); if (it!=s.lineage.end()) root=it->second;
    }
    s.lineage[id]=root; ++s.n_tracks; if (root==-1) ++s.unknown_ancestry;
    if (nuclear || benchmark) {
      w_.track_id=id; w_.parent_id=parent; w_.root_id=root; w_.pdg=pdg;
      w_.energy=t->GetKineticEnergy()/MeV; w_.time=t->GetGlobalTime()/ns;
      auto dir=t->GetMomentumDirection(),pos=t->GetPosition();
      w_.dx=dir.x(); w_.dy=dir.y(); w_.dz=dir.z(); w_.x=pos.x()/mm; w_.y=pos.y()/mm; w_.z=pos.z()/mm;
      auto* ion=dynamic_cast<const G4Ions*>(t->GetDefinition());
      w_.excitation=ion?ion->GetExcitationEnergy()/MeV:0;
      w_.process=benchmark?"correlated_pair_benchmark":"RadioactiveDecay"; w_.emissions->Fill();
      if (std::abs(pdg)==12) s.neutrino_emitted_MeV+=w_.energy;
      if (pdg==11 && parent==s.source_ion_track_id) { ++s.n_beta; s.emitted_beta_MeV+=w_.energy; }
      if (ion && parent==s.source_ion_track_id) s.beta_daughter_excitation_MeV=w_.excitation;
    }
  }
 private: Writer& w_;
};
class Stepping final:public G4UserSteppingAction {
 public: explicit Stepping(Writer& w):w_(w) {}
  void UserSteppingAction(const G4Step* step) override {
    auto& s=w_.state; auto* t=step->GetTrack(); auto* pre=step->GetPreStepPoint(); auto* post=step->GetPostStepPoint();
    auto* pv=pre->GetPhysicalVolume(); auto* next=post->GetPhysicalVolume();
    int id=t->GetTrackID(),root=s.lineage.at(id),vid=w_.volumes.at(pv),reg=w_.volumeRegions.at(vid);
    double edep=step->GetTotalEnergyDeposit()/MeV;
    if (edep>0) {
      s.edep_regions[reg]+=edep; auto& d=s.deposits[{root,vid}];
      if (d.steps==0) d.first=pre->GetGlobalTime()/ns;
      ++d.steps; d.energy+=edep; auto p=(pre->GetPosition()+post->GetPosition())/2;
      d.x+=edep*p.x()/mm; d.y+=edep*p.y()/mm; d.z+=edep*p.z()/mm;
      if (root>=0) s.gammas.at(root).edep_regions[reg]+=edep;
    }
    auto* proc=post->GetProcessDefinedStep(); std::string process=proc?std::string(proc->GetProcessName()):"none";
    if (process=="RadioactiveDecay") {
      double released=edep-pre->GetKineticEnergy()/MeV;
      for (const auto* secondary:*step->GetSecondaryInCurrentStep()) released+=secondary->GetKineticEnergy()/MeV;
      s.released_total_MeV+=released; ++s.n_nuclear_vertices;
      w_.track_id=id; w_.pdg=t->GetDefinition()->GetPDGEncoding();
      auto* ion=dynamic_cast<const G4Ions*>(t->GetDefinition());
      w_.excitation=ion?ion->GetExcitationEnergy()/MeV:0.; w_.energy=released; w_.time=post->GetGlobalTime()/ns;
      w_.vertices->Fill();
    }
    if (root==id && t->GetDefinition()==G4Gamma::Gamma() && process!="Transportation" && process!="none") {
      auto& gamma=s.gammas.at(root);
      if (reg==0) ++gamma.bgo_interactions;
      if (!gamma.first_valid) {
        gamma.first_valid=true; gamma.first_volume_id=vid; gamma.first_process=process;
        auto pos=post->GetPosition(); gamma.first_x_mm=pos.x()/mm; gamma.first_y_mm=pos.y()/mm;
        gamma.first_z_mm=pos.z()/mm; gamma.first_time_ns=post->GetGlobalTime()/ns;
        gamma.first_path_mm=t->GetTrackLength()/mm;
      }
    }
    if (!next || (next->GetName()=="detector_frame_world" && pv->GetName()!="detector_frame_world")) {
      double energy=post->GetKineticEnergy()/MeV;
      // EM available-energy convention: escaped e+ retains the annihilation reservoir.
      if (t->GetDefinition()->GetPDGEncoding()==-11) energy+=2*electron_mass_c2/MeV;
      if (std::abs(t->GetDefinition()->GetPDGEncoding())==12) s.escaped_neutrino_MeV+=energy;
      else { s.escaped_visible_MeV+=energy; if (root>=0) s.gammas.at(root).escaped_MeV+=energy; }
      w_.root_id=root; w_.track_id=id; w_.pdg=t->GetDefinition()->GetPDGEncoding();
      w_.energy=energy; w_.time=post->GetGlobalTime()/ns; w_.escapes->Fill();
      t->SetTrackStatus(fStopAndKill);
    }
  }
 private: Writer& w_;
};

int main(int argc,char** argv) {
  try {
    if (argc!=8) throw std::runtime_error("usage: co60_transport INPUT_DIR full|pair|natural_clock_test|ni60_timing_test N SEED CUT_MM OUTPUT.root CORRELATED_0_OR_1");
    fs::path input=argv[1],output=argv[6]; std::string mode=argv[2];
    int n=std::stoi(argv[3]); long seed=std::stol(argv[4]); double cut=std::stod(argv[5]); bool correlated=std::stoi(argv[7]);
    if (n<1 || n>110000 || cut<=0 || (mode!="full" && mode!="pair" && mode!="natural_clock_test" && mode!="ni60_timing_test"))
      throw std::runtime_error("invalid bounded run configuration");
    if (!correlated) throw std::runtime_error("correlations must be enabled for Co-60 production/benchmark");
    if (mode=="natural_clock_test" && n>100) throw std::runtime_error("natural-clock diagnostic limited to 100 events");
    if (mode=="ni60_timing_test" && n>5000) throw std::runtime_error("Ni-60 diagnostic limited to 5000 events");
    std::cout<<std::setprecision(12)<<"Geant4 "<<G4Version<<" mode="<<mode<<" n="<<n<<" seed="<<seed<<" cut_mm="<<cut<<'\n';
    CLHEP::HepRandom::setTheSeed(seed);
    Writer writer(output,mode,n,seed,cut,input,correlated);
    auto manager=std::make_unique<G4RunManager>();
    auto* detector=new Detector(input); detector->sourceCut=std::min(.001,cut)*mm;
    manager->SetUserInitialization(detector);
    manager->SetUserInitialization(new Physics(cut*mm,mode!="natural_clock_test",correlated));
    manager->SetUserAction(new Generator(writer,input));
    manager->SetUserAction(new EventAction(writer)); manager->SetUserAction(new Tracking(writer));
    manager->SetUserAction(new Stepping(writer)); manager->Initialize(); detector->Probes(); writer.SetupVolumes(*detector);
    auto* co=G4IonTable::GetIonTable()->GetIon(27,60,0.);
    std::cout<<"Co60 physical mean_lifetime_s="<<co->GetPDGLifeTime()/s
             <<" correlated="<<G4NuclearLevelData::GetInstance()->GetParameters()->CorrelatedGamma()<<'\n';
    manager->BeamOn(n); writer.Finish();
    std::cout<<"SUCCESS ROOT complete: "<<output<<'\n'; return 0;
  } catch (const std::exception& e) { std::cerr<<"FAIL "<<e.what()<<'\n'; return 1; }
}
