#include "AnalyzerCore.h"
#include "FLUKANtuple.h"
#include "GEANT4Ntuple.h"

AnalyzerCore::AnalyzerCore(){
  MaxEvent = -1;
  NSkipEvent = 0;
  LogEvery = 1000;
  MCSample = "";
  Simulator = "";
  Userflags.clear();
  outfile = NULL;
  smear = new SmearParticles();
}

AnalyzerCore::~AnalyzerCore(){

  //=== hist maps
  
  for(std::map< TString, TH1D* >::iterator mapit = maphist_TH1D.begin(); mapit!=maphist_TH1D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH1D.clear();

  for(std::map< TString, TH2D* >::iterator mapit = maphist_TH2D.begin(); mapit!=maphist_TH2D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH2D.clear();

  for(std::map< TString, TH3D* >::iterator mapit = maphist_TH3D.begin(); mapit!=maphist_TH3D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH3D.clear();

  //==== output rootfile
  
  if(outfile) outfile->Close();
  delete outfile;

  //==== Tools
  delete smear;

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  cout << "[AnalyzerCore::~AnalyzerCore] JOB FINISHED " << printcurrunttime() << endl;

}
//==================
//==== Read Trees
//==================
Int_t AnalyzerCore::GetEntry(Long64_t entry)
{
  // Read contents of entry 
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void AnalyzerCore::Loop(){

  Long64_t nentries = fChain->GetEntries();
  if(MaxEvent>0){
    nentries = std::min(nentries,MaxEvent);
  }

  cout << "[AnalyzerCore::Loop] MaxEvent = " << MaxEvent << endl;
  cout << "[AnalyzerCore::Loop] NSkipEvent = " << NSkipEvent << endl;
  cout << "[AnalyzerCore::Loop] LogEvery = " << LogEvery << endl;
  cout << "[SKFlatNtuple::Loop] MCSample = " << MCSample << endl;
  cout << "[AnalyzerCore::Loop] Simulator = " << Simulator << endl;
  cout << "[AnalyzerCore::Loop] Userflags = {" << endl;
  for(unsigned int i=0; i<Userflags.size(); i++){
    cout << "[AnalyzerCore::Loop]   \"" << Userflags.at(i) << "\"," << endl;
  }
  cout << "[AnalyzerCore::Loop] }" << endl;

  cout << "[AnalyzerCore::Loop] Event Loop Started " << printcurrunttime() << endl;

  for(Long64_t jentry=0; jentry<nentries;jentry++){

    if(jentry<NSkipEvent){
      continue;
    }

    if(jentry%LogEvery==0){
      cout << "[AnalyzerCore::Loop RUNNING] " << jentry << "/" << nentries << " ("<<100.*jentry/nentries<<" %) @ " << printcurrunttime() << endl;
    }

    fChain->GetEntry(jentry);

    executeEvent();
  }

  cout << "[AnalyzerCore::Loop] LOOP END " << printcurrunttime() << endl;

}


//==== Attach the historams to ai different direcotry, not outfile
//==== We will write these histograms in WriteHist() to outfile
void AnalyzerCore::SwitchToTempDir(){

  gROOT->cd();
  TDirectory *tempDir = NULL;
  int counter = 0;
  while (!tempDir) {
    //==== First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << "AnalyzerCore" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    //==== Let's try to make this directory
    tempDir = gROOT->mkdir((dirname.str()).c_str());
  }
  tempDir->cd();

}

void AnalyzerCore::SetOutfilePath(TString outname){
  outfile = new TFile(outname,"RECREATE");
};


Event AnalyzerCore::GetEvent(){

  Event ev;
  ev.SetSimulator(Simulator);

  return ev;

}

std::vector<Gen> AnalyzerCore::GetAllParticles(){
  
  std::vector<Gen> out;
  // -- Save particle vector with using proper class for each Simulator (GEANT4, FLUKA ...)
  return out;
  
}

std::vector<Gen> AnalyzerCore::GetAllParticles_GEANT4(){

  std::vector<Gen> out;
  for(unsigned int i=0; i< this_GEANT4Ntuple.PDGcode->size(); i++){
    Gen current_particle;
    current_particle.SetPIDPosition(this_GEANT4Ntuple.PDGcode->at(i), this_GEANT4Ntuple.X->at(i), this_GEANT4Ntuple.Y->at(i), this_GEANT4Ntuple.Z->at(i));
    current_particle.SetPxPyPzE(this_GEANT4Ntuple.Px->at(i) / 1000., this_GEANT4Ntuple.Py->at(i) / 1000., this_GEANT4Ntuple.Pz->at(i) / 1000., this_GEANT4Ntuple.E->at(i) / 1000.); // MeV to GeV
    current_particle.SetInterType(this_GEANT4Ntuple.interType->at(i));

    out.push_back(current_particle);
  }

  return out;
}

std::vector<Gen> AnalyzerCore::GetAllParticles_FLUKA(){
  
  std::vector<Gen> out;
  // -- Push beam first
  if(this_FLUKANtuple.NIneHits >0){
    Gen current_beam;
    current_beam.SetPIDPosition(this_FLUKANtuple.IdIne[0], this_FLUKANtuple.PosIne[0][0], this_FLUKANtuple.PosIne[0][1], this_FLUKANtuple.PosIne[0][2]);
    current_beam.SetPxPyPzE(this_FLUKANtuple.PIne[0][0], this_FLUKANtuple.PIne[0][1], this_FLUKANtuple.PIne[0][2], this_FLUKANtuple.PIne[0][3]);
    current_beam.SetInterType(this_FLUKANtuple.TypeIne[0]);
    out.push_back(current_beam);
  }

  // -- Push secondary particles
  const unsigned int N_secondary_particles = this_FLUKANtuple.NSecIne[0];
  for(unsigned int i=0; i < N_secondary_particles; i++){
    Gen current_particle;
    current_particle.SetPIDPosition(this_FLUKANtuple.IdSecIne[i], this_FLUKANtuple.PosIne[0][0], this_FLUKANtuple.PosIne[0][1], this_FLUKANtuple.PosIne[0][2]);
    current_particle.SetPxPyPzE(this_FLUKANtuple.PSec[i][0], this_FLUKANtuple.PSec[i][1], this_FLUKANtuple.PSec[i][2], this_FLUKANtuple.PSec[i][3]);
    current_particle.SetInterType(this_FLUKANtuple.TypeIne[0]);

    out.push_back(current_particle);
  }

  return out;

}

std::vector<Gen> AnalyzerCore::GetPiplus(const std::vector<Gen>& particles, double min_P){

  std::vector<Gen> out;
  for(unsigned int i=0; i<particles.size(); i++){
    if(particles.at(i).PID() == 211 && particles.at(i).P() > min_P) out.push_back(particles.at(i));
  }

  return out;
}

std::vector<Gen> AnalyzerCore::GetPiminus(const std::vector<Gen>& particles, double min_P){

  std::vector<Gen> out;
  for(unsigned int i=0; i<particles.size(); i++){
    if(particles.at(i).PID() == -211 && particles.at(i).P() > min_P) out.push_back(particles.at(i));
  }

  return out;
}

std::vector<Gen> AnalyzerCore::GetProtons(const std::vector<Gen>& particles, double min_P){

  std::vector<Gen> out;
  for(unsigned int i=0; i<particles.size(); i++){
    if(particles.at(i).PID() == 2212 && particles.at(i).P() > min_P) out.push_back(particles.at(i));
  }

  std::sort(out.begin(), out.end(), PComparing);
  return out;
}

std::vector<Gen> AnalyzerCore::GetNeutrons(const std::vector<Gen>& particles, double min_P){

  std::vector<Gen> out;
  for(unsigned int i=0; i<particles.size(); i++){
    if(particles.at(i).PID() == 2112 && particles.at(i).P() > min_P) out.push_back(particles.at(i));
  }

  std::sort(out.begin(), out.end(), PComparing);
  return out;
}

std::vector<Gen> AnalyzerCore::GetPizeros(const std::vector<Gen>& particles, double min_P){

  std::vector<Gen> out;
  for(unsigned int i=0; i<particles.size(); i++){
    if(particles.at(i).PID() == 111 && particles.at(i).P() > min_P) out.push_back(particles.at(i));
  }

  std::sort(out.begin(), out.end(), PComparing);
  return out;
}

std::vector<Gen> AnalyzerCore::GetBkgParticles(const std::vector<Gen>& particles, double min_P){
  // == Collect : pi- (-111), Kaon +- (+- 321), K0 (321), muon (+- 13)

  std::vector<Gen> out;
  for(unsigned int i=0; i<particles.size(); i++){
    if(particles.at(i).PID() == -211 || abs(particles.at(i).PID()) == 321 || abs(particles.at(i).PID()) == 13){ 
      if(particles.at(i).P() > min_P) out.push_back(particles.at(i));
    }
  }

  return out;
}

std::vector<Gen> AnalyzerCore::GetNuclei(const std::vector<Gen>& particles){

  std::vector<Gen> out;
  if(Simulator.Contains("GEANT")){
    for(unsigned int i=0; i<particles.size(); i++){
      if(particles.at(i).PID() > 1000000000) out.push_back(particles.at(i));
    }
  }
  if(Simulator.Contains("FLUKA")){
    for(unsigned int i=0; i<particles.size(); i++){
      if(particles.at(i).PID() > 9000) out.push_back(particles.at(i));
    }
  }

  return out;
}
//==================
// Initialize
//==================
void AnalyzerCore::initializeAnalyzerTools(){

}

void AnalyzerCore::Init(){
  
  cout << "Let initiallize!" << endl;
  if(Simulator.Contains("GEANT")){
    this_GEANT4Ntuple.Init_GEANT4(fChain);
  }
  else if(Simulator.Contains("FLUKA")){
    this_FLUKANtuple.Init_FLUKA(fChain);
  }
  else return;

}

//==================
// Functions for Particles
//==================
int AnalyzerCore::GetAtomicNumber(int pid){

  int out;

  if(pid < 1000000000) return 0;
  
  out = (pid - 1000000000) / 10000;

  return out;

}

int AnalyzerCore::GetAtomicMass(int pid){

  int out;
  
  if(pid < 1000000000) return 0;

  int current_atomic_number = GetAtomicNumber(pid);
  
  //cout << "[[AnalyzerCore::GetAtomicMass]] pid - 1000000000 : " << pid - 1000000000 << endl;
  //cout << "[[AnalyzerCore::GetAtomicMass]] pid - 1000000000 current_atomic_number * 10000 : "<< pid - 1000000000 - current_atomic_number * 10000<< endl;

  out = (pid - 1000000000 - current_atomic_number * 10000) / 10;

  return out;

}

//==================
//==== Plotting
//==================
TH1D* AnalyzerCore::GetHist1D(TString histname){

  TH1D *h = NULL;
  std::map<TString, TH1D*>::iterator mapit = maphist_TH1D.find(histname);
  if(mapit != maphist_TH1D.end()) return mapit->second;

  return h;

}

TH2D* AnalyzerCore::GetHist2D(TString histname){

  TH2D *h = NULL;
  std::map<TString, TH2D*>::iterator mapit = maphist_TH2D.find(histname);
  if(mapit != maphist_TH2D.end()) return mapit->second;

  return h;

}

TH3D* AnalyzerCore::GetHist3D(TString histname){

  TH3D *h = NULL;
  std::map<TString, TH3D*>::iterator mapit = maphist_TH3D.find(histname);
  if(mapit != maphist_TH3D.end()) return mapit->second;

  return h;

}

void AnalyzerCore::FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

void AnalyzerCore::FillHist(TString histname, double value, double weight, int n_bin, double *xbins){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, xbins);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y,
			    double weight,
			    int n_binx, double x_min, double x_max,
			    int n_biny, double y_min, double y_max){

  TH2D *this_hist = GetHist2D(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    this_hist->SetDirectory(NULL);
    maphist_TH2D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y,
			    double weight,
			    int n_binx, double *xbins,
			    int n_biny, double *ybins){

  TH2D *this_hist = GetHist2D(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
    this_hist->SetDirectory(NULL);
    maphist_TH2D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y, double value_z,
			    double weight,
			    int n_binx, double x_min, double x_max,
			    int n_biny, double y_min, double y_max,
			    int n_binz, double z_min, double z_max){

  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    this_hist = new TH3D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y, double value_z,
			    double weight,
			    int n_binx, double *xbins,
			    int n_biny, double *ybins,
			    int n_binz, double *zbins){

  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    this_hist = new TH3D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, weight);

}

TH1D* AnalyzerCore::JSGetHist1D(TString suffix, TString histname){

  TH1D *h = NULL;

  std::map< TString, std::map<TString, TH1D*> >::iterator mapit = JSmaphist_TH1D.find(suffix);
  if(mapit==JSmaphist_TH1D.end()){
    return h;
  }
  else{

    std::map<TString, TH1D*> this_maphist = mapit->second;
    std::map<TString, TH1D*>::iterator mapitit = this_maphist.find(histname);
    if(mapitit != this_maphist.end()) return mapitit->second;

  }

  return h;

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = JSGetHist1D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    (JSmaphist_TH1D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value, weight);

}

TH2D* AnalyzerCore::JSGetHist2D(TString suffix, TString histname){

  TH2D *h = NULL;

  std::map< TString, std::map<TString, TH2D*> >::iterator mapit = JSmaphist_TH2D.find(suffix);
  if(mapit==JSmaphist_TH2D.end()){
    return h;
  }
  else{

    std::map<TString, TH2D*> this_maphist = mapit->second;
    std::map<TString, TH2D*>::iterator mapitit = this_maphist.find(histname);
    if(mapitit != this_maphist.end()) return mapitit->second;

  }

  return h;

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname,
			      double value_x, double value_y,
			      double weight,
			      int n_binx, double x_min, double x_max,
			      int n_biny, double y_min, double y_max){

  TH2D *this_hist = JSGetHist2D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH2D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    (JSmaphist_TH2D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname,
			      double value_x, double value_y,
			      double weight,
			      int n_binx, double *xbins,
			      int n_biny, double *ybins){

  TH2D *this_hist = JSGetHist2D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
    (JSmaphist_TH2D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::WriteHist(){

  outfile->cd();
  for(std::map< TString, TH1D* >::iterator mapit = maphist_TH1D.begin(); mapit!=maphist_TH1D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
  for(std::map< TString, TH2D* >::iterator mapit = maphist_TH2D.begin(); mapit!=maphist_TH2D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
  for(std::map< TString, TH3D* >::iterator mapit = maphist_TH3D.begin(); mapit!=maphist_TH3D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }

  outfile->cd();
  for(std::map< TString, std::map<TString, TH1D*> >::iterator mapit=JSmaphist_TH1D.begin(); mapit!=JSmaphist_TH1D.end(); mapit++){

    TString this_suffix = mapit->first;
    std::map< TString, TH1D* > this_maphist = mapit->second;


    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);

    for(std::map< TString, TH1D* >::iterator mapit = this_maphist.begin(); mapit!=this_maphist.end(); mapit++){
      mapit->second->Write();
    }

    outfile->cd();

  }

  for(std::map< TString, std::map<TString, TH2D*> >::iterator mapit=JSmaphist_TH2D.begin(); mapit!=JSmaphist_TH2D.end(); mapit++){

    TString this_suffix = mapit->first;
    std::map< TString, TH2D* > this_maphist = mapit->second;

    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);

    for(std::map< TString, TH2D* >::iterator mapit = this_maphist.begin(); mapit!=this_maphist.end(); mapit++){
      mapit->second->Write();
    }

    outfile->cd();

  }
}
