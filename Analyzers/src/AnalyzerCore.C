#include "AnalyzerCore.h"
#include "PDSPTree.h"

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

std::vector<Daughter> AnalyzerCore::GetAllDaughters(){


  // ==== Beam related variables                                                             
  TVector3 reco_beam_end(evt.reco_beam_endX, evt.reco_beam_endY, evt.reco_beam_endZ);
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;

  vector<Daughter> out;
  //cout << "[PionXsec::GetAllRecoDaughters] evt.reco_daughter_allTrack_ID->size() : " << evt.reco_daughter_allTrack_ID->size() << endl;
  for (size_t i = 0; i < evt.reco_daughter_allTrack_ID->size(); i++){
    if((*evt.reco_daughter_allTrack_dQdX_SCE)[i].empty()) continue;

    Daughter this_Daughter;
    //cout << "[PionXsec::GetAllDaughters] i : " << i << endl;
    this_Daughter.SetIsEmpty(false);
    if(evt.MC){
      this_Daughter.Set_PFP_true_byHits_PDG((*evt.reco_daughter_PFP_true_byHits_PDG).at(i));
      this_Daughter.Set_PFP_true_byHits_ID((*evt.reco_daughter_PFP_true_byHits_ID).at(i));
      this_Daughter.Set_PFP_true_byHits_origin((*evt.reco_daughter_PFP_true_byHits_origin).at(i));
      this_Daughter.Set_PFP_true_byHits_parID((*evt.reco_daughter_PFP_true_byHits_parID).at(i));
      this_Daughter.Set_PFP_true_byHits_parPDG((*evt.reco_daughter_PFP_true_byHits_parPDG).at(i));
      this_Daughter.Set_PFP_true_byHits_process((*evt.reco_daughter_PFP_true_byHits_process).at(i));
      this_Daughter.Set_PFP_true_byHits_sharedHits((*evt.reco_daughter_PFP_true_byHits_sharedHits).at(i));
      this_Daughter.Set_PFP_true_byHits_emHits((*evt.reco_daughter_PFP_true_byHits_emHits).at(i));
      this_Daughter.Set_PFP_true_byHits_len((*evt.reco_daughter_PFP_true_byHits_len).at(i));
      this_Daughter.Set_PFP_true_byHits_startX((*evt.reco_daughter_PFP_true_byHits_startX).at(i));
      this_Daughter.Set_PFP_true_byHits_startY((*evt.reco_daughter_PFP_true_byHits_startY).at(i));
      this_Daughter.Set_PFP_true_byHits_startZ((*evt.reco_daughter_PFP_true_byHits_startZ).at(i));
      this_Daughter.Set_PFP_true_byHits_endX((*evt.reco_daughter_PFP_true_byHits_endX).at(i));
      this_Daughter.Set_PFP_true_byHits_endY((*evt.reco_daughter_PFP_true_byHits_endY).at(i));
      this_Daughter.Set_PFP_true_byHits_endZ((*evt.reco_daughter_PFP_true_byHits_endZ).at(i));
      this_Daughter.Set_PFP_true_byHits_startPx((*evt.reco_daughter_PFP_true_byHits_startPx).at(i));
      this_Daughter.Set_PFP_true_byHits_startPy((*evt.reco_daughter_PFP_true_byHits_startPy).at(i));
      this_Daughter.Set_PFP_true_byHits_startPz((*evt.reco_daughter_PFP_true_byHits_startPz).at(i));
      this_Daughter.Set_PFP_true_byHits_startP((*evt.reco_daughter_PFP_true_byHits_startP).at(i));
      this_Daughter.Set_PFP_true_byHits_startE((*evt.reco_daughter_PFP_true_byHits_startE).at(i));
      this_Daughter.Set_PFP_true_byHits_endProcess((*evt.reco_daughter_PFP_true_byHits_endProcess).at(i));
      this_Daughter.Set_PFP_true_byHits_purity((*evt.reco_daughter_PFP_true_byHits_purity).at(i));
      this_Daughter.Set_PFP_true_byHits_completeness((*evt.reco_daughter_PFP_true_byHits_completeness).at(i));
      this_Daughter.Set_PFP_true_byE_PDG((*evt.reco_daughter_PFP_true_byE_PDG).at(i));
      this_Daughter.Set_PFP_true_byE_len((*evt.reco_daughter_PFP_true_byE_len).at(i));
    }

    this_Daughter.Set_PFP_ID((*evt.reco_daughter_PFP_ID).at(i));
    this_Daughter.Set_PFP_nHits((*evt.reco_daughter_PFP_nHits).at(i));
    this_Daughter.Set_PFP_nHits_collection((*evt.reco_daughter_PFP_nHits_collection).at(i));
    this_Daughter.Set_PFP_trackScore((*evt.reco_daughter_PFP_trackScore).at(i));
    this_Daughter.Set_PFP_emScore((*evt.reco_daughter_PFP_emScore).at(i));
    this_Daughter.Set_PFP_michelScore((*evt.reco_daughter_PFP_michelScore).at(i));
    this_Daughter.Set_PFP_trackScore_collection((*evt.reco_daughter_PFP_trackScore_collection).at(i));
    this_Daughter.Set_PFP_emScore_collection((*evt.reco_daughter_PFP_emScore_collection).at(i));
    this_Daughter.Set_PFP_michelScore_collection((*evt.reco_daughter_PFP_michelScore_collection).at(i));
    this_Daughter.Set_allTrack_ID((*evt.reco_daughter_allTrack_ID).at(i));
    this_Daughter.Set_allTrack_EField_SCE((*evt.reco_daughter_allTrack_EField_SCE).at(i));
    this_Daughter.Set_allTrack_resRange_SCE((*evt.reco_daughter_allTrack_resRange_SCE).at(i));
    this_Daughter.Set_allTrack_resRange_SCE_plane0((*evt.reco_daughter_allTrack_resRange_plane0).at(i)); // == FIXME, to SCE
    this_Daughter.Set_allTrack_resRange_SCE_plane1((*evt.reco_daughter_allTrack_resRange_plane1).at(i)); // == FIXME, to SCE
    this_Daughter.Set_allTrack_calibrated_dEdX_SCE((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE).at(i));
    this_Daughter.Set_allTrack_calibrated_dEdX_SCE_plane0((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE_plane0).at(i));
    this_Daughter.Set_allTrack_calibrated_dEdX_SCE_plane1((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE_plane1).at(i));
    this_Daughter.Set_allTrack_Chi2_proton((*evt.reco_daughter_allTrack_Chi2_proton).at(i));
    this_Daughter.Set_allTrack_Chi2_pion((*evt.reco_daughter_allTrack_Chi2_pion).at(i));
    this_Daughter.Set_allTrack_Chi2_muon((*evt.reco_daughter_allTrack_Chi2_muon).at(i));
    this_Daughter.Set_allTrack_Chi2_ndof((*evt.reco_daughter_allTrack_Chi2_ndof).at(i));
    this_Daughter.Set_allTrack_Chi2_ndof_pion((*evt.reco_daughter_allTrack_Chi2_ndof_pion).at(i));
    this_Daughter.Set_allTrack_Chi2_ndof_muon((*evt.reco_daughter_allTrack_Chi2_ndof_muon).at(i));
    this_Daughter.Set_allTrack_Theta((*evt.reco_daughter_allTrack_Theta).at(i));
    this_Daughter.Set_allTrack_Phi((*evt.reco_daughter_allTrack_Phi).at(i));
    // == FIXME
    //this_Daughter.Set_allTrack_startDirX((*evt.reco_daughter_allTrack_startDirX).at(i));
    //this_Daughter.Set_allTrack_startDirY((*evt.reco_daughter_allTrack_startDirY).at(i));
    //this_Daughter.Set_allTrack_startDirZ((*evt.reco_daughter_allTrack_startDirZ).at(i));
    this_Daughter.Set_allTrack_alt_len((*evt.reco_daughter_allTrack_alt_len).at(i));
    this_Daughter.Set_allTrack_startX((*evt.reco_daughter_allTrack_startX).at(i));
    this_Daughter.Set_allTrack_startY((*evt.reco_daughter_allTrack_startY).at(i));
    this_Daughter.Set_allTrack_startZ((*evt.reco_daughter_allTrack_startZ).at(i));
    this_Daughter.Set_allTrack_endX((*evt.reco_daughter_allTrack_endX).at(i));
    this_Daughter.Set_allTrack_endY((*evt.reco_daughter_allTrack_endY).at(i));
    this_Daughter.Set_allTrack_endZ((*evt.reco_daughter_allTrack_endZ).at(i));
    this_Daughter.Set_allTrack_vertex_michel_score((*evt.reco_daughter_allTrack_vertex_michel_score).at(i));
    this_Daughter.Set_allTrack_vertex_nHits((*evt.reco_daughter_allTrack_vertex_nHits).at(i));
    this_Daughter.Set_pandora_type((*evt.reco_daughter_pandora_type).at(i));

    TVector3 unit_daughter((*evt.reco_daughter_allTrack_endX).at(i) - (*evt.reco_daughter_allTrack_startX).at(i),
                           (*evt.reco_daughter_allTrack_endY).at(i) - (*evt.reco_daughter_allTrack_startY).at(i),
                           (*evt.reco_daughter_allTrack_endZ).at(i) - (*evt.reco_daughter_allTrack_startZ).at(i) );
    unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;
    double cos_theta = cos(unit_daughter.Angle(reco_unit_beam));
    this_Daughter.Set_Beam_Cos(cos_theta);

    TVector3 reco_daughter_start((*evt.reco_daughter_allTrack_startX).at(i), (*evt.reco_daughter_allTrack_startY).at(i), (*evt.reco_daughter_allTrack_startZ).at(i) );
    double dist_beam_end = (reco_daughter_start - reco_beam_end).Mag();
    this_Daughter.Set_Beam_Dist(dist_beam_end);

    out.push_back(this_Daughter);
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetPions(const vector<Daughter> in){

  vector<Daughter> out;

  double cut_cos_beam = 0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 60.;
  double cut_startZ = 220.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 > cut_chi2_proton
       && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist && this_in.allTrack_startZ() < cut_startZ){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetProtons(const vector<Daughter> in){

  vector<Daughter> out;

  double cut_cos_beam =0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 50.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() /this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 < cut_chi2_proton && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetTruePions(const vector<Daughter> in){

  vector<Daughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(abs(this_true_PID) == 211){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetTrueProtons(const vector<Daughter> in){

  vector<Daughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(this_true_PID == 2212){
      out.push_back(this_in);
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
