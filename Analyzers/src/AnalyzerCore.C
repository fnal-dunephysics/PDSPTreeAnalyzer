#include "AnalyzerCore.h"
#include "PDSPTree.h"

AnalyzerCore::AnalyzerCore(){
  MaxEvent = -1;
  NSkipEvent = 0;
  LogEvery = 1000;
  MCSample = "";
  Beam_Momentum = -1.;
  Userflags.clear();
  outfile = NULL;
  MCCorr = new MCCorrection();
  G4Xsec = new GEANT4_XS();
  map_BB[13] = new BetheBloch(13);
  map_BB[211] = new BetheBloch(211);
  map_BB[321] = new BetheBloch(321);
  map_BB[2212] = new BetheBloch(2212);
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
  delete MCCorr;
  delete G4Xsec;
  //delete map_BB;

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
  cout << "[AnalyzerCore::Loop] MCSample = " << MCSample << endl;
  cout << "[AnalyzerCore::Loop] Beam_Momentum = " << Beam_Momentum << endl;
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
    Init_evt();
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
  ev.SetBeam_Momentum(Beam_Momentum);

  return ev;

}

//==================
// Get Daughters
//==================
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

std::vector<Daughter> AnalyzerCore::GetPions(const vector<Daughter>& in){

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

std::vector<Daughter> AnalyzerCore::GetProtons(const vector<Daughter>& in){

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

std::vector<Daughter> AnalyzerCore::GetTruePions(const vector<Daughter>& in){

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

std::vector<Daughter> AnalyzerCore::GetTrueProtons(const vector<Daughter>& in){

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
// Set Beam Variables
//==================
void AnalyzerCore::SetPandoraSlicePDG(int pdg){
  pandora_slice_pdg = pdg;
};

int AnalyzerCore::GetPiParType(){

  if (!evt.MC){
    return pi::kData;
  }
  else if (evt.event%2){ // divide half of MC as fake data
    return pi::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){ // the true beam track is not selected
    if (evt.reco_beam_true_byE_origin == 2) {
      return pi::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){ // the selected track is a pion (but not true beam pion, so it is a secondary pion)
      return pi::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return pi::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return pi::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return pi::kMIDeg;
    }
    else {
      return pi::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == -13){
    return pi::kMuon;
  }
  else if (evt.true_beam_PDG == 211){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      return pi::kPiInel;
    }
    else return pi::kPiElas;
  }

  return pi::kMIDother;
}

//==================
// Event Selections
//==================
bool AnalyzerCore::Pass_Beam_PID(int PID){

  // == Study only proton and pion/muon beams
  if(PID != 2212 && PID != 211 && abs(PID) != 13) return false;

  vector<int> PDGID_candiate_vec;
  PDGID_candiate_vec.clear();
  if(PID == 2212){
    PDGID_candiate_vec.push_back(2212);
  }

  if(PID == 211 || PID == 13){
    PDGID_candiate_vec.push_back(211);
    PDGID_candiate_vec.push_back(13);
    PDGID_candiate_vec.push_back(-13);
  }

  if (evt.reco_reconstructable_beam_event == 0) return false; // First, remove empty events
  if (evt.MC){
    for (size_t i = 0; i < PDGID_candiate_vec.size(); ++i){
      if (evt.true_beam_PDG == PDGID_candiate_vec[i]) return true; // Truth matched
    }
    return false;
  }
  else{ // data
    if (evt.beam_inst_trigger == 8) return false; // Is cosmics
    if (evt.beam_inst_nMomenta != 1 || evt.beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<PDGID_candiate_vec.size(); ++i){
      for (size_t j = 0; j<evt.beam_inst_PDG_candidates->size(); ++j){
        if ((*evt.beam_inst_PDG_candidates)[j] == PDGID_candiate_vec[i]) return true; // Matching data PID result and PDGID candidates
      }
    }
    return false;
  }
}

bool AnalyzerCore::PassBeamScraperCut() const{

  // == Define beam plug circle in [cm] unit
  double center_x = 0.;
  double center_y = 0.;
  double radius = 0.;
  double N_sigma = 0.;
  if(evt.MC){
    center_x = -29.6;
    center_y = 422.;
    radius = 4.8;
    N_sigma = 1.4;
  }
  else{
    center_x = -32.16;
    center_y = 422.7;
    radius = 4.8;
    N_sigma = 1.2;
  }

  double this_distance = sqrt( pow(evt.beam_inst_X - center_x , 2.) + pow(evt.beam_inst_Y - center_y, 2.) );
  bool out = this_distance < radius * N_sigma;

  return out;
}

bool AnalyzerCore::PassBeamMomentumWindowCut() const{

  bool out = false;
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  out = ((P_beam_inst > beam_momentum_low) && (P_beam_inst < beam_momentum_high));
  return out;
}

bool AnalyzerCore::PassPandoraSliceCut() const{
  bool out = false;
  out = (evt.reco_beam_type == pandora_slice_pdg);
  return out;
}

bool AnalyzerCore::PassCaloSizeCut() const{ // Require hits information in the collection plane
  bool out = false;
  out = !(evt.reco_beam_calo_wire->empty());
  return out;
}

bool AnalyzerCore::PassAPA3Cut(double cut){ // only use track in the first TPC
  bool out = false;
  out = evt.reco_beam_calo_endZ < cut;
  return out;
}

bool AnalyzerCore::PassMichelScoreCut(double cut){ // further veto muon tracks according to Michel score
  bool out = false;
  out = daughter_michel_score < cut;
  return out;
}

bool AnalyzerCore::PassBeamCosCut(double cut){
  bool out = false;
  out = beam_costh > cut;
  return out;
}

bool AnalyzerCore::PassBeamStartZCut(double cut){
  bool out = false;
  out = evt.reco_beam_calo_startZ > cut;
  return out;
}

bool AnalyzerCore::PassProtonVetoCut(double cut){ // to remove proton background
  bool out = false;
  out = chi2_proton > cut;
  return out;
}

bool AnalyzerCore::PassMuonVetoCut(double cut){ // to remove muon background
  bool out = false;
  out = chi2_muon > cut;
  return out;
}

bool AnalyzerCore::PassStoppedPionVetoCut(double cut){ // to remove stopped pion background
  bool out = false;
  out = chi2_pion > cut;
  return out;
}

//==================
// Initialize
//==================
void AnalyzerCore::initializeAnalyzerTools(){

}

void AnalyzerCore::Init(){

  cout << "Let initiallize!" << endl;
  evt.Init_PDSPTree(fChain);

  // == Additional Root files
  TString datapath = getenv("DATA_DIR");
  TFile *file_profile = TFile::Open(datapath + "/dEdx_profiles/dEdxrestemplates.root");
  map_profile[13] = (TProfile *)file_profile -> Get("dedx_range_mu");
  map_profile[211] = (TProfile *)file_profile -> Get("dedx_range_pi");
  map_profile[321] = (TProfile *)file_profile -> Get("dedx_range_ka");
  map_profile[2212] = (TProfile *)file_profile -> Get("dedx_range_pro");
  cout << "[[AnalyzerCore::Init]] Called Profiles" << endl;

  // == Beam Window cut
  if(evt.MC) P_beam_inst_scale = Beam_Momentum;
  beam_momentum_low = Beam_Momentum * 1000. * 0.8;
  beam_momentum_high = Beam_Momentum * 1000. * 1.2;
  cout << "[[AnalyzerCore::Init]] Called beam window cuts" << endl;

  // == Pandora
  SetPandoraSlicePDG(13);
  cout << "[[AnalyzerCore::Init]] Set Pandora sllice PDG" << endl;

}

void AnalyzerCore::Init_evt(){
  daughter_michel_score = -999.;
  beam_costh = -999;
  chi2_proton = -1.;
  if (!evt.reco_beam_calo_wire->empty()){
    if (evt.reco_beam_vertex_nHits) daughter_michel_score = evt.reco_beam_vertex_michel_score_weight_by_charge;

    TVector3 pt0(evt.reco_beam_calo_startX,
                 evt.reco_beam_calo_startY,
                 evt.reco_beam_calo_startZ);
    TVector3 pt1(evt.reco_beam_calo_endX,
                 evt.reco_beam_calo_endY,
                 evt.reco_beam_calo_endZ);
    TVector3 dir = pt1 - pt0;
    dir = dir.Unit();
    if (evt.MC){
      TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180),
                       cos(beam_angleY_mc*TMath::Pi()/180),
                       cos(beam_angleZ_mc*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }

    else{
      TVector3 beamdir(cos(beam_angleX_data*TMath::Pi()/180),
                       cos(beam_angleY_data*TMath::Pi()/180),
                       cos(beam_angleZ_data*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }

    chi2_proton = evt.reco_beam_Chi2_proton/evt.reco_beam_Chi2_ndof;
    chi2_pion = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 211);
    chi2_muon = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), true, 13);
  }
}

//==================
// Additional Functions
//==================
double AnalyzerCore::Convert_P_Spectrometer_to_P_ff(double P_beam_inst, TString particle, TString key, int syst){

  double out = P_beam_inst;
  double delta_P = 0.;
  double p0 = 0., p1 = 0., p2 = 0.;

  if(particle.Contains("pion")){
    if(key == "AllTrue"){
      if(syst == 0){
        p0 = 182.7;
        p1 = -0.4893;
        p2 = 0.0003186;
      }
      else if(syst == -1){
        p0 = 140.2;
        p1 = -0.3979;
        p2 = 0.0002678;
      }
      else if(syst == 1){
        p0 = 235.2;
        p1 = -0.6022;
        p2 = 0.0003807;
      }
      else{
        return out;
      }
    }
  }
  else{
    return out;
  }

  delta_P = p0 + p1 * out + p2 * out * out;

  return out - delta_P;
}

double AnalyzerCore::Particle_chi2(const vector<double> & dEdx, const vector<double> & ResRange, int PID, double dEdx_res_frac){

  //cout << "[AnalyzerCore::Particle_chi2] Start" << endl; 
  if(PID != 2212 && PID != 13 && PID != 211){
    return 99999.;
  }
  if( dEdx.size() < 1 || ResRange.size() < 1 ) return 88888.;

  int N_max_hits = 1000;
  int this_N_calo = dEdx.size();
  int this_N_hits = min(N_max_hits, this_N_calo);
  int N_skip = 1;
  double dEdx_truncate_upper = 1000.;
  double dEdx_truncate_bellow = 0.;
  double this_chi2 = 0.;
  int npt = 0;
  for(int j = N_skip; j < this_N_hits - N_skip; j++){

    double dEdx_measured = dEdx.at(j);
    if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue;                                                                                              

    double this_res_length = ResRange.at(j);
    int bin = map_profile[PID] -> FindBin( this_res_length );
    if( bin >= 1 && bin <= map_profile[PID]->GetNbinsX() ){
      double template_dedx = map_profile[PID]->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( map_profile[PID]->GetBinContent( bin - 1 ) + map_profile[PID]->GetBinContent( bin + 1 ) ) / 2.;
      }

      double template_dedx_err = map_profile[PID]->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( map_profile[PID]->GetBinError( bin - 1 ) + map_profile[PID]->GetBinError( bin + 1 ) ) / 2.;
      }

      double dedx_res = 0.04231 + 0.0001783 * dEdx_measured * dEdx_measured;
      dedx_res *= dEdx_measured * dEdx_res_frac;

      this_chi2 += ( pow( (dEdx_measured - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }
  if(npt == 0) return 77777.;

  return this_chi2 / npt;
}

double AnalyzerCore::Fit_HypTrkLength_Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam){

  //cout << "[HadAna::Fit_dEdx_Residual_Length] Start" << endl;
  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }

  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  int N_skip = 3;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 15){
    //cout << "[HadAna::Fit_dEdx_Residual_Length] Too small number of hits!" << endl;
    return -9999.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  vector<double> chi2_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits
      double this_res_length = ResRange.at(j) + this_additional_res_length;

      double this_KE = map_BB[PID]->KEFromRangeSpline(this_res_length);
      //double dEdx_theory = dEdx_Bethe_Bloch(this_KE, this_mass);
      double dEdx_theory = map_BB[PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(save_graph){
      // == Save vectors for graphes
      chi2_vector.push_back(this_chi2);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  // == Save graph
  if(save_graph){
    // == Vectors for graphes
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = N_skip; i < this_N_hits - N_skip; i++){
      double this_range_original = ResRange.at(i);
      range_original.push_back(this_range_original);

      double this_range_bestfit = ResRange.at(i) + best_additional_res_length;
      range_bestfit.push_back(this_range_bestfit);
      range_reco.push_back(ResRange.at(i));
      dEdx_ordered.push_back(dEdx.at(i));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits - 1));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &chi2_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    chi2_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
  }

  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  if(this_is_beam) original_res_length = ResRange.at(0); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    best_chi2 = 99999.;
    return -9999.;
  }
  else if(best_chi2 > 99990.){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    best_chi2 = 99999.;
    return -9999.;
  }
  else if(best_chi2 < 1.0e-11){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    best_chi2 = 99999.;
    return -9999.;
  }

  return best_total_res_length;
}

double AnalyzerCore::Fit_HypTrkLength_Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam){
  // == only for charged pions and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }

  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  int N_skip = 3;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow =dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length =max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialization
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 15){
    return -9999.; // == Too small number of hits
  }
  int this_N_hits = this_N_calo; // == Use how many hits
  int i_bestfit = -1;

  // == Fit
  vector<double> m2lnL_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_BB[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      
      // == Likelihood
      double this_pitch = fabs(ResRange.at(j - 1) - ResRange.at(j + 1)) / 2.0;
      double this_likelihood = map_BB[PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      this_m2lnL += (-2.0) * log(this_likelihood);
    }
    //this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(save_graph){
      // == Save vectors for graphes
      m2lnL_vector.push_back(this_m2lnL);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  // == Save graph
  if(save_graph){
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = N_skip; i < this_N_hits - N_skip; i++){
      double this_range_original = ResRange.at(i);
      range_original.push_back(this_range_original);

      double this_range_bestfit = ResRange.at(i) + best_additional_res_length;
      range_bestfit.push_back(this_range_bestfit);
      range_reco.push_back(ResRange.at(i));
      dEdx_ordered.push_back(dEdx.at(i));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits - 1));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &m2lnL_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
  }

  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  if(this_is_beam) original_res_length = ResRange.at(0); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    return -9999.;
  }
  else if(best_m2lnL > 99990.){
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    return -9999.;
  }
  else if(best_m2lnL < 1.0e-11){
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    return -9999.;
  }

  // == Return
  m2lnL_vector.clear();
  additional_res_legnth_vector.clear();
  return best_total_res_length;
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
