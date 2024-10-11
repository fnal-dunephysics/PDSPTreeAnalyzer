#include "PionSkimTreeMaker.h"

void PionSkimTreeMaker::initializeAnalyzer(){
  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

  outfile->cd();
  fEventTree = new TTree("PionBeamSkim", "PionBeamSkim");
  fEventTree->Branch("run", &_run);
  fEventTree->Branch("subrun", &_subrun);
  fEventTree->Branch("event", &_event);
  fEventTree->Branch("beam_reco_KE_ff", &_beam_reco_KE_ff);
  fEventTree->Branch("beam_reco_KE_end", &_beam_reco_KE_end);
  fEventTree->Branch("beam_reco_dir_X", &_beam_reco_dir_X);
  fEventTree->Branch("beam_reco_dir_Y", &_beam_reco_dir_Y);
  fEventTree->Branch("beam_reco_dir_Z", &_beam_reco_dir_Z);
  fEventTree->Branch("beam_true_KE_ff", &_beam_true_KE_ff);
  fEventTree->Branch("beam_true_KE_end", &_beam_true_KE_end);
  fEventTree->Branch("beam_true_dir_X", &_beam_true_dir_X);
  fEventTree->Branch("beam_true_dir_Y", &_beam_true_dir_Y);
  fEventTree->Branch("beam_true_dir_Z", &_beam_true_dir_Z);
  fEventTree->Branch("N_daughter_reco", &_N_daughter_reco);
  fEventTree->Branch("daughter_reco_dir_X", "std::vector<double>", &_daughter_reco_dir_X);
  fEventTree->Branch("daughter_reco_dir_Y", "std::vector<double>", &_daughter_reco_dir_Y);
  fEventTree->Branch("daughter_reco_dir_Z", "std::vector<double>", &_daughter_reco_dir_Z);
  fEventTree->Branch("daughter_reco_KE_range_muon", "std::vector<double>", &_daughter_reco_KE_range_muon);
  fEventTree->Branch("daughter_reco_KE_range_pion", "std::vector<double>", &_daughter_reco_KE_range_pion);
  fEventTree->Branch("daughter_reco_KE_range_proton", "std::vector<double>", &_daughter_reco_KE_range_proton);
  //fEventTree->Branch("daughter_reco_KE_hypfit_pion_gaus", "std::vector<double>", &_daughter_reco_KE_hypfit_pion_gaus);  // Commented out
  fEventTree->Branch("daughter_reco_KE_hypfit_pion_likelihood", "std::vector<double>", &_daughter_reco_KE_hypfit_pion_likelihood);
  fEventTree->Branch("daughter_reco_KE_hypfit_proton_gaus", "std::vector<double>", &_daughter_reco_KE_hypfit_proton_gaus);
  //fEventTree->Branch("daughter_reco_KE_hypfit_proton_likelihood", "std::vector<double>", &_daughter_reco_KE_hypfit_proton_likelihood);  // Commented out
  fEventTree->Branch("daughter_reco_cos_beam", "std::vector<double>", &_daughter_reco_cos_beam);
  fEventTree->Branch("daughter_reco_chi2_muon", "std::vector<double>", &_daughter_reco_chi2_muon);
  fEventTree->Branch("daughter_reco_chi2_pion", "std::vector<double>", &_daughter_reco_chi2_pion);
  fEventTree->Branch("daughter_reco_chi2_proton", "std::vector<double>", &_daughter_reco_chi2_proton);
  fEventTree->Branch("daughter_reco_dist_beam", "std::vector<double>", &_daughter_reco_dist_beam);
  fEventTree->Branch("daughter_reco_trackscore", "std::vector<double>", &_daughter_reco_trackscore);
  fEventTree->Branch("daughter_reco_startZ", "std::vector<double>", &_daughter_reco_startZ);
  fEventTree->Branch("daughter_reco_true_PDG", "std::vector<int>", &_daughter_reco_true_PDG);
  fEventTree->Branch("daughter_reco_true_dir_X", "std::vector<double>", &_daughter_reco_true_dir_X);
  fEventTree->Branch("daughter_reco_true_dir_Y", "std::vector<double>", &_daughter_reco_true_dir_Y);
  fEventTree->Branch("daughter_reco_true_dir_Z", "std::vector<double>", &_daughter_reco_true_dir_Z);
  fEventTree->Branch("daughter_reco_true_ID", "std::vector<int>", &_daughter_reco_true_ID);
  fEventTree->Branch("daughter_reco_true_KE", "std::vector<double>", &_daughter_reco_true_KE);
  fEventTree->Branch("daughter_reco_true_P", "std::vector<double>", &_daughter_reco_true_P);
  fEventTree->Branch("daughter_reco_true_mass", "std::vector<double>", &_daughter_reco_true_mass);
  fEventTree->Branch("daughter_reco_true_purity", "std::vector<double>", &_daughter_reco_true_purity);
  fEventTree->Branch("daughter_reco_true_completeness", "std::vector<double>", &_daughter_reco_true_completeness);
  fEventTree->Branch("daughter_reco_true_cos_true_beam", "std::vector<double>", &_daughter_reco_true_cos_true_beam);
  fEventTree->Branch("daughter_reco_true_endProcess", "std::vector<std::string>", &_daughter_reco_true_endProcess);
  fEventTree->Branch("_N_daughter_all_true", &_N_daughter_all_true);
  fEventTree->Branch("daughter_all_true_PDG", "std::vector<int>", &_daughter_all_true_PDG);
  fEventTree->Branch("daughter_all_true_ID", "std::vector<int>", &_daughter_all_true_ID);
  fEventTree->Branch("daughter_all_true_startX", "std::vector<double>", &_daughter_all_true_startX);
  fEventTree->Branch("daughter_all_true_startY", "std::vector<double>", &_daughter_all_true_startY);
  fEventTree->Branch("daughter_all_true_startZ", "std::vector<double>", &_daughter_all_true_startZ);
  fEventTree->Branch("daughter_all_true_startPx", "std::vector<double>", &_daughter_all_true_startPx);
  fEventTree->Branch("daughter_all_true_startPy", "std::vector<double>", &_daughter_all_true_startPy);
  fEventTree->Branch("daughter_all_true_startPz", "std::vector<double>", &_daughter_all_true_startPz);
  fEventTree->Branch("daughter_all_true_startP", "std::vector<double>", &_daughter_all_true_startP);
  fEventTree->Branch("daughter_all_true_endX", "std::vector<double>", &_daughter_all_true_endX);
  fEventTree->Branch("daughter_all_true_endY", "std::vector<double>", &_daughter_all_true_endY);
  fEventTree->Branch("daughter_all_true_endZ", "std::vector<double>", &_daughter_all_true_endZ);
  fEventTree->Branch("daughter_all_true_Process", "std::vector<std::string>", &_daughter_all_true_Process);
  fEventTree->Branch("daughter_all_true_endProcess", "std::vector<std::string>", &_daughter_all_true_endProcess);
}

void PionSkimTreeMaker::executeEvent(){

  pi_type = GetPiParType();
  pi_type_str = Form("%d", pi_type);

  ////////////////////////////////////////
  // == List of selections for pion beam
  ////////////////////////////////////////
  
  // -- 1. Beam instruments
  if(!PassBeamScraperCut()) return;
  if(!PassBeamMomentumWindowCut()) return;

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  KE_beam_inst = map_BB[211] -> MomentumtoKE(P_beam_inst);
  exp_trk_len_beam_inst = map_BB[211] -> RangeFromKESpline(KE_beam_inst);
  trk_len_ratio = evt.reco_beam_alt_len / exp_trk_len_beam_inst;

  mass_beam = 139.57;
  P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;
  
  double P_reweight = 1.;
  if(!IsData){
    KE_ff_true = Get_true_ffKE();
    P_ff_true = KE_to_P(KE_ff_true, evt.true_beam_PDG);
  }

  if(!Pass_Beam_PID(211)) return;
  
  // -- 2. TPC info
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;
  if(!IsData) P_reweight = MCCorr -> MomentumReweight_SF("TrkLength", P_beam_inst, 0.);
  if(!PassAPA3Cut()) return;
  if(!PassMichelScoreCut()) return;
  if(!Pass_beam_start_Z_cut(2.0)) return;
  if(!Pass_beam_delta_X_cut(2.0)) return;
  if(!Pass_beam_delta_Y_cut(2.0)) return;
  //if(chi2_proton < 170 || chi2_proton > 270) return;

  FillSkimTree();

  /*
  if(!PassBeamCosCut()) return;
  if(!PassBeamStartZCut()) return;
  if(!PassProtonVetoCut()) return;
  if(!PassMuonVetoCut()) return;
  if(!PassStoppedPionVetoCut()) return;
  */
  //cout << "[PionSkimTreeMaker::executeEvent] Passed all beam selections" << endl;

  //Study_with_daughters(P_reweight);
}

void PionSkimTreeMaker::FillSkimTree(){

  // == Initialize variables
  _run = -999;
  _subrun = -999;
  _event = -999;
  _beam_reco_KE_ff = -999.0;
  _beam_reco_KE_end = -999.0;
  _beam_reco_dir_X = -999.0;
  _beam_reco_dir_Y = -999.0;
  _beam_reco_dir_Z = -999.0;
  _beam_true_KE_ff = -999.0;
  _beam_true_KE_end = -999.0;
  _beam_true_dir_X = -999.0;
  _beam_true_dir_Y = -999.0;
  _beam_true_dir_Z = -999.0;
  _N_daughter_reco = -999;
  _N_daughter_all_true = -999;
  _daughter_reco_dir_X.clear();
  _daughter_reco_dir_Y.clear();
  _daughter_reco_dir_Z.clear();
  _daughter_reco_KE_range_muon.clear();
  _daughter_reco_KE_range_pion.clear();
  _daughter_reco_KE_range_proton.clear();
  //_daughter_reco_KE_hypfit_pion_gaus.clear();
  _daughter_reco_KE_hypfit_pion_likelihood.clear();
  _daughter_reco_KE_hypfit_proton_gaus.clear();
  //_daughter_reco_KE_hypfit_proton_likelihood.clear();
  _daughter_reco_cos_beam.clear();
  _daughter_reco_chi2_muon.clear();
  _daughter_reco_chi2_pion.clear();
  _daughter_reco_chi2_proton.clear();
  _daughter_reco_dist_beam.clear();
  _daughter_reco_trackscore.clear();
  _daughter_reco_startZ.clear();
  _daughter_reco_true_PDG.clear();
  _daughter_reco_true_ID.clear();
  _daughter_reco_true_dir_X.clear();
  _daughter_reco_true_dir_Y.clear();
  _daughter_reco_true_dir_Z.clear();
  _daughter_reco_true_KE.clear();
  _daughter_reco_true_P.clear();
  _daughter_reco_true_mass.clear();
  _daughter_reco_true_purity.clear();
  _daughter_reco_true_completeness.clear();
  _daughter_reco_true_cos_true_beam.clear();
  _daughter_reco_true_endProcess.clear();
  _daughter_all_true_PDG.clear();
  _daughter_all_true_ID.clear();
  _daughter_all_true_startX.clear();
  _daughter_all_true_startY.clear();
  _daughter_all_true_startZ.clear();
  _daughter_all_true_startPx.clear();
  _daughter_all_true_startPy.clear();
  _daughter_all_true_startPz.clear();
  _daughter_all_true_startP.clear();
  _daughter_all_true_endX.clear();
  _daughter_all_true_endY.clear();
  _daughter_all_true_endZ.clear();
  _daughter_all_true_Process.clear();
  _daughter_all_true_endProcess.clear();

  
  double P_beam_true =  -999.;
  double KE_beam_true = -999.;
  vector<Daughter> daughter_all = GetAllDaughters();

  TVector3 p_vec_beam(evt.true_beam_endPx, evt.true_beam_endPy, evt.true_beam_endPz);
  p_vec_beam = 1000. * p_vec_beam;

  TVector3 pt0(evt.reco_beam_calo_startX,
               evt.reco_beam_calo_startY,
               evt.reco_beam_calo_startZ);
  TVector3 pt1(evt.reco_beam_calo_endX,
               evt.reco_beam_calo_endY,
               evt.reco_beam_calo_endZ);
  TVector3 dir = pt1 - pt0;
  TVector3 reco_unit_beam = (1. / dir.Mag() ) * dir;

  // == Set reco beam info
  _run = evt.run;
  _subrun = evt.subrun;
  _event = evt.event;
  _beam_reco_KE_ff = KE_ff_reco;
  _beam_reco_KE_end = KE_end_reco;
  _beam_reco_dir_X = reco_unit_beam.X();
  _beam_reco_dir_Y = reco_unit_beam.Y();
  _beam_reco_dir_Z = reco_unit_beam.Z();
  _N_daughter_reco = daughter_all.size();
  
  for(size_t i = 0; i < daughter_all.size(); i++){
    // == Set reco daughter info
    Daughter this_daughter = daughter_all.at(i);
    double this_KE_range_muon =	map_BB[13] -> KEFromRangeSpline(this_daughter.allTrack_resRange_SCE().back());
    double this_KE_range_pion = map_BB[211] -> KEFromRangeSpline(this_daughter.allTrack_resRange_SCE().back());
    double this_KE_range_proton = map_BB[2212] -> KEFromRangeSpline(this_daughter.allTrack_resRange_SCE().back());
    double this_hypfit_pion_likelihood = KE_Hypfit_Likelihood(this_daughter, 211);
    double this_hypfit_proton_gaus = KE_Hypfit_Gaussian(this_daughter, 2212);
    double this_reco_cos_beam = this_daughter.Beam_Cos();
    double this_chi2_muon = this_daughter.allTrack_Chi2_muon() / this_daughter.allTrack_Chi2_ndof();
    double this_chi2_pion = this_daughter.allTrack_Chi2_pion() / this_daughter.allTrack_Chi2_ndof();
    double this_chi2_proton = this_daughter.allTrack_Chi2_proton() / this_daughter.allTrack_Chi2_ndof();
    double this_dist_beam = this_daughter.Beam_Dist();
    double this_trackscore = this_daughter.PFP_trackScore();
    double this_startZ = this_daughter.allTrack_startZ();

    TVector3 this_unit_daughter(this_daughter.allTrack_endX() - this_daughter.allTrack_startX(),
			   this_daughter.allTrack_endY() - this_daughter.allTrack_startY(),
			   this_daughter.allTrack_endZ() - this_daughter.allTrack_startZ());
    this_unit_daughter = (1./ this_unit_daughter.Mag() ) * this_unit_daughter;
    _daughter_reco_dir_X.push_back(this_unit_daughter.X());
    _daughter_reco_dir_Y.push_back(this_unit_daughter.Y());
    _daughter_reco_dir_Z.push_back(this_unit_daughter.Z());
    _daughter_reco_KE_range_muon.push_back(this_KE_range_muon);
    _daughter_reco_KE_range_pion.push_back(this_KE_range_pion);
    _daughter_reco_KE_range_proton.push_back(this_KE_range_proton);
    _daughter_reco_KE_hypfit_pion_likelihood.push_back(this_hypfit_pion_likelihood);
    _daughter_reco_KE_hypfit_proton_gaus.push_back(this_hypfit_proton_gaus);
    _daughter_reco_cos_beam.push_back(this_reco_cos_beam);
    _daughter_reco_chi2_muon.push_back(this_chi2_muon);
    _daughter_reco_chi2_pion.push_back(this_chi2_pion);
    _daughter_reco_chi2_proton.push_back(this_chi2_proton);
    _daughter_reco_dist_beam.push_back(this_dist_beam);
    _daughter_reco_trackscore.push_back(this_trackscore);
    _daughter_reco_startZ.push_back(this_startZ);
    
    // == Set true daughter info
    if(!IsData){
      double this_true_PDG = this_daughter.PFP_true_byHits_PDG(); 
      double this_true_ID = this_daughter.PFP_true_byHits_ID();

      double this_true_P = this_daughter.PFP_true_byHits_startP() * 1000.;
      double this_true_E = this_daughter.PFP_true_byHits_startE() * 1000.;
      double this_true_mass = pow(pow(this_true_E, 2.) - pow(this_true_P, 2.), 0.5);
      double this_true_KE = this_true_E - this_true_mass;

      double this_true_purity = this_daughter.PFP_true_byHits_purity();
      double this_true_completeness = this_daughter.PFP_true_byHits_completeness();

      TVector3 this_p_vec(this_daughter.PFP_true_byHits_startPx(), this_daughter.PFP_true_byHits_startPy(), this_daughter.PFP_true_byHits_startPz());
      this_p_vec = 1000. * this_p_vec;
      double this_true_cos_true_beam = p_vec_beam.Dot(this_p_vec) / (p_vec_beam.Mag() * this_p_vec.Mag());

      string this_true_endProcess = this_daughter.PFP_true_byHits_endProcess();

      _daughter_reco_true_PDG.push_back(this_true_PDG);
      _daughter_reco_true_ID.push_back(this_true_ID);
      _daughter_reco_true_dir_X.push_back(this_p_vec.X() / this_p_vec.Mag());
      _daughter_reco_true_dir_Y.push_back(this_p_vec.Y() / this_p_vec.Mag());
      _daughter_reco_true_dir_Z.push_back(this_p_vec.Z() / this_p_vec.Mag());
      _daughter_reco_true_KE.push_back(this_true_KE);
      _daughter_reco_true_P.push_back(this_true_P);
      _daughter_reco_true_mass.push_back(this_true_mass);
      _daughter_reco_true_purity.push_back(this_true_purity);
      _daughter_reco_true_completeness.push_back(this_true_completeness);
      _daughter_reco_true_cos_true_beam.push_back(this_true_cos_true_beam);
      _daughter_reco_true_endProcess.push_back(this_true_endProcess);
    }
  }
  
  if(!IsData){
    double true_P_end = evt.true_beam_endP;
    true_P_end = true_P_end * 1000.;
    double true_KE_end = -999.;
    
    if(evt.true_beam_PDG == 211){
      true_KE_end = sqrt(pow(true_P_end, 2) + pow(M_pion, 2)) - M_pion;
    }
    else if(evt.true_beam_PDG == 2212){
      true_KE_end = sqrt(pow(true_P_end, 2) + pow(M_proton, 2)) - M_proton;
    }

    // == Set true beam info
    _beam_true_KE_ff = KE_ff_true;
    _beam_true_KE_end = true_KE_end;
    _beam_true_dir_X = p_vec_beam.X() / p_vec_beam.Mag();
    _beam_true_dir_Y = p_vec_beam.Y() / p_vec_beam.Mag();
    _beam_true_dir_Z = p_vec_beam.Z() / p_vec_beam.Mag();

    vector<TrueDaughter> true_daughter_all = GetAllTrueDaughters();
    _N_daughter_all_true = true_daughter_all.size();
    
    for(size_t i = 0; i < true_daughter_all.size(); i++){
      TrueDaughter this_daughter = true_daughter_all.at(i);
       _daughter_all_true_PDG.push_back(this_daughter.PDG());
      _daughter_all_true_ID.push_back(this_daughter.ID());
      _daughter_all_true_startX.push_back(this_daughter.startX());
      _daughter_all_true_startY.push_back(this_daughter.startY());
      _daughter_all_true_startZ.push_back(this_daughter.startZ());
      _daughter_all_true_startPx.push_back(this_daughter.startPx() * 1000);
      _daughter_all_true_startPy.push_back(this_daughter.startPy() * 1000);
      _daughter_all_true_startPz.push_back(this_daughter.startPz() * 1000);
      _daughter_all_true_startP.push_back(this_daughter.startP() * 1000);
      _daughter_all_true_endX.push_back(this_daughter.endX());
      _daughter_all_true_endY.push_back(this_daughter.endY());
      _daughter_all_true_endZ.push_back(this_daughter.endZ());
      _daughter_all_true_Process.push_back(this_daughter.Process());
      _daughter_all_true_endProcess.push_back(this_daughter.endProcess());
    }

  }  

  fEventTree->Fill();
}

double PionSkimTreeMaker::KE_to_P(double KE, int PID){

  if(KE < 0.) return -999.;
  double mass = -999.;
  if(abs(PID) == 13) mass = M_mu;
  else if(abs(PID) == 211) mass = M_pion;
  else if(abs(PID) == 2212) mass = M_proton;
  else return -999.;

  double this_P = sqrt(KE*KE + 2.0 * mass * KE);
  return this_P;
}

bool PionSkimTreeMaker::Pass_beam_start_Z_cut(double N_sigma){

  bool out = false;
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_data) / Beam_startZ_sigma_data;
  if(!IsData) Beam_startZ_over_sigma= (evt.reco_beam_calo_startZ - Beam_startZ_mu_mc) / Beam_startZ_sigma_mc;

  if(fabs(Beam_startZ_over_sigma) < N_sigma) out = true;

  return out;
}

bool PionSkimTreeMaker::Pass_beam_delta_X_cut(double N_sigma){

  bool out = false;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  if(!IsData) Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;

  if(fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

  return out;
}

bool PionSkimTreeMaker::Pass_beam_delta_Y_cut(double N_sigma){

  bool out = false;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData) Beam_delta_Y_over_sigma =(delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;

  if(fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

  return out;
}

bool PionSkimTreeMaker::Pass_beam_TPC_theta_cut(double N_sigma){
  bool out = false;
  double Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_data) / Beam_TPC_theta_sigma_data;
  if(!IsData) Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_mc) / Beam_TPC_theta_sigma_mc;

  if(fabs(Beam_TPC_theta_over_sigma) < N_sigma) out = true;

  return out;
}

bool PionSkimTreeMaker::Pass_beam_TPC_phi_cut(double N_sigma){
  bool out = false;
  double Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_data) / Beam_TPC_phi_sigma_data;
  if(!IsData) Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_mc) / Beam_TPC_phi_sigma_mc;

  if(fabs(Beam_TPC_phi_over_sigma) < N_sigma) out = true;

  return out;
}

PionSkimTreeMaker::PionSkimTreeMaker(){

}

PionSkimTreeMaker::~PionSkimTreeMaker(){

}
