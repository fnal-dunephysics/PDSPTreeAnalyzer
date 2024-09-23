#include "Pion2Proton.h"

void Pion2Proton::initializeAnalyzer(){
  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;
}

void Pion2Proton::executeEvent(){

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
    Beam_true_Eloss(); // == Compare E-loss due to beam plug between proton and pion beam particles
  }

  //if(!IsData) True_KI_study("");
  
  if(!Pass_Beam_PID(211)) return;
  
  // -- 2. TPC info
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;
  if(!IsData) P_reweight = MCCorr -> MomentumReweight_SF("TrkLength", P_beam_inst, 0.);
  FillBeamPlots("Beam_CaloSize", P_reweight);

  if(!PassAPA3Cut()) return;
  FillBeamPlots("Beam_APA3", P_reweight);

  if(!PassMichelScoreCut()) return;
  FillBeamPlots("Beam_MichelScor", P_reweight);

  if(!Pass_beam_start_Z_cut(2.0)) return;
  FillBeamPlots("Beam_BeamStartZ", P_reweight);

  if(!Pass_beam_delta_X_cut(2.0)) return;
  FillBeamPlots("Beam_deltaX", P_reweight);

  if(!Pass_beam_delta_Y_cut(2.0)) return;
  FillBeamPlots("Beam_deltaY", P_reweight);
  
  if(chi2_proton < 170 || chi2_proton > 270) return;
  FillBeamPlots("Beam_chi2p", P_reweight);

  Study_with_daughters(P_reweight);
}

void Pion2Proton::True_KI_study(TString prefix){

  double true_ke_ff = Get_true_ffKE();
  if(true_ke_ff < 0.) return;
  vector<TrueDaughter> true_daughter_all = GetAllTrueDaughters();
  vector<TrueDaughter> true_daughter_pion = GetPionTrueDaughters(true_daughter_all);
  vector<TrueDaughter> true_daughter_proton = GetProtonTrueDaughters(true_daughter_all);
  
  if(abs(evt.true_beam_PDG) == 211){
    JSFillHist(prefix + "true_pion", "N_true_pion", true_daughter_pion.size(), 1., 10., -0.5, 9.5);
    JSFillHist(prefix + "true_pion", "N_true_proton", true_daughter_proton.size(), 1., 10., -0.5, 9.5);
    JSFillHist(prefix + "true_pion", "N_true_pion_vs_N_true_proton", true_daughter_pion.size(), true_daughter_proton.size(), 1., 5., -0.5, 4.5, 5., -0.5, 4.5);

    if(true_daughter_pion.size() == 1 && true_daughter_proton.size() > 0){
    
      TrueDaughter this_pion = true_daughter_pion.at(0);
      TrueDaughter this_proton = true_daughter_proton.at(0);
    
      TVector3 p_vec_beam(evt.true_beam_endPx, evt.true_beam_endPy, evt.true_beam_endPz);
      TVector3 p_vec_pion(this_pion.startPx(), this_pion.startPy(), this_pion.startPz());
      TVector3 p_vec_proton(this_proton.startPx(), this_proton.startPy(), this_proton.startPz());
      p_vec_beam = 1000. * p_vec_beam;
      p_vec_pion = 1000. * p_vec_pion;
      p_vec_proton = 1000. * p_vec_proton;
      double E_proton = pow( pow(p_vec_proton.Mag(), 2.) + pow(M_proton, 2.), 0.5);
      double KE_proton = E_proton - M_proton;
      JSFillHist(prefix + "true_pion_tki", "KE_proton", KE_proton, 1., 2000., 0., 2000.);

      double cos_theta_beam_sec_pion = p_vec_beam.Dot(p_vec_pion) / (p_vec_beam.Mag() * p_vec_pion.Mag());
      double this_EQE_pion = Get_EQE_NC_Pion(p_vec_pion.Mag(), cos_theta_beam_sec_pion, 4., -1.);
      JSFillHist(prefix + "true_pion_tki", "EQEm_pion",  this_EQE_pion, 1., 2000., 0., 2000.);

      double this_deltaEQE_pion = sqrt(pow(p_vec_beam.Mag(), 2) + pow(M_pion,2)) - this_EQE_pion;
      JSFillHist(prefix + "true_pion_tki", "deltaEQEm_pion",  this_deltaEQE_pion, 1., 2000., -1000., 1000.);

      double cos_theta_beam_sec_proton = p_vec_beam.Dot(p_vec_proton) / (p_vec_beam.Mag() * p_vec_proton.Mag());
      double this_EQEm_proton = Get_EQE_NC_Proton(p_vec_proton.Mag(), cos_theta_beam_sec_proton, 4., -1);
      double this_deltaEQEm_proton = sqrt(pow(p_vec_beam.Mag(), 2) + pow(M_pion,2)) - this_EQEm_proton;
      double this_EQEp_proton = Get_EQE_NC_Proton(p_vec_proton.Mag(), cos_theta_beam_sec_proton, 4., 1);
      double this_deltaEQEp_proton = sqrt(pow(p_vec_beam.Mag(), 2) + pow(M_pion,2)) - this_EQEp_proton;
      JSFillHist(prefix + "true_pion_tki", "EQEm_proton",  this_EQEm_proton, 1., 2000., 0., 2000.);
      JSFillHist(prefix + "true_pion_tki", "deltaEQEm_proton",  this_deltaEQEm_proton, 1., 2000., -1000., 1000.);
      
      TString n_proton_str = Form("%zu_protons", true_daughter_proton.size());
      FillTKIVar(prefix + "true_pion_tki", p_vec_beam, p_vec_pion, p_vec_proton, 1.);
      FillTKIVar(prefix + "true_pion_tki_" + n_proton_str, p_vec_beam, p_vec_pion, p_vec_proton, 1.);
      if(fabs(this_deltaEQE_pion) < 50.){
	FillTKIVar(prefix + "true_pion_tki_deltaEQEcut", p_vec_beam, p_vec_pion, p_vec_proton, 1.);
      }

      if(KE_proton > 40.){
	JSFillHist(prefix + "true_pion_tki_protonKE40cut", "EQEm_proton",  this_EQEm_proton, 1., 2000., 0., 2000.);
	JSFillHist(prefix + "true_pion_tki_protonKE40cut", "deltaEQEm_proton",  this_deltaEQEm_proton, 1., 2000., -1000., 1000.);
	FillTKIVar(prefix + "true_pion_tki_protonKE40cut", p_vec_beam, p_vec_pion, p_vec_proton, 1.);
      }
      if(KE_proton > 100.){
	JSFillHist(prefix + "true_pion_tki_protonKE100cut", "EQEm_proton",  this_EQEm_proton, 1., 2000., 0., 2000.);
        JSFillHist(prefix + "true_pion_tki_protonKE100cut", "deltaEQEm_proton",  this_deltaEQEm_proton, 1., 2000., -1000., 1000.);
        FillTKIVar(prefix + "true_pion_tki_protonKE100cut", p_vec_beam, p_vec_pion, p_vec_proton, 1.);
      }
      
      //if(true_daughter_proton.size() == 1) cout << Form("this_EQE (pion, proton+, proton-) : (%f, %f, %f), deltaEQE : (%f, %f, %f)", this_EQE_pion, this_EQEp_proton, this_EQEm_proton, this_deltaEQE_pion, this_deltaEQEp_proton, this_deltaEQEm_proton) << endl;
    }
  }
}

void Pion2Proton::Study_with_daughters(double weight){
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> pions = GetPions(daughters_all);
  vector<Daughter> stopping_pions = GetStoppingPions(daughters_all, 6.);
  vector<Daughter> stopping_protons = GetProtons(daughters_all);

  JSFillHist("count", "N_2nd_stopping_pions", stopping_pions.size(), weight, 10, -0.5, 9.5);
  JSFillHist("count", "N_2nd_stopping_protons", stopping_protons.size(), weight, 10, -0.5, 9.5);

  if(stopping_protons.size() > 0){
    Daughter leading_p = stopping_protons.at(0);
    double max_len = stopping_protons.at(0).allTrack_alt_len();
    if(stopping_protons.size() > 1){
      for(unsigned int i = 1; i < stopping_protons.size(); i++){
	double this_len = stopping_protons.at(i).allTrack_alt_len();
	if(this_len > max_len){
	  max_len = this_len;
	  leading_p = stopping_protons.at(i);
	}
      }
    }

    Study_leading_proton(leading_p, weight);
  }

  if(stopping_protons.size() == 1){
    //True_KI_study("Reco_1p_");
  }
}

void Pion2Proton::Study_leading_proton(Daughter leading_p, double weight){

  FillLeadingProtonPlot("nocut", leading_p, weight);
  
}

void Pion2Proton::FillLeadingProtonPlot(TString suffix, Daughter leading_p, double weight){
  TString this_dir = "leading_p_" + suffix;
  double proton_KE = map_BB[2212] -> KEFromRangeSpline(leading_p.allTrack_resRange_SCE().back());

  JSFillHist(this_dir, "beam_KE_ff_" + this_dir + "_" + pi_type_str, KE_ff_reco, weight, 2000., 0., 2000.);
  JSFillHist(this_dir, "beam_KE_end_" + this_dir + "_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);
  JSFillHist(this_dir, "leading_p_KE_" + this_dir + "_" + pi_type_str, proton_KE, weight, 2000., 0., 2000.);
  JSFillHist(this_dir, "cosine_beam_p_" + this_dir + "_" + pi_type_str, leading_p.Beam_Cos(), weight, 2000., -1., 1.);
  JSFillHist(this_dir, "beam_KE_end_vs_proton_KE_" + this_dir + "_" + pi_type_str, KE_end_reco, proton_KE, weight, 1500., 0., 1500., 1500., 0., 1500.);

  if(KE_end_reco > 650. && KE_end_reco < 750.){
    JSFillHist(this_dir, "cosine_beam_p_vs_proton_KE_at_KE_end_reco_650_to_750_" + this_dir + "_" + pi_type_str, leading_p.Beam_Cos(), proton_KE, weight, 200., -1., 1., 150., 0., 1500.);
  }
  else if(KE_end_reco > 750. && KE_end_reco < 850.){
    JSFillHist(this_dir, "cosine_beam_p_vs_proton_KE_at_KE_end_reco_750_to_850_" + this_dir + "_" + pi_type_str, leading_p.Beam_Cos(), proton_KE, weight, 200., -1., 1., 150., 0., 1500.);
  }
}

void Pion2Proton::Beam_true_Eloss(){

  int true_beam_pdg = evt.true_beam_PDG;
  TString true_beam_pdg_str = "";
  double P_beam_true = evt.true_beam_startP * 1000.;
  double KE_beam_true = -999.;
  if(true_beam_pdg == 211){
    true_beam_pdg_str = "pion";
    KE_beam_true = sqrt(pow(P_beam_true, 2) + pow(M_pion, 2)) - M_pion;
  }
  else if(true_beam_pdg == 2212){
    true_beam_pdg_str = "proton";
    KE_beam_true = sqrt(pow(P_beam_true, 2) + pow(M_proton, 2)) - M_proton;
  }
  else return;

  if(KE_ff_true > 0.){
    JSFillHist("true_beam", true_beam_pdg_str + "_delta_KE_true_ff", KE_beam_true - KE_ff_true, 1., 500., 0., 50.);
    JSFillHist("true_beam", true_beam_pdg_str + "_KE_beam_true_vs_delta_KE_true_ff", KE_beam_true, KE_beam_true - KE_ff_true, 1., 2000., 0., 2000., 500., 0., 50.);
  }
}

void Pion2Proton::FillBeamPlots(TString beam_selec_str, double weight){

  // == Fit results after calo-size cut
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_data) / Beam_startZ_sigma_data;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  //double Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_data) / Beam_TPC_theta_sigma_data;
  //double Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_data) / Beam_TPC_phi_sigma_data;
  if(!IsData){
    Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_mc) / Beam_startZ_sigma_mc;
    Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;
    Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;
    //Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_mc) / Beam_TPC_theta_sigma_mc;
    //Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_mc) / Beam_TPC_phi_sigma_mc;
  }

  // == Comparison between beam spectrometer track and TPC reco track
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX", evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY", evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ", evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_over_sigma", Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ", evt.reco_beam_calo_endZ, weight, 1000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst", P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh", beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton", chi2_proton, weight, 10000., 0., 1000.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX_" + pi_type_str, evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY_" + pi_type_str, evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_" + pi_type_str, evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_over_sigma_" + pi_type_str, Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_" + pi_type_str, delta_X_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_" + pi_type_str, delta_Y_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_X_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_Y_over_sigma, weight, 2000., -10., 10.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_cos_delta_spec_TPC_" + pi_type_str, cos_delta_spec_TPC, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh_" + pi_type_str, beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_theta_" + pi_type_str, beam_TPC_theta, weight, 5000., -1., 4.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_phi_" + pi_type_str, beam_TPC_phi, weight, 8000., -4., 4.);
  //JSFillHist(beam_selec_str, "Beam_TPC_theta_over_sigma_" + pi_type_str,Beam_TPC_theta_over_sigma, weight, 2000., -10., 10.);
  //JSFillHist(beam_selec_str, "Beam_TPC_phi_over_sigma_" + pi_type_str,Beam_TPC_phi_over_sigma, weight, 2000., -10., 10.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ_" + pi_type_str, evt.reco_beam_calo_endZ, weight, 1000., 0., 1000.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_alt_len_" + pi_type_str, evt.reco_beam_alt_len, weight, 1000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst_" + pi_type_str, P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton_" + pi_type_str, chi2_proton, weight, 10000., 0., 1000.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_trk_len_ratio_" + pi_type_str, trk_len_ratio, weight, 1000., 0., 10.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_ff_" + pi_type_str, KE_ff_reco, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_end_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);
}

void Pion2Proton::drawKIVar(Daughter pion, Daughter proton, double weight){

  TVector3 p_vec_beam(evt.reco_beam_trackDirX, evt.reco_beam_trackDirY, evt.reco_beam_trackDirZ);
  p_vec_beam = p_vec_beam * P_ff_reco;

  TVector3 p_vec_pion(pion.allTrack_endX() - pion.allTrack_startX(), pion.allTrack_endY() - pion.allTrack_startY(), pion.allTrack_endZ() - pion.allTrack_startZ()); 
  p_vec_pion = p_vec_pion.Unit();
  double pion_KE_BB = map_BB[211] -> KEFromRangeSpline(pion.allTrack_resRange_SCE().back());
  double pion_p = map_BB[211] -> KEtoMomentum(pion_KE_BB);
  p_vec_pion = p_vec_pion * pion_p;

  TVector3 p_vec_proton(proton.allTrack_endX() - proton.allTrack_startX(), proton.allTrack_endY() - proton.allTrack_startY(), proton.allTrack_endZ() - proton.allTrack_startZ());
  p_vec_proton = p_vec_proton.Unit();
  double proton_KE_BB = map_BB[2212] -> KEFromRangeSpline(proton.allTrack_resRange_SCE().back());
  double proton_p = map_BB[2212] -> KEtoMomentum(proton_KE_BB);
  p_vec_proton = p_vec_proton * proton_p;

  FillTKIVar("pion_tki", p_vec_beam, p_vec_pion, p_vec_proton, 1.);
}

double Pion2Proton::KE_to_P(double KE, int PID){

  if(KE < 0.) return -999.;
  double mass = -999.;
  if(abs(PID) == 13) mass = M_mu;
  else if(abs(PID) == 211) mass = M_pion;
  else if(abs(PID) == 2212) mass = M_proton;
  else return -999.;

  double this_P = sqrt(KE*KE + 2.0 * mass * KE);
  return this_P;
}

bool Pion2Proton::Pass_beam_start_Z_cut(double N_sigma){

  bool out = false;
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_data) / Beam_startZ_sigma_data;
  if(!IsData) Beam_startZ_over_sigma= (evt.reco_beam_calo_startZ - Beam_startZ_mu_mc) / Beam_startZ_sigma_mc;

  if(fabs(Beam_startZ_over_sigma) < N_sigma) out = true;

  return out;
}

bool Pion2Proton::Pass_beam_delta_X_cut(double N_sigma){

  bool out = false;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  if(!IsData) Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;

  if(fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

  return out;
}

bool Pion2Proton::Pass_beam_delta_Y_cut(double N_sigma){

  bool out = false;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData) Beam_delta_Y_over_sigma =(delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;

  if(fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

  return out;
}

bool Pion2Proton::Pass_beam_TPC_theta_cut(double N_sigma){
  bool out = false;
  double Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_data) / Beam_TPC_theta_sigma_data;
  if(!IsData) Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_mc) / Beam_TPC_theta_sigma_mc;

  if(fabs(Beam_TPC_theta_over_sigma) < N_sigma) out = true;

  return out;
}

bool Pion2Proton::Pass_beam_TPC_phi_cut(double N_sigma){
  bool out = false;
  double Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_data) / Beam_TPC_phi_sigma_data;
  if(!IsData) Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_mc) / Beam_TPC_phi_sigma_mc;

  if(fabs(Beam_TPC_phi_over_sigma) < N_sigma) out = true;

  return out;
}

void Pion2Proton::FillTKIVar(TString dir, TVector3 p_vec_beam, TVector3 p_vec_f1, TVector3 p_vec_f2, double weight){

  double this_beam_phi = p_vec_beam.Phi();

  p_vec_beam.RotateZ(-1. * this_beam_phi);
  double this_beam_theta = p_vec_beam.Theta();
  p_vec_beam.RotateY(-1. * this_beam_theta);
  p_vec_f1.RotateZ(-1. * this_beam_phi);
  p_vec_f1.RotateY(-1. * this_beam_theta);
  p_vec_f2.RotateZ(-1. * this_beam_phi);
  p_vec_f2.RotateY(-1. * this_beam_theta);

  TVector3 p_vec_f1_T(p_vec_f1.X(), p_vec_f1.Y(), 0);
  TVector3 p_vec_f2_T(p_vec_f2.X(), p_vec_f2.Y(), 0);
  TVector3 delta_p_T_vec = p_vec_f1_T + p_vec_f2_T;
  double delta_p_T = delta_p_T_vec.Mag();
  double delta_alpha_T = TMath::ACos((-1. * p_vec_f1_T.Dot(delta_p_T_vec)/(p_vec_f1_T.Mag() * delta_p_T_vec.Mag())));
  delta_alpha_T = 360. * delta_alpha_T / (2. * TMath::Pi());
  double delta_phi_T = TMath::ACos((-1. * p_vec_f1_T.Dot(p_vec_f2_T)/(p_vec_f1_T.Mag() * p_vec_f2_T.Mag())));
  delta_phi_T = 360. * delta_phi_T / (2. * TMath::Pi());

  TVector3 delta_p_vec = p_vec_f1 + p_vec_f2 - p_vec_beam;
  double delta_p = delta_p_vec.Mag();

  JSFillHist(dir, dir + "_beam_p_" + pi_type_str, p_vec_beam.Mag(), weight, 2000., 0., 2000.); 

  JSFillHist(dir, dir + "_delta_p_" + pi_type_str, delta_p, weight, 2000., 0., 2000.);
  JSFillHist(dir, dir + "_delta_p_T_" + pi_type_str, delta_p_T, weight, 2000., 0., 2000.);
  JSFillHist(dir, dir + "_delta_alpha_T_" + pi_type_str, delta_alpha_T, weight, 400., 0., 400.);
  JSFillHist(dir, dir + "_delta_phi_T_" + pi_type_str, delta_phi_T, weight, 400., 0., 400.);

  JSFillHist(dir, dir + "_delta_p_T_vs_delta_alpha_T_" + pi_type_str, delta_p_T, delta_alpha_T, weight, 2000., 0., 2000., 400., 0., 400.);

  if(delta_alpha_T < 45.){
    JSFillHist(dir, dir + "_delta_p_T_with_delta_alpha_T_0_to_45_" + pi_type_str, delta_p_T, weight, 2000., 0., 2000.);
  }
  else if(delta_alpha_T < 90.){
    JSFillHist(dir, dir + "_delta_p_T_with_delta_alpha_T_45_to_90_" + pi_type_str, delta_p_T, weight, 2000., 0., 2000.);
  }
  else if(delta_alpha_T < 135.){
    JSFillHist(dir, dir + "_delta_p_T_with_delta_alpha_T_90_to_135_" + pi_type_str, delta_p_T, weight, 2000., 0., 2000.);
  }
  else{
    JSFillHist(dir, dir + "_delta_p_T_with_delta_alpha_T_135_to_180_" + pi_type_str, delta_p_T, weight, 2000., 0., 2000.);
  }

  // -- mass of outgoing nucleus

  double f1_mass = M_pion;
  if(dir.Contains("proton_tki")) f1_mass = M_proton;
  ROOT::Math::PxPyPzEVector beam_v4, f1_v4, f2_v4, Ar_target_v4;
  beam_v4.SetXYZT(p_vec_beam.X(), p_vec_beam.Y(), p_vec_beam.Z(), pow( pow(f1_mass, 2.) + pow(p_vec_beam.Mag(), 2), 0.5 ));
  f1_v4.SetXYZT(p_vec_f1.X(), p_vec_f1.Y(), p_vec_f1.Z(), pow( pow(f1_mass, 2.) + pow(p_vec_f1.Mag(), 2), 0.5 ));
  f2_v4.SetXYZT(p_vec_f2.X(), p_vec_f2.Y(), p_vec_f2.Z(), pow( pow(M_proton, 2.) + pow(p_vec_f2.Mag(), 2), 0.5 ));
  Ar_target_v4.SetXYZT(0., 0., 0., 37225.0);
  
  ROOT::Math::PxPyPzEVector Xp = beam_v4 + Ar_target_v4 - f1_v4 - f2_v4;
  double mXp = Xp.M();
  mXp = mXp * 0.001;
  if(dir.Contains("true")){
    //cout << dir << ", mXp : " << mXp << ", beam_v4.M() : " << beam_v4.M() << ", f1_v4.M() : " << f1_v4.M() << ", f2_v4.M() : " << f2_v4.M() << endl;
  }
  JSFillHist(dir, dir + "_mXp_" + pi_type_str, mXp, weight, 10000., 30., 40.);

}

Pion2Proton::Pion2Proton(){

}

Pion2Proton::~Pion2Proton(){

}
