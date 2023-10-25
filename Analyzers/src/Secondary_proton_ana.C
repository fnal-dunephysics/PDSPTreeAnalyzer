#include "Secondary_proton_ana.h"

void Secondary_proton_ana::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void Secondary_proton_ana::executeEvent(){

  p_type = GetPParType();
  p_type_str = Form("%d", p_type);

  ////////////////////////////////////
  // === Event selections for beam
  ////////////////////////////////////

  // == Event selections using beam instruments
  if(!PassBeamScraperCut()) return;
  if(!PassBeamMomentumWindowCut()) return;
  if(!IsData) Beam_true_Eloss(); // == Compare E-loss due to beam plug between proton and pion beam particles
  if(!Pass_Beam_PID(2212)) return;

  // == Event selections using TPC info
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  if(!IsData){
    KE_ff_true = Get_true_ffKE();
    P_ff_true = KE_to_P(KE_ff_true, evt.true_beam_PDG);
  }

  FillBeamPlots("Beam_CaloSize");

  //if(!PassAPA3Cut(100.)) return;
  //FillBeamPlots("Beam_APA3");

  //if(!PassBeamCosCut()) return;
  //FillBeamPlots("Beam_BeamCos");

  if(!Pass_beam_start_Z_cut(2.0)) return;
  //if(!PassBeamStartZCut()) return;
  FillBeamPlots("Beam_BeamStartZ");

  //if(chi2_proton < 10.) return;
  //FillBeamPlots("Beam_chi2_proton");

  if(!Pass_beam_delta_X_cut(2.0)) return;
  FillBeamPlots("Beam_deltaX");

  if(!Pass_beam_delta_Y_cut(2.0)) return;
  FillBeamPlots("Beam_deltaY");

  if(!Pass_beam_TPC_theta_cut(2.0)) return;
  FillBeamPlots("Beam_TPCtheta");

  if(!Pass_beam_TPC_phi_cut(2.0)) return;
  FillBeamPlots("Beam_TPCphi");

  if(evt.reco_beam_alt_len < 5.) return;
  FillBeamPlots("Beam_AltLen");


  // == Study relation between P_beam_inst and P_ff
  //Study_Beam_Proton_Eloss("Beam_BeamStartZ");

  mass_beam = M_proton;
  P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;
 
  JSFillHist("reco_beam", "htrack_P_beam", P_beam_inst, 1., 5000., 0., 5000.);

  ////////////////////////////////////
  // == Event selections for daughters
  ////////////////////////////////////

  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  JSFillHist("Beam_BeamStartZ", "Beam_N_nocut_daughters", daughters_all.size(), 1., 10, -0.5, 9.5);
  JSFillHist("Beam_BeamStartZ", "Beam_N_nocut_daughters_" + p_type_str, daughters_all.size(), 1., 10, -0.5, 9.5);

  vector<Daughter> daughters_cut = GetDaughters(daughters_all);
  JSFillHist("Beam_BeamStartZ", "Beam_N_cut_daughters", daughters_cut.size(), 1., 10, -0.5, 9.5);
  JSFillHist("Beam_BeamStartZ", "Beam_N_cut_daughters_" + p_type_str, daughters_cut.size(), 1., 10, -0.5, 9.5);

  vector<Daughter> protons = GetProtons(daughters_all);

  // == MC true daughters
  if(!IsData){
    vector<Daughter> true_protons = GetTrueProtons(daughters_all);
    vector<Daughter> true_pions = GetTruePions(daughters_all);
    True_Daughter_study(true_protons, true_pions);
  }

  // == Sideband : no basic cut daughter
  if(daughters_cut.size() == 0){
    FillBeamPlots("Beam_SB_N_cut_daughters");
  }
  
  if(daughters_cut.size() < 1) return;
  FillBeamPlots("Beam_N_cut_daughters");
  JSFillHist("Beam_N_cut_daughters", "Beam_N_cut_daughters", daughters_cut.size(), 1., 10, -0.5, 9.5);
  JSFillHist("Beam_N_cut_daughters", "Beam_N_cut_daughters_" + p_type_str, daughters_cut.size(), 1., 10, -0.5, 9.5);

  Daughter_study(protons);

  return;
}

void Secondary_proton_ana::FillBeamPlots(TString beam_selec_str){

  // == Fit results after calo-size cut
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_data) / Beam_startZ_sigma_data;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  double Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_data) / Beam_TPC_theta_sigma_data;
  double Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_data) / Beam_TPC_phi_sigma_data;
  if(!IsData){
    Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_mc) / Beam_startZ_sigma_mc;
    Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;
    Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;
    Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_mc) / Beam_TPC_theta_sigma_mc;
    Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_mc) / Beam_TPC_phi_sigma_mc;
  }

  // == Comparison between beam spectrometer track and TPC reco track
  JSFillHist(beam_selec_str, "Beam_startX", evt.reco_beam_calo_startX, 1., 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startY", evt.reco_beam_calo_startY, 1., 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ", evt.reco_beam_calo_startZ, 1., 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ_over_sigma", Beam_startZ_over_sigma, 1., 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_endZ", evt.reco_beam_calo_endZ, 1., 1000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_P_beam_inst", P_beam_inst, 1., 2000., 0., 2000.);
  JSFillHist(beam_selec_str, "Beam_costh", beam_costh, 1., 2000., -1., 1.);
  JSFillHist(beam_selec_str, "Beam_chi2_proton", chi2_proton, 1., 10000., 0., 1000.);

  JSFillHist(beam_selec_str, "Beam_startX_" + p_type_str, evt.reco_beam_calo_startX, 1., 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startY_" + p_type_str, evt.reco_beam_calo_startY, 1., 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ_" + p_type_str, evt.reco_beam_calo_startZ, 1., 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ_over_sigma_" + p_type_str, Beam_startZ_over_sigma, 1., 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_delta_X_spec_TPC_" + p_type_str, delta_X_spec_TPC, 1., 2000., -100., 100.);
  JSFillHist(beam_selec_str, "Beam_delta_Y_spec_TPC_" + p_type_str, delta_Y_spec_TPC, 1., 2000., -100., 100.);
  JSFillHist(beam_selec_str, "Beam_delta_X_spec_TPC_over_sigma_" + p_type_str, Beam_delta_X_over_sigma, 1., 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_delta_Y_spec_TPC_over_sigma_" + p_type_str, Beam_delta_Y_over_sigma, 1., 2000., -10., 10.);

  JSFillHist(beam_selec_str, "Beam_cos_delta_spec_TPC_" + p_type_str, cos_delta_spec_TPC, 1., 2000., -1., 1.);
  JSFillHist(beam_selec_str, "Beam_costh_" + p_type_str, beam_costh, 1., 2000., -1., 1.);
  JSFillHist(beam_selec_str, "Beam_TPC_theta_" + p_type_str, beam_TPC_theta, 1., 5000., -1., 4.);
  JSFillHist(beam_selec_str, "Beam_TPC_phi_" + p_type_str, beam_TPC_phi, 1., 8000., -4., 4.);
  JSFillHist(beam_selec_str, "Beam_TPC_theta_over_sigma_" + p_type_str,Beam_TPC_theta_over_sigma, 1., 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_TPC_phi_over_sigma_" + p_type_str,Beam_TPC_phi_over_sigma, 1., 2000., -10., 10.);

  JSFillHist(beam_selec_str, "Beam_endZ_" + p_type_str, evt.reco_beam_calo_endZ, 1., 1000., 0., 1000.);

  JSFillHist(beam_selec_str, "Beam_alt_len_" + p_type_str, evt.reco_beam_alt_len, 1., 1000., 0., 1000.);
  JSFillHist(beam_selec_str, "Beam_P_beam_inst_" + p_type_str, P_beam_inst, 1., 2000., 0., 2000.);
  JSFillHist(beam_selec_str, "Beam_chi2_proton_" + p_type_str, chi2_proton, 1., 10000., 0., 1000.);

}

void Secondary_proton_ana::Beam_true_Eloss(){

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
    //cout << "[Secondary_proton_ana::Beam_true_Elos] KE_beam_true : " << KE_beam_true << ", KE_ff_true : " << KE_ff_true << endl;
    JSFillHist("true_beam", true_beam_pdg_str + "_delta_KE_true_ff", KE_beam_true - KE_ff_true, 1., 500., 0., 50.);
  }
}

void Secondary_proton_ana::Study_Beam_Proton_Eloss(TString beam_cut_str){

  // == Study using true info
  if(!IsData){
    int true_beam_pdg = evt.true_beam_PDG;
    TString true_beam_pdg_str = "";
    if(true_beam_pdg == 211){
      true_beam_pdg_str = "pion";
    }
    else if(true_beam_pdg == 2212){
      true_beam_pdg_str = "proton";
    }
    if(KE_ff_true > 0. && P_ff_true > 0.){
      double delta_P_reco_beam_true_ff = P_beam_inst - P_ff_true;
      JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_P_ff", P_ff_true, 1., 2000., 0., 2000.);
      JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_P_beam_reco_vs_P_beam_reco_minus_P_ff_true", P_beam_inst, delta_P_reco_beam_true_ff, 1., 2000., 0., 2000., 400., -200., 200.);
    }
  }

  // == Study using reco stopping beam protons
  if(chi2_proton > 10. || (*evt.reco_beam_resRange_SCE).at(0) < 40.) return;
  double this_beam_hyp_length = Fit_HypTrkLength_Gaussian((*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 2212, false, true);
  double this_P_Hypfit = -999.;
  if(this_beam_hyp_length > 0.){
    double this_KE_HypFit = map_BB[2212] -> KEFromRangeSpline(this_beam_hyp_length);
    this_P_Hypfit = KE_to_P(this_KE_HypFit, 2212);
  }
  if(this_P_Hypfit < 0.) return;

  double delta_P_reco_beam_hypfit_ff = P_beam_inst - this_P_Hypfit;
  JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_hypfit", this_P_Hypfit, 1., 2000., 0., 2000.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_reco_vs_P_beam_reco_minus_P_ff_hypfit", P_beam_inst, delta_P_reco_beam_hypfit_ff, 1., 2000., 0., 2000., 400., -200., 200.);

  double beam_res_range = (*evt.reco_beam_resRange_SCE).at(0);
  double this_KE_range = map_BB[2212] -> KEFromRangeSpline(beam_res_range);
  double this_P_range = KE_to_P(this_KE_range, 2212);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_range", this_P_range, 1., 2000., 0., 2000.);

  if(!IsData && P_ff_true > 0.){
    JSFillHist(beam_cut_str, "true_proton_P_ff_chi2_cut", P_ff_true, 1., 2000., 0., 2000.);
    if(p_type == 2){
      JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_hypfit_true_Elas", this_P_Hypfit, 1., 2000., 0., 2000.);
      JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_range_true_Elas", this_P_range, 1., 2000., 0., 2000.);
    } 
  }
}
void Secondary_proton_ana::True_Daughter_study(const vector<Daughter>& protons, const vector<Daughter>& pions){

  int N_protons = protons.size();
  int N_pions = pions.size();
  JSFillHist("Daughter", "N_true_protons_by_Pandora", N_protons, 1., 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_pions_by_Pandora", N_pions, 1., 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_protons", evt.true_daughter_nProton, 1., 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_pions", evt.true_daughter_nPiPlus + evt.true_daughter_nPiMinus, 1., 10., -0.5, 9.5);
}

void Secondary_proton_ana::Daughter_study(const vector<Daughter>& protons){

  int N_protons = protons.size();
  JSFillHist("Daughter", "N_protons", N_protons, 1., 10., -0.5, 9.5);
  if(protons.size() < 1) return;

  for(unsigned int i = 0; i < protons.size(); i++){
    Daughter this_proton = protons.at(i);
    double this_chi2 = this_proton.allTrack_Chi2_proton() / this_proton.allTrack_Chi2_ndof();
    JSFillHist("Daughter", "Proton_chi2", this_chi2, 1., 200., 0., 200.);
  }

}

double Secondary_proton_ana::KE_to_P(double KE, int PID){

  if(KE < 0.) return -999.;
  double mass = -999.;
  if(abs(PID) == 13) mass = M_mu;
  else if(abs(PID) == 211) mass = M_pion;
  else if(abs(PID) == 2212) mass = M_proton;
  else return -999.;

  double this_P = sqrt(KE*KE + 2.0 * mass * KE);
  return this_P;
}

bool Secondary_proton_ana::Pass_beam_start_Z_cut(double N_sigma){

  bool out = false;
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - Beam_startZ_mu_data) / Beam_startZ_sigma_data;
  if(!IsData) Beam_startZ_over_sigma= (evt.reco_beam_calo_startZ - Beam_startZ_mu_mc) / Beam_startZ_sigma_mc;

  if(fabs(Beam_startZ_over_sigma) < N_sigma) out = true;

  return out;
}

bool Secondary_proton_ana::Pass_beam_delta_X_cut(double N_sigma){

  bool out = false;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  if(!IsData) Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;

  if(fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

  return out;
}

bool Secondary_proton_ana::Pass_beam_delta_Y_cut(double N_sigma){

  bool out = false;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData) Beam_delta_Y_over_sigma =(delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;

  if(fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

  return out;
}

bool Secondary_proton_ana::Pass_beam_TPC_theta_cut(double N_sigma){
  bool out = false;
  double Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_data) / Beam_TPC_theta_sigma_data;
  if(!IsData) Beam_TPC_theta_over_sigma = (beam_TPC_theta - Beam_TPC_theta_mu_mc) / Beam_TPC_theta_sigma_mc;

  if(fabs(Beam_TPC_theta_over_sigma) < N_sigma) out = true;

  return out;
}

bool Secondary_proton_ana::Pass_beam_TPC_phi_cut(double N_sigma){
  bool out = false;
  double Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_data) / Beam_TPC_phi_sigma_data;
  if(!IsData) Beam_TPC_phi_over_sigma = (beam_TPC_phi - Beam_TPC_phi_mu_mc) / Beam_TPC_phi_sigma_mc;

  if(fabs(Beam_TPC_phi_over_sigma) < N_sigma) out = true;

  return out;
}

Secondary_proton_ana::Secondary_proton_ana(){

}

Secondary_proton_ana::~Secondary_proton_ana(){

}
