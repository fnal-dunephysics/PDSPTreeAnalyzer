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
  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  KE_beam_inst = map_BB[2212] -> MomentumtoKE(P_beam_inst);
  exp_trk_len_beam_inst = map_BB[2212] -> RangeFromKESpline(KE_beam_inst);
  trk_len_ratio = evt.reco_beam_alt_len / exp_trk_len_beam_inst;

  double P_reweight = 1.;
  if(!IsData){
    KE_ff_true = Get_true_ffKE();
    P_ff_true = KE_to_P(KE_ff_true, evt.true_beam_PDG);
    Beam_true_Eloss(); // == Compare E-loss due to beam plug between proton and pion beam particles
  }

  if(!Pass_Beam_PID(2212)) return;

  if(evt.beam_inst_TOF->size() != 1) return;

  //cout << "evt.beam_inst_PDG_candidates->size() : " <<evt.beam_inst_PDG_candidates->size() << endl;

  //cout << "[Secondary_proton_ana::executeEvent] evt.beam_inst_TOF->size() : " << evt.beam_inst_TOF->size() << endl;
  //cout << "[Secondary_proton_ana::executeEvent] (*evt.beam_inst_TOF).at(0) : " << (*evt.beam_inst_TOF).at(0) << endl;
  beam_TOF = (*evt.beam_inst_TOF).at(0);
  P_beam_TOF = TOF_to_P((*evt.beam_inst_TOF).at(0), 2212);
  //cout << "[Secondary_proton_ana::executeEvent] P_beam_TOF : " << P_beam_TOF << ", P_beam_inst : " << P_beam_inst << endl;

  /*
  double P_TOF_129ns = TOF_to_P(129., 2212);
  double P_TOF_130ns = TOF_to_P(130., 2212);
  double P_TOF_131ns = TOF_to_P(131., 2212);
  cout << Form("(P_TOF_129ns, P_TOF_130ns, P_TOF_131ns) = (%.1f, %.1f, %.1f)", P_TOF_129ns, P_TOF_130ns, P_TOF_131ns) << endl; 
  */

  // == Event selections using TPC info
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;

  if(!IsData) P_reweight = MCCorr -> MomentumReweight_SF("Proton_CaloSize", P_beam_inst, 0.);
  //cout << "[Secondary_proton_ana::executeEvent] P_reweight : " << P_reweight << endl;

  FillBeamPlots("Beam_CaloSize", P_reweight);
  Study_Beam_Proton_Eloss("Beam_CaloSize");

  //if(!PassAPA3Cut(100.)) return;
  //FillBeamPlots("Beam_APA3");

  //if(!PassBeamCosCut()) return;
  //FillBeamPlots("Beam_BeamCos");

  if(!Pass_beam_start_Z_cut(2.0)) return;
  //if(!PassBeamStartZCut()) return;
  FillBeamPlots("Beam_BeamStartZ", P_reweight);

  //if(chi2_proton < 10.) return;
  //FillBeamPlots("Beam_chi2_proton");

  if(!Pass_beam_delta_X_cut(2.0)) return;
  FillBeamPlots("Beam_deltaX", P_reweight);

  if(!Pass_beam_delta_Y_cut(2.0)) return;
  FillBeamPlots("Beam_deltaY", P_reweight);

  if(!Pass_beam_TPC_theta_cut(2.0)) return;
  FillBeamPlots("Beam_TPCtheta", P_reweight);

  if(!Pass_beam_TPC_phi_cut(2.0)) return;
  FillBeamPlots("Beam_TPCphi", P_reweight);

  if(evt.reco_beam_alt_len < 5.) return;
  FillBeamPlots("Beam_AltLen", P_reweight);

  // == Study relation between P_beam_inst and P_ff
  Study_Beam_Proton_Eloss("Beam_AltLen");

  if(trk_len_ratio > 0.7) return;
  FillBeamPlots("Beam_TrkLenRatio", P_reweight);
  
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
  JSFillHist("Beam_BeamStartZ", "Beam_N_nocut_daughters", daughters_all.size(), P_reweight, 10, -0.5, 9.5);
  JSFillHist("Beam_BeamStartZ", "Beam_N_nocut_daughters_" + p_type_str, daughters_all.size(), P_reweight, 10, -0.5, 9.5);

  vector<Daughter> daughters_cut = GetDaughters(daughters_all);
  JSFillHist("Beam_BeamStartZ", "Beam_N_cut_daughters", daughters_cut.size(), P_reweight, 10, -0.5, 9.5);
  JSFillHist("Beam_BeamStartZ", "Beam_N_cut_daughters_" + p_type_str, daughters_cut.size(), P_reweight, 10, -0.5, 9.5);

  vector<Daughter> protons = GetProtons(daughters_all);

  // == MC true daughters
  if(!IsData){
    vector<Daughter> true_protons = GetTrueProtons(daughters_all);
    vector<Daughter> true_pions = GetTruePions(daughters_all);
    True_Daughter_study(true_protons, true_pions, P_reweight);
  }

  // == Sideband : no basic cut daughter
  if(daughters_cut.size() == 0){
    //FillBeamPlots("Beam_SB_N_cut_daughters");
    FillBeamPlots("Beam_SB_N_cut_daughters", P_reweight);
  }
  
  if(daughters_cut.size() < 1) return;
  FillBeamPlots("Beam_N_cut_daughters", P_reweight);
  //FillBeamPlots("Beam_N_cut_daughters");
  JSFillHist("Beam_N_cut_daughters", "Beam_N_cut_daughters", daughters_cut.size(), P_reweight, 10, -0.5, 9.5);
  JSFillHist("Beam_N_cut_daughters", "Beam_N_cut_daughters_" + p_type_str, daughters_cut.size(), P_reweight, 10, -0.5, 9.5);

  Daughter_study(protons);

  return;
}

void Secondary_proton_ana::FillBeamPlots(TString beam_selec_str, double weight){

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
  JSFillHist(beam_selec_str, "Beam_startX", evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startY", evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ", evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ_over_sigma", Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_endZ", evt.reco_beam_calo_endZ, weight, 1000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_P_beam_inst", P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, "Beam_costh", beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, "Beam_chi2_proton", chi2_proton, weight, 10000., 0., 1000.);

  JSFillHist(beam_selec_str, "Beam_startX_" + p_type_str, evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startY_" + p_type_str, evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ_" + p_type_str, evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, "Beam_startZ_over_sigma_" + p_type_str, Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_delta_X_spec_TPC_" + p_type_str, delta_X_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, "Beam_delta_Y_spec_TPC_" + p_type_str, delta_Y_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, "Beam_delta_X_spec_TPC_over_sigma_" + p_type_str, Beam_delta_X_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_delta_Y_spec_TPC_over_sigma_" + p_type_str, Beam_delta_Y_over_sigma, weight, 2000., -10., 10.);

  JSFillHist(beam_selec_str, "Beam_cos_delta_spec_TPC_" + p_type_str, cos_delta_spec_TPC, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, "Beam_costh_" + p_type_str, beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, "Beam_TPC_theta_" + p_type_str, beam_TPC_theta, weight, 5000., -1., 4.);
  JSFillHist(beam_selec_str, "Beam_TPC_phi_" + p_type_str, beam_TPC_phi, weight, 8000., -4., 4.);
  JSFillHist(beam_selec_str, "Beam_TPC_theta_over_sigma_" + p_type_str,Beam_TPC_theta_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, "Beam_TPC_phi_over_sigma_" + p_type_str,Beam_TPC_phi_over_sigma, weight, 2000., -10., 10.);

  JSFillHist(beam_selec_str, "Beam_endZ_" + p_type_str, evt.reco_beam_calo_endZ, weight, 1000., 0., 1000.);

  JSFillHist(beam_selec_str, "Beam_alt_len_" + p_type_str, evt.reco_beam_alt_len, weight, 1000., 0., 1000.);
  JSFillHist(beam_selec_str, "Beam_P_beam_inst_" + p_type_str, P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, "Beam_chi2_proton_" + p_type_str, chi2_proton, weight, 10000., 0., 1000.);

  JSFillHist(beam_selec_str, "Beam_trk_len_ratio_" + p_type_str, trk_len_ratio, weight, 1000., 0., 10.);
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
    JSFillHist("true_beam", true_beam_pdg_str + "_KE_beam_true_vs_delta_KE_true_ff", KE_beam_true, KE_beam_true - KE_ff_true, 1., 2000., 0., 2000., 500., 0., 50.);
  }
}

void Secondary_proton_ana::Study_Beam_Proton_Eloss(TString beam_cut_str, double weight){

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

      double delta_P_reco_beam_true_ff = P_beam_inst - P_ff_true;
      JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_P_ff", P_ff_true, 1., 2000., 0., 2000.);
      JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_P_beam_reco_vs_P_beam_reco_minus_P_ff_true", P_beam_inst, delta_P_reco_beam_true_ff, weight, 2000., 0., 2000., 400., -200., 200.);
      JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_delta_KE_true_ff", KE_beam_true - KE_ff_true, weight, 500., 0., 50.);
      JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_KE_beam_true_vs_delta_KE_true_ff", KE_beam_true, KE_beam_true - KE_ff_true, weight, 2000., 0., 2000., 500., 0., 50.);
    }
    
    double P_beam_true = evt.true_beam_startP * 1000.;
    double delta_P_reco_beam_true = P_beam_inst - P_beam_true;
    JSFillHist(beam_cut_str, "true_" + true_beam_pdg_str + "_P_beam_reco_vs_P_beam_reco_minus_P_beam_true", P_beam_inst, delta_P_reco_beam_true, weight, 2000., 0., 2000., 400., -200., 200.);
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
  JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_hypfit", this_P_Hypfit, weight, 2000., 0., 2000.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_reco_vs_P_beam_reco_minus_P_ff_hypfit", P_beam_inst, delta_P_reco_beam_hypfit_ff, weight, 2000., 0., 2000., 400., -200., 200.);

  double beam_res_range = (*evt.reco_beam_resRange_SCE).at(0);
  double this_KE_range = map_BB[2212] -> KEFromRangeSpline(beam_res_range);
  double this_P_range = KE_to_P(this_KE_range, 2212);
  double delta_P_reco_beam_hypfit_range = P_beam_inst - this_P_range;
  JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_range", this_P_range, weight, 2000., 0., 2000.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_reco_vs_P_beam_reco_minus_P_ff_range", P_beam_inst, delta_P_reco_beam_hypfit_range, weight, 2000., 0., 2000., 400., -200., 200.);

  double delta_P_TOF_beam_reco = P_beam_inst - P_beam_TOF;
  double delta_P_TOF_beam_range_ff = P_beam_TOF - this_P_range;
  double delta_P_TOF_beam_hypfit_ff = P_beam_TOF - this_P_Hypfit;
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_TOF_vs_P_beam_TOF_minus_P_spec", P_beam_TOF, delta_P_TOF_beam_reco, weight, 2000., 0., 2000., 400., -200., 200.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_TOF_vs_P_beam_TOF_minus_P_ff_range", P_beam_TOF, delta_P_TOF_beam_range_ff, weight, 2000., 0., 2000., 400., -200., 200.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_TOF_vs_P_beam_TOF_minus_P_ff_hypfit", P_beam_TOF, delta_P_TOF_beam_hypfit_ff, weight, 2000., 0., 2000., 400., -200., 200.);

  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_inst_vs_P_beam_TOF_minus_P_spec", P_beam_inst, delta_P_TOF_beam_reco, weight, 2000., 0., 2000., 400., -200., 200.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_inst_vs_P_beam_TOF_minus_P_ff_range", P_beam_inst, delta_P_TOF_beam_range_ff, weight, 2000., 0., 2000., 400., -200., 200.);
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_inst_vs_P_beam_TOF_minus_P_ff_hypfit", P_beam_inst, delta_P_TOF_beam_hypfit_ff, weight, 2000., 0., 2000., 400., -200., 200.);

  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_inst_vs_reco_stop_proton_P_beam_TOF", P_beam_inst, P_beam_TOF, weight, 2000., 0., 2000., 2000., 0., 2000.);

  //beam_TOF;
  JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_inst_vs_reco_stop_proton_beam_TOF", P_beam_inst, beam_TOF, weight, 2000., 0., 2000., 2000., 0., 200.);

 
  if(!IsData){
    //double delta_P_true_beam_TOF = P_beam_true - P_beam_TOF;
    //JSFillHist(beam_cut_str, "reco_stop_proton_P_beam_TOF_vs_P_beam_TOF_minus_P_spec", P_beam_TOF, delta_P_TOF_beam_reco, weight, 2000., 0., 2000., 400., -200., 200.);

    // == Using P_ff_true
    if(P_ff_true > 0.){
      JSFillHist(beam_cut_str, "true_proton_P_ff_chi2_cut", P_ff_true, weight, 2000., 0., 2000.);
      if(p_type == 2){
	JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_hypfit_true_Elas", this_P_Hypfit, weight, 2000., 0., 2000.);
	JSFillHist(beam_cut_str, "reco_stop_proton_P_ff_range_true_Elas", this_P_range, weight, 2000., 0., 2000.);
      } 
    }
  }
}

void Secondary_proton_ana::True_Daughter_study(const vector<Daughter>& protons, const vector<Daughter>& pions, double weight){

  int N_protons = protons.size();
  int N_pions = pions.size();
  JSFillHist("Daughter", "N_true_protons_by_Pandora", N_protons, weight, 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_pions_by_Pandora", N_pions, weight, 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_protons", evt.true_daughter_nProton, weight, 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_pions", evt.true_daughter_nPiPlus + evt.true_daughter_nPiMinus, weight, 10., -0.5, 9.5);

  JSFillHist("Daughter", "N_true_protons_by_Pandora_" + p_type_str, N_protons, weight, 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_pions_by_Pandora_" + p_type_str, N_pions, weight, 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_protons_" + p_type_str, evt.true_daughter_nProton, weight, 10., -0.5, 9.5);
  JSFillHist("Daughter", "N_true_pions_" + p_type_str, evt.true_daughter_nPiPlus + evt.true_daughter_nPiMinus, weight, 10., -0.5, 9.5);

  // == Loop for reco daughter matched with true protons
  for(unsigned int i = 0; i < protons.size(); i++){
    Daughter this_proton = protons.at(i);

    double this_chi2_proton = this_proton.allTrack_Chi2_proton() / this_proton.allTrack_Chi2_ndof();
    unsigned int N_hit_colletion = this_proton.allTrack_calibrated_dEdX_SCE().size();
    TString this_true_end_process = this_proton.PFP_true_byHits_endProcess();
    double true_length = this_proton.PFP_true_byHits_len();
    double KE_true_length = map_BB[2212] -> KEFromRangeSpline(true_length);
    double true_KE = this_proton.PFP_true_byHits_startE() * 1000. - M_proton;

    TVector3 unit_yz(0.,
		     this_proton.allTrack_endY() - this_proton.allTrack_startY(),
		     this_proton.allTrack_endZ() - this_proton.allTrack_startZ());
    unit_yz = (1. / unit_yz.Mag()) * unit_yz;
    TVector3 unit_z(0., 0., 1.);
    double cos_theta_yz = cos(unit_yz.Angle(unit_z));
    if(fabs(KE_true_length - true_KE) < 10. && this_true_end_process.Contains("hIoni")){
      //cout << "this_true_end_process : " << this_true_end_process << ", true_KE : " << true_KE << ", KE_true_length : " << KE_true_length << ", this_chi2_proton : " << this_chi2_proton << ", N_hit_colletion : " << N_hit_colletion << endl;
      JSFillHist("Daughter", "Reco_daughter_matched_proton_chi2_proton", this_chi2_proton, weight, 200., 0., 200.);
      JSFillHist("Daughter", "Reco_daughter_matched_proton_chi2_proton_vs_nhit_collection", this_chi2_proton, N_hit_colletion, weight, 200., 0., 200., 100., -0.5, 99.5);
      JSFillHist("Daughter", "Reco_daughter_matched_proton_chi2_proton_vs_cos_theta_yz", this_chi2_proton, cos_theta_yz, weight, 200., 0., 200., 100., -1., 1.);
    }
  }

  


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
