#include "dEdx_res.h"
TRandom3 gRan(1800);

double dEdx_from_dqdx(double dqdx, double Efield);
double Ccal_from_dqdx_dedx(double dqdx,double dedx, double Efield);

void dEdx_res::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void dEdx_res::executeEvent(){

  //cout << "test evt.beam_inst_P : " << evt.beam_inst_P * 1000. << endl;
  
  if(!PassBeamScraperCut()) return;
  if(!PassBeamMomentumWindowCut()) return;
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;
  //if(!PassAPA3Cut()) return;
  //if(!PassMichelScoreCut()) return;
  if(!PassBeamCosCut()) return;
  if(!PassBeamStartZCut()) return;
  //if(!PassProtonVetoCut()) return;
  //if(!PassMuonVetoCut()) return;
  //if(!PassStoppedPionVetoCut()) return;

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  mass_beam = 139.57;
  P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;

  JSFillHist("test", "htrack_P_beam", P_beam_inst, 1., 5000., 0., 5000.);

  // == Functions to study beam
  //Run_beam_dEdx_vector();
 
  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> pions = GetPions(daughters_all);
  vector<Daughter> protons = GetProtons(daughters_all);
  //if(pions.size() > 0) Run_Daughter(pions);

  //Run_Beam(2212);
  Run_Beam(13);
}

void dEdx_res::Run_Beam(int PID){

  if(PID == 2212 && !Pass_Beam_PID(PID)) return;
  if(PID == 13 && (evt.reco_beam_calo_endZ < 260. || evt.reco_beam_calo_endZ > 420. || !Pass_Beam_PID(PID) || daughter_michel_score < 0.6)) return;

  //if(PID == 13) cout << "muons!" << endl;
  //cout << "[dEdx_res::Run_Beam] Beam ID cut" << endl;

  double this_chi2 = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), PID, 1.);
  //cout << "this_chi2 : " << this_chi2 << endl;
  if(PID == 2212 && this_chi2 > 10.) return;
  if(PID == 13){
    vector<Daughter> daughters_all = GetAllDaughters();
    //if(daughters_all.size() > 0) return;
  }

  //cout << "[dEdx_res::Run_Beam] this_chi2 cut" << endl;
  TString beam_particle = "";
  if(PID == 2212) beam_particle = "proton";
  else if(PID == 211) beam_particle = "pion";
  else if(PID == 13) beam_particle = "muon";
  else return;
  //cout << "[dEdx_res::Run_Beam] PID cut" << endl;

  double this_chi2_muon = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 13, 1.);
  double this_chi2_pion = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 211, 1.);
  double this_chi2_proton = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 2212, 1.);
  double this_chi2_muon_skip4 = Particle_chi2_skip( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 13, 1.);
  double this_chi2_pion_skip4 = Particle_chi2_skip( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 211, 1.);
  double this_chi2_proton_skip4 = Particle_chi2_skip( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 2212, 1.);

  vector<double> corrected_dEdx_vec;
  for(unsigned int i = 0; i < (*evt.reco_beam_calibrated_dEdX_SCE).size(); i++){
    double this_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    if(!IsData) this_dEdx = dEdx_scaled(13, this_dEdx);
    corrected_dEdx_vec.push_back(this_dEdx);
  }
  double this_chi2_muon_dEdx_corr = Particle_chi2(corrected_dEdx_vec, (*evt.reco_beam_resRange_SCE), 13, 1.);
  double this_chi2_pion_dEdx_corr = Particle_chi2(corrected_dEdx_vec, (*evt.reco_beam_resRange_SCE), 211, 1.);
  double this_chi2_proton_dEdx_corr = Particle_chi2(corrected_dEdx_vec, (*evt.reco_beam_resRange_SCE), 2212, 1.);
  double this_chi2_muon_dEdx_corr_skip4 = Particle_chi2_skip(corrected_dEdx_vec, (*evt.reco_beam_resRange_SCE), 13, 1.);
  double this_chi2_pion_dEdx_corr_skip4 = Particle_chi2_skip(corrected_dEdx_vec, (*evt.reco_beam_resRange_SCE), 211, 1.);
  double this_chi2_proton_dEdx_corr_skip4 = Particle_chi2_skip(corrected_dEdx_vec, (*evt.reco_beam_resRange_SCE), 2212, 1.);

  vector<double> Abbey_recom_dEdx_vec;
  for(unsigned int i = 0; i < (*evt.reco_beam_calibrated_dEdX_SCE).size(); i++){
    //double cal_Efield = 0.5;
    double cal_Efield = (*evt.reco_beam_EField_SCE).at(i);
    //if(IsData) cal_Efield = MCCorr -> SCE_Corrected_E((*evt.reco_beam_calo_X_allTrack).at(i), (*evt.reco_beam_calo_Y_allTrack).at(i), (*evt.reco_beam_calo_Z_allTrack).at(i));
    double this_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    if(IsData) this_dEdx = Use_Abbey_Recom_Params(this_dEdx, cal_Efield, 0.9488);
    Abbey_recom_dEdx_vec.push_back(this_dEdx);
  }
  double this_chi2_muon_Abbey = Particle_chi2(Abbey_recom_dEdx_vec, (*evt.reco_beam_resRange_SCE), 13, 1.);
  double this_chi2_pion_Abbey = Particle_chi2(Abbey_recom_dEdx_vec, (*evt.reco_beam_resRange_SCE), 211, 1.);
  double this_chi2_proton_Abbey = Particle_chi2(Abbey_recom_dEdx_vec, (*evt.reco_beam_resRange_SCE), 2212, 1.);
  double this_chi2_muon_Abbey_skip4 = Particle_chi2_skip(Abbey_recom_dEdx_vec, (*evt.reco_beam_resRange_SCE), 13, 1.);
  double this_chi2_pion_Abbey_skip4 = Particle_chi2_skip(Abbey_recom_dEdx_vec, (*evt.reco_beam_resRange_SCE), 211, 1.);
  double this_chi2_proton_Abbey_skip4 = Particle_chi2_skip(Abbey_recom_dEdx_vec, (*evt.reco_beam_resRange_SCE), 2212, 1.);

  JSFillHist(beam_particle, "Chi2_muon", this_chi2_muon, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_pion", this_chi2_pion, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_proton", this_chi2_proton, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_muon_dEdx_corr", this_chi2_muon_dEdx_corr, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_pion_dEdx_corr", this_chi2_pion_dEdx_corr, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_proton_dEdx_corr", this_chi2_proton_dEdx_corr, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_muon_Abbey", this_chi2_muon_Abbey, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_pion_Abbey", this_chi2_pion_Abbey, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_proton_Abbey", this_chi2_proton_Abbey, 1., 1000., 0., 1000.);

  JSFillHist(beam_particle, "Chi2_muon_skip4", this_chi2_muon_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_pion_skip4", this_chi2_pion_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_proton_skip4", this_chi2_proton_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_muon_dEdx_corr_skip4", this_chi2_muon_dEdx_corr_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_pion_dEdx_corr_skip4", this_chi2_pion_dEdx_corr_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_proton_dEdx_corr_skip4", this_chi2_proton_dEdx_corr_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_muon_Abbey_skip4", this_chi2_muon_Abbey_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_pion_Abbey_skip4", this_chi2_pion_Abbey_skip4, 1., 1000., 0., 1000.);
  JSFillHist(beam_particle, "Chi2_proton_Abbey_skip4", this_chi2_proton_Abbey_skip4, 1., 1000., 0., 1000.);

  bool is_2nd_peak = false;
  int total_N_hits = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  for(int i = 0; i < total_N_hits; i++){
    JSFillHist(beam_particle, "ResRange_vs_dEdx", (*evt.reco_beam_resRange_SCE).at(i), (*evt.reco_beam_calibrated_dEdX_SCE).at(i), 1., 200., 0., 200., 1000., 0., 50.);

    double this_dqdx = (*evt.reco_beam_dQdX_SCE).at(i);
    double this_dqdx_calib = (*evt.reco_beam_calibrated_dQdX_SCE).at(i);
    //cout << "this_dqdx : " << this_dqdx << ", this_dqdx_calib : " << this_dqdx_calib << ", ratio : " << this_dqdx_calib / this_dqdx << endl;
    double this_cal_dEdx = dEdx_from_dqdx(this_dqdx_calib, (*evt.reco_beam_EField_SCE).at(i));

    double corr_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    double smear_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    double Abbey_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    if(!IsData){
      corr_dEdx = dEdx_scaled(PID, corr_dEdx);
      smear_dEdx = dEdx_smeared(PID, corr_dEdx);
    }

    double this_Ccal = Ccal_from_dqdx_dedx(this_dqdx, Abbey_dEdx, (*evt.reco_beam_EField_SCE).at(i));
    //cout << "this_cal_dEdx : " << this_cal_dEdx << ", PDSP dEdx : " << Abbey_dEdx << ", this_Ccal : " << this_Ccal << endl;

    double dEdx_central = Abbey_dEdx;
    double dEdx_alpha_down_beta_down = Abbey_dEdx;
    double dEdx_alpha_down_beta_central = Abbey_dEdx;
    double dEdx_alpha_down_beta_up = Abbey_dEdx;
    double dEdx_alpha_central_beta_down = Abbey_dEdx;
    double dEdx_alpha_central_beta_up = Abbey_dEdx;
    double dEdx_alpha_up_beta_down = Abbey_dEdx;
    double dEdx_alpha_up_beta_central = Abbey_dEdx;
    double dEdx_alpha_up_beta_up = Abbey_dEdx;

    if(IsData){
      double cal_Efield = MCCorr -> SCE_Corrected_E((*evt.reco_beam_calo_X_allTrack).at(i), (*evt.reco_beam_calo_Y_allTrack).at(i), (*evt.reco_beam_calo_Z_allTrack).at(i));
      //double cal_Efield = (*evt.reco_beam_EField_SCE).at(i);
      //Abbey_dEdx = Use_Abbey_Recom_Params(Abbey_dEdx, (*evt.reco_beam_EField_SCE).at(i), 0.9488);
      //cout << "Before, Abbey_dEdx : " << Abbey_dEdx << endl;
      //Abbey_dEdx = Use_Abbey_Recom_Params(Abbey_dEdx, cal_Efield, 0.9488);
      Abbey_dEdx = MCCorr->Use_Abbey_Recom_Params(Abbey_dEdx, cal_Efield, 1.);
      //cout << "After, Abbey_dEdx : " << Abbey_dEdx << endl;
      //cout << "------------" << endl;
      dEdx_alpha_down_beta_down    = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.91, 0.210, 0.957);
      dEdx_alpha_down_beta_central = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.91, 0.212, 0.958);
      dEdx_alpha_down_beta_up      = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.91, 0.214, 0.958);
      dEdx_alpha_central_beta_down = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.93, 0.210, 0.995);
      dEdx_alpha_central_beta_up   = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.93, 0.214, 0.995);
      dEdx_alpha_up_beta_down      = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.95, 0.210, 1.030);
      dEdx_alpha_up_beta_central   = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.95, 0.212, 1.034);
      dEdx_alpha_up_beta_up        = MCCorr->Use_Other_Mod_Box_Params(dEdx_central, cal_Efield, 0.95, 0.214, 1.033);

      /*
      cout << "dEdx_central : " << dEdx_central << endl;
      cout << "dEdx_alpha_down_beta_down : " << dEdx_alpha_down_beta_down << endl;
      cout << "dEdx_alpha_down_beta_central : " << dEdx_alpha_down_beta_central << endl;
      cout << "dEdx_alpha_down_beta_up : " << dEdx_alpha_down_beta_up << endl;
      cout << "dEdx_alpha_central_beta_down : " << dEdx_alpha_central_beta_down << endl;
      cout << "dEdx_alpha_central_beta_up : " << dEdx_alpha_central_beta_up << endl;
      cout << "dEdx_alpha_up_beta_down : " << dEdx_alpha_up_beta_down << endl;
      cout << "dEdx_alpha_up_beta_central : " << dEdx_alpha_up_beta_central << endl;
      cout << "dEdx_alpha_up_beta_up : " << dEdx_alpha_up_beta_up << endl;
      */
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_down_beta_down", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_down_beta_down, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_down_beta_central", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_down_beta_central, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_down_beta_up", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_down_beta_up, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_central_beta_down", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_central_beta_down, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_central_beta_up", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_central_beta_up, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_up_beta_down", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_up_beta_down, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_up_beta_central", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_up_beta_central, 1., 200., 0., 200., 1000., 0., 50.);
      JSFillHist(beam_particle, "ResRange_vs_dEdx_alpha_up_beta_up", (*evt.reco_beam_resRange_SCE).at(i), dEdx_alpha_up_beta_up, 1., 200., 0., 200., 1000., 0., 50.);
    }

    //if(fabs(Abbey_dEdx - (*evt.reco_beam_calibrated_dEdX_SCE).at(i)) > 0.0001) cout << "Case" << endl;
    JSFillHist(beam_particle, "ResRange_vs_dEdx_corr", (*evt.reco_beam_resRange_SCE).at(i), corr_dEdx, 1., 200., 0., 200., 1000., 0., 50.);
    JSFillHist(beam_particle, "ResRange_vs_dEdx_smeared", (*evt.reco_beam_resRange_SCE).at(i), smear_dEdx, 1., 200., 0., 200., 1000., 0., 50.);
    JSFillHist(beam_particle, "ResRange_vs_dEdx_Abbey", (*evt.reco_beam_resRange_SCE).at(i), Abbey_dEdx, 1., 200., 0., 200., 1000., 0., 50.);

    if(i != 0 && i != total_N_hits - 1){
      double this_pitch = 0.5 * fabs((*evt.reco_beam_resRange_SCE).at(i + 1) - (*evt.reco_beam_resRange_SCE).at(i - 1));
      JSFillHist(beam_particle, "ResRange_vs_pitch", (*evt.reco_beam_resRange_SCE).at(i), this_pitch, 1., 200., 0., 200., 200., 0., 2.);
    }
    if(PID == 13 && (*evt.reco_beam_resRange_SCE).at(i) < 5. && (*evt.reco_beam_calibrated_dEdX_SCE).at(i) > 1.5 && (*evt.reco_beam_calibrated_dEdX_SCE).at(i) < 2.5){
      is_2nd_peak = true;
    }
  }

  /*
  if(PID == 13){
    JSFillHist(beam_particle, "Beam_end_z", evt.reco_beam_calo_endZ, 1., 1000., 0., 1000.);
    if(is_2nd_peak){
      vector<Daughter> daughters_all = GetAllDaughters();
      cout << "==============" << endl;
      cout << "daughters_all.size() : " << daughters_all.size() << endl;
      cout << "evt.reco_beam_calo_endZ : " << evt.reco_beam_calo_endZ << endl;
      for(unsigned int i_d = 0; i_d < daughters_all.size(); i_d++){
	cout << Form("daughters_all.at(%d).Beam_Cos() : %f", i_d, daughters_all.at(i_d).Beam_Cos()) << endl;
	cout << Form("daughters_all.at(%d).allTrack_resRange_SCE().back() : %f", i_d, daughters_all.at(i_d).allTrack_resRange_SCE().back()) << endl;
	cout << Form("daughters_all.at(%d).allTrack_startZ() : : %f", i_d, daughters_all.at(i_d).allTrack_startZ()) << endl;
      }
    }
  }
  */

  return;
}

double dEdx_res::dEdx_scaled(int PID, double MC_dEdx){
  double f_const_p0 = 0.997;
  double f_pol1_p0 = 0.945;
  double f_pol1_p1 = 0.021;

  double scale = 1.;
  if(PID = 13){
    if(MC_dEdx < 2.44) scale = f_const_p0;
    else scale = f_pol1_p0 + MC_dEdx * f_pol1_p1;
    //if(scale > 1.2) scale = 1.2;
  }

  //cout << "[dEdx_res::dEdx_scaled] scale : " << scale << endl;
  return scale * MC_dEdx;
}

double dEdx_res::dEdx_smeared(int PID, double MC_dEdx){
  double func_p0 = 6.99163;
  double func_p1 = -1.71138;
  double func_p2 = 2.02863;

  double smear_sigma = 0.;
  if(PID = 13){
    if(MC_dEdx > 2.27) smear_sigma = fabs(func_p0 + func_p1 /(MC_dEdx - func_p2));
  }

  smear_sigma = 0.06;
  double this_gaus = gRan.Gaus(0,smear_sigma * 0.01);
  double out = MC_dEdx * (1. + this_gaus); 

  //cout << "[dEdx_res::dEdx_smeared] MC_dEdx : " << MC_dEdx << ", smear_sigma : " << smear_sigma << ", this_gaus : " << this_gaus << ", out : " << out << endl;

  return out;
}


double dEdx_from_dqdx(double dqdx, double Efield){

  double alpha = 0.93;
  double beta = 0.212;
  double Rho = 1.40;
  double Wion = 23.6e-6;
  double Ccal = 1.038e-3;

  return (exp((dqdx/Ccal)*(beta/(Rho*Efield)*Wion))-alpha)/(beta/(Rho*Efield));

}

double Ccal_from_dqdx_dedx(double dqdx, double dedx, double Efield){
  double alpha = 0.93;
  double beta = 0.212;
  double Rho = 1.40;
  double Wion = 23.6e-6;

  double Ccal = ((beta * Wion * dqdx) / (Rho * Efield)) * (1 / (log(alpha + (beta * dedx) / (Rho * Efield)) ) );

  return Ccal;
}

double dEdx_res::Use_Abbey_Recom_Params(double dEdx, double Efield, double calib_const_ratio){

  double alpha_default = 0.93;
  double beta_default = 0.212;
  //double rho = 1.396;
  double rho = 1.39;

  double alpha_Abbey = 0.905;
  double beta_Abbey = 0.220;

  double exp_term = exp( calib_const_ratio * (beta_Abbey / beta_default) * log(alpha_default + beta_default * dEdx / (rho * Efield)) );
  //exp_term = exp_term;
  double new_dEdx = (exp_term - alpha_Abbey) * rho * Efield / beta_Abbey;

  //cout << "[dEdx_res::Use_Abbey_Recom_Params] dEdx : " << dEdx << ", new_dEdx : " << new_dEdx << endl;
  /*
  unsigned int N_hits = (*evt.reco_beam_resRange_SCE).size();
  cout << "(*evt.reco_beam_resRange_SCE).size() : " << (*evt.reco_beam_resRange_SCE).size() << ", (*evt.reco_beam_calo_X_allTrack).size() : " << (*evt.reco_beam_calo_X_allTrack).size() << endl;
  double cal_res_range = (*evt.reco_beam_resRange_SCE).back();
  if(N_hits < 1) return 0.;
  for(int i = N_hits - 1; i >= 0; i--){
    double this_pitch = 0.;
    if(i != N_hits - 1) this_pitch = pow( pow((*evt.reco_beam_calo_X_allTrack).at(i) - (*evt.reco_beam_calo_X_allTrack).at(i + 1),2)
					  + pow((*evt.reco_beam_calo_Y_allTrack).at(i) - (*evt.reco_beam_calo_Y_allTrack).at(i + 1),2)
					  + pow((*evt.reco_beam_calo_Z_allTrack).at(i) - (*evt.reco_beam_calo_Z_allTrack).at(i + 1),2) , 0.5);
    cal_res_range = cal_res_range + this_pitch;

    double PDSP_Efield = (*evt.reco_beam_EField_SCE).at(i);
    double cal_Efield = MCCorr -> SCE_Corrected_E((*evt.reco_beam_calo_X_allTrack).at(i), (*evt.reco_beam_calo_Y_allTrack).at(i), (*evt.reco_beam_calo_Z_allTrack).at(i));
    cout << "[dEdx_res::Use_Abbey_Recom_Params] " << i << "'th res_range : " << (*evt.reco_beam_resRange_SCE).at(i) << ", cal_res_range : " << cal_res_range
	 << ", PDSP_Efield : " << PDSP_Efield << ", cal_Efield : " << cal_Efield << endl;
  }
  */

  return new_dEdx;
}

dEdx_res::dEdx_res(){

}

dEdx_res::~dEdx_res(){

}
