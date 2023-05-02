#include "dEdx_res.h"
TRandom3 gRan(1800);

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

  Run_Beam(2212);
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

  bool is_2nd_peak = false;
  int total_N_hits = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  for(int i = 0; i < total_N_hits; i++){
    JSFillHist(beam_particle, "ResRange_vs_dEdx", (*evt.reco_beam_resRange_SCE).at(i), (*evt.reco_beam_calibrated_dEdX_SCE).at(i), 1., 200., 0., 200., 1000., 0., 50.);

    double corr_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    double smear_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
    if(!IsData){
      corr_dEdx = dEdx_scaled(PID, corr_dEdx);
      smear_dEdx = dEdx_smeared(PID, corr_dEdx);
    }
    JSFillHist(beam_particle, "ResRange_vs_dEdx_corr", (*evt.reco_beam_resRange_SCE).at(i), corr_dEdx, 1., 200., 0., 200., 1000., 0., 50.);
    JSFillHist(beam_particle, "ResRange_vs_dEdx_smeared", (*evt.reco_beam_resRange_SCE).at(i), smear_dEdx, 1., 200., 0., 200., 1000., 0., 50.);

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
    if(scale > 1.2) scale = 1.2;
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

dEdx_res::dEdx_res(){

}

dEdx_res::~dEdx_res(){

}
