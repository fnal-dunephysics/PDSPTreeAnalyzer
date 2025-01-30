#include "pionbeamsel.h"
#include "TLorentzVector.h"

void pionbeamsel::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void pionbeamsel::executeEvent(){

  pi_type = GetPi2ParType();
  pi_type_str = Form("%d", pi_type);
  FillHist("beam_cut_flow", 0.5, 1., 20, 0., 20.);

  pi_truetype = GetPiTrueType();
  //if (pi_truetype == pitrue::kQE) FillQEPlots("QE_All");

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  KE_beam_inst = map_BB[211] -> MomentumtoKE(P_beam_inst);
  exp_trk_len_beam_inst = map_BB[211] -> RangeFromKESpline(KE_beam_inst);
  trk_len_ratio = evt.reco_beam_alt_len / exp_trk_len_beam_inst;
  mass_beam = 139.57;
  KE_ff_reco = KE_beam_inst - 30.; // -- checked upstream KE loss using stopping muons, central value is 30 MeV but we need to assign syst. uncert. on it
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;

  double KE_ff_reweight = 1.;
  if(!IsData) KE_ff_reweight = MCCorr -> MomentumReweight_SF("Pion_chip_KE_0p5", KE_ff_reco, 0.);
  
  // -- 1. Beam instruments
  if(P_beam_inst < 420. || P_beam_inst > 580.) return;
  //if(!PassBeamMomentumWindowCut()) return;

  if(!Pass_Beam_PID(211)) return;
  FillBeamPlots("Beam_PID", 1.);

  if(!PassBeamScraperCut()) return;
  FillHist("beam_cut_flow", 1.5, 1., 20, 0., 20.);
  FillBeamPlots("Beam_scraper", 1.);

  if(!Pass_BeamStartZ(2.)) return;
  FillBeamPlots("Beam_startZ", 1.);

  if(!Pass_beam_delta_X_cut(2.)) return;
  if(!Pass_beam_delta_Y_cut(2.)) return;
  FillBeamPlots("Beam_deltaXY", 1.);

  if(chi2_proton > 300. || chi2_proton < 50.) return;
  FillBeamPlots("Beam_chi2proton", KE_ff_reweight);
  
  if(trk_len_ratio > 1.){
    //MuonKELoss("Beam_chi2proton", 1.);
  }
}

void pionbeamsel::FillBeamPlots(TString beam_selec_str, double weight){

  // == Fit results after calo-size cut
  double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - beam_start_z_mu_data) / beam_start_z_sigma_data;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData){
     Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - beam_start_z_mu_mc) / beam_start_z_sigma_mc;
     Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;
     Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;
  }
  
  // == Comparison between beam spectrometer track and TPC reco track
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX", evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY", evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ", evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ", evt.reco_beam_calo_endZ, weight, 1000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_beam_inst_X", evt.beam_inst_X, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_beam_inst_Y", evt.beam_inst_Y, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst", P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh", beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton", chi2_proton, weight, 10000., 0., 1000.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX_" + pi_type_str, evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY_" + pi_type_str, evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_" + pi_type_str, evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_inst_X_" + pi_type_str, evt.beam_inst_X, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_inst_Y_" + pi_type_str, evt.beam_inst_Y, weight, 10000., -100., 900.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_" + pi_type_str, delta_X_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_" + pi_type_str, delta_Y_spec_TPC, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_cos_delta_spec_TPC_" + pi_type_str, cos_delta_spec_TPC, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh_" + pi_type_str, beam_costh, weight, 2000., -1., 1.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_theta_" + pi_type_str, beam_TPC_theta, weight, 5000., -1., 4.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_phi_" + pi_type_str, beam_TPC_phi, weight, 8000., -4., 4.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ_" + pi_type_str, evt.reco_beam_calo_endZ, weight, 1000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_alt_len_" + pi_type_str, evt.reco_beam_alt_len, weight, 1000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst_" + pi_type_str, P_beam_inst, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton_" + pi_type_str, chi2_proton, weight, 10000., 0., 1000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_trk_len_ratio_" + pi_type_str, trk_len_ratio, weight, 1000., 0., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_ff_" + pi_type_str, KE_ff_reco, weight, 2000., 0., 2000.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_end_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);

  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_over_sigma_" + pi_type_str, Beam_startZ_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_X_over_sigma, weight, 2000., -10., 10.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_Y_over_sigma, weight, 2000., -10., 10.);
}

void pionbeamsel::MuonKELoss(TString beam_selec_str, double weight){

  double this_muon_beam_inst_KE = map_BB[13] -> MomentumtoKE(P_beam_inst);

  double this_rr_at_z10cm = GetBeamRRatZ10cm((*evt.reco_beam_resRange_SCE), (*evt.reco_beam_calo_Z));
  double this_muon_KE_z10cm = map_BB[13] -> KEFromRangeSpline(this_rr_at_z10cm);

  double this_KELoss = this_muon_beam_inst_KE - this_muon_KE_z10cm;
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KELoss_" + pi_type_str, this_KELoss, weight, 2000., -100., 100.);
  JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_beam_inst_vs_Beam_KELoss_" + pi_type_str, this_muon_beam_inst_KE, this_KELoss, weight, 100., 0., 1000., 200., -100., 100.);
}

double pionbeamsel::GetBeamRRatZ10cm(const vector<double> & ResRange, const vector<double> & calo_Z){

  double out = -1.;
  if(ResRange.size() < 1 || calo_Z.size() < 1) return out;

  
  int this_N_calo = calo_Z.size();
  //cout << "ResRange.size(): " << ResRange.size() << ", calo_Z.size(): " << calo_Z.size() << endl; // -- confirmed that the two vectors have exactly the same sizes
  
  for(int i = 0; i < this_N_calo; i++){
    //cout << "ResRange " << i << ": " << ResRange.at(i) << ", calo_Z: " << calo_Z.at(i) << endl; // -- confirmed that rr is in descending order
    double this_calo_Z = calo_Z.at(i);
    if(this_calo_Z > 10.){
      double this_rr = ResRange.at(i);
      double prev_rr = ResRange.at(i - 1);
      double prev_calo_Z = calo_Z.at(i - 1);

      double this_ratio = (10. - prev_calo_Z) / (this_calo_Z - prev_calo_Z);
      double this_delta_rr = prev_rr - this_rr;
      out = prev_rr - this_delta_rr * this_ratio;
      //cout << "this_rr_at_z10cm: " << this_rr_at_z10cm << endl; // -- confirmed that the output value is reasonable
      break;
    }
  }
  
  return out;
}

bool pionbeamsel::Pass_BeamStartZ(double N_sigma){
  double this_beam_start_Z = evt.reco_beam_calo_startZ;
  bool out = false;
  if(IsData){
    out = fabs(beam_start_z_mu_data - this_beam_start_Z) < (N_sigma * beam_start_z_sigma_data);
  }
  else{
    out	= fabs(beam_start_z_mu_mc - this_beam_start_Z) < (N_sigma * beam_start_z_sigma_mc);
  }

  return out;
}

bool pionbeamsel::Pass_beam_delta_X_cut(double N_sigma){

  bool out = false;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  if(!IsData) Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;

  if(fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

  return out;
}

bool pionbeamsel::Pass_beam_delta_Y_cut(double N_sigma){

  bool out = false;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData) Beam_delta_Y_over_sigma =(delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;

  if(fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

  return out;
}

pionbeamsel::pionbeamsel(){

}

pionbeamsel::~pionbeamsel(){

}
