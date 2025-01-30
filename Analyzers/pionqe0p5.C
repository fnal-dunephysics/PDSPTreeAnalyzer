#include "pionqe0p5.h"
#include "TLorentzVector.h"

void pionqe0p5::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void pionqe0p5::executeEvent(){

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
  if(!PassBeamScraperCut()) return;
  FillHist("beam_cut_flow", 1.5, 1., 20, 0., 20.);

  if(!Pass_BeamStartZ(2.)) return;
  if(!Pass_beam_delta_X_cut(2.)) return;
  if(!Pass_beam_delta_Y_cut(2.)) return;
  if(chi2_proton > 300. || chi2_proton < 50.) return;

  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> loose_pions = SelectLoosePions(daughters_all);

  FillRecoPionPlots("loose_pion", loose_pions, KE_ff_reweight);
}

void pionqe0p5::FillRecoPionPlots(TString daughter_sec_str, const vector<Daughter> pions, double weight){

  JSFillHist(daughter_sec_str, daughter_sec_str + "_N_reco_pions_" + pi_type_str, pions.size(), weight, 10., -0.5, 9.5);
  if(pions.size() == 0){
    JSFillHist(daughter_sec_str, daughter_sec_str + "_Beam_KE_end_0pi_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);
  }
  else{
    JSFillHist(daughter_sec_str, daughter_sec_str + "_Beam_KE_end_least1pi_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);
  }
}

double pionqe0p5::GetBeamRRatZ10cm(const vector<double> & ResRange, const vector<double> & calo_Z){

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

bool pionqe0p5::Pass_BeamStartZ(double N_sigma){
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

bool pionqe0p5::Pass_beam_delta_X_cut(double N_sigma){

  bool out = false;
  double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
  if(!IsData) Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;

  if(fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

  return out;
}

bool pionqe0p5::Pass_beam_delta_Y_cut(double N_sigma){

  bool out = false;
  double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
  if(!IsData) Beam_delta_Y_over_sigma =(delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;

  if(fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

  return out;
}

std::vector<Daughter> pionqe0p5::SelectLoosePions(const vector<Daughter>& in){

  vector<Daughter> out;
  double cut_trackScore = 0.5;
  double cut_chi2_proton = 60.;
  double cut_trk_len_upper = 180.;
  double cut_trk_len_lower = 10.;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_chi2 > cut_chi2_proton && this_in.allTrack_alt_len() < cut_trk_len_upper && this_in.allTrack_alt_len() > cut_trk_len_lower) out.push_back(this_in);
  }
  return out;
}


pionqe0p5::pionqe0p5(){

}

pionqe0p5::~pionqe0p5(){

}
