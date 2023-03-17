#include "MCS_Performance.h"

void MCS_Performance::initializeAnalyzer(){

  cout << "[[MCS_Performance::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void MCS_Performance::executeEvent(){

  if(!PassBeamScraperCut()) return;
  if(!PassBeamMomentumWindowCut()) return;
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;
  if(!PassAPA3Cut()) return;
  if(!PassMichelScoreCut()) return;
  if(!PassBeamCosCut()) return;
  if(!PassBeamStartZCut()) return;
  if(!PassProtonVetoCut()) return;
  if(!PassMuonVetoCut()) return;
  if(!PassStoppedPionVetoCut()) return;

  pi_type = GetPiParType();
  pi_type_str = Form("%d", pi_type);

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  mass_beam = 139.57;
  P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;

  JSFillHist("test", "htrack_P_beam", P_beam_inst, 1., 5000., 0., 5000.);

  // == Functions to study beam
  Run_Beam();

  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> pions = GetPions(daughters_all);
  if(pions.size() > 0){
    Run_Daughter(pions);
  }

  return;
}

void MCS_Performance::Run_Beam(){

  int total_N_hits = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  int this_N_hits = total_N_hits;
  double this_KE_BB = KE_ff_reco;
  int this_PDG = 211;
  if(pi_type == 0 || pi_type == 1 || pi_type == 2) this_PDG = 211;
  else if(pi_type == 3) this_PDG = 13;
  else return;

  vector<TVector3> true_position_vec;
  vector<TVector3> true_Pvec_vec;
  vector<double> true_P_vec;
  for(unsigned int i_true_hit = 0; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++){
    TVector3 this_position((*evt.true_beam_traj_X).at(i_true_hit), (*evt.true_beam_traj_Y).at(i_true_hit), (*evt.true_beam_traj_Z).at(i_true_hit));
    TVector3 this_Pvec((*evt.true_beam_traj_Px).at(i_true_hit), (*evt.true_beam_traj_Py).at(i_true_hit), (*evt.true_beam_traj_Pz).at(i_true_hit));
    this_Pvec = this_Pvec * 1000.;// * P_beam_inst_scale;
    double this_P = this_Pvec.Mag();

    true_position_vec.push_back(this_position);
    true_Pvec_vec.push_back(this_position);
    true_P_vec.push_back(this_P);
    //cout << Form("[PionKEScale::Run_beam_MCS] %d this_P : %.2f", i_true_hit, this_P) << endl;
  }

  vector<double> reco_X = (*evt.reco_beam_calo_X);
  vector<double> reco_Y = (*evt.reco_beam_calo_Y);
  vector<double> reco_Z = (*evt.reco_beam_calo_Z);
  vector<TVector3> reco_position_vec;
  for(unsigned int i_reco_hit = 0; i_reco_hit < reco_Z.size(); i_reco_hit++){
    TVector3 this_position(reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit));
    reco_position_vec.push_back(this_position);
    //cout << Form("[PionKEScale::Run_beam_MCS] reco %d : (%.2f, %.2f, %.2f)", i_reco_hit, reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit)) << endl; 
  }
  
  
  TVector3 first_reco_hit = reco_position_vec.at(0);
  double true_start_P = -999.;
  for(unsigned int i = 0; i < true_position_vec.size() - 1; i++){
    TVector3 this_true_hit = true_position_vec.at(i);
    TVector3 next_true_hit = true_position_vec.at(i + 1);
    if(first_reco_hit.Z() > this_true_hit.Z() && first_reco_hit.Z() < next_true_hit.Z()){
      double distance = (this_true_hit - first_reco_hit).Mag();
      if(distance < 10.) true_start_P = map_BB[this_PDG] -> KEtoMomentum(map_BB[this_PDG] -> KEAtLength(map_BB[this_PDG] -> MomentumtoKE(true_P_vec.at(i)), distance));
      break;
    }
  }

  /*
  Run_Beam_Segments(reco_position_vec, true_start_P, this_PDG, 4., "4cm");
  Run_Beam_Segments(reco_position_vec, true_start_P, this_PDG, 5., "5cm");
  Run_Beam_Segments(reco_position_vec, true_start_P, this_PDG, 8., "8cm");
  Run_Beam_Segments(reco_position_vec, true_start_P, this_PDG, 10., "10cm"); 
  */  
  Run_Beam_Segments(reco_position_vec, true_start_P, this_PDG, 14., "14cm");

  
  return;
}

void MCS_Performance::Run_Beam_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name){

  unsigned int min_N_segment = 3;
  vector<MCSSegment> segments = SplitIntoSegments(reco_position_vec, segment_size);
  if(segments.size() < min_N_segment) return;

  TString this_PdgID_str = Form("%d", this_PdgID);
  double total_range = 0.;
  for(unsigned int i = 0; i < segments.size(); i++){
    double this_segment_range = segments.at(i).Range();
    total_range = total_range + this_segment_range;
  }
  double KE_from_range = map_BB[this_PdgID] -> KEFromRangeSpline(total_range);
  double P_from_range = map_BB[this_PdgID] -> KEtoMomentum(KE_from_range);

  vector<TH1D*> this_likelihood_hist = Calculate_Likelihoods_for_Performance(segments, segment_size, this_PdgID);
  //cout << "[MCS_Performance::Run_Beam_Segments] segments.size() : " << segments.size() << ", this_likelihood_hist.size() : " << this_likelihood_hist.size() << endl;
  for(unsigned int i = 0; i < this_likelihood_hist.size(); i++){
    TString this_hist_name = Form("Run%d_Evt%d_Segment%d_likelihood", evt.run, evt.event, i);
    this_likelihood_hist.at(i) -> SetName(this_hist_name);
    //maphist_TH1D[this_hist_name] = (TH1D*)this_likelihood_hist.at(i) -> Clone();
  }
  TH1D *this_likelihood_sum = (TH1D*)this_likelihood_hist.at(0) -> Clone();
  for(unsigned int i = 1; i < this_likelihood_hist.size(); i++){
    this_likelihood_sum -> Add(this_likelihood_hist.at(i));
    if(i > 1){
      int this_N_segment = i + 1;
      TString this_N_segments_str = Get_N_segment_str(segments.size());
      int this_bin_max = this_likelihood_sum -> GetMaximumBin();
      double P_fitted = this_likelihood_sum -> GetBinCenter(this_bin_max);
      double res = (P_fitted - true_P) / true_P;
      double inv_res = (1./P_fitted - 1./true_P) / (1./true_P);
      JSFillHist("Beam", "Beam_MCS_Res_" + name + "_" + this_N_segments_str + "_" + pi_type_str, res, 1., 400., -2., 2.);
      JSFillHist("Beam", "Beam_MCS_InvRes_" + name + "_" + this_N_segments_str + "_" + pi_type_str, inv_res, 1., 400., -2., 2.);
      JSFillHist("Beam", "Beam_MCS_P_true_vs_P_MCS_" + name + "_" + this_N_segments_str + "_" + pi_type_str, true_P, P_fitted, 1., 400., 0., 2000., 400., 0., 2000.);
    }
  }

  TString this_hist_name = Form("Run%d_Evt%d_SegmentSuumed_likelihood", evt.run, evt.event);
  this_likelihood_sum -> SetName(this_hist_name);
  //maphist_TH1D[this_hist_name] = (TH1D*)this_likelihood_sum -> Clone();

  int bin_max_y = this_likelihood_sum -> GetMaximumBin();
  double P_fitted = this_likelihood_sum -> GetBinCenter(bin_max_y);
  double res = (P_fitted - true_P) / true_P;
  //cout << Form("[MCS_Performance::Run_Beam_Segments] " + name + " Run%d_Evt%d, PID : %d, true_P : %.2f, Start P : %.2f, P_fitted : %.2f, res : %.2f", evt.run, evt.event, this_PdgID, true_P, P_from_range, P_fitted, res) << endl;
  JSFillHist("Beam", "Beam_MCS_Res_" + name + "_" + pi_type_str, res, 1., 400., -2., 2.);

  // -- Clean up
  for(unsigned i = 0; i < this_likelihood_hist.size(); i++){
    delete this_likelihood_hist.at(i);
  }
  this_likelihood_hist.clear();

  return;
}

void MCS_Performance::Run_Daughter(const vector<Daughter>& pions){

  int true_beam_ID = evt.true_beam_ID;
  double beam_last_X = (*evt.reco_beam_calo_X).back();
  double beam_last_Y = (*evt.reco_beam_calo_Y).back();
  double beam_last_Z = (*evt.reco_beam_calo_Z).back();
  TVector3 beam_end(beam_last_X, beam_last_Y, beam_last_Z);

  for(unsigned int i_pion = 0; i_pion < pions.size(); i_pion++){
    Daughter this_daughter = pions.at(i_pion);

    int this_true_ID = this_daughter.PFP_true_byHits_ID();
    if(true_beam_ID == this_true_ID) continue;

    int this_PdgID = this_daughter.PFP_true_byHits_PDG();

    vector<double> reco_X = this_daughter.allTrack_calo_X();
    vector<double> reco_Y = this_daughter.allTrack_calo_Y();
    vector<double> reco_Z = this_daughter.allTrack_calo_Z();
    vector<TVector3> reco_position_vec;
    for(int i_reco_hit = reco_Z.size() - 1; i_reco_hit > -1; i_reco_hit--){
      TVector3 this_position(reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit));
      reco_position_vec.push_back(this_position);
    }
    double true_P = this_daughter.PFP_true_byHits_startP() * 1000.;

    if(abs(this_PdgID) == 13 || abs(this_PdgID) == 211 || abs(this_PdgID) == 2212){
      /*
      Run_Daughter_Segments(reco_position_vec, true_P, this_PdgID, 4., "4cm");
      Run_Daughter_Segments(reco_position_vec, true_P, this_PdgID, 5., "5cm");
      Run_Daughter_Segments(reco_position_vec, true_P, this_PdgID, 8., "8cm");
      Run_Daughter_Segments(reco_position_vec, true_P, this_PdgID, 10., "10cm");
      */
      Run_Daughter_Segments(reco_position_vec, true_P, this_PdgID, 14., "14cm");
    }
  }

  return;
}

void MCS_Performance::Run_Daughter_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name){

  vector<MCSSegment> segments = SplitIntoSegments(reco_position_vec, segment_size);
  if(segments.size() < 3) return;

  TString PdgID_sign = "";
  if(this_PdgID < 0){
    PdgID_sign = "m";
  }
  else if(this_PdgID > 0){
    PdgID_sign = "p";
  }
  TString this_PdgID_str = Form(PdgID_sign + "%d", abs(this_PdgID));
  double total_range = 0.;
  for(unsigned int i = 0; i < segments.size(); i++){
    double this_segment_range = segments.at(i).Range();
    total_range = total_range + this_segment_range;
  }
  double KE_from_range = map_BB[abs(this_PdgID)] -> KEFromRangeSpline(total_range);
  double P_from_range = map_BB[abs(this_PdgID)] -> KEtoMomentum(KE_from_range);

  vector<TH1D*> this_likelihood_hist = Calculate_Likelihoods_for_Performance(segments, segment_size, this_PdgID);
  //cout << "[MCS_Performance::Run_Daughter_Segments] segments.size() : " << segments.size() << ", this_likelihood_hist.size() : " << this_likelihood_hist.size() << endl;
  TH1D *this_likelihood_sum = (TH1D*)this_likelihood_hist.at(0) -> Clone();
  for(unsigned int i = 1; i < this_likelihood_hist.size(); i++){
    this_likelihood_sum -> Add(this_likelihood_hist.at(i));
    if(i > 1){
      int this_N_segment = i + 1;
      TString this_N_segments_str = Get_N_segment_str(segments.size());
      int this_bin_max = this_likelihood_sum -> GetMaximumBin();
      double P_fitted = this_likelihood_sum -> GetBinCenter(this_bin_max);
      double res = (P_fitted - true_P) / true_P;
      double inv_res = (1./P_fitted - 1./true_P) / (1./true_P);
      JSFillHist("Daughter", "Daughter_MCS_Res_" + name + "_" + this_N_segments_str + "_" + this_PdgID_str, res, 1., 400., -2., 2.);
      JSFillHist("Daughter", "Daughter_MCS_InvRes_" + name + "_" + this_N_segments_str + "_" + this_PdgID_str, inv_res, 1., 400., -2., 2.);
      JSFillHist("Daughter", "Daughter_MCS_P_true_vs_P_MCS_" + name + "_" + this_N_segments_str + "_" + this_PdgID_str, true_P, P_fitted, 1., 400., 0., 2000., 400., 0., 2000.);
      JSFillHist("Daughter", "Daughter_MCS_P_true_vs_MCS_Res_" + name + "_" + this_N_segments_str + "_" + this_PdgID_str, true_P, res, 1., 400., 0., 2000., 400., -2., 2.);
      JSFillHist("Daughter", "Daughter_MCS_P_true_vs_MCS_InvRes_" + name + "_" + this_N_segments_str + "_" + this_PdgID_str, true_P, inv_res, 1., 400., 0., 2000., 400., -2., 2.);
    }
  }
  
  int bin_max_y = this_likelihood_sum -> GetMaximumBin();
  double P_fitted = this_likelihood_sum -> GetBinCenter(bin_max_y);
  double res = (P_fitted - true_P) / true_P;
  //cout << Form("[MCS_Performance::Run_Daughter_Segments] " + name + " Run%d_Evt%d, PID : %d, true_P : %.2f, Start P : %.2f, P_fitted : %.2f, res : %.2f", evt.run, evt.event, this_PdgID, true_P, P_from_range, P_fitted, res) << endl;
  JSFillHist("Daughter", "Daughter_MCS_Res_" + this_PdgID_str, res, 1., 400., -2., 2.);

  // -- Clean up
  for(unsigned i = 0; i < this_likelihood_hist.size(); i++){
    delete this_likelihood_hist.at(i);
    //this_likelihood_hist.at(i) = nullptr;
  }
  this_likelihood_hist.clear();

  return;
}

vector<TH1D*> MCS_Performance::Calculate_Likelihoods_for_Performance(const vector<MCSSegment> & segments, double segment_size, int PID){

  vector<TH1D*> out;
  unsigned int min_N_segment = 3;
  if(segments.size() < min_N_segment) return out;

  double this_mass = M_pion;
  if(PID == 13) this_mass = M_mu;

  double total_range = 0.;
  for(unsigned int i = 0; i < segments.size(); i++){
    double this_segment_range = segments.at(i).Range();
    total_range = total_range + this_segment_range;
  }

  double KE_from_range = map_BB[abs(PID)] -> KEFromRangeSpline(total_range);
  double P_from_range = map_BB[abs(PID)] -> KEtoMomentum(KE_from_range);
  double KE_step = 10.;
  int N_step = 200;

  vector<double> range_vec;
  range_vec.push_back(0.);
  vector<double> theta_xz_vec;
  vector<double> theta_yz_vec;
  for(unsigned int i = 0; i < segments.size() - 1; i++){
    double this_partial_range = segments.at(i).Range();
    range_vec.push_back(this_partial_range);

    TVector3 this_vec = segments.at(i).FittedVec();
    TVector3 next_vec = segments.at(i + 1).FittedVec();
    TVector3 rotated_this_vec = RotateToZaxis(this_vec, this_vec);
    TVector3 rotated_next_vec = RotateToZaxis(this_vec, next_vec);
    double this_theta_xz = TMath::ATan(rotated_next_vec.X() / rotated_next_vec.Z());
    double this_theta_yz = TMath::ATan(rotated_next_vec.Y() / rotated_next_vec.Z());
    theta_xz_vec.push_back(this_theta_xz);
    theta_yz_vec.push_back(this_theta_yz);
 
    TString this_segment_str = Form("segment%d", i);
    TH1D* h_check = (TH1D*)gROOT -> FindObject(this_segment_str);
    delete h_check;
    TH1D* this_hist = new TH1D(this_segment_str, this_segment_str, 1500., 0., 3000.);
    out.push_back(this_hist);
  }

  for(int i = 0; i < N_step; i++){
    double this_start_KE = KE_from_range + KE_step * (i + 0.);
    double this_start_P = map_BB[abs(PID)] -> KEtoMomentum(this_start_KE);
    double this_KE = this_start_KE;
    //cout << "[MCS_Performance::Calculate_Likelihoods_for_Performance] " << i << ", this_KE : " << this_KE << endl;
    //double this_start_P = P_from_range + P_step * (i + 0.);
    for(unsigned int j = 0; j < out.size(); j++){
      this_KE = map_BB[abs(PID)] -> KEAtLength(this_KE, range_vec.at(j));
      double this_P = map_BB[abs(PID)] -> KEtoMomentum(this_KE);
      double this_HL_sigma = MCS_Get_HL_Sigma(segment_size, this_P, this_mass);
      double N_sigma_xz = theta_xz_vec.at(j) / this_HL_sigma;
      double N_sigma_yz = theta_yz_vec.at(j) / this_HL_sigma;
      //double this_xz_likelihood = MCS_Get_Likelihood(this_HL_sigma, theta_xz_vec.at(j));
      //double this_yz_likelihood = MCS_Get_Likelihood(this_HL_sigma, theta_yz_vec.at(j));
      double this_xz_likelihood = MCS_Get_Likelihood_tail_model(this_P, this_HL_sigma, theta_xz_vec.at(j), (int) segment_size);
      double this_yz_likelihood = MCS_Get_Likelihood_tail_model(this_P, this_HL_sigma, theta_yz_vec.at(j), (int) segment_size);
      double this_likelihood = 0.;
      if(N_sigma_xz < 3.) this_likelihood = this_likelihood + this_xz_likelihood;
      if(N_sigma_yz < 3.) this_likelihood = this_likelihood + this_yz_likelihood;
      //cout << "[MCS_Performance::Calculate_Likelihoods_for_Performance] P : " << this_P << ", xz : " << theta_xz_vec.at(j) << ", yz : " << theta_yz_vec.at(j) << ", HL : " << this_HL_sigma << endl;
      int this_bin_number = out.at(j) -> GetXaxis() -> FindBin(this_start_P);
      
      out.at(j) -> SetBinContent(this_bin_number, (-1.) * this_likelihood + 10000.);
    }
  }

  return out;

}

double MCS_Performance::MCS_Get_Likelihood_tail_model(double this_P, double HL_sigma, double delta_angle, int segment_size){
  double out = TMath::Log(HL_sigma) + 0.5 * pow(delta_angle / HL_sigma, 2.0);

  double ratio_sigma = 2.9;
  double ratio_p0 = 0.0616;
  double ratio_p1 = 0.217;
  double ratio_p2 = 0.00421;
  double ratio_p3 = 252.;
  if(segment_size == 14){
    ratio_sigma = 2.9;
    ratio_p0 = 0.0616;
    ratio_p1 = 0.217;
    ratio_p2 = 0.00421;
    ratio_p3 = 252.;
  }

  double ratio_constant = (ratio_p0 + ratio_p1 * exp(-1. * ratio_p2 * (this_P - ratio_p3)));
  out = 2.0 * TMath::Log(HL_sigma * (1. + ratio_constant * ratio_sigma)) - TMath::Log(ratio_constant) + 0.5 * pow(delta_angle / HL_sigma, 2.0) + 0.5 * pow(delta_angle / (ratio_sigma * HL_sigma), 2.0);
  //cout << "[MCS_Performance::MCS_Get_Likelihood_tail_model] this_P : " << this_P << ", HL_sigma : " << HL_sigma << ", delta_angle : " << delta_angle << ", segment_size : " << segment_size << ", out : " << out << endl;

  return out;
}


TString MCS_Performance::Get_N_segment_str(int N_segment){
  TString out = "";
  if(N_segment >= 3 && N_segment <= 5) out = "Nseg3to5";
  else if(N_segment >= 6 && N_segment <= 9) out = "Nseg6to9";
  else if(N_segment >= 10 && N_segment <= 15) out= "Nseg10to15";
  else if(N_segment >= 16 && N_segment <= 30) out= "Nseg16to30";
  else if(N_segment > 30) out= "NsegOver30";

  return out;
}

MCS_Performance::MCS_Performance(){

}

MCS_Performance::~MCS_Performance(){

}
