#include "PionKEScale.h"

void PionKEScale::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void PionKEScale::executeEvent(){

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
  //Run_beam_dEdx_vector();
  //Run_beam_MCS();

  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> pions = GetPions(daughters_all);
  if(pions.size() > 0){
    Run_Daughter_HypFit(pions);
    //Run_Daughter_MCS(pions);
  }

  return;
}

void PionKEScale::Run_beam_dEdx_vector(){

  int total_N_hits = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  if(total_N_hits < 11) return;
  
  // == Checked that 0th index is first hit with the largest resRange
  /*
  cout << "P_beam_inst : " << P_beam_inst << "\ti\tdEdx\tRange" << endl;
  for(int i_hit = 0; i_hit < total_N_hits; i_hit++){
    cout << i_hit << "\t" << (*evt.reco_beam_calibrated_dEdX_SCE).at(i_hit) << "\t" << (*evt.reco_beam_resRange_SCE).at(i_hit) << endl;
  }
  */

  int this_N_hits = total_N_hits;
  int skip_N_hits = 0;
  double this_KE_BB = KE_ff_reco;
  while(this_N_hits > 10 && this_KE_BB > 0.1){
    vector<double> this_dEdx_vec;
    vector<double> this_range_vec;
    this_N_hits = 0;
    for(int i_hit = skip_N_hits; i_hit < total_N_hits; i_hit++){
      this_dEdx_vec.push_back((*evt.reco_beam_calibrated_dEdX_SCE).at(i_hit));
      this_range_vec.push_back((*evt.reco_beam_resRange_SCE).at(i_hit));
      this_N_hits++;
    }

    double this_alt_length = (*evt.reco_beam_resRange_SCE).at(0) - (*evt.reco_beam_resRange_SCE).at(skip_N_hits);
    //cout << "this_alt_length : " << this_alt_length << endl;
    this_KE_BB = map_BB[211]->KEAtLength(KE_ff_reco, this_alt_length);
    
    double Length_HypFit_Gaussian = Fit_HypTrkLength_Gaussian(this_dEdx_vec, this_range_vec, 211, false, true);
    double Length_HypFit_Likelihood = Fit_HypTrkLength_Likelihood(this_dEdx_vec, this_range_vec, 211, false, true);
    double KE_HypFit_Gaussian = -9999.;
    double KE_HypFit_Likelihood = -9999.;
    //cout << "this_KE_BB : " << this_KE_BB << "\tthis_N_hits\t" << this_N_hits << endl;
    if(Length_HypFit_Gaussian > 0.){
      KE_HypFit_Gaussian = map_BB[211] -> KEFromRangeSpline(Length_HypFit_Gaussian);
      //double this_res = (this_KE_BB - KE_HypFit_Gaussian) / KE_HypFit_Gaussian;
      //cout << "KE_HypFit_Gaussian : " << KE_HypFit_Gaussian << "\tRes\t" << this_res << endl;
    }
    if(Length_HypFit_Likelihood > 0.){
      KE_HypFit_Likelihood = map_BB[211] -> KEFromRangeSpline(Length_HypFit_Likelihood);
      //double this_res= (this_KE_BB - KE_HypFit_Likelihood) / KE_HypFit_Likelihood;
      //cout << "KE_HypFit_Likelihood : " << KE_HypFit_Likelihood << "\tRes\t" << this_res <<endl;
    }

    int N_hits_step = 30;
    int N_hits_down = N_hits_step * (this_N_hits / N_hits_step);
    int N_hits_up = N_hits_step * (this_N_hits / N_hits_step + 1);
    TString N_hits_str = Form("Nhits%dto%d", N_hits_down, N_hits_up);

    JSFillHist("Denom", "KE_beam_" + N_hits_str, this_KE_BB, 1., 1500., 0., 1500.);
    if(Length_HypFit_Likelihood > 0.){
      JSFillHist("Gaussian", "Gaussian_KE_fit_vs_KE_beam_" + N_hits_str, KE_HypFit_Gaussian, this_KE_BB, 1., 300., 0., 1500., 300., 0., 1500.);
    }
    if(Length_HypFit_Likelihood > 0.){
      JSFillHist("Likelihood", "Likelihood_KE_fit_vs_KE_beam_" + N_hits_str, KE_HypFit_Likelihood, this_KE_BB, 1., 300., 0., 1500., 300., 0., 1500.);
    }

    this_dEdx_vec.clear();
    this_range_vec.clear();
    skip_N_hits++;
    this_N_hits = total_N_hits - skip_N_hits;
  }

  return;
}

void PionKEScale::Run_beam_MCS(){

  if(pi_type != 1 && pi_type != 3) return;

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
  vector<double> reco_range = (*evt.reco_beam_resRange_SCE);
  vector<TVector3> reco_position_vec;
  for(unsigned int i_reco_hit = 0; i_reco_hit < reco_Z.size(); i_reco_hit++){
    TVector3 this_position(reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit));
    reco_position_vec.push_back(this_position);
    //cout << Form("[PionKEScale::Run_beam_MCS] reco %d : (%.2f, %.2f, %.2f)", i_reco_hit, reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit)) << endl;
  }


  // == Test 3D linear fit function
  TString evt_run_str = Form("Run%d_Evt%d", evt.run, evt.event);
  //TVector3 test_fit = Fitter->line3Dfit(reco_position_vec, true, evt_run_str);
  //cout << "[PionKEScale::Run_beam_MCS] " << evt_run_str << ", test_fit.Mag() : " << test_fit.Mag() << endl;


  MCS_Plot_Angles_for_Segment_Size_True(reco_position_vec, true_position_vec, true_P_vec, 4., "4cm");
  MCS_Plot_Angles_for_Segment_Size_True(reco_position_vec, true_position_vec, true_P_vec, 5., "5cm");
  MCS_Plot_Angles_for_Segment_Size_True(reco_position_vec, true_position_vec, true_P_vec, 8., "8cm");
  MCS_Plot_Angles_for_Segment_Size_True(reco_position_vec, true_position_vec, true_P_vec, 10., "10cm");
  MCS_Plot_Angles_for_Segment_Size_True(reco_position_vec, true_position_vec, true_P_vec, 14., "14cm");

  true_position_vec.clear();
  true_Pvec_vec.clear();
  true_P_vec.clear();
  reco_X.clear();
  reco_Y.clear();
  reco_Z.clear();
  reco_range.clear();
  reco_position_vec.clear();

  return;
}

void PionKEScale::MCS_Plot_Angles_for_Segment_Size_True(const vector<TVector3> & reco_position_vec, const vector<TVector3> & true_position_vec, const vector<double> & true_P_vec, double segment_size, TString name){

  int this_PDG = 211;
  if(pi_type == 3) this_PDG = 13;
  vector<MCSSegment> this_segments = SplitIntoSegments(reco_position_vec, segment_size);
  if(this_segments.size() < 2) return;
  vector<double> segment_true_P = GetSegmentTrueP(this_segments, true_position_vec, true_P_vec, this_PDG);
  for(unsigned int i_seg = 0; i_seg < segment_true_P.size(); i_seg++){
    if(segment_true_P.at(i_seg) < 0.) return;
  }

  for(unsigned int i_seg = 0; i_seg < segment_true_P.size() - 1; i_seg++){
    TVector3 this_vec = this_segments.at(i_seg).FittedVec();
    TVector3 next_vec = this_segments.at(i_seg + 1).FittedVec();
    TVector3 rotated_this_vec = RotateToZaxis(this_vec, this_vec);
    TVector3 rotated_next_vec = RotateToZaxis(this_vec, next_vec);

    double theta_yz = TMath::ATan(rotated_next_vec.Y() / rotated_next_vec.Z());
    double theta_xz = TMath::ATan(rotated_next_vec.X() / rotated_next_vec.Z());

    double theta_3D = rotated_next_vec.Theta();

    JSFillHist("Beam_MCS", "Beam_MCS_true_P_vs_theta_yz_" + name + "_" + pi_type_str, segment_true_P.at(i_seg), theta_yz, 1., 3000., 0., 3000., 1000., -0.5, 0.5);
    JSFillHist("Beam_MCS", "Beam_MCS_true_P_vs_theta_xz_" + name + "_" + pi_type_str, segment_true_P.at(i_seg), theta_xz, 1., 3000., 0., 3000., 1000., -0.5, 0.5);
    JSFillHist("Beam_MCS", "Beam_MCS_true_P_vs_theta_3D_" + name + "_" + pi_type_str, segment_true_P.at(i_seg), theta_3D, 1., 3000., 0., 3000., 1000., -0.5, 0.5);
  }

  this_segments.clear();
  segment_true_P.clear();
  return;
}

void PionKEScale::Run_Daughter_HypFit(const vector<Daughter>& pions){

  for(unsigned int i_pion = 0; i_pion < pions.size(); i_pion++){

    Daughter this_daughter = pions.at(i_pion);

    double this_chi2_pion = Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), 211);

    TString particle_str = "Data";
    double true_KE = -9999.;
    double this_KE_BB = map_BB[211] -> KEFromRangeSpline(this_daughter.allTrack_resRange_SCE().back());
    if(evt.MC){
      int this_PdgID = this_daughter.PFP_true_byHits_PDG();
      if(this_PdgID == 2212) particle_str = "proton";
      else if(abs(this_PdgID) == 211) particle_str = "pion";
      else if(abs(this_PdgID) == 13) particle_str = "muon";
      else particle_str = "other";

      JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_chi2_pion_" + particle_str, this_daughter.PFP_true_byHits_startE() * 1000. - M_pion, this_chi2_pion, 1., 400., 0., 2000., 100.,0., 100.);
      double true_P = this_daughter.PFP_true_byHits_startP() * 1000.;
      double true_E = this_daughter.PFP_true_byHits_startE() * 1000.;
      true_KE = true_E - pow(true_E * true_E - true_P * true_P , 0.5);
      JSFillHist("Daughter_pion", "Daughter_pion_true_start_KEv2_vs_chi2_pion_" + particle_str, true_KE, this_chi2_pion, 1., 400., 0., 2000., 100.,0., 100.);
      JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_KE_BB_" + particle_str, true_KE, this_KE_BB, 1., 400., 0., 2000., 400., 0., 2000.);

      double this_range = this_daughter.allTrack_resRange_SCE().back();
      JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_Range_" + particle_str, true_KE, this_range, 1., 400., 0., 2000., 200., 0., 1000.);
    }

    JSFillHist("Daughter_pion", "Daughter_pion_chi2_pion_" + particle_str, this_chi2_pion, 1., 1000., 0., 1000.);

    // == Select stopped pion
    if(this_chi2_pion < 6.){
      if(evt.MC){
	JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_KE_BB_chi2_pion_cut_" + particle_str, true_KE, this_KE_BB, 1., 400., 0., 2000., 400., 0., 2000.);
	double this_range_res = (this_KE_BB - true_KE) / true_KE;
	double this_ragne_invres = (1./this_KE_BB - 1./true_KE) / (1./true_KE);
	JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_KE_BB_Res_chi2_pion_cut_" + particle_str, true_KE, this_range_res, 1., 400., 0., 2000., 400., -2., 2.);
        JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_KE_BB_InvRes_chi2_pion_cut_" + particle_str, true_KE, this_ragne_invres, 1., 400., 0., 2000., 400., -2., 2.);
      }
      FitWithVectors(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), particle_str, true_KE);
    }
  }

  return;
}

void PionKEScale::FitWithVectors(const vector<double>& dEdx, const vector<double>& range, TString particle_str, double true_KE){

  int total_N_hits = dEdx.size();
  if(total_N_hits < 11) return;

  int this_N_hits = total_N_hits;
  int skip_N_hits = 0;

  //cout << "range.at(0) : " << range.at(0) << ", range.at(total_N_hits - 1) : " << range.at(total_N_hits - 1) << endl;
  double this_KE_BB = map_BB[211] -> KEFromRangeSpline(range.at(this_N_hits - 1));
  //cout << "[PionKEScale::FitWithVectors] range.at(this_N_hits - 1) : " << range.at(this_N_hits - 1) << ", this_KE_BB : " << this_KE_BB << endl;

  while(this_N_hits > 10 && this_KE_BB > 0.1){
    vector<double> this_dEdx_vec;
    vector<double> this_range_vec;
    this_N_hits = 0;
    for(int i_hit = skip_N_hits; i_hit < total_N_hits; i_hit++){
      this_dEdx_vec.push_back(dEdx.at(i_hit));
      this_range_vec.push_back(range.at(i_hit) - range.at(skip_N_hits));
      this_N_hits++;
    }

    //double this_alt_length = range.at(total_N_hits - 1) - range.at(skip_N_hits);
    //cout << "this_alt_length : " << this_alt_length << endl;
    //this_KE_BB = map_BB[211]->KEAtLength(this_KE_BB, this_alt_length);

    double Length_HypFit_Gaussian = Fit_HypTrkLength_Gaussian(this_dEdx_vec, this_range_vec, 211, false, false);
    double Length_HypFit_Likelihood = Fit_HypTrkLength_Likelihood(this_dEdx_vec, this_range_vec, 211, false, false);
    double KE_HypFit_Gaussian = -9999.;
    double KE_HypFit_Likelihood = -9999.;
    //cout << "this_KE_BB : " << this_KE_BB << "\tthis_N_hits\t" << this_N_hits << endl;
    if(Length_HypFit_Gaussian > 0.){
      KE_HypFit_Gaussian = map_BB[211] -> KEFromRangeSpline(Length_HypFit_Gaussian);
    }
    if(Length_HypFit_Likelihood > 0.){
      KE_HypFit_Likelihood = map_BB[211] -> KEFromRangeSpline(Length_HypFit_Likelihood);
    }

    int N_hits_step = 30;
    int N_hits_down = N_hits_step * (this_N_hits / N_hits_step);
    int N_hits_up = N_hits_step * (this_N_hits / N_hits_step + 1);
    TString N_hits_str = Form("Nhits%dto%d", N_hits_down, N_hits_up);

    JSFillHist("Denom", "KE_beam_" + N_hits_str + "_" + particle_str, this_KE_BB, 1., 1500., 0., 1500.);
    if(Length_HypFit_Likelihood > 0.){
      JSFillHist("Gaussian", "Gaussian_KE_fit_vs_KE_BB_" + N_hits_str + "_" + particle_str, KE_HypFit_Gaussian, this_KE_BB, 1., 300., 0., 1500., 300., 0., 1500.);
      double this_Gaus_fitted_res = (KE_HypFit_Gaussian - this_KE_BB) / this_KE_BB;
      double this_Gaus_fitted_invres = (1./KE_HypFit_Gaussian - 1./this_KE_BB) / (1./this_KE_BB);
      JSFillHist("Gaussian", "Gaussian_KE_BB_vs_KE_fit_Res_" + N_hits_str + "_" + particle_str, this_KE_BB, this_Gaus_fitted_res, 1., 300., 0., 1500., 400., -2., 2.);
      JSFillHist("Gaussian", "Gaussian_KE_BB_vs_KE_fit_InvRes_" + N_hits_str + "_" + particle_str,this_KE_BB, this_Gaus_fitted_invres, 1., 300., 0., 1500., 400., -2., 2.);

      if(evt.MC){
	double this_Gaus_true_res = (KE_HypFit_Gaussian - true_KE) / true_KE;
	double this_Gaus_true_invres = (1./KE_HypFit_Gaussian - 1./true_KE) / (1./true_KE);
	JSFillHist("Gaussian", "Gaussian_KE_true_vs_KE_fit_Res_" + N_hits_str + "_" + particle_str, true_KE, this_Gaus_true_res, 1., 300., 0., 1500., 400., -2., 2.);
	JSFillHist("Gaussian", "Gaussian_KE_true_vs_KE_fit_InvRes_" + N_hits_str + "_" + particle_str, true_KE, this_Gaus_true_invres, 1., 300., 0., 1500., 400., -2., 2.);
      }
    }
    if(Length_HypFit_Likelihood > 0.){
      JSFillHist("Likelihood", "Likelihood_KE_fit_vs_KE_BB_" + N_hits_str + "_" + particle_str, KE_HypFit_Likelihood, this_KE_BB, 1., 300., 0., 1500., 300., 0., 1500.);
      double this_Likelihood_fitted_res = (KE_HypFit_Likelihood - this_KE_BB) /this_KE_BB;
      double this_Likelihood_fitted_invres = (1./KE_HypFit_Likelihood -1./this_KE_BB) / (1./this_KE_BB);
      JSFillHist("Likelihood", "Likelihood_KE_BB_vs_KE_fit_Res_" + N_hits_str + "_" + particle_str,this_KE_BB, this_Likelihood_fitted_res, 1., 300., 0., 1500., 400., -2., 2.);
      JSFillHist("Likelihood", "Likelihood_KE_BB_vs_KE_fit_InvRes_" + N_hits_str + "_" + particle_str,this_KE_BB, this_Likelihood_fitted_invres, 1., 300., 0., 1500., 400., -2., 2.);
      if(evt.MC){
	double this_Likelihood_true_res = (KE_HypFit_Likelihood - true_KE) / true_KE;
	double this_Likelihood_true_invres = (1./KE_HypFit_Likelihood - 1./true_KE) / (1./true_KE);
	JSFillHist("Likelihood", "Likelihood_KE_true_vs_KE_fit_Res_" + N_hits_str + "_" + particle_str, true_KE, this_Likelihood_true_res, 1., 300., 0., 1500., 400., -2., 2.);
	JSFillHist("Likelihood", "Likelihood_KE_true_vs_KE_fit_InvRes_" + N_hits_str + "_" + particle_str, true_KE, this_Likelihood_true_invres, 1., 300., 0., 1500., 400., -2., 2.);
      }
    }

    this_dEdx_vec.clear();
    this_range_vec.clear();
    skip_N_hits++;
    this_N_hits = total_N_hits - skip_N_hits;
  }

  return;
}

void PionKEScale::Run_Daughter_MCS(const vector<Daughter>& pions){
  
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

      //double this_distance = (beam_end - this_position).Mag();
      //cout << Form("%d (reco_X, reco_Y, reco_Z) = (%.2f, %.2f, %.2f), where %.2f from the beam end : (%.2f, %.2f, %.2f)", i_reco_hit, reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit), this_distance, beam_last_X, beam_last_Y, beam_last_Z) << endl;

      
      reco_position_vec.push_back(this_position);
    }
    double true_P = this_daughter.PFP_true_byHits_startP() * 1000.;

    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 4., "4cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 5., "5cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 8., "8cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 10., "10cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 14., "14cm");
  }

  return;
}

void PionKEScale::Run_Daughter_MCS_for_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name){

  vector<MCSSegment> this_segments = SplitIntoSegments(reco_position_vec, segment_size);
  if(this_segments.size() < 3) return;

  TString PdgID_sign = "";
  if(this_PdgID < 0){
    PdgID_sign = "m";
  }
  else if(this_PdgID > 0){
    PdgID_sign = "p";
  }
  TString this_PdgID_str = Form(PdgID_sign + "%d", abs(this_PdgID));

  TVector3 this_vec = this_segments.at(0).FittedVec();
  TVector3 next_vec = this_segments.at(1).FittedVec();
  TVector3 rotated_this_vec = RotateToZaxis(this_vec, this_vec);
  TVector3 rotated_next_vec = RotateToZaxis(this_vec, next_vec);
  /*
  cout << Form("this_vec (%.5f, %.5f, %.5f) %.2f, next_vec (%.5f, %.5f, %.5f) %.2f",
	       this_vec.X(), this_vec.Y(), this_vec.Z(), this_vec.Mag(), next_vec.X(), next_vec.Y(), next_vec.Z(), next_vec.Mag()) << endl;
  cout << Form("rotated_this_vec (%.5f, %.5f, %.5f) %.2f, rotated_next_vec (%.5f, %.5f, %.5f) %.2f",
	       rotated_this_vec.X(), rotated_this_vec.Y(), rotated_this_vec.Z(), rotated_this_vec.Mag(), rotated_next_vec.X(), rotated_next_vec.Y(), rotated_next_vec.Z(), rotated_next_vec.Mag()) << endl;
  */
  double this_vec_theta = this_vec.Theta();
  double theta_yz = TMath::ATan(rotated_next_vec.Y() / rotated_next_vec.Z());
  double theta_xz = TMath::ATan(rotated_next_vec.X() / rotated_next_vec.Z());
  double theta_3D = rotated_next_vec.Theta();

  if(this_vec_theta > 1.55 && this_vec_theta < 1.59) return;

  if(name == "14cm" && fabs(theta_xz) < 0.001){
    cout << Form("this_vec (%.5f, %.5f, %.5f) %.2f, next_vec (%.5f, %.5f, %.5f) %.2f",
		 this_vec.X(), this_vec.Y(), this_vec.Z(), this_vec.Mag(), next_vec.X(), next_vec.Y(), next_vec.Z(), next_vec.Mag()) << endl;
    cout << Form("rotated_this_vec (%.5f, %.5f, %.5f) %.2f, rotated_next_vec (%.5f, %.5f, %.5f) %.2f",
		 rotated_this_vec.X(), rotated_this_vec.Y(), rotated_this_vec.Z(), rotated_this_vec.Mag(), rotated_next_vec.X(), rotated_next_vec.Y(), rotated_next_vec.Z(), rotated_next_vec.Mag()) << endl;
  }

  //cout << "[PionKEScale::Run_Daughter_MCS_for_Segments] " << name << ", theta_yz : " << theta_yz << ", theta_xz : " << theta_xz << ", theta_3D : " << theta_3D << endl;
  JSFillHist("Daughter_MCS", "Daughter_MCS_true_P_vs_theta_yz_" + name + "_" + this_PdgID_str, true_P, theta_yz, 1., 3000., 0., 3000., 1000., -0.5, 0.5);
  JSFillHist("Daughter_MCS", "Daughter_MCS_true_P_vs_theta_xz_" + name + "_" + this_PdgID_str, true_P, theta_xz, 1., 3000., 0., 3000., 1000., -0.5, 0.5);
  JSFillHist("Daughter_MCS", "Daughter_MCS_true_P_vs_theta_3D_" + name + "_" + this_PdgID_str, true_P, theta_3D, 1., 3000., 0., 3000., 1000., -0.5, 0.5);
  JSFillHist("Daughter_MCS", "Daughter_MCS_this_vec_theta_vs_theta_xz_" + name + "_" + this_PdgID_str, this_vec_theta, theta_xz, 1., 1000., 0., 4., 1000., -0.5, 0.5);
  JSFillHist("Daughter_MCS", "Daughter_MCS_this_vec_theta_vs_theta_yz_" + name + "_" + this_PdgID_str, this_vec_theta, theta_yz, 1., 1000., 0., 4., 1000., -0.5, 0.5);

  return;
}

vector<double> PionKEScale::GetSegmentTrueP(const vector<MCSSegment> & segments, const vector<TVector3> & true_position_vec, const vector<double> true_P_vec, int PDG){
  vector<double> out;

  if(segments.size() == 0) return out;

  //cout << "[GetSegmentTrueP] Start" << endl;

  for(unsigned int i = 0; i < segments.size(); i++){

    MCSSegment this_segment = segments.at(i);
    TVector3 this_start_hit = this_segment.RecoStart();
    double true_P = -999.;
    for(unsigned int j = 0; j < true_position_vec.size() - 1; j++){
      TVector3 this_true_hit = true_position_vec.at(j);
      TVector3 next_true_hit = true_position_vec.at(j + 1);
      if(this_start_hit.Z() > this_true_hit.Z() && this_start_hit.Z() < next_true_hit.Z()){
	double distance = (this_true_hit - this_start_hit).Mag();
	if(distance < 10.) true_P = map_BB[PDG] -> KEtoMomentum(map_BB[PDG] -> KEAtLength(map_BB[PDG] -> MomentumtoKE(true_P_vec.at(j)), distance));
	break;
      }
    }
    //cout << i << ", true_P : " << true_P << endl;
    out.push_back(true_P);
  }


  return out;
}

PionKEScale::PionKEScale(){

}

PionKEScale::~PionKEScale(){

}
