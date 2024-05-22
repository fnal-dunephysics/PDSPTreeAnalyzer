#include "PionQE.h"

void PionQE::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void PionQE::executeEvent(){

  //cout << "test evt.beam_inst_P : " << evt.beam_inst_P * 1000. << endl;
  
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
  if(pions.size() > 0) Run_Daughter(pions);
}

void PionQE::Run_beam_dEdx_vector(){

  int total_N_hits = (*evt.reco_beam_calibrated_dEdX_SCE).size();
  if(total_N_hits < 15) return;
  
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
  while(this_N_hits > 15 && this_KE_BB > 0.1){
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

}

void PionQE::Run_Daughter(const vector<Daughter>& pions){

  for(unsigned int i_pion = 0; i_pion < pions.size(); i_pion++){
    Daughter this_daughter = pions.at(i_pion);
    double this_chi2_pion = Particle_chi2(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE(), 211);
    JSFillHist("Daughter_pion", "Daughter_pion_chi2_pion", this_chi2_pion, 1., 1000., 0., 1000.);
    //cout << "this_daughter.PFP_true_byHits_startE() : " << this_daughter.PFP_true_byHits_startE() << endl;
    JSFillHist("Daughter_pion", "Daughter_pion_true_start_KE_vs_chi2_pion", this_daughter.PFP_true_byHits_startE() * 1000. - M_pion, this_chi2_pion, 1., 400., 0., 2000., 100., 0., 100.); 

    // == Select stopped pion
    if(this_chi2_pion < 6.){
      FitWithVectors(this_daughter.allTrack_calibrated_dEdX_SCE(), this_daughter.allTrack_resRange_SCE());
    }
  }

}

void PionQE::FitWithVectors(const vector<double>& dEdx, const vector<double>& range){

  int total_N_hits = dEdx.size();
  if(total_N_hits < 15) return;

  int this_N_hits = total_N_hits;
  int skip_N_hits = 0;

  //cout << "range.at(0) : " << range.at(0) << ", range.at(total_N_hits - 1) : " << range.at(total_N_hits - 1) << endl;
  double this_KE_BB = map_BB[211] -> KEFromRangeSpline(range.at(total_N_hits - 1));

  while(this_N_hits > 15 && this_KE_BB > 0.1){
    vector<double> this_dEdx_vec;
    vector<double> this_range_vec;
    this_N_hits = 0;
    for(int i_hit = skip_N_hits; i_hit < total_N_hits; i_hit++){
      this_dEdx_vec.push_back(dEdx.at(i_hit));
      this_range_vec.push_back(range.at(i_hit));
      this_N_hits++;
    }

    //double this_alt_length = range.at(total_N_hits - 1) - range.at(skip_N_hits);
    //cout << "this_alt_length : " << this_alt_length << endl;
    //this_KE_BB = map_BB[211]->KEAtLength(this_KE_BB, this_alt_length);

    double Length_HypFit_Gaussian = Fit_HypTrkLength_Gaussian(this_dEdx_vec, this_range_vec, 211, false, true);
    double Length_HypFit_Likelihood = Fit_HypTrkLength_Likelihood(this_dEdx_vec, this_range_vec, 211, false, true);
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

    JSFillHist("Denom", "KE_beam_" + N_hits_str, this_KE_BB, 1., 1500., 0., 1500.);
    if(Length_HypFit_Likelihood > 0.){
      JSFillHist("Gaussian", "Gaussian_KE_fit_vs_KE_BB_" + N_hits_str, KE_HypFit_Gaussian, this_KE_BB, 1., 300., 0., 1500., 300., 0., 1500.);
    }
    if(Length_HypFit_Likelihood > 0.){
      JSFillHist("Likelihood", "Likelihood_KE_fit_vs_KE_BB_" + N_hits_str, KE_HypFit_Likelihood, this_KE_BB, 1., 300., 0., 1500., 300., 0., 1500.);
    }

    this_dEdx_vec.clear();
    this_range_vec.clear();
    skip_N_hits++;
    this_N_hits = total_N_hits - skip_N_hits;
  }
}

PionQE::PionQE(){

}

PionQE::~PionQE(){

}
