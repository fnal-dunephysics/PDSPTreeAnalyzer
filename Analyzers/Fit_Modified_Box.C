#include "Fit_Modified_Box.h"

void Fit_Modified_Box::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void Fit_Modified_Box::executeEvent(){

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

  //cout << "[Fit_Modified_Box::executeEvent] alpha_new : " << alpha_new << ", beta_new : " << beta_new << endl;

  //Run_Beam(2212);
  Run_Beam(13);
}

void Fit_Modified_Box::Run_Beam(int PID){

  if(PID == 2212 && !Pass_Beam_PID(PID)) return;
  if(PID == 13 && (evt.reco_beam_calo_endZ < 260. || evt.reco_beam_calo_endZ > 420. || !Pass_Beam_PID(PID) || daughter_michel_score < 0.6)) return;

  double this_chi2 = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), PID, 1.);
  if(PID == 2212 && this_chi2 > 10.) return;

  TString beam_particle = "";

  if(PID == 2212) beam_particle = "proton";
  else if(PID == 211) beam_particle = "pion";
  else if(PID == 13) beam_particle = "muon";
  else return;


  double calib_ratio_low = 0.80;
  double calib_ratio_high = 1.20;
  double calib_ratio_step = 0.01;
  int N_calib_ratio_steps = (calib_ratio_high - calib_ratio_low) / calib_ratio_step;
  //cout << "N_calib_ratio_steps : " << N_calib_ratio_steps << endl;
  for(int i_cali = 0; i_cali < N_calib_ratio_steps + 2; i_cali++){
    double this_calib_ratio = calib_ratio_low + calib_ratio_step * (i_cali + 0.);
    TString this_calib_ratio_str = Form("Calib_Ratio_%.2f", this_calib_ratio);
    int total_N_hits = (*evt.reco_beam_calibrated_dEdX_SCE).size();
    for(int i = 0; i < total_N_hits; i++){
      //JSFillHist(beam_particle, "ResRange_vs_dEdx", (*evt.reco_beam_resRange_SCE).at(i), (*evt.reco_beam_calibrated_dEdX_SCE).at(i), 1., 200., 0., 200., 1000., 0., 50.);

      double new_recom_dEdx = (*evt.reco_beam_calibrated_dEdX_SCE).at(i);
      if(IsData){
	double cal_Efield = MCCorr -> SCE_Corrected_E((*evt.reco_beam_calo_X_allTrack).at(i), (*evt.reco_beam_calo_Y_allTrack).at(i), (*evt.reco_beam_calo_Z_allTrack).at(i));
	//new_recom_dEdx = Use_Abbey_Recom_Params(Abbey_dEdx, (*evt.reco_beam_EField_SCE).at(i), 0.9488);
	new_recom_dEdx = dEdx_different_recom(new_recom_dEdx, cal_Efield, this_calib_ratio);
      }

      JSFillHist(beam_particle, "ResRange_vs_dEdx_" + this_calib_ratio_str, (*evt.reco_beam_resRange_SCE).at(i), new_recom_dEdx, 1., 200., 0., 200., 1000., 0., 50.);
    }
  }

  return;
}

void Fit_Modified_Box::Extract_MPVs(){

  outfile->cd();

  double calib_ratio_low = 0.80;
  double calib_ratio_high = 1.20;
  double calib_ratio_step = 0.01;
  int N_calib_ratio_steps = (calib_ratio_high - calib_ratio_low) / calib_ratio_step;
  double res_range_cut = 100.; // == Use hits with residual range < 100 cm

  for(int i_cali = 0; i_cali < N_calib_ratio_steps + 2; i_cali++){
    double this_calib_ratio = calib_ratio_low + calib_ratio_step * (i_cali + 0.);
    TString this_calib_ratio_str = Form("Calib_Ratio_%.2f", this_calib_ratio);
    TString this_2D_histname = "ResRange_vs_dEdx_" + this_calib_ratio_str;
    TH2D * this_hist_2D = (TH2D*)(JSmaphist_TH2D["muon"])[this_2D_histname] -> Clone();
    
    int N_binsX = this_hist_2D -> GetNbinsX();
    int N_binsY = this_hist_2D -> GetNbinsY();

    vector<double> MPVs;
    vector<double> MPV_errs;
    vector<double> Res_ranges;
    vector<double> Res_range_errs;

    for(int i = 1; i < N_binsX + 1; i++){
      TString i_str = Form("%d", i);
      double this_ResRange = this_hist_2D -> GetXaxis() -> GetBinCenter(i);
      if(this_ResRange > res_range_cut) break;
      double this_ResRange_err = 0.5 * this_hist_2D -> GetXaxis() -> GetBinWidth(i);
      TString ResRange_range_str = Form("ResRange%.1fto%.1fcm", this_ResRange - this_ResRange_err, this_ResRange + this_ResRange_err);
      TString ResRange_range_latex = Form("Residal range : %.1f - %.1f cm", this_ResRange -this_ResRange_err, this_ResRange + this_ResRange_err);
      TString this_hist_name = ResRange_range_str;

      // == Make 1D dE/dx distribution for a residual range
      TH1D * this_hist_1D = new TH1D(this_2D_histname + this_hist_name, this_2D_histname + this_hist_name, N_binsY, 0., 50.);
      for(int j = 1; j < N_binsY + 1; j++){
	double this_hist_content = this_hist_2D -> GetBinContent(i, j);
	double this_hist_error = this_hist_2D -> GetBinError(i, j);
	this_hist_1D -> SetBinContent(j, this_hist_content);
	this_hist_1D -> SetBinError(j, this_hist_error);
      }
      double this_integ = this_hist_1D -> Integral();
      if(this_integ < 50.) continue;
      
      // == Fit using Landau+Gaussian
      double max_y = this_hist_1D -> GetMaximum();
      double max_x = this_hist_1D -> GetBinCenter(this_hist_1D -> GetMaximumBin());
      Double_t fitting_range[2];
      fitting_range[0] = 0.;
      fitting_range[1] = 15.;
      Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
      sv[0] = 0.1;
      sv[1] = max_x;
      sv[2] = this_hist_1D -> Integral() * 0.05;
      sv[3] = 0.2;
      for(int j=0; j<4; ++j){
	pllo[j] = 0.01*sv[j];
	plhi[j] = 100*sv[j];
      }

      Double_t chisqr;
      Int_t    ndf;
      Int_t    status;

      TF1 *this_Langau_fit = langaufit(this_hist_1D, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_" + this_2D_histname + this_hist_name);
      double this_MPV = this_Langau_fit -> GetParameter(1);
      double this_MPV_err = this_Langau_fit -> GetParError(1);

      MPVs.push_back(this_MPV);
      MPV_errs.push_back(this_MPV_err);
      Res_ranges.push_back(this_ResRange);
      Res_range_errs.push_back(this_ResRange_err);
    }

    cout << "Res_vs_MPV_" + this_calib_ratio_str << ", MPVs.size() : " << MPVs.size() << endl;
    TGraphErrors *this_gr = new TGraphErrors(MPVs.size(), &Res_ranges[0], &MPVs[0], &Res_range_errs[0], &MPV_errs[0]);
    this_gr -> SetName("Res_vs_MPV_" + this_calib_ratio_str);
    this_gr -> Write();
  }
}

double Fit_Modified_Box::dEdx_different_recom(double dEdx, double Efield, double calib_const_ratio){

  double alpha_default = 0.93;
  double beta_default = 0.212;
  double rho = 1.39;

  double exp_term = exp( calib_const_ratio * (beta_new / beta_default) * log(alpha_default + beta_default * dEdx / (rho * Efield)) );
  double new_dEdx = (exp_term - alpha_new) * rho * Efield / beta_new;

  return new_dEdx;
}

Fit_Modified_Box::Fit_Modified_Box(){

}

Fit_Modified_Box::~Fit_Modified_Box(){

}
