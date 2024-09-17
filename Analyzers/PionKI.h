#ifndef PionKI_h
#define PionKI_h

#include "AnalyzerCore.h"

class PionKI : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void True_KI_study(TString prefix = "");
  void Study_with_daughters(double weight);
  void Beam_true_Eloss();
  void FillBeamPlots(TString beam_selec_str, double weight);
  void drawKIVar(Daughter pion, Daughter proton, double weight);

  double KE_to_P(double KE, int PID);

  // == Beam selections
  double Beam_startZ_mu_data = 3.243;
  double Beam_startZ_sigma_data = 1.283;
  double Beam_startZ_mu_mc = 0.122;
  double Beam_startZ_sigma_mc = 0.207;

  double Beam_delta_X_mu_data = -2.632;
  double Beam_delta_X_sigma_data = 0.972;
  double Beam_delta_X_mu_mc = 1.573;
  double Beam_delta_X_sigma_mc = 0.766;

  double Beam_delta_Y_mu_data = -2.312;
  double Beam_delta_Y_sigma_data = 1.547;
  double Beam_delta_Y_mu_mc = -0.676;
  double Beam_delta_Y_sigma_mc = 0.796;

  double Beam_TPC_theta_mu_data = 0.317;
  double Beam_TPC_theta_sigma_data = 0.045;
  double Beam_TPC_theta_mu_mc = 0.288;
  double Beam_TPC_theta_sigma_mc = 0.038;

  double Beam_TPC_phi_mu_data = -2.223;
  double Beam_TPC_phi_sigma_data = 0.129;
  double Beam_TPC_phi_mu_mc = -2.378;
  double Beam_TPC_phi_sigma_mc = 0.130;

  bool Pass_beam_start_Z_cut(double N_sigma);
  bool Pass_beam_delta_X_cut(double N_sigma);
  bool Pass_beam_delta_Y_cut(double N_sigma);
  bool Pass_beam_TPC_theta_cut(double N_sigma);
  bool Pass_beam_TPC_phi_cut(double N_sigma);

  void FillTKIVar(TString dir, TVector3 p_vec_beam, TVector3 p_vec_f1, TVector3 p_vec_f2, double weight);
   
  PionKI();
  ~PionKI();

};

#endif
