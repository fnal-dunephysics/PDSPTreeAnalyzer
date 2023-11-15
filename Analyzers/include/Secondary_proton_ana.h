#ifndef Secondary_proton_ana_h
#define Secondary_proton_ana_h

#include "AnalyzerCore.h"

class Secondary_proton_ana : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void FillBeamPlots(TString beam_selec_str, double weight = 1.);
  void Beam_true_Eloss();
  void Study_Beam_Proton_Eloss(TString beam_cut_str, double weight = 1.);
  void True_Daughter_study(const vector<Daughter>& protons, const vector<Daughter>& pions, double weight = 1.);
  void Daughter_study(const vector<Daughter>& protons);

  double KE_to_P(double KE, int PID);

  // == Beam selections
  double Beam_startZ_mu_data = 3.175;
  double Beam_startZ_sigma_data = 1.192;
  double Beam_startZ_mu_mc = 0.056;
  double Beam_startZ_sigma_mc = 0.196;

  double Beam_delta_X_mu_data = -2.649;
  double Beam_delta_X_sigma_data = 1.245;
  double Beam_delta_X_mu_mc = 1.543;
  double Beam_delta_X_sigma_mc = 1.006;

  double Beam_delta_Y_mu_data = -2.305;
  double Beam_delta_Y_sigma_data = 1.773;
  double Beam_delta_Y_mu_mc = -0.695;
  double Beam_delta_Y_sigma_mc = 1.023;

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

  Secondary_proton_ana();
  ~Secondary_proton_ana();

};

#endif
