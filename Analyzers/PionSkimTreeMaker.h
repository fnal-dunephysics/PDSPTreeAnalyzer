#ifndef PionSkimTreeMaker_h
#define PionSkimTreeMaker_h

#include "AnalyzerCore.h"

class PionSkimTreeMaker : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void FillSkimTree();

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

  PionSkimTreeMaker();
  ~PionSkimTreeMaker();

  // == output tree branches
  //TTree            *fEventTree;
  int              _run;
  int              _subrun;
  int              _event;
  double           _beam_reco_KE_ff;
  double      	   _beam_reco_KE_end;
  double           _beam_reco_dir_X;
  double           _beam_reco_dir_Y;
  double           _beam_reco_dir_Z;
  double      	   _beam_true_KE_ff;
  double           _beam_true_KE_end;
  double           _beam_true_dir_X;
  double      	   _beam_true_dir_Y;
  double      	   _beam_true_dir_Z;
  int              _N_daughter_reco;
  vector<double>   _daughter_reco_dir_X;
  vector<double>   _daughter_reco_dir_Y;
  vector<double>   _daughter_reco_dir_Z;
  vector<double>   _daughter_reco_KE_range_muon;
  vector<double>   _daughter_reco_KE_range_pion;
  vector<double>   _daughter_reco_KE_range_proton;
  //vector<double>   _daughter_reco_KE_hypfit_pion_gaus;
  vector<double>   _daughter_reco_KE_hypfit_pion_likelihood;
  vector<double>   _daughter_reco_KE_hypfit_proton_gaus;
  //vector<double>   _daughter_reco_KE_hypfit_proton_likelihood;
  vector<double>   _daughter_reco_cos_beam;
  vector<double>   _daughter_reco_chi2_muon;
  vector<double>   _daughter_reco_chi2_pion;
  vector<double>   _daughter_reco_chi2_proton;
  vector<double>   _daughter_reco_dist_beam;
  vector<double>   _daughter_reco_trackscore;
  vector<double>   _daughter_reco_startZ;
  vector<int>      _daughter_reco_true_PDG;
  vector<int>      _daughter_reco_true_ID;
  vector<double>   _daughter_reco_true_dir_X;
  vector<double>   _daughter_reco_true_dir_Y;
  vector<double>   _daughter_reco_true_dir_Z;
  vector<double>   _daughter_reco_true_KE;
  vector<double>   _daughter_reco_true_P;
  vector<double>   _daughter_reco_true_mass;
  vector<double>   _daughter_reco_true_purity;
  vector<double>   _daughter_reco_true_completeness;
  vector<double>   _daughter_reco_true_cos_true_beam;
  vector<string>   _daughter_reco_true_endProcess;
  int              _N_daughter_all_true;
  vector<int>      _daughter_all_true_PDG;
  vector<int>      _daughter_all_true_ID;
  vector<double>   _daughter_all_true_startX;
  vector<double>   _daughter_all_true_startY;
  vector<double>   _daughter_all_true_startZ;
  vector<double>   _daughter_all_true_startPx;
  vector<double>   _daughter_all_true_startPy;
  vector<double>   _daughter_all_true_startPz;
  vector<double>   _daughter_all_true_startP;
  vector<double>   _daughter_all_true_endX;
  vector<double>   _daughter_all_true_endY;
  vector<double>   _daughter_all_true_endZ;
  vector<string>   _daughter_all_true_Process;
  vector<string>   _daughter_all_true_endProcess;
};

#endif
