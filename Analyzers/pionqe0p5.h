#ifndef pionqe0p5_h
#define pionqe0p5_h

#include "AnalyzerCore.h"

class pionqe0p5 : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  //==================
  // QE variables
  //==================
  double QE_Q2 = -999;
  double QE_KEPi0 = -999;
  double QE_KEPi1 = -999;
  double QE_AngPi = -999;
  double QE_nu = -999;
  double QE_EQE = -999;
  double recoQE_Q2 = -999;
  double recoQE_KEPi0 = -999;
  double recoQE_KEPi1 = -999;
  double recoQE_AngPi = -999;
  double recoQE_nu = -999;
  double recoQE_EQE = -999;

  //==================
  // beam sel
  //==================
  bool Pass_BeamStartZ(double N_sigma = 2.);
  double beam_start_z_mu_mc = 0.143;
  double beam_start_z_sigma_mc = 0.219;
  double beam_start_z_mu_data = 3.252;
  double beam_start_z_sigma_data = 1.076;

  bool Pass_beam_delta_X_cut(double N_sigma = 2.);
  bool Pass_beam_delta_Y_cut(double N_sigma = 2.);
  double Beam_delta_X_mu_mc = 1.573;
  double Beam_delta_X_sigma_mc = 1.439;
  double Beam_delta_X_mu_data = -2.647;
  double Beam_delta_X_sigma_data = 1.802;

  double Beam_delta_Y_mu_mc = -0.644;
  double Beam_delta_Y_sigma_mc = 1.473;
  double Beam_delta_Y_mu_data =	-2.269;
  double Beam_delta_Y_sigma_data = 2.314;

  //==================
  // beam reco
  //==================
  void MuonKELoss(TString beam_selec_str, double weight);
  double GetBeamRRatZ10cm(const vector<double> & ResRange, const vector<double> & calo_Z);

  //==================
  // daughter selection
  //==================
  std::vector<Daughter> SelectLoosePions(const vector<Daughter>& in);
  void FillRecoPionPlots(TString daughter_sec_str, const vector<Daughter> pions, double weight);
  
  pionqe0p5();
  ~pionqe0p5();
};

#endif
