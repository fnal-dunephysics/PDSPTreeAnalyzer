#ifndef MCS_Performance_h
#define MCS_Performance_h

#include "AnalyzerCore.h"

class MCS_Performance : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_Beam();
  void Run_Beam_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name);
  void Run_Daughter(const vector<Daughter>& pions);
  void Run_Daughter_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name);

  vector<TH1D*> Calculate_Likelihoods_for_Performance(const vector<MCSSegment> & segments, double segment_size, int PID);
  double MCS_Get_Likelihood_tail_model(double this_P, double HL_sigma, double delta_angle, int segment_size);
  TString Get_N_segment_str(int N_segment);

  MCS_Performance();
  ~MCS_Performance();

};

#endif
