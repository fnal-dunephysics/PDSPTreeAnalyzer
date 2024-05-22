#ifndef PionKEScale_h
#define PionKEScale_h

#include "AnalyzerCore.h"

class PionKEScale : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_beam_dEdx_vector();
  void Run_beam_MCS();
  void MCS_Plot_Angles_for_Segment_Size_True(const vector<TVector3> & reco_position_vec, const vector<TVector3> & true_position_vec, const vector<double> & true_P_vec, double segment_size, TString name);
  void Run_Daughter_HypFit(const vector<Daughter>& pions);
  void FitWithVectors(const vector<double>& dEdx, const vector<double>& range, const vector<double>& E_field, TString particle_str, double true_KE);
  void Run_Daughter_MCS(const vector<Daughter>& pions);
  void Run_Daughter_MCS_for_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name);
  vector<double> GetSegmentTrueP(const vector<MCSSegment> & segments, const vector<TVector3> & true_position_vec, const vector<double> true_P_vec, int PDG);
  double Smear_dedx(double this_dEdx, double smear_with);

  PionKEScale();
  ~PionKEScale();

};

#endif
