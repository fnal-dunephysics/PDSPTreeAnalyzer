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
  void Run_Daughter(const vector<Daughter>& pions);
  void FitWithVectors(const vector<double>& dEdx, const vector<double>& range, TString particle_str);
  void MCSAnglesWithVectors();
  vector<double> GetSegmentTrueP(const vector<MCSSegment> & segments, const vector<TVector3> & true_position_vec, const vector<double> true_P_vec, int PDG);

  PionKEScale();
  ~PionKEScale();

};

#endif
