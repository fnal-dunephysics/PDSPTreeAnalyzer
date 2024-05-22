#ifndef MCS_Angle_Shapes_h
#define MCS_Angle_Shapes_h

#include "AnalyzerCore.h"

class MCS_Angle_Shapes : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_beam_MCS();
  void Run_beam_MCS_for_Segments(const vector<TVector3> & reco_position_vec, const vector<TVector3> & true_position_vec, const vector<double> & true_P_vec, double segment_size, TString name);
  void Run_Daughter_MCS(const vector<Daughter>& pions);
  void Run_Daughter_MCS_for_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name);
  vector<double> GetSegmentTrueP(const vector<MCSSegment> & segments, const vector<TVector3> & true_position_vec, const vector<double> true_P_vec, int PDG);

  MCS_Angle_Shapes();
  ~MCS_Angle_Shapes();

};

#endif
