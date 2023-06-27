#ifndef MCS_Tree_Maker_h
#define MCS_Tree_Maker_h

#include "AnalyzerCore.h"

class MCS_Tree_Maker : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_beam_MCS();
  void Run_beam_MCS_for_Segments(const vector<TVector3> & reco_position_vec, const vector<TVector3> & true_position_vec, const vector<double> & true_P_vec, double segment_size, TString name);
  void Run_Daughter_MCS(const vector<Daughter>& pions);
  void Run_Daughter_MCS_for_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name);
  vector<double> GetSegmentTrueP(const vector<MCSSegment> & segments, const vector<TVector3> & true_position_vec, const vector<double> true_P_vec, int PDG);

  MCS_Tree_Maker();
  ~MCS_Tree_Maker();

  // == For output Tree
  TTree *MCS_tree;
  vector<int> PID;
  vector<int> Segment_size;
  vector<vector<double>> P_true;
  vector<vector<double>> Theta_xz;
  vector<vector<double>> Theta_yz;
  vector<vector<double>> Theta_3D;
  void Write_Tree();
};

#endif
