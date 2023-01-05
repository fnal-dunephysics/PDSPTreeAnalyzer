#ifndef PionKEScale_h
#define PionKEScale_h

#include "AnalyzerCore.h"

class PionKEScale : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_beam_dEdx_vector();
  void Run_Daughter(const vector<Daughter>& pions);
  void FitWithVectors(const vector<double>& dEdx, const vector<double>& range);

  PionKEScale();
  ~PionKEScale();

};

#endif
