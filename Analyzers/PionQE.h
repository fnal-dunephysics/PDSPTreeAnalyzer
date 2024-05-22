#ifndef PionQE_h
#define PionQE_h

#include "AnalyzerCore.h"

class PionQE : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_beam_dEdx_vector();
  void Run_Daughter(const vector<Daughter>& pions);
  void FitWithVectors(const vector<double>& dEdx, const vector<double>& range);

  PionQE();
  ~PionQE();

};

#endif
