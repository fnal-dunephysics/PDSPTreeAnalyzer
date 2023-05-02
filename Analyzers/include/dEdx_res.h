#ifndef dEdx_res_h
#define dEdx_res_h

#include "AnalyzerCore.h"

class dEdx_res : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_Beam(int PID);
  double dEdx_scaled(int PID, double MC_dEdx);
  double dEdx_smeared(int PID, double MC_dEdx);

  dEdx_res();
  ~dEdx_res();

};

#endif
