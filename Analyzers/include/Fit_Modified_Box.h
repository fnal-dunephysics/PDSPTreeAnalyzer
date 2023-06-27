#ifndef Fit_Modified_Box_h
#define Fit_Modified_Box_h

#include "AnalyzerCore.h"

class Fit_Modified_Box : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void Run_Beam(int PID);

  double alpha_new = 0.93;
  double beta_new = 0.212;
  double dEdx_different_recom(double dEdx, double Efield, double calib_const);

  void Extract_MPVs();
  
  Fit_Modified_Box();
  ~Fit_Modified_Box();

};

#endif
