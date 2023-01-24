#ifndef MCSSegment_h
#define MCSSegment_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include "TVector3.h"
#include <iostream>

using namespace std;

class MCSSegment{
  
 public:

  MCSSegment();
  ~MCSSegment();

  void SetMCSSegment(TVector3 i_reco_hit_start, TVector3 i_reco_hit_end, TVector3 i_fitted_start, TVector3 i_fitted_end, TVector3 i_fitted_vec);
  void TestClass();

  inline TVector3 RecoStart() const { return j_reco_hit_start; }
  inline TVector3 RecoEnd() const { return j_reco_hit_end; }
  inline TVector3 FittedStart() const { return j_fitted_start; }
  inline TVector3 FittedEnd() const { return j_fitted_end; }
  inline TVector3 FittedVec() const { return j_fitted_vec; }

 private:
  
  TVector3 j_reco_hit_start;
  TVector3 j_reco_hit_end;
  TVector3 j_fitted_start;
  TVector3 j_fitted_end;
  TVector3 j_fitted_vec;
};

#endif
