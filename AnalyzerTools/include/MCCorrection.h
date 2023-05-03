#ifndef MCCorrection_h
#define MCCorrection_h

#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TString.h"
#include "TRegexp.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

#include "TDirectoryHelper.h"
#include "TRandom3.h"

using namespace std;

class MCCorrection{

public:

  MCCorrection();
  ~MCCorrection();

  TDirectory *histDir;
  void ReadHistograms();

  bool IsData;
  void SetIsData(bool b);

  bool IgnoreNoHist;

  double MomentumReweight_SF(TString flag, double P_beam_inst, int sys);
  std::map< TString, TH1D* > map_hist_MomentumReweight;

  double dEdx_scaled(double MC_dEdx);

};

#endif
