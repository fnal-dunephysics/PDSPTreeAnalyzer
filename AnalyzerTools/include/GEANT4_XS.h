#ifndef GEANT4_XS_h
#define GEANT4_XS_h

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

class GEANT4_XS{

public:

  GEANT4_XS();
  ~GEANT4_XS();

  TDirectory *histDir;
  void ReadHistograms();

  bool IsData;
  void SetIsData(bool b);

  bool IgnoreNoHist;

  std::map< TString, TGraph* > map_graph_Xsec;

};

#endif
