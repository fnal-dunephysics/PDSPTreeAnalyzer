#ifndef PionQE_h
#define PionQE_h

#include "AnalyzerCore.h"

class PionQE : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  void FillBeamPlots(TString beam_selec_str, double weight);
  void Run_Daughter(TString daughter_sec_str, const vector<Daughter>& daughters);
  void TrueDaughterStudy(const vector<TrueDaughter>& true_daughters);

  std::vector<Daughter> SelectLoosePions(const vector<Daughter>& in);
  std::vector<Daughter> SelectPions_trkscore(const vector<Daughter>& in, double cut_trackScore);
  std::vector<Daughter> SelectPions_chi2(const vector<Daughter>& in, double cut_chi2_proton);
  std::vector<Daughter> SelectPions_trklen_upper(const vector<Daughter>& in, double cut_trk_len_upper);
  std::vector<Daughter> SelectPions_trklen_lower(const vector<Daughter>& in, double cut_trk_len_lower);

  PionQE();
  ~PionQE();
};

#endif
