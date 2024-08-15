#ifndef KaonTOF_h
#define KaonTOF_h

#include "AnalyzerCore.h"

class KaonTOF : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  KaonTOF();
  ~KaonTOF();

};

#endif
