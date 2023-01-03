#ifndef PionKEScale_h
#define PionKEScale_h

#include "AnalyzerCore.h"

class PionKEScale : public AnalyzerCore {

 public:

  void initializeAnalyzer();
  void executeEvent();

  PionKEScale();
  ~PionKEScale();

};

#endif
