#include "PionKEScale.h"

void PionKEScale::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void PionKEScale::executeEvent(){

  cout << "test evt.beam_inst_P : " << evt.beam_inst_P * 1000. << endl;

}

PionKEScale::PionKEScale(){

}

PionKEScale::~PionKEScale(){

}
