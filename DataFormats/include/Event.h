#ifndef Event_h
#define Event_h

#include "TString.h"
#include "TObject.h"

class Event : public TObject {
public:

  Event();
  ~Event();

  void SetSimulator(TString Beam_Momentum);
  inline TString Simulator() const {return j_Beam_Momentum;}

private:
  TString j_Beam_Momentum;

  ClassDef(Event,1)
};

#endif
