#ifndef Event_h
#define Event_h

#include "TString.h"
#include "TObject.h"

class Event : public TObject {
public:

  Event();
  ~Event();

  void SetBeam_Momentum(TString Beam_Momentum);
  inline TString Beam_Momentum() const {return j_Beam_Momentum;}

private:
  TString j_Beam_Momentum;

  ClassDef(Event,1)
};

#endif
