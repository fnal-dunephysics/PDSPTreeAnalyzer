#include "Event.h"

ClassImp(Event)

Event::Event(){
  j_Beam_Momentum = "";
}

Event::~Event(){

}

void Event::SetBeam_Momentum(TString Beam_Momentum){
  j_Beam_Momentum = Beam_Momentum;
}
