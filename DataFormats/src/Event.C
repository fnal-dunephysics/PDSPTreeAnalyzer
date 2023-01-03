#include "Event.h"

ClassImp(Event)

Event::Event(){
  j_Beam_Momentum = "";
}

Event::~Event(){

}

void Event::SetSimulator(TString Beam_Momentum){
  j_Beam_Momentum = Beam_Momentum;
}
