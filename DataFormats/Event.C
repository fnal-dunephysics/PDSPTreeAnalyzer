#include "Event.h"

Event::Event(){
  j_Beam_Momentum = -1.;
}

Event::~Event(){

}

void Event::SetBeam_Momentum(double Beam_Momentum){
  j_Beam_Momentum = Beam_Momentum;
}
