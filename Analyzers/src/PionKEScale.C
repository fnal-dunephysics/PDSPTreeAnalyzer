#include "PionKEScale.h"

void PionKEScale::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void PionKEScale::executeEvent(){

  //cout << "test evt.beam_inst_P : " << evt.beam_inst_P * 1000. << endl;
  
  if(!PassBeamScraperCut()) return;
  if(!PassBeamMomentumWindowCut()) return;
  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;
  if(!PassAPA3Cut()) return;
  if(!PassMichelScoreCut()) return;
  if(!PassBeamCosCut()) return;
  if(!PassBeamStartZCut()) return;
  if(!PassProtonVetoCut()) return;
  if(!PassMuonVetoCut()) return;
  if(!PassStoppedPionVetoCut()) return;

  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  double mass_beam = 139.57;
  double KE_beam_inst = sqrt(pow(P_beam_inst, 2) + pow(mass_beam, 2)) - mass_beam;

  double P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  double KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  double KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  double E_end_reco = KE_end_reco + mass_beam;

  cout << "E_end_reco : " << E_end_reco << endl;
  
  JSFillHist("test", "htrack_P_beam", P_beam_inst, 1., 5000., 0., 5000.);
  
}

PionKEScale::PionKEScale(){

}

PionKEScale::~PionKEScale(){

}
