#include "PionKI.h"

void PionKI::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void PionKI::executeEvent(){

  if(!Pass_Beam_PID(211)) return;
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

  //cout << "[PionKI::executeEvent] Passed all beam selections" << endl;

  pi_type = GetPiParType();
  pi_type_str = Form("%d", pi_type);

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  mass_beam = 139.57;
  P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;

  JSFillHist("test", "htrack_P_beam", P_beam_inst, 1., 5000., 0., 5000.);

  // == Functions to study beam

  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> pions = GetPions(daughters_all);
  vector<Daughter> stopping_pions = GetStoppingPions(daughters_all, 6.);
  vector<Daughter> stopping_protons = GetProtons(daughters_all);
  
  JSFillHist("count", "N_2nd_stopping_pions", stopping_pions.size(), 1., 10, -0.5, 9.5);
  JSFillHist("count", "N_2nd_stopping_protons", stopping_protons.size(), 1., 10, -0.5, 9.5);

  if(stopping_pions.size() > 0 && stopping_protons.size() > 0){
    JSFillHist("count", "N_1pi_1p", 1., 1., 2., 0.5, 2.5);
  }
  return;
}

PionKI::PionKI(){

}

PionKI::~PionKI(){

}
