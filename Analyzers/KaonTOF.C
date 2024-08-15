#include "KaonTOF.h"

void KaonTOF::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

}

void KaonTOF::executeEvent(){

  //cout << "test evt.beam_inst_P : " << evt.beam_inst_P * 1000. << endl;
  
  if(!PassBeamScraperCut()) return;
  if(!PassBeamMomentumWindowCut()) return;

  int this_true_beam_PDG = evt.true_beam_PDG;
  JSFillHist("true_beam", "true_beam_PDGID", this_true_beam_PDG, 1., 10000., -4999.5, 5000.5);

  if(!PassPandoraSliceCut()) return;
  if(!PassCaloSizeCut()) return;
  if(!PassAPA3Cut()) return;
  if(!PassMichelScoreCut()) return;
  if(!PassBeamCosCut()) return;
  if(!PassBeamStartZCut()) return;
  if(!PassProtonVetoCut()) return;
  if(!PassMuonVetoCut()) return;
  if(!PassStoppedPionVetoCut()) return;

  P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  mass_beam = 139.57;
  P_ff_reco = Convert_P_Spectrometer_to_P_ff(P_beam_inst, "pion", "AllTrue", 0);
  KE_ff_reco = sqrt(pow(P_ff_reco, 2) + pow(mass_beam, 2)) - mass_beam;
  KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, evt.reco_beam_alt_len);
  E_end_reco = KE_end_reco + mass_beam;

  JSFillHist("test", "htrack_P_beam", P_beam_inst, 1., 5000., 0., 5000.);

}

KaonTOF::KaonTOF(){

}

KaonTOF::~KaonTOF(){

}
