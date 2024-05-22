#include "MCS_Tree_Maker.h"

void MCS_Tree_Maker::initializeAnalyzer(){

  cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
  debug_mode = true;
  debug_mode = false;

  outfile->cd();
  MCS_tree = new TTree("MCS_tree", "MCS_tree");
  MCS_tree -> Branch("PID", "vector<int>", &PID);
  MCS_tree -> Branch("Segment_size", "vector<int>", &Segment_size);
  MCS_tree -> Branch("P_true", "vector<vector<double>>", &P_true);
  MCS_tree -> Branch("Theta_xz", "vector<vector<double>>", &Theta_xz);
  MCS_tree -> Branch("Theta_yz", "vector<vector<double>>", &Theta_yz);
  MCS_tree -> Branch("Theta_3D", "vector<vector<double>>", &Theta_3D);
}

void MCS_Tree_Maker::executeEvent(){

  PID.clear();
  Segment_size.clear();
  P_true.clear();
  Theta_xz.clear();
  Theta_yz.clear();
  Theta_3D.clear();

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
  Run_beam_MCS();

  // == Functions to study daughters
  vector<Daughter> daughters_all = GetAllDaughters();
  vector<Daughter> pions = GetPions(daughters_all);
  if(pions.size() > 0){
    Run_Daughter_MCS(pions);
  }

  MCS_tree -> Fill();

  return;
}

void MCS_Tree_Maker::Run_beam_MCS(){

  if(pi_type != 1 && pi_type != 3) return;

  vector<TVector3> true_position_vec;
  vector<TVector3> true_Pvec_vec;
  vector<double> true_P_vec;
  for(unsigned int i_true_hit = 0; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++){
    TVector3 this_position((*evt.true_beam_traj_X).at(i_true_hit), (*evt.true_beam_traj_Y).at(i_true_hit), (*evt.true_beam_traj_Z).at(i_true_hit));
    TVector3 this_Pvec((*evt.true_beam_traj_Px).at(i_true_hit), (*evt.true_beam_traj_Py).at(i_true_hit), (*evt.true_beam_traj_Pz).at(i_true_hit));
    this_Pvec = this_Pvec * 1000.;// * P_beam_inst_scale;
    double this_P = this_Pvec.Mag();

    true_position_vec.push_back(this_position);
    true_Pvec_vec.push_back(this_position);
    true_P_vec.push_back(this_P);
    //cout << Form("[MCS_Tree_Maker::Run_beam_MCS] %d this_P : %.2f", i_true_hit, this_P) << endl;
  }

  vector<double> reco_X = (*evt.reco_beam_calo_X);
  vector<double> reco_Y = (*evt.reco_beam_calo_Y);
  vector<double> reco_Z = (*evt.reco_beam_calo_Z);
  vector<double> reco_range = (*evt.reco_beam_resRange_SCE);
  vector<TVector3> reco_position_vec;
  for(unsigned int i_reco_hit = 0; i_reco_hit < reco_Z.size(); i_reco_hit++){
    TVector3 this_position(reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit));
    reco_position_vec.push_back(this_position);
    //cout << Form("[MCS_Tree_Maker::Run_beam_MCS] reco %d : (%.2f, %.2f, %.2f)", i_reco_hit, reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit)) << endl;
  }


  // == Test 3D linear fit function
  TString evt_run_str = Form("Run%d_Evt%d", evt.run, evt.event);
  //TVector3 test_fit = Fitter->line3Dfit(reco_position_vec, true, evt_run_str);
  //cout << "[MCS_Tree_Maker::Run_beam_MCS] " << evt_run_str << ", test_fit.Mag() : " << test_fit.Mag() << endl;

  Run_beam_MCS_for_Segments(reco_position_vec, true_position_vec, true_P_vec, 4., "4cm");
  Run_beam_MCS_for_Segments(reco_position_vec, true_position_vec, true_P_vec, 5., "5cm");
  Run_beam_MCS_for_Segments(reco_position_vec, true_position_vec, true_P_vec, 8., "8cm");
  Run_beam_MCS_for_Segments(reco_position_vec, true_position_vec, true_P_vec, 10., "10cm");
  Run_beam_MCS_for_Segments(reco_position_vec, true_position_vec, true_P_vec, 14., "14cm");

  true_position_vec.clear();
  true_Pvec_vec.clear();
  true_P_vec.clear();
  reco_X.clear();
  reco_Y.clear();
  reco_Z.clear();
  reco_range.clear();
  reco_position_vec.clear();

  return;
}

void MCS_Tree_Maker::Run_beam_MCS_for_Segments(const vector<TVector3> & reco_position_vec, const vector<TVector3> & true_position_vec, const vector<double> & true_P_vec, double segment_size, TString name){

  int this_PDG = 211;
  if(pi_type == 3) this_PDG = 13;
  vector<MCSSegment> this_segments = SplitIntoSegments(reco_position_vec, segment_size);
  if(this_segments.size() < 2) return;
  vector<double> segment_true_P = GetSegmentTrueP(this_segments, true_position_vec, true_P_vec, this_PDG);
  for(unsigned int i_seg = 0; i_seg < segment_true_P.size(); i_seg++){
    if(segment_true_P.at(i_seg) < 0.) return;
  }

  vector<double> this_theta_xz_vec;
  vector<double> this_theta_yz_vec;
  vector<double> this_theta_3D_vec;

  for(unsigned int i_seg = 0; i_seg < segment_true_P.size() - 1; i_seg++){
    TVector3 this_vec = this_segments.at(i_seg).FittedVec();
    TVector3 next_vec = this_segments.at(i_seg + 1).FittedVec();
    TVector3 rotated_this_vec = RotateToZaxis(this_vec, this_vec);
    TVector3 rotated_next_vec = RotateToZaxis(this_vec, next_vec);

    double this_theta_yz = TMath::ATan(rotated_next_vec.Y() / rotated_next_vec.Z());
    double this_theta_xz = TMath::ATan(rotated_next_vec.X() / rotated_next_vec.Z());

    double this_theta_3D = rotated_next_vec.Theta();

    this_theta_xz_vec.push_back(this_theta_xz);
    this_theta_yz_vec.push_back(this_theta_yz);
    this_theta_3D_vec.push_back(this_theta_3D);
  }

  int segment_size_int = segment_size + 0.;
  PID.push_back(this_PDG);
  Segment_size.push_back(segment_size_int);
  P_true.push_back(segment_true_P);
  Theta_xz.push_back(this_theta_xz_vec);
  Theta_yz.push_back(this_theta_yz_vec);
  Theta_3D.push_back(this_theta_3D_vec);

  this_segments.clear();
  segment_true_P.clear();
  return;
}

void MCS_Tree_Maker::Run_Daughter_MCS(const vector<Daughter>& pions){
  
  int true_beam_ID = evt.true_beam_ID;
  double beam_last_X = (*evt.reco_beam_calo_X).back();
  double beam_last_Y = (*evt.reco_beam_calo_Y).back();
  double beam_last_Z = (*evt.reco_beam_calo_Z).back();
  TVector3 beam_end(beam_last_X, beam_last_Y, beam_last_Z);

  for(unsigned int i_pion = 0; i_pion < pions.size(); i_pion++){
    Daughter this_daughter = pions.at(i_pion);

    int this_true_ID = this_daughter.PFP_true_byHits_ID();
    if(true_beam_ID == this_true_ID) continue;

    int this_PdgID = this_daughter.PFP_true_byHits_PDG();

    vector<double> reco_X = this_daughter.allTrack_calo_X();
    vector<double> reco_Y = this_daughter.allTrack_calo_Y();
    vector<double> reco_Z = this_daughter.allTrack_calo_Z();
    vector<TVector3> reco_position_vec;
    for(int i_reco_hit = reco_Z.size() - 1; i_reco_hit > -1; i_reco_hit--){

      TVector3 this_position(reco_X.at(i_reco_hit), reco_Y.at(i_reco_hit), reco_Z.at(i_reco_hit));
      reco_position_vec.push_back(this_position);
    }
    double true_P = this_daughter.PFP_true_byHits_startP() * 1000.;

    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 4., "4cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 5., "5cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 8., "8cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 10., "10cm");
    Run_Daughter_MCS_for_Segments(reco_position_vec, true_P, this_PdgID, 14., "14cm");
  }

  return;
}

void MCS_Tree_Maker::Run_Daughter_MCS_for_Segments(const vector<TVector3> & reco_position_vec, double true_P, int this_PdgID, double segment_size, TString name){

  vector<MCSSegment> this_segments = SplitIntoSegments(reco_position_vec, segment_size);
  if(this_segments.size() < 3) return;

  TString PdgID_sign = "";
  if(this_PdgID < 0){
    PdgID_sign = "m";
  }
  else if(this_PdgID > 0){
    PdgID_sign = "p";
  }
  TString this_PdgID_str = Form(PdgID_sign + "%d", abs(this_PdgID));

  vector<double> this_P_vec;
  vector<double> this_theta_xz_vec;
  vector<double> this_theta_yz_vec;
  vector<double> this_theta_3D_vec;
  for(unsigned int i_seg = 0; i_seg < this_segments.size() - 1; i_seg++){
    TVector3 this_vec = this_segments.at(i_seg).FittedVec();
    TVector3 next_vec = this_segments.at(i_seg + 1).FittedVec();
    TVector3 rotated_this_vec = RotateToZaxis(this_vec, this_vec);
    TVector3 rotated_next_vec = RotateToZaxis(this_vec, next_vec);

    double this_theta_yz = TMath::ATan(rotated_next_vec.Y() / rotated_next_vec.Z());
    double this_theta_xz = TMath::ATan(rotated_next_vec.X() / rotated_next_vec.Z());
    double this_theta_3D = rotated_next_vec.Theta();

    this_theta_xz_vec.push_back(this_theta_xz);
    this_theta_yz_vec.push_back(this_theta_yz);
    this_theta_3D_vec.push_back(this_theta_3D);
  }

  this_P_vec.push_back(true_P);
  int segment_size_int = segment_size + 0.;

  PID.push_back(abs(this_PdgID));
  Segment_size.push_back(segment_size_int);
  P_true.push_back(this_P_vec);
  Theta_xz.push_back(this_theta_xz_vec);
  Theta_yz.push_back(this_theta_yz_vec);
  Theta_3D.push_back(this_theta_3D_vec);

  return;
}

vector<double> MCS_Tree_Maker::GetSegmentTrueP(const vector<MCSSegment> & segments, const vector<TVector3> & true_position_vec, const vector<double> true_P_vec, int PDG){
  vector<double> out;

  if(segments.size() == 0) return out;

  //cout << "[GetSegmentTrueP] Start" << endl;

  for(unsigned int i = 0; i < segments.size(); i++){

    MCSSegment this_segment = segments.at(i);
    TVector3 this_start_hit = this_segment.RecoStart();
    double true_P = -999.;
    for(unsigned int j = 0; j < true_position_vec.size() - 1; j++){
      TVector3 this_true_hit = true_position_vec.at(j);
      TVector3 next_true_hit = true_position_vec.at(j + 1);
      if(this_start_hit.Z() > this_true_hit.Z() && this_start_hit.Z() < next_true_hit.Z()){
	double distance = (this_true_hit - this_start_hit).Mag();
	if(distance < 10.) true_P = map_BB[PDG] -> KEtoMomentum(map_BB[PDG] -> KEAtLength(map_BB[PDG] -> MomentumtoKE(true_P_vec.at(j)), distance));
	break;
      }
    }
    //cout << i << ", true_P : " << true_P << endl;
    out.push_back(true_P);
  }


  return out;
}

void MCS_Tree_Maker::Write_Tree(){
  outfile->cd();
  MCS_tree -> Write(); 
}

MCS_Tree_Maker::MCS_Tree_Maker(){
  MCS_tree = NULL;
}

MCS_Tree_Maker::~MCS_Tree_Maker(){

}
