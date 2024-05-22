#include "MCSSegment.h"

MCSSegment::MCSSegment(){

}

void MCSSegment::TestClass(){
  cout << "[[MCSSegment::TestClass]] Test " << endl;
}

void MCSSegment::SetMCSSegment(TVector3 i_reco_hit_start, TVector3 i_reco_hit_end, TVector3 i_fitted_start, TVector3 i_fitted_end, TVector3 i_fitted_vec, double i_range){
  j_reco_hit_start = i_reco_hit_start;
  j_reco_hit_end = i_reco_hit_end;
  j_fitted_start = i_fitted_start;
  j_fitted_end = i_fitted_end;
  j_fitted_vec = i_fitted_vec;
  j_range = i_range;

  /*
  cout << Form("[MCSSegment::SetMCSSegment] Reco hit start (%.2f, %.2f, %.2f), Reco hit end (%.2f, %.2f, %.2f), fitted vec : (%.2f, %.2f, %.2f)",
	       i_reco_hit_start.X(), i_reco_hit_start.Y(), i_reco_hit_start.Z(),
	       i_reco_hit_end.X(), i_reco_hit_end.Y(), i_reco_hit_end.Z(),
	       i_fitted_vec.X(), i_fitted_vec.Y(), i_fitted_vec.Z()) << endl;
  */
}

MCSSegment::~MCSSegment(){

}
