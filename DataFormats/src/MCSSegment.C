#include "MCSSegment.h"

MCSSegment::MCSSegment(){

}

void MCSSegment::TestClass(){
  cout << "[[MCSSegment::TestClass]] Test " << endl;
}

void MCSSegment::SetMCSSegment(TVector3 i_reco_hit_start, TVector3 i_reco_hit_end, TVector3 i_fitted_start, TVector3 i_fitted_end, TVector3 i_fitted_vec){
  j_reco_hit_start = i_reco_hit_start;
  j_reco_hit_end = i_reco_hit_end;
  j_fitted_start = i_fitted_start;
  j_fitted_end = i_fitted_end;
  j_fitted_vec = i_fitted_vec;
}

MCSSegment::~MCSSegment(){

}
