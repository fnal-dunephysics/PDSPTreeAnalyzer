#include "Fittings.h"
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

vector<TVector3> Fittings::Linear3Dfit(const vector<TVector3> & hits, bool save_gr, TString outname){

  // == Refer https://root.cern.ch/doc/master/line3Dfit_8C.html
  //cout << "[Fittings::line3Dfit] Start" << endl;

  vector<TVector3> out; // at(0) starting point, at(1) direction vector

  // == Make a TGraph2D for fitting using hit vectors
  int N_hits = 0;
  vector<double> this_X_vec;
  vector<double> this_Y_vec;
  vector<double> this_Z_vec;
  for(unsigned int i_hit = 0; i_hit < hits.size(); i_hit++){
    this_X_vec.push_back(hits.at(i_hit).X());
    this_Y_vec.push_back(hits.at(i_hit).Y());
    this_Z_vec.push_back(hits.at(i_hit).Z());
    N_hits++;
  }
  //cout << "[Fittings::line3Dfit] Saved vectors for N_hits : " << N_hits << ", (hits.size() : " << hits.size() << ")" << ", this_X_vec.size() : " << this_X_vec.size() << endl;
  TGraph2D * this_gr = new TGraph2D(N_hits, &this_X_vec[0], &this_Y_vec[0], &this_Z_vec[0]);
  //cout << "[Fittings::line3Dfit] Saved graph" << endl;

  // == Fit
  ROOT::Fit::Fitter fitter;
  SumDistance2 chi2_def(this_gr);
  ROOT::Math::Functor fitted_obj(chi2_def,4);
  double param_start[4] = {1,1,1,1};
  fitter.SetFCN(fitted_obj, param_start);
  for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

  bool ok = fitter.FitFCN();
  if (!ok) {
    //Error("[Fittings::line3Dfit]","Fit failed");
    this_gr -> Clear();   
    return out;
  }
  const ROOT::Fit::FitResult & result = fitter.Result();  

  // == Summarize result
  const double * parFit = result.GetParams();
  TVector3 dir_vec(parFit[1], parFit[3], 1.);
  double first_hit_Z = hits.at(0).Z();
  double last_hit_Z = hits.back().Z();
  TVector3 fitted_start(parFit[0] + parFit[1] * first_hit_Z, parFit[2] + parFit[3] * first_hit_Z, first_hit_Z);
  TVector3 fitted_end(parFit[0] + parFit[1] * last_hit_Z, parFit[2] + parFit[3] * last_hit_Z, last_hit_Z);
  out.push_back(fitted_start);
  out.push_back(fitted_end - fitted_start);

  //cout << "[Fittings::line3Dfit] first_hit_Z : " << first_hit_Z << ", last_hit_Z : " << last_hit_Z  << endl;

  this_X_vec.clear();
  this_Y_vec.clear();
  this_Z_vec.clear();

  
  if(save_gr){
    this_gr -> Draw("p0");
    int n = 1000;
    double t0 = first_hit_Z;
    double t1 = last_hit_Z;
    double dt = (t1 - t0) / (n + 0.);
    TGraph2D *l = new TGraph2D();
    for (int i = 0; i <n; ++i) {
      double t = t0+ dt*(i + 0.);
      double x,y,z;
      line_3D(t,parFit,x,y,z);
      l->SetPoint(i,x,y,z);
    }
    l->SetLineColor(kRed);
    //l->Draw("same");
    this_gr -> SetName(outname + "_Original");
    l -> SetName(outname + "_Fitted");
    map_TGraph2D[outname + "_Original"] = this_gr;
    map_TGraph2D[outname + "_Fitted"] = l;
  }
  this_gr -> Clear();
  return out;
}

Fittings::Fittings(){

}

Fittings::~Fittings(){

}
