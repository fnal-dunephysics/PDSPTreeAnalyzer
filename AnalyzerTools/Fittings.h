#ifndef FITTINGS_H
#define FITTINGS_H

#include <map>
#include "TMath.h"
#include "TGraph2D.h"
#include "TRandom2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "Math/Functor.h"
#include "TPolyLine3D.h"
#include "TVector3.h"
#include "Fit/Fitter.h"
#include "TF1.h"

using namespace std;

class Fittings {

 public:

  Fittings();
  ~Fittings();

  map< TString, TGraph2D* > map_TGraph2D;
  map< TString, TPolyLine3D* > map_TPolyLine3D;
  // == Funtions for linear fitting on 3D points 
  void line_3D(double t, const double *p, double &x, double &y, double &z) {
    // line = (x0, y0, z0) + (x1, y1, z1)*t
    // If line is not parallel to x-y plane, we can select a point with z0 = 0
    // Since magnitude of the direction vector is not required to be 1, we can set z1 = 1 without loosing generality
    x = p[0] + p[1]*t;
    y = p[2] + p[3]*t;
    z = p[4] + t;
  }

  struct SumDistance2 {
    // the TGraph is a data member of the object
    TGraph2D *fGraph;

    SumDistance2(TGraph2D *g) : fGraph(g) {}

    // calculate distance line-point
    double distance2(double x,double y,double z, const double *p) {
      // distance line point is D= | (xp-x0) cross  ux |
      // where ux is direction of line and x0 is a point in the line (like t = 0)
      TVector3 xp(x,y,z);
      TVector3 x0(p[0], p[2], p[4] );
      TVector3 x1(p[0] + p[1], p[2] + p[3], p[4] + 1. );
      TVector3 u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      return d2;
    }

    // implementation of the function to be minimized
    double operator() (const double *par) {
      assert(fGraph != 0);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i  = 0; i < npoints; ++i) {
	double d = distance2(x[i],y[i],z[i],par);
	sum += d;
      }
      return sum;
    }
  };

  vector<TVector3> Linear3Dfit(const vector<TVector3> & hits, bool save_gr = false, TString outname = "");

};

#endif
