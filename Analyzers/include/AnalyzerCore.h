#ifndef AnalyzerCore_h
#define AnalyzerCore_h

#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"
#include "TH3.h"
#include "TH2D.h"
#include <sstream>
#include "TRandom.h"

#include "PDSPTree.h"
#include "Event.h"
#include "Daughter.h"
#include "BetheBloch.h"
#include "GEANT4_XS.h"
#include "MCCorrection.h"

#define M_Z 91.1876
#define M_W 80.379
#define M_mu 105.65837/1e3
#define M_neutron 939.565/1e3
#define M_proton 938.272/1e3
#define M_pion 139.570/1e3
#define M_Kaon 493.677/1e3
#define M_e 0.510998/1e3
#define M_pizero 134.976/1e3

class AnalyzerCore {

public:

  AnalyzerCore();
  ~AnalyzerCore();

  PDSPTree evt;
  Long64_t MaxEvent, NSkipEvent;
  int LogEvery;
  TString MCSample;
  TString Beam_Momentum;
  vector<TString> Userflags;

  bool debug_mode = false;

  virtual void initializeAnalyzer(){

  };

  virtual void executeEvent(){

  };

  //==================
  // Read Tree
  //==================
  virtual void SetTreeName(){
    TString tname = "/pduneana/beamana";
    fChain = new TChain(tname);
  }

  virtual void AddFile(TString filename){
    fChain->Add(filename);
  }

  Int_t GetEntry(Long64_t entry);

  TChain *fChain;

  void Init();

  std::string printcurrunttime(){

    std::stringstream out;
    TDatime datime;
    out << datime.GetYear()<<"-"<<AddZeroToTime(datime.GetMonth())<<"-"<<AddZeroToTime(datime.GetDay())<<" "<<AddZeroToTime(datime.GetHour())<<":"<<AddZeroToTime(datime.GetMinute())<<":"<<AddZeroToTime(datime.GetSecond());
    return out.str();

  }

  std::string AddZeroToTime(int twodigit){
    if(twodigit<10){
      return "0"+std::to_string(twodigit);
    }
    else{
      return std::to_string(twodigit);
    }
  }

  void Loop();

  //==================
  // Get Particles
  //==================
  inline static bool PComparing(const Particle& p1, const Particle& p2){ return (p1.P() > p2.P()); }

  Event GetEvent();

  std::vector<Gen> GetAllParticles();
  std::vector<Gen> GetNeutrons(const std::vector<Gen>& particles, double min_P);

  //==================
  // Tools
  //==================
  std::map< int, BetheBloch* > map_BB;
  MCCorrection MCCorr;
  GEANT4_XS G4Xsec;
  void initializeAnalyzerTools();

  //==================
  // Functions
  //==================
  

  //===Plotting
  std::map< TString, TH1D* > maphist_TH1D;
  std::map< TString, TH2D* > maphist_TH2D;
  std::map< TString, TH3D* > maphist_TH3D;

  TH1D* GetHist1D(TString histname);
  TH2D* GetHist2D(TString histname);
  TH3D* GetHist3D(TString histname);

  void FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  void FillHist(TString histname, double value, double weight, int n_bin, double *xbins);
  void FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double x_min, double x_max,
                int n_biny, double y_min, double y_max);
  void FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double *xbins,
                int n_biny, double *ybins);
  void FillHist(TString histname,
                double value_x, double value_y, double value_z,
                double weight,
                int n_binx, double x_min, double x_max,
                int n_biny, double y_min, double y_max,
                int n_binz, double z_min, double z_max);
  void FillHist(TString histname,
                double value_x, double value_y, double value_z,
                double weight,
                int n_binx, double *xbins,
                int n_biny, double *ybins,
                int n_binz, double *zbins);

  //==== JSFillHist : 1D
  std::map< TString, std::map<TString, TH1D*> > JSmaphist_TH1D;
  TH1D* JSGetHist1D(TString suffix, TString histname);
  void JSFillHist(TString suffix, TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  //==== JSFillHist : 2D
  std::map< TString, std::map<TString, TH2D*> > JSmaphist_TH2D;
  TH2D* JSGetHist2D(TString suffix, TString histname);
  void JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double x_min, double x_max,
                  int n_biny, double y_min, double y_max);
  void JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double *xbins,
                  int n_biny, double *ybins);

  virtual void WriteHist();
  
  //==== Output rootfile 
  void SwitchToTempDir();
  TFile *outfile;
  void SetOutfilePath(TString outname);

};

#endif
