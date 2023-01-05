#ifndef AnalyzerCore_h
#define AnalyzerCore_h

#include "TLorentzVector.h"
#include "TString.h"
#include "TMath.h"
#include "TH3.h"
#include "TH2D.h"
#include <sstream>
#include "TRandom.h"
#include "TProfile.h"

#include "PDSPTree.h"
#include "Event.h"
#include "Daughter.h"
#include "BetheBloch.h"
#include "GEANT4_XS.h"
#include "MCCorrection.h"

//#define M_Z 91.1876
//#define M_W 80.379
#define M_mu 105.65837
#define M_neutron 939.565
#define M_proton 938.272
#define M_pion 139.570
#define M_Kaon 493.677
#define M_e 0.510998
#define M_pizero 134.976

class AnalyzerCore {

public:

  AnalyzerCore();
  virtual ~AnalyzerCore();

  PDSPTree evt;
  Long64_t MaxEvent, NSkipEvent;
  int LogEvery;
  TString MCSample;
  double Beam_Momentum;
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
  void Init_evt();

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
  Event GetEvent();
  std::vector<Daughter> GetAllDaughters();
  std::vector<Daughter> GetPions(const vector<Daughter>& in);
  std::vector<Daughter> GetProtons(const vector<Daughter>& in);
  std::vector<Daughter> GetTruePions(const vector<Daughter>& in);
  std::vector<Daughter> GetTrueProtons(const vector<Daughter>& in);

  //==================
  // AnalyzerTools
  //==================
  std::map< int, BetheBloch* > map_BB;
  MCCorrection *MCCorr;
  GEANT4_XS *G4Xsec;
  void initializeAnalyzerTools();

  //==================
  // Set Beam Variables
  //==================
  double P_beam_inst = 1000.;
  double mass_beam = 139.57;
  double P_ff_reco = -999.;
  double KE_ff_reco = -999.;
  double KE_end_reco = -999.;
  double E_end_reco = -999.;

  int pandora_slice_pdg;
  void SetPandoraSlicePDG(int pdg);
  double daughter_michel_score;
  double beam_costh;
  double chi2_proton;
  double chi2_pion;
  double chi2_muon;
  std::map< int, TProfile* > map_profile;
  int GetPiParType();

  //==================
  // Event Selections
  //==================
  bool Pass_Beam_PID(int PID);
  bool PassBeamScraperCut() const;
  double P_beam_inst_scale = 1.;
  double beam_momentum_low = 0.;
  double beam_momentum_high = 10000.;
  bool PassBeamMomentumWindowCut() const;
  bool PassPandoraSliceCut() const;
  bool PassCaloSizeCut() const;  
  bool PassAPA3Cut(double cut = 215.);
  bool PassMichelScoreCut(double cut = 0.55);
  const double beam_angleX_data = 100.464;
  const double beam_angleY_data = 103.442;
  const double beam_angleZ_data = 17.6633;
  const double beam_angleX_mc = 101.579;
  const double beam_angleY_mc = 101.212;
  const double beam_angleZ_mc = 16.5822;
  bool PassBeamCosCut(double cut = 0.95);
  bool PassBeamStartZCut(double cut = -5.);
  bool PassProtonVetoCut(double cut = 80.);
  bool PassMuonVetoCut(double cut = 8.);
  bool PassStoppedPionVetoCut(double cut = 13.);

  //==================
  // Additional Functions
  //==================
  double Convert_P_Spectrometer_to_P_ff(double P_beam_inst, TString particle, TString key, int syst);
  double Particle_chi2(const vector<double> & dEdx, const vector<double> & ResRange, int PID, double dEdx_res_frac = 1.);
  double max_additional_res_length_pion = 450.; // == [cm] 
  double max_additional_res_length_proton = 120.; // == [cm]
  double res_length_step_pion = 1.0; // == [cm]
  double res_length_step_proton = 0.2; // == [cm]
  double dEdx_truncate_upper_pion = 5.;
  double dEdx_truncate_bellow_pion = 0.5;
  double dEdx_truncate_upper_proton = 20.;
  double dEdx_truncate_bellow_proton = 0.2;
  double Fit_HypTrkLength_Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam);
  double Fit_HypTrkLength_Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam);

  //==================
  //===Plotting
  //==================
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

namespace pi{
  const unsigned int nIntTypes = 9;
  const char intTypeName[nIntTypes+1][100] = {"Data",
					      "PiInel",
					      "PiElas",
					      "Muon",
					      "misID:cosmic",
					      "misID:p",
					      "misID:pi",
					      "misID:mu",
					      "misID:e/#gamma",
					      "misID:other"};
  enum intType{
    kData,
    kPiInel,
    kPiElas,
    kMuon,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}

namespace p{
  const unsigned int nIntTypes = 8;
  const char intTypeName[nIntTypes+1][100] = {"Data",
					      "PInel",
					      "PElas",
					      "misID:cosmic",
					      "misID:p",
					      "misID:pi",
					      "misID:mu",
					      "misID:e/#gamma",
					      "misID:other"};
  enum intType{
    kData,
    kPInel,
    kPElas,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}
#endif
