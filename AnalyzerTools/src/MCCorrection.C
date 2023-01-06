#include "MCCorrection.h"

MCCorrection::MCCorrection() :
  IgnoreNoHist(false)
{

  histDir = TDirectoryHelper::GetTempDirectory("MCCorrection");

}

void MCCorrection::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");

  TDirectory* origDir = gDirectory;

  //==== Momentum reweight
  TString MomentumReweight_path = datapath + "/Momentum_reweight/";
  TString MomentumReweight_path_xrootd = "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/sungbino/PDSP_data/Momentum_reweight/"; 
  //cout << "[MCCorrection::ReadHistograms] MomentumReweight_path  : " << MomentumReweight_path << endl;
  string elline;
  ifstream in(MomentumReweight_path + "histmap.txt");
  while(getline(in,elline)){
    std::istringstream is( elline );

    TString tstring_elline = elline;
    if(tstring_elline.Contains("#")) continue;
    //cout << "[MCCorrection::ReadHistograms] tstring_elline : " << tstring_elline << endl;
    TString a,b,c,d,e,f;
    is >> a; // Beam,
    is >> b; // SF, Eff
    is >> c; // <Flag>
    is >> d; // <rootfilename>
    is >> e; // <histname>
    is >> f; // Class
    TFile *file = TFile::Open(MomentumReweight_path_xrootd + d);
    //cout << "MomentumReweight_path + d : " << MomentumReweight_path + d << endl;
    //cout << "a : " << a << ", b : " << b << ", c : << " << c << ", d : " << d << ", e : " << e << ", f : " << f << endl;

    if(f=="TH1D"){
      histDir->cd();
      map_hist_MomentumReweight[a + "_" + b + "_" + c] = (TH1D *)file -> Get(e) -> Clone();
    }
  }

  cout << "[MCCorrection::ReadHistograms] map_hist_MomentumReweight :" << endl;
  for(std::map< TString, TH1D* >::iterator it = map_hist_MomentumReweight.begin(); it != map_hist_MomentumReweight.end(); it++){
    cout << "[MCCorrection::ReadHistograms] key = " << it -> first << endl;
  }

  cout << "[MCCorrection::ReadHistograms] END" << endl;
}

MCCorrection::~MCCorrection(){

}

void MCCorrection::SetIsData(bool b){
  IsData = b;
}

double MCCorrection::MomentumReweight_SF(TString flag, double P_beam_inst, int sys){

  double value = 1.;
  double error = 0.;

  TH1D *this_hist = map_hist_MomentumReweight["Beam_SF_" + flag];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::MuonReco_SF] No "<< "Beam_SF_" + flag << endl;
      exit(EXIT_FAILURE);
    }
  }

  int this_bin = this_hist->FindBin(P_beam_inst);
  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  return value+double(sys)*error;

}
