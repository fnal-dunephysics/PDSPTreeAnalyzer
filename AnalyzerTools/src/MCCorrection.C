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
  TString MomentumReweight_path_xrootd = "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/scratch/users/sungbino/PDSP_data/Momentum_reweight/"; 
  cout << "[MCCorrection::ReadHistograms] MomentumReweight_path  : " << MomentumReweight_path << endl;
  string elline;
  ifstream in(MomentumReweight_path + "histmap.txt");
  while(getline(in,elline)){
    std::istringstream is( elline );

    TString tstring_elline = elline;
    if(tstring_elline.Contains("#")) continue;
    cout << "[MCCorrection::ReadHistograms] tstring_elline : " << tstring_elline << endl;
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

double MCCorrection::dEdx_scaled(double MC_dEdx){
  double f_const_p0 = 0.997;
  double f_pol1_p0 = 0.945;
  double f_pol1_p1 = 0.021;

  double scale = 1.;
  if(MC_dEdx < 2.44) scale = f_const_p0;
  else scale = f_pol1_p0 + MC_dEdx * f_pol1_p1;

  if(scale > 1.2) scale = 1.2;

  //cout << "[dEdx_res::dEdx_scaled] scale : " << scale << endl; 
  return scale * MC_dEdx;
}

double MCCorrection::Use_Abbey_Recom_Params(double dEdx, double Efield, double calib_const_ratio){

  double alpha_default = 0.93;
  double beta_default = 0.212;
  double rho = 1.396;

  double alpha_Abbey = 0.905;
  double beta_Abbey = 0.220;

  double exp_term = exp( calib_const_ratio * (beta_Abbey / beta_default) * log(alpha_default + beta_default * dEdx / (rho * Efield)) );
  double new_dEdx = (exp_term - alpha_Abbey) * rho * Efield / beta_Abbey;

  cout << "[MCCorrection::Use_Abbey_Recom_Params] dEdx : " << dEdx << ", new_dEdx : " << new_dEdx << endl;

  return new_dEdx;
}
