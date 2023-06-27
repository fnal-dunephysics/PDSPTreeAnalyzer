R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(libGpad.so)
R__LOAD_LIBRARY(libGraf3d.so)
R__LOAD_LIBRARY(libAnalyzerTools.so)
R__LOAD_LIBRARY(libDataFormats.so)
R__LOAD_LIBRARY(libAnalyzers.so)

void run_recom_fit(){

  Fit_Modified_Box m;
  m.MaxEvent = -1;
  //m.MaxEvent = 1000;
  //m.NSkipEvent = 0;
  //m.MaxEvent = m.NSkipEvent + 1000;
  m.LogEvery = 1000;
  m.MCSample = "Data";
  //m.MCSample = "MC";
  m.Beam_Momentum = 1.0;
  m.alpha_new = 0.905;
  m.beta_new = 0.220;
  m.SetTreeName();
  //m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/20/91/32/85/PDSPProd4a_MC_2GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"); // 2 GeV MC
  //m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/18/80/01/67/PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root"); // 1 GeV MC
  //m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/22/59/77/18/PDSPProd4a_MC_0.5GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_04.root"); // 0.5 GeV MC
  //m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/54/29/PDSPProd4_data_2GeV_reco2_ntuple_v09_42_03_01.root"); // 2 GeV data
  m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/52/19/PDSPProd4_data_1GeV_reco2_ntuple_v09_41_00_04.root"); // 1 GeV data
  //m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/58/25/PDSPProd4_data_0.5GeV_reco2_ntuple_v09_41_00_04.root"); // 0.5 GeV data
  m.SetOutfilePath("hists_Data_Fit_Modified_Box_1.0GeV.root");
  //m.SetOutfilePath("hists_MC_dEdx_res_1.0GeV_no_trun.root");
  m.Init();
  m.initializeAnalyzer();  
  m.initializeAnalyzerTools();  
  m.SwitchToTempDir();
  m.Loop();
  m.Extract_MPVs();
  //m.WriteHist();
}
