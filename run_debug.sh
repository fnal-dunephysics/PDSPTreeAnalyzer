root -l -b<<EOF
cout << "aaaaa" << endl
PionKEScale m;
m.MaxEvent = 1000;
m.LogEvery = 1000;
m.MCSample = "Piplus_1GeV_Ar"
m.Beam_Momentum = "1GeV"
m.SetTreeName()
m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/mc/physics/PDSPProd4a/18/80/01/67/PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1_ntuple_v09_41_00_03.root") // 1GeV MC
m.AddFile("xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2022/detector/physics/PDSPProd4/00/00/52/19/PDSPProd4_data_1GeV_reco2_ntuple_v09_41_00_04.root") // 1 GeV data
m.SetOutfilePath("hists.root");
m.Init();
cout << "Running" << endl;
m.initializeAnalyzer();
m.initializeAnalyzerTools();
m.SwitchToTempDir();
m.Loop();

m.WriteHist();

EOF
