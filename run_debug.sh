root -l -b<<EOF
SmearParticles a
a.TestClass()
//cout << "aaaaa" << endl
PionAnalyzer m
m.MaxEvent = 2000000
m.LogEvery = 100000
m.MCSample = "Piplus_1GeV_Ar"
m.Simulator = "GEANT4";
//m.Simulator = "FLUKA"
m.SetTreeName()
//m.AddFile("/data6/Users/suoh/FNAL/pion_study/samples/GEANT4/XS_Pi+_G4_lAr_1GeV_10b_0.root") // == GEANT4
//m.AddFile("/data6/Users/suoh/FNAL/pion_study/samples/FLUKA/18files/pion1GeVa001_Out.root") // == FLUKA 
m.AddFile("/dune/data/users/tjyang/G4HadStudies/XS_Pi+_G4_lAr_1GeV_10b_0.root")
m.SetOutfilePath("hists.root");
m.Init();
cout << "Running" << endl;
m.initializeAnalyzer();
m.initializeAnalyzerTools();
m.SwitchToTempDir();
m.Loop();

m.WriteHist();

EOF
