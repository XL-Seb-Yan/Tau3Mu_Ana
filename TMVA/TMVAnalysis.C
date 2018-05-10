void TMVAnalysis( )
{
  TFile* outFile = TFile::Open("TMVA.root", "RECREATE");

  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile *infile_s = TFile::Open("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/MCFlat/ntuples/dstau_Tau3Mu_select_MRC3.root");
  TFile *infile_b = TFile::Open("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/DataFlat/ntuples/data_background_select_MRC3.root");

  factory->AddSignalTree((TTree*)infile_s->Get("Events"));
  factory->AddBackgroundTree((TTree*)infile_b->Get("Events"));
  factory->AddVariable("BDT_ptMin", 'F');
  factory->AddVariable("BDT_etaMax", 'F');
  factory->AddVariable("BDT_trkKinkMax", 'F');
  factory->AddVariable("BDT_d0Min", 'F');
  factory->AddVariable("BDT_muNchi2Max", 'F');
  factory->AddVariable("BDT_nMatchStnMin", 'F');
  factory->AddVariable("BDT_nTkLayersMin", 'F');
  factory->AddVariable("BDT_nC", 'F');
  factory->AddVariable("BDT_deltaRMin", 'F');
  factory->AddVariable("BDT_deltaRMax", 'F');
  factory->AddVariable("BDT_vfIp", 'F');
  factory->AddVariable("BDT_triIsoNtrk", 'F');

  //factory->PrepareTrainingAndTestTree("","nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=-1" );

  //factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outFile->Close();
  delete factory;
}
