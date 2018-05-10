void TMVApplicationD( )
{

  // Plots
  TH1F* h = new TH1F("1","Invariant Mass",24,1.65,1.89);
  TH1F* h2 = new TH1F("2","BDT response",100,-1,1);
  // Local variables to store to outfile
  Float_t mass;

  // Create output file
  TFile *outFile = TFile::Open("InvMassB.root", "RECREATE");
  TTree *outTree = new TTree("Events","Events");
  outTree->Branch("mass", &mass,"mass/F");
  
  TMVA::Reader *reader = new TMVA::Reader("Color");
  // Link local values and weight file values
  Float_t var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12, BDT_mass; 
  reader->AddVariable( "BDT_ptMin", &var1 );
  reader->AddVariable( "BDT_etaMax", &var2 );
  reader->AddVariable( "BDT_trkKinkMax", &var3 );
  reader->AddVariable( "BDT_d0Min", &var4 );
  reader->AddVariable( "BDT_muNchi2Max", &var5 );
  reader->AddVariable( "BDT_nMatchStnMin", &var6 );
  reader->AddVariable( "BDT_nTkLayersMin", &var7 );
  reader->AddVariable( "BDT_nC", &var8 );
  reader->AddVariable( "BDT_deltaRMin", &var9 );
  reader->AddVariable( "BDT_deltaRMax", &var10 );
  reader->AddVariable( "BDT_vfIp", &var11 );
  reader->AddVariable( "BDT_triIsoNtrk", &var12 );
  
  reader->BookMVA( "BDT classifier", "weights/MVAnalysis_BDT.weights.xml" );
  TFile *input = TFile::Open("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/DataFlat/ntuples/data_evaluation_select_MRC1.root");
  TTree* theTree = (TTree*)input->Get("Events");
  // Improt variables for BDT evaluation
  theTree->SetBranchAddress("BDT_ptMin", &var1);
  theTree->SetBranchAddress("BDT_etaMax", &var2);
  theTree->SetBranchAddress("BDT_trkKinkMax", &var3);
  theTree->SetBranchAddress("BDT_d0Min", &var4);
  theTree->SetBranchAddress("BDT_muNchi2Max", &var5);
  theTree->SetBranchAddress("BDT_nMatchStnMin", &var6);
  theTree->SetBranchAddress("BDT_nTkLayersMin", &var7);
  theTree->SetBranchAddress("BDT_nC", &var8);
  theTree->SetBranchAddress("BDT_deltaRMin", &var9);
  theTree->SetBranchAddress("BDT_deltaRMax", &var10);
  theTree->SetBranchAddress("BDT_vfIp", &var11);
  theTree->SetBranchAddress("BDT_triIsoNtrk", &var12);

  // Import variables for output
  theTree->SetBranchAddress("BDT_mass", &BDT_mass);
  
  for (int ievt = 0; ievt<theTree->GetEntries();ievt++) {
    theTree->GetEntry(ievt);
    // BDT evaluation
    Float_t response = reader->EvaluateMVA( "BDT classifier" );

    h2->Fill(response);
    // Cut on BDT response
    if(response < 0) continue;
    // Blind signal region
    if(BDT_mass > 1.72 && BDT_mass < 1.82) continue;
    mass  = BDT_mass;
    h->Fill(mass);
    outTree->Fill();
  }
  outFile->Write();
  outFile->Close();
  delete reader;

  TCanvas *c1 = new TCanvas("1","1",1200,900);
  c1->cd();
  h->Draw();
  c1->Print("out.png");

  TCanvas *c2 = new TCanvas("2","2",1200,900);
  c2->cd();
  h2->Draw();
  c2->Print("res.png");
}
