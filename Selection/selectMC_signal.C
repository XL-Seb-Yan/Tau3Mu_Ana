//================================================================================================
//
// Select tau->mumumu candidates RECO level
//
//  * outputs ROOT files of events passing selection
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include "TLorentzVector.h"         // 4-vector class
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include "TH1D.h"
#include "TRandom.h"
#include "ConfParse.hh"             // input conf file parser
#include "../Utils/CSample.hh"      // helper class to handle samples
// define structures to read in ntuple
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
// lumi section selection with JSON files
#include "BaconAna/Utils/interface/RunLumiRangeMap.hh"
#include "../Utils/LeptonIDCuts.hh" // helper functions for lepton ID selection
#include "../Utils/MyTools.hh"      // various helper functions
#endif

void selectMC(const TString conf="samples.conf", // input file
	       const TString outputDir=".",  // output directory
	       const Bool_t  doScaleCorr=0   // apply energy scale corrections?
	       ) {
  gBenchmark->Start("select3Mu");
  gSystem->Load("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/Tau3Mu_Ana/Selection/lib/libMylib.so");
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  //TH1F* hist0 = new TH1F("tau -> 3 muon 0","MC RECO m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}",60,1.5,2.1);
  TH1F* hist0 = new TH1F("tau -> 3 muon 0","MC RECO m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}",80,1.4,2.2);
  TH1F* hist1a = new TH1F("tau -> 3 muon 1a","MC RECO m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}",80,1.4,2.2);
  TH1F* hist1b = new TH1F("tau -> 3 muon 1b","MC RECO m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}",80,1.4,2.2);
  TH1F* hist1c = new TH1F("tau -> 3 muon 1c","MC RECO m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}",80,1.4,2.2);
  TH1F* hist1d = new TH1F("tau -> 3 muon 1d","MC RECO m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}",50,0,50);
  
  TH1F* hist2 = new TH1F("tau -> 3 muon 2","MC RECO p_{T}",125,0,25);
  TH1F* hist3 = new TH1F("tau -> 3 muon 3","MC RECO #eta",50,-3,3);
  TH1F* hist4 = new TH1F("tau -> 3 muon 4","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist5 = new TH1F("tau -> 3 muon 5","MC RECO p_{T}",125,0,25);
  TH1F* hist6 = new TH1F("tau -> 3 muon 6","MC RECO #eta",50,-3,3);
  TH1F* hist7 = new TH1F("tau -> 3 muon 7","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist8 = new TH1F("tau -> 3 muon 8","MC RECO p_{T}",125,0,25);
  TH1F* hist9 = new TH1F("tau -> 3 muon 9","MC RECO #eta",50,-3,3);
  TH1F* hist10 = new TH1F("tau -> 3 muon 10","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist2b = new TH1F("tau -> 3 muon 2b","MC RECO p_{T}",125,0,255);
  TH1F* hist3b = new TH1F("tau -> 3 muon 3b","MC RECO #eta",50,-3,3);
  TH1F* hist4b = new TH1F("tau -> 3 muon 4b","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist5b = new TH1F("tau -> 3 muon 5b","MC RECO p_{T}",125,0,25);
  TH1F* hist6b = new TH1F("tau -> 3 muon 6b","MC RECO #eta",50,-3,3);
  TH1F* hist7b = new TH1F("tau -> 3 muon 7b","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist8b = new TH1F("tau -> 3 muon 8b","MC RECO p_{T}",125,0,25);
  TH1F* hist9b = new TH1F("tau -> 3 muon 9b","MC RECO #eta",50,-3,3);
  TH1F* hist10b = new TH1F("tau -> 3 muon 10b","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist2c = new TH1F("tau -> 3 muon 2c","MC RECO p_{T}",125,0,25);
  TH1F* hist3c = new TH1F("tau -> 3 muon 3c","MC RECO #eta",50,-3,3);
  TH1F* hist4c = new TH1F("tau -> 3 muon 4c","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist5c = new TH1F("tau -> 3 muon 5c","MC RECO p_{T}",125,0,25);
  TH1F* hist6c = new TH1F("tau -> 3 muon 6c","MC RECO #eta",50,-3,3);
  TH1F* hist7c = new TH1F("tau -> 3 muon 7c","MC RECO #varphi",25,-3.5,3.5);
  TH1F* hist8c = new TH1F("tau -> 3 muon 8c","MC RECO p_{T}",125,0,25);
  TH1F* hist9c = new TH1F("tau -> 3 muon 9c","MC RECO #eta",50,-3,3);
  TH1F* hist10c = new TH1F("tau -> 3 muon 10c","MC RECO #varphi",25,-3.5,3.5);
  UInt_t count1=0, count2=0, count3=0, count4=0, count5=0, count6=0;
  gStyle->SetOptStat(0);

  
  // load trigger menu
  const baconhep::TTrigger triggerMenu("../../BaconAna/DataFormats/data/HLT_50nsGRun");

  // load pileup reweighting file
  TFile *f_rw = TFile::Open("../Utils/data/puWeights_76x.root", "read");
  TH1D *h_rw = (TH1D*) f_rw->Get("puWeights");
  TH1D *h_rw_up = (TH1D*) f_rw->Get("puWeightsUp");
  TH1D *h_rw_down = (TH1D*) f_rw->Get("puWeightsDown");


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  

  vector<TString>  snamev;      // sample name (for output files)  
  vector<CSample*> samplev;     // data/MC samples

  // parse .conf file
  confParse(conf, snamev, samplev);
  const Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // Create output directory
  gSystem->mkdir(outputDir,kTRUE);
  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  
  // Data structures to store info from TTrees
  baconhep::TEventInfo *info  = new baconhep::TEventInfo();
  baconhep::TGenEventInfo *geninfo  = new baconhep::TGenEventInfo();
  TClonesArray *genPartArr = new TClonesArray("baconhep::TGenParticle");
  TClonesArray *vertexArr  = new TClonesArray("baconhep::TVertex");
  UInt_t* runnum = new UInt_t();
  UInt_t* evtnum = new UInt_t();
  std::vector<float> *muon_pt = new std::vector<float>();
  std::vector<float> *muon_eta = new std::vector<float>();
  std::vector<float> *muon_phi = new std::vector<float>();
  std::vector<float> *muon_ptErr = new std::vector<float>();
  std::vector<float> *muon_staPt = new std::vector<float>();
  std::vector<float> *muon_staEta = new std::vector<float>();
  std::vector<float> *muon_staPhi = new std::vector<float>();
  std::vector<float> *muon_pfPt = new std::vector<float>();
  std::vector<float> *muon_pfEta = new std::vector<float>();
  std::vector<float> *muon_pfPhi = new std::vector<float>();
  std::vector<int>   *muon_q = new std::vector<int>();
  std::vector<float> *muon_trkIso = new std::vector<float>();
  std::vector<float> *muon_ecalIso = new std::vector<float>();
  std::vector<float> *muon_hcalIso = new std::vector<float>();
  std::vector<float> *muon_chHadIso = new std::vector<float>();
  std::vector<float> *muon_gammaIso = new std::vector<float>();
  std::vector<float> *muon_neuHadIso = new std::vector<float>();
  std::vector<float> *muon_puIso = new std::vector<float>();
  std::vector<float> *muon_d0 = new std::vector<float>();
  std::vector<float> *muon_dz = new std::vector<float>();
  std::vector<float> *muon_sip3d = new std::vector<float>();
  std::vector<float> *muon_tkNchi2 = new std::vector<float>();
  std::vector<float> *muon_muNchi2 = new std::vector<float>();
  std::vector<float> *muon_trkKink = new std::vector<float>();
  std::vector<float> *muon_glbKink = new std::vector<float>();
  std::vector<int>   *muon_nValidHits = new std::vector<int>();
  std::vector<unsigned int>   *muon_typeBits = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_selectorBits = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_pogIDBits = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_nTkHits = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_nPixHits = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_nTkLayers = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_nPixLayers = new std::vector<unsigned int>();
  std::vector<unsigned int>   *muon_nMatchStn = new std::vector<unsigned int>();
  std::vector<int>   *muon_trkID = new std::vector<int>();
  std::vector<TriggerObjects>   *muon_hltMatchBits = new std::vector<std::bitset<256> >();
  std::vector<float> *vf_tC = new std::vector<float>();
  std::vector<float> *vf_dOF = new std::vector<float>();
  std::vector<float> *vf_nC = new std::vector<float>();
  std::vector<float> *vf_Prob = new std::vector<float>();
  std::vector<int>   *tri_category = new std::vector<int>();
  std::vector<int>   *vf_valid = new std::vector<int>();
  std::vector<float> *tri_invmass = new std::vector<float>();
  TFile *infile=0;
  TTree *eventTree=0;

  // Data structure for output files.
  float sysinvmass;

  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    Bool_t isData=kFALSE;
    if(isam==0 && !hasData) continue;
    else if (isam==0) isData=kTRUE;
    Bool_t isSignal = (snamev[isam].CompareTo("dstau",TString::kIgnoreCase)==0);
    CSample* samp = samplev[isam];

    //
    // Set up output ntuple
    TString outfilename = ntupDir + TString("/") + snamev[isam] + TString("_Tau3Mu_select.root");
    //if(isam!=0 && !doScaleCorr) outfilename = ntupDir + TString("/") + snamev[isam] + TString("_select.raw.root");
    TFile *outFile = new TFile(outfilename,"RECREATE"); 
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("sysinvmass", &sysinvmass,"sysinvmass/F");
    cout<<"begin loop over files"<<endl;

    // loop through files
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {  
      
      // Read input file and get the TTrees
      cout << "Processing " << samp->fnamev[ifile] << " [xsec = " << samp->xsecv[ifile] << " pb] ... "; cout.flush();
      infile = TFile::Open(samp->fnamev[ifile]); 
      assert(infile);

      // Access Event Tree
      eventTree = (TTree*)infile->Get("Events");
      assert(eventTree);
      eventTree->SetBranchAddress("runNum", &runnum);               TBranch *runNumBr = eventTree->GetBranch("runNum");
      eventTree->SetBranchAddress("evtNum", &evtnum);               TBranch *evtNumBr = eventTree->GetBranch("evtNum");
      eventTree->SetBranchAddress("Info", &info);                   TBranch *infoBr = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("GenEvtInfo", &geninfo);          TBranch *geninfoBr = eventTree->GetBranch("GenEvtInfo");
      eventTree->SetBranchAddress("PV", &vertexArr);                TBranch *vertexBr = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("GenParticle", &genPartArr);      TBranch *genPartBr = eventTree->GetBranch("GenParticle");
      eventTree->SetBranchAddress("MuonPt", &muon_pt);              TBranch *muonPtBr = eventTree->GetBranch("MuonPt");
      eventTree->SetBranchAddress("MuonEta", &muon_eta);            TBranch *muonEtaBr = eventTree->GetBranch("MuonEta");
      eventTree->SetBranchAddress("MuonPhi", &muon_phi);            TBranch *muonPhiBr = eventTree->GetBranch("MuonPhi");
      eventTree->SetBranchAddress("MuonPtErr", &muon_ptErr);        TBranch *muonPtErrBr = eventTree->GetBranch("MuonPtErr");
      eventTree->SetBranchAddress("MuonStaPt", &muon_staPt);        TBranch *muonStaPtBr = eventTree->GetBranch("MuonStaPt");
      eventTree->SetBranchAddress("MuonStaEta", &muon_staEta);      TBranch *muonStaEtaBr = eventTree->GetBranch("MuonStaEta");
      eventTree->SetBranchAddress("MuonStaPhi", &muon_staPhi);      TBranch *muonStaPhiBr = eventTree->GetBranch("MuonStaPhi");
      eventTree->SetBranchAddress("MuonPfPt", &muon_pfPt);          TBranch *muonPfPtBr = eventTree->GetBranch("MuonPfPt");
      eventTree->SetBranchAddress("MuonPfEta", &muon_pfEta);        TBranch *muonPfEtaBr = eventTree->GetBranch("MuonPfEta");
      eventTree->SetBranchAddress("MuonPfPhi", &muon_pfPhi);        TBranch *muonPfPhiBr = eventTree->GetBranch("MuonPfPhi");
      eventTree->SetBranchAddress("MuonQ", &muon_q);                TBranch *muonQBr = eventTree->GetBranch("MuonQ");
      eventTree->SetBranchAddress("MuonTrkIso", &muon_trkIso);      TBranch *muonTrkIsoBr = eventTree->GetBranch("MuonTrkIso");
      eventTree->SetBranchAddress("MuonEcalIso", &muon_ecalIso);    TBranch *muonEcalIsoBr = eventTree->GetBranch("MuonEcalIso");
      eventTree->SetBranchAddress("MuonHcalIso", &muon_hcalIso);    TBranch *muonHcalIsoBr = eventTree->GetBranch("MuonHcalIso");
      eventTree->SetBranchAddress("MuonChHadIso", &muon_chHadIso);  TBranch *muonChHadIsoBr = eventTree->GetBranch("MuonChHadIso");
      eventTree->SetBranchAddress("MuonGammaIso", &muon_gammaIso);  TBranch *muonGammaIsoBr = eventTree->GetBranch("MuonGammaIso");
      eventTree->SetBranchAddress("MuonNeuHadIso", &muon_neuHadIso);   TBranch *muonNeuHadIsoBr = eventTree->GetBranch("MuonNeuHadIso");
      eventTree->SetBranchAddress("MuonPuIso", &muon_puIso);           TBranch *muonPuIsoBr = eventTree->GetBranch("MuonPuIso");
      eventTree->SetBranchAddress("MuonD0", &muon_d0);                 TBranch *muonD0Br = eventTree->GetBranch("MuonD0");
      eventTree->SetBranchAddress("MuonDz", &muon_dz);                 TBranch *muonDzBr = eventTree->GetBranch("MuonDz");
      eventTree->SetBranchAddress("MuonSip3d", &muon_sip3d);           TBranch *muonSip3dBr = eventTree->GetBranch("MuonSip3d");
      eventTree->SetBranchAddress("MuonTkNchi2", &muon_tkNchi2);       TBranch *muonTkNchi2Br = eventTree->GetBranch("MuonTkNchi2");
      eventTree->SetBranchAddress("MuonMuNchi2", &muon_muNchi2);       TBranch *muonMuNchi2Br = eventTree->GetBranch("MuonMuNchi2");
      eventTree->SetBranchAddress("MuonTrkKink", &muon_trkKink);       TBranch *muonTrkKinkBr = eventTree->GetBranch("MuonTrkKink");
      eventTree->SetBranchAddress("MuonGlbKink", &muon_glbKink);       TBranch *muonGlbKinkBr = eventTree->GetBranch("MuonGlbKink");
      eventTree->SetBranchAddress("MuonNValidHits", &muon_nValidHits); TBranch *muonNValidHitsBr = eventTree->GetBranch("MuonNValidHits");
      eventTree->SetBranchAddress("MuonTypeBits", &muon_typeBits);     TBranch *muonTypeBitsBr = eventTree->GetBranch("MuonTypeBits");
      eventTree->SetBranchAddress("MuonSelectorBits", &muon_selectorBits);    TBranch *muonSelectorBitsBr = eventTree->GetBranch("MuonSelectorBits");
      eventTree->SetBranchAddress("MuonPogIDBits", &muon_pogIDBits);          TBranch *muonPogIDBitsBr = eventTree->GetBranch("MuonPogIDBits");
      eventTree->SetBranchAddress("MuonNTkHits", &muon_nTkHits);              TBranch *muonNTkHitsBr = eventTree->GetBranch("MuonNTkHits");
      eventTree->SetBranchAddress("MuonNPixHits", &muon_nPixHits);            TBranch *muonNPixHitsBr = eventTree->GetBranch("MuonNPixHits");
      eventTree->SetBranchAddress("MuonNTkLayers", &muon_nTkLayers);          TBranch *muonNTkLayersBr = eventTree->GetBranch("MuonNTkLayers");
      eventTree->SetBranchAddress("MuonNPixLayers", &muon_nPixLayers);        TBranch *muonNPixLayersBr = eventTree->GetBranch("MuonNPixLayers");
      eventTree->SetBranchAddress("MuonNMatchStn", &muon_nMatchStn);          TBranch *muonNMatchStnBr = eventTree->GetBranch("MuonNMatchStn");
      eventTree->SetBranchAddress("MuonTrkID", &muon_trkID);                  TBranch *muonTrkIDBr = eventTree->GetBranch("MuonTrkID");
      eventTree->SetBranchAddress("MuonHltMatchBits", &muon_hltMatchBits);    TBranch *muonHltMatchBitsBr = eventTree->GetBranch("MuonHltMatchBits");
      eventTree->SetBranchAddress("VfTc", &vf_tC);                            TBranch *vfTcBr = eventTree->GetBranch("VfTc");
      eventTree->SetBranchAddress("VfDof", &vf_dOF);                          TBranch *vfDofBr = eventTree->GetBranch("VfDof");
      eventTree->SetBranchAddress("VfNc", &vf_nC);                            TBranch *vfNcBr = eventTree->GetBranch("VfNc");
      eventTree->SetBranchAddress("VfProb", &vf_Prob);                        TBranch *vfProbBr = eventTree->GetBranch("VfProb");
      eventTree->SetBranchAddress("Category", &tri_category);                 TBranch *triCategoryBr = eventTree->GetBranch("Category");
      eventTree->SetBranchAddress("VfValid", &vf_valid);                      TBranch *vfValidBr = eventTree->GetBranch("VfValid");
      eventTree->SetBranchAddress("InvMass", &tri_invmass);                   TBranch *triInvMassBr = eventTree->GetBranch("InvMass");
      Bool_t hasGen = eventTree->GetBranchStatus("GenEvtInfo");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<10000; ientry++) {
        if(ientry%500==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

	// Get Events
	infoBr->GetEntry(ientry);        geninfoBr->GetEntry(ientry);
	genPartArr->Clear();             genPartBr->GetEntry(ientry);
	vertexArr->Clear();              vertexBr->GetEntry(ientry);
	runNumBr->GetEntry(ientry);
	evtNumBr->GetEntry(ientry);
	muon_pt->clear();                muonPtBr->GetEntry(ientry);
	muon_eta->clear();               muonEtaBr->GetEntry(ientry);
	muon_phi->clear();               muonPhiBr->GetEntry(ientry);
	muon_ptErr->clear();             muonPtErrBr->GetEntry(ientry);
	muon_staPt->clear();             muonStaPtBr->GetEntry(ientry);
	muon_staEta->clear();            muonStaEtaBr->GetEntry(ientry);
	muon_staPhi->clear();            muonStaPhiBr->GetEntry(ientry);
	muon_pfPt->clear();              muonPfPtBr->GetEntry(ientry);
	muon_pfEta->clear();             muonPfEtaBr->GetEntry(ientry);
	muon_pfPhi->clear();             muonPfPhiBr->GetEntry(ientry);
	muon_q->clear();                 muonQBr->GetEntry(ientry);
	muon_trkIso->clear();            muonTrkIsoBr->GetEntry(ientry);
	muon_ecalIso->clear();           muonEcalIsoBr->GetEntry(ientry);
	muon_hcalIso->clear();           muonHcalIsoBr->GetEntry(ientry);
	muon_chHadIso->clear();          muonChHadIsoBr->GetEntry(ientry);
	muon_gammaIso->clear();          muonGammaIsoBr->GetEntry(ientry);
	muon_neuHadIso->clear();         muonNeuHadIsoBr->GetEntry(ientry);
	muon_puIso->clear();             muonPuIsoBr->GetEntry(ientry);
	muon_d0->clear();                muonD0Br->GetEntry(ientry);
	muon_dz->clear();                muonDzBr->GetEntry(ientry);
	muon_sip3d->clear();             muonSip3dBr->GetEntry(ientry);
	muon_tkNchi2->clear();           muonTkNchi2Br->GetEntry(ientry);
	muon_muNchi2->clear();           muonMuNchi2Br->GetEntry(ientry);
	muon_trkKink->clear();           muonTrkKinkBr->GetEntry(ientry);
	muon_glbKink->clear();           muonGlbKinkBr->GetEntry(ientry);
	muon_nValidHits->clear();        muonNValidHitsBr->GetEntry(ientry);
	muon_typeBits->clear();          muonTypeBitsBr->GetEntry(ientry);
	muon_selectorBits->clear();      muonSelectorBitsBr->GetEntry(ientry);
	muon_pogIDBits->clear();         muonPogIDBitsBr->GetEntry(ientry);
	muon_nTkHits->clear();           muonNTkHitsBr->GetEntry(ientry);
	muon_nPixHits->clear();          muonNPixHitsBr->GetEntry(ientry);
	muon_nTkLayers->clear();         muonNTkLayersBr->GetEntry(ientry);
	muon_nPixLayers->clear();        muonNPixLayersBr->GetEntry(ientry);
	muon_nMatchStn->clear();         muonNMatchStnBr->GetEntry(ientry);
	muon_trkID->clear();             muonTrkIDBr->GetEntry(ientry);
	muon_hltMatchBits->clear();      muonHltMatchBitsBr->GetEntry(ientry);
	vf_tC->clear();                  vfTcBr->GetEntry(ientry);
	vf_dOF->clear();                 vfDofBr->GetEntry(ientry);
	vf_nC->clear();                  vfNcBr->GetEntry(ientry);
	vf_Prob->clear();                vfProbBr->GetEntry(ientry);
	vf_valid->clear();               vfValidBr->GetEntry(ientry);
	tri_category->clear();           triCategoryBr->GetEntry(ientry);
	tri_invmass->clear();            triInvMassBr->GetEntry(ientry);

	std::vector<float> pt, eta, phi, ptErr, staPt, staEta, staPhi, pfPt, pfEta, pfPhi;
	std::vector<float> trkIso, ecalIso, hcalIso, chHadIso, gammaIso, neuHadIso, puIso, d0, dz, sip3d;
	std::vector<float> tkNchi2, muNchi2, trkKink, glbKink;
	std::vector<int> q, nValidHits, trkID, category, Valid;
	std::vector<unsigned int> typeBits, selectorBits, pogIDBits, nTkHits, nPixHits, nTkLayers, nPixLayers, nMatchStn;
	std::vector<float> tC, dOF, nC, Prob, invMass;
	std::vector<TriggerObjects> hltMatchBits;

        for(vector<float>::iterator it = muon_pt->begin(); it != muon_pt->end(); it++)
	  pt.push_back(*it);
	for(vector<float>::iterator it = muon_eta->begin(); it != muon_eta->end(); it++)
	  eta.push_back(*it);
	for(vector<float>::iterator it = muon_phi->begin(); it != muon_phi->end(); it++)
	  phi.push_back(*it);
	for(vector<float>::iterator it = muon_ptErr->begin(); it != muon_ptErr->end(); it++)
	  ptErr.push_back(*it);
	for(vector<float>::iterator it = muon_staPt->begin(); it != muon_staPt->end(); it++)
	  staPt.push_back(*it);
	for(vector<float>::iterator it = muon_staEta->begin(); it != muon_staEta->end(); it++)
	  staEta.push_back(*it);
	for(vector<float>::iterator it = muon_staPhi->begin(); it != muon_staPhi->end(); it++)
	  staPhi.push_back(*it);
	for(vector<float>::iterator it = muon_pfPt->begin(); it != muon_pfPt->end(); it++)
	  pfPt.push_back(*it);
	for(vector<float>::iterator it = muon_pfEta->begin(); it != muon_pfEta->end(); it++)
	  pfEta.push_back(*it);
	for(vector<float>::iterator it = muon_pfPhi->begin(); it != muon_pfPhi->end(); it++)
	  pfPhi.push_back(*it);
	for(vector<float>::iterator it = muon_trkIso->begin(); it != muon_trkIso->end(); it++)
	  trkIso.push_back(*it);
	for(vector<float>::iterator it = muon_ecalIso->begin(); it != muon_ecalIso->end(); it++)
	  ecalIso.push_back(*it);
	for(vector<float>::iterator it = muon_hcalIso->begin(); it != muon_hcalIso->end(); it++)
	  hcalIso.push_back(*it);
	for(vector<float>::iterator it = muon_chHadIso->begin(); it != muon_chHadIso->end(); it++)
	  chHadIso.push_back(*it);
	for(vector<float>::iterator it = muon_gammaIso->begin(); it != muon_gammaIso->end(); it++)
	  gammaIso.push_back(*it);
	for(vector<float>::iterator it = muon_neuHadIso->begin(); it != muon_neuHadIso->end(); it++)
	  neuHadIso.push_back(*it);
	for(vector<float>::iterator it = muon_puIso->begin(); it != muon_puIso->end(); it++)
	  puIso.push_back(*it);
	for(vector<float>::iterator it = muon_d0->begin(); it != muon_d0->end(); it++)
	  d0.push_back(*it);
	for(vector<float>::iterator it = muon_dz->begin(); it != muon_dz->end(); it++)
	  dz.push_back(*it);
	for(vector<float>::iterator it = muon_sip3d->begin(); it != muon_sip3d->end(); it++)
	  sip3d.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_typeBits->begin(); it != muon_typeBits->end(); it++)
	  typeBits.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_selectorBits->begin(); it != muon_selectorBits->end(); it++)
	  selectorBits.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_pogIDBits->begin(); it != muon_pogIDBits->end(); it++)
	  pogIDBits.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_nTkHits->begin(); it != muon_nTkHits->end(); it++)
	  nTkHits.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_nPixHits->begin(); it != muon_nPixHits->end(); it++)
	  nPixHits.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_nTkLayers->begin(); it != muon_nTkLayers->end(); it++)
	  nTkLayers.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_nPixLayers->begin(); it != muon_nPixLayers->end(); it++)
	  nPixLayers.push_back(*it);
	for(vector<unsigned int>::iterator it = muon_nMatchStn->begin(); it != muon_nMatchStn->end(); it++)
	  nMatchStn.push_back(*it);
	for(vector<int>::iterator it = muon_q->begin(); it != muon_q->end(); it++)
	  q.push_back(*it);
	for(vector<int>::iterator it = muon_nValidHits->begin(); it != muon_nValidHits->end(); it++)
	  nValidHits.push_back(*it);
	for(vector<int>::iterator it = muon_trkID->begin(); it != muon_trkID->end(); it++)
	  trkID.push_back(*it);
	for(vector<TriggerObjects>::iterator it = muon_hltMatchBits->begin(); it != muon_hltMatchBits->end(); it++)
	  hltMatchBits.push_back(*it);	
	for(vector<float>::iterator it = muon_tkNchi2->begin(); it != muon_tkNchi2->end(); it++)
	  tkNchi2.push_back(*it);
	for(vector<float>::iterator it = muon_muNchi2->begin(); it != muon_muNchi2->end(); it++)
	  muNchi2.push_back(*it);
	for(vector<float>::iterator it = muon_trkKink->begin(); it != muon_trkKink->end(); it++)
	  trkKink.push_back(*it);
	for(vector<float>::iterator it = muon_glbKink->begin(); it != muon_glbKink->end(); it++)
	  glbKink.push_back(*it);
	for(vector<float>::iterator it = vf_tC->begin(); it != vf_tC->end(); it++)
	  tC.push_back(*it);
	for(vector<float>::iterator it = vf_dOF->begin(); it != vf_dOF->end(); it++)
	  dOF.push_back(*it);
	for(vector<float>::iterator it = vf_nC->begin(); it != vf_nC->end(); it++)
	  nC.push_back(*it);
	for(vector<float>::iterator it = vf_Prob->begin(); it != vf_Prob->end(); it++)
	  Prob.push_back(*it);
	for(vector<int>::iterator it = tri_category->begin(); it != tri_category->end(); it++)
	  category.push_back(*it);
	for(vector<float>::iterator it = tri_invmass->begin(); it != tri_invmass->end(); it++)
	  invMass.push_back(*it);
	for(vector<int>::iterator it = vf_valid->begin(); it != vf_valid->end(); it++)
	  Valid.push_back(*it);

	//Store GEN level tau muon
	vector<baconhep::TGenParticle*> genmuonArr;
	for(int i=0; i<genPartArr->GetEntries(); i++){
	  baconhep::TGenParticle *genpar = (baconhep::TGenParticle*)((*genPartArr)[i]);
	  if(genpar->pdgId != 13 && genpar->pdgId != -13) continue;
	  if(genpar->status != 1) continue;
	  Int_t parentid1=dynamic_cast<baconhep::TGenParticle *>(genPartArr->At(genpar->parent>-1 ? genpar->parent : 0))->pdgId;
	  if(parentid1 != 15 && parentid1 != -15) continue;
	  genmuonArr.push_back(genpar);
	}

	//Only study events contain single decay 229041
	if(genmuonArr.size() > 3) continue;
	count1++;
	
	TLorentzVector temp1,temp2,temp3;
	float minmass = 0;
	int maxmassnum = -99;
	for(int i=0; i<category.size(); i++){
	  temp1.SetPtEtaPhiM(pt[i*3],eta[i*3],phi[i*3],0.105658369);
	  temp2.SetPtEtaPhiM(pt[i*3+1],eta[i*3+1],phi[i*3+1],0.105658369);
	  temp3.SetPtEtaPhiM(pt[i*3+2],eta[i*3+2],phi[i*3+2],0.105658369);//0.13957
	  float invmass_temp = (temp1+temp2+temp3).M();

	  if(category[i] != 1) continue;

	  // Mass cuts
	  if(q[i*3] * q[i*3+1] < 0){
	    double mass_temp = (temp1+temp2).M();
	    if(mass_temp < 0.45) continue;
	    if(abs(mass_temp - 0.782) < 0.02) continue;
	  }
	  if(q[i*3] * q[i*3+2] < 0){
	    double mass_temp = (temp1+temp3).M();
	    if(mass_temp < 0.45) continue;
	    if(abs(mass_temp - 0.782) < 0.02) continue;
	  }
	  if(q[i*3+1] * q[i*3+2] < 0){
	    double mass_temp = (temp2+temp3).M();
	    if(mass_temp < 0.45) continue;
	    if(abs(mass_temp - 0.782) < 0.02) continue;
	  }

	  if(!Valid[i]) continue;
	  
	  bool MuonTriObj1 = (triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltDoubleMu3TrkTau3muL3Filtered",hltMatchBits[i*3]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltL1fL1sL1DoubleMuorTripleMuL1Filtered0",hltMatchBits[i*3]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltL2fL1sL1DoubleMuorTripleMuL1f0L2PreFiltered0",hltMatchBits[i*3]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltTau3muTkVertexFilter",hltMatchBits[i*3]));
	  bool MuonTriObj2 = (triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltDoubleMu3TrkTau3muL3Filtered",hltMatchBits[i*3+1]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltL1fL1sL1DoubleMuorTripleMuL1Filtered0",hltMatchBits[i*3+1]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltL2fL1sL1DoubleMuorTripleMuL1f0L2PreFiltered0",hltMatchBits[i*3+1]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltTau3muTkVertexFilter",hltMatchBits[i*3+1]));
	  bool MuonTriObj3 = (triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltDoubleMu3TrkTau3muL3Filtered",hltMatchBits[i*3+2]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltL1fL1sL1DoubleMuorTripleMuL1Filtered0",hltMatchBits[i*3+2]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltL2fL1sL1DoubleMuorTripleMuL1f0L2PreFiltered0",hltMatchBits[i*3+2]) ||
			      triggerMenu.passObj("HLT_DoubleMu3_Trk_Tau3mu_v*","hltTau3muTkVertexFilter",hltMatchBits[i*3+2]));
	  if(ientry<50) cout<<MuonTriObj1<<MuonTriObj2<<MuonTriObj3<<endl;
	  if(!MuonTriObj1) continue;
	  if(!MuonTriObj2) continue;
	  if(!MuonTriObj3) continue;
	  
	  if(invmass_temp > minmass){
	    minmass = invmass_temp;
	    maxmassnum = i;
	  }
	}
	if(maxmassnum < 0) continue;

	// Fill histogram
	// Mass
	TLorentzVector v1,v2,v3;
	v1.SetPtEtaPhiM(pt[maxmassnum*3],eta[maxmassnum*3],phi[maxmassnum*3],0.105658369);
	v2.SetPtEtaPhiM(pt[maxmassnum*3+1],eta[maxmassnum*3+1],phi[maxmassnum*3+1],0.105658369);
	v3.SetPtEtaPhiM(pt[maxmassnum*3+2],eta[maxmassnum*3+2],phi[maxmassnum*3+2],0.13957);
	float invmass = (v1+v2+v3).M();
	hist0->Fill(invmass);
	count2++;
	Int_t B = 0;
	Int_t M = 0;
	Int_t E = 0;
	if (fabs(eta[maxmassnum*3])<0.8) B++;
	if (fabs(eta[maxmassnum*3+1])<0.8) B++;
	if (fabs(eta[maxmassnum*3+2])<0.8) B++;
	if (fabs(eta[maxmassnum*3])>=0.8 && fabs(eta[maxmassnum*3])<1.6) M++;
	if (fabs(eta[maxmassnum*3+1])>=0.8 && fabs(eta[maxmassnum*3+1])<1.6) M++;
	if (fabs(eta[maxmassnum*3+2])>=0.8 && fabs(eta[maxmassnum*3+2])<1.6) M++;
	if (fabs(eta[maxmassnum*3])>=1.6) E++;
	if (fabs(eta[maxmassnum*3+1])>=1.6) E++;
	if (fabs(eta[maxmassnum*3+2])>=1.6) E++;
	if (B==3){
	  count3++;
	  hist1a->Fill(invmass);}//store events come from tau
	else if(B==2&&M==1){
	  count3++;
	  hist1a->Fill(invmass);}
	else if(B==2&&E==1){
	  count3++;
	  hist1a->Fill(invmass);}
	else if(M==3){
	  count4++;
	  hist1b->Fill(invmass);}
	else if(M==2&&B==1){
	  count4++;
	  hist1b->Fill(invmass);}
	else if(M==2&&E==1){
	  count4++;
	  hist1b->Fill(invmass);}
	else if(E==3){
	  count5++;
	  hist1c->Fill(invmass);}
	else if(E==2&&B==1){
	  count5++;
	  hist1c->Fill(invmass);}
	else if(E==2&&M==1){
	  count5++;
	  hist1c->Fill(invmass);}
	else if(B==1&&M==1&&E==1){
	  hist1d->Fill(invmass);
	  count6++;}	      
	else{
	  cout<<B<<M<<E<<endl;
	}

	// Kinematics
	//Sort muon array
	int sort[3]={-99,-99,-99};
	Double_t maxpt=-99, submaxpt=-99;
	for(int l=0; l<3; l++){
	  if(pt[maxmassnum*3+l] >= maxpt){
	    submaxpt = maxpt;
	    maxpt = pt[maxmassnum*3+l];
	    sort[2] = sort[1];
	    sort[1] = sort[0];
	    sort[0] = maxmassnum*3+l;
	  }
	  else if(pt[maxmassnum*3+l] < maxpt && pt[maxmassnum*3+l] >= submaxpt){
	    submaxpt = pt[maxmassnum*3+l];
	    sort[2] = sort[1];
	    sort[1] = maxmassnum*3+l;
	  }
	  else{
	    sort[2] = maxmassnum*3+l;
	  }
	}
	hist2->Fill(pt[sort[0]]); hist5->Fill(pt[sort[1]]); hist8->Fill(pt[sort[2]]);
	hist3->Fill(eta[sort[0]]); hist6->Fill(eta[sort[1]]); hist9->Fill(eta[sort[2]]);
	hist4->Fill(phi[sort[0]]); hist7->Fill(phi[sort[1]]); hist10->Fill(phi[sort[2]]);\

	// Fill output tree
	sysinvmass = invmass;
	outTree->Fill();
      }//end of event loop
      
      hist0->Scale(1/hist0->GetEntries());
      hist1a->Scale(1/hist1a->GetEntries());
      hist1b->Scale(1/hist1b->GetEntries());
      hist1c->Scale(1/hist1c->GetEntries());
      hist1d->Scale(1/hist1d->GetEntries());
      hist2->Scale(1/hist2->GetEntries());
      hist3->Scale(1/hist3->GetEntries());
      hist4->Scale(1/hist4->GetEntries());
      hist5->Scale(1/hist5->GetEntries());
      hist6->Scale(1/hist6->GetEntries());
      hist7->Scale(1/hist7->GetEntries());
      hist8->Scale(1/hist8->GetEntries());
      hist9->Scale(1/hist9->GetEntries());
      hist10->Scale(1/hist10->GetEntries());
      
      delete infile;
      infile=0, eventTree=0;    
    }
    outFile->Write();
    outFile->Close();
  }
  delete h_rw;
  delete h_rw_up;
  delete h_rw_down;
  delete f_rw;
  delete info;
  delete geninfo;
  delete genPartArr;
  delete vertexArr;
  
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
  cout<<count1<<" "<<count2<<" "<<count3<<" "<<count4<<" "<<count5<<" "<<count6<<endl;

  //THStack *all = new THStack("tau -> 3 muons","invariant mass");
  TCanvas *c0 = new TCanvas("c0","invariant mass",1200,900);
  TAxis *xaxis = hist0->GetXaxis();
  TAxis *yaxis = hist0->GetYaxis();
  xaxis->SetTitle("m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}");
  yaxis->SetTitle("a.u. / 10 MeV");
  yaxis->SetTitleOffset(1.3);
  c0->cd();

  hist0->SetFillColor(40);
  //hist0->SetFillStyle(0);
  //hist1->SetFillColor(7);
  //hist3->SetFillColor(7);
  //hist4->SetFillColor(8);


  hist0->Draw();
  //hist1->Draw("SAME");
  
  // all->Add(hist1);
  //all->Add(hist2);
  //all->Add(hist3);
  //all->Add(hist4);
  //c0->cd();
  //c0->SetLogy();
  /*
    all->Draw();
    TAxis *xaxis = all->GetXaxis();
    TAxis *yaxis = all->GetYaxis();
    xaxis->SetTitle("Invariant mass (GeV)");
    yaxis->SetTitle("Entries");
    yaxis->SetTitleOffset(1.2);
    all->SetMinimum(8.);
    all->SetMaximum(120000.);
  */
  
  TF1 *total = new TF1("mstotal","gaus(0)+gaus(3)",1.7,1.85);
  Double_t par[6]={0.02,1.75,0.02,0.01,1.75,0.1};
  //hist0->Fit(f1,"R0");
  //hist2->Fit(f2,"R0+");
  //f1->GetParameters(&par[0]);
  //f2->GetParameters(&par[3]);
  total->SetParameters(par);
  hist0->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  total->Draw("SAME");
  //f1->SetLineColor(3);
  //f1->Draw("SAME");

  //auto legend = new TLegend(0.5,0.7,0.7,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  //legend->AddEntry(hist0,"Tau -> 3 Mu MC","f");
  //legend->Draw();
  

  c0->Print("invariant mass.png");

  TCanvas *c0a = new TCanvas("c0a","invariant mass",1200,900);
  xaxis = hist1a->GetXaxis();
  yaxis = hist1a->GetYaxis();
  xaxis->SetTitle("m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}");
  yaxis->SetTitle("a.u. / 10 MeV");
  yaxis->SetTitleOffset(1.3);
  c0a->cd();
  hist1a->SetFillColor(30);
  //hist1a->SetFillStyle(0);
  hist1a->Draw();
  Double_t par1[6]={10,1.7,0.02,1,1.7,0.1};
  total->SetParameters(par1);
  hist1a->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  //total->Draw("SAME");
  c0a->Print("invariant massa.png");

  TCanvas *c0b = new TCanvas("c0b","invariant mass",1200,900);
  xaxis = hist1b->GetXaxis();
  yaxis = hist1b->GetYaxis();
  xaxis->SetTitle("m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}");
  yaxis->SetTitle("a.u. / 10 MeV");
  yaxis->SetTitleOffset(1.3);
  c0b->cd();
  hist1b->SetFillColor(40);
  //hist1b->SetFillStyle(0);
  hist1b->Draw();
  Double_t par2[6]={8,1.7,0.02,1,1.7,0.1};
  total->SetParameters(par2);
  hist1b->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  //total->Draw("SAME");
  c0b->Print("invariant massb.png");

  TCanvas *c0c = new TCanvas("c0c","invariant mass",1200,900);
  xaxis = hist1c->GetXaxis();
  yaxis = hist1c->GetYaxis();
  xaxis->SetTitle("m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}");
  yaxis->SetTitle("a.u. / 10 MeV");
  yaxis->SetTitleOffset(1.3);
  c0c->cd();
  hist1c->SetFillColor(41);
  //hist1c->SetFillStyle(0);
  hist1c->Draw();
  Double_t par3[6]={5,1.7,0.02,1,1.7,0.1};
  total->SetParameters(par3);
  hist1c->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  //total->Draw("SAME");
  c0c->Print("invariant massc.png");

  TCanvas *c0d = new TCanvas("c0d","invariant mass",1200,900);
  cout<<hist1d->GetEntries()<<endl;
  xaxis = hist1d->GetXaxis();
  yaxis = hist1d->GetYaxis();
  xaxis->SetTitle("m_{#mu^{#pm}#mu^{#pm}#mu^{#pm}}");
  yaxis->SetTitle("a.u. / 1 GeV");
  yaxis->SetTitleOffset(1.3);
  c0d->cd();
  hist1d->SetFillColor(46);
  //hist1d->SetFillStyle(0);
  hist1d->Draw();
  total->SetParameters(par);
  hist1d->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  //total->Draw("SAME");
  c0d->Print("invariant massd.png");


  TCanvas *c1 = new TCanvas("1","muon pT",1200,900);
  xaxis = hist5->GetXaxis();
  yaxis = hist5->GetYaxis();
  xaxis->SetTitle("p_{T} (GeV)");
  yaxis->SetTitle("a.u. / 0.2 GeV");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0,1);
  c1->cd();
  hist2->SetLineColor(4);
  hist5->SetLineColor(4);
  hist8->SetLineColor(4);
  hist2->SetLineStyle(1);
  hist5->SetLineStyle(2);
  hist8->SetLineStyle(5);
  hist2->SetFillStyle(0);
  hist5->SetFillStyle(0);
  hist8->SetFillStyle(0);
  hist2->Draw();
  hist5->Draw("SAME");
  hist8->Draw("SAME");
  //legend = new TLegend(0.35,0.75,0.6,0.85);
  //legend->AddEntry(hist2,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist5,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist8,"Tau -> 3 Mu MC","f");
  //legend->Draw();
  c1->Print("pt.png");

  TCanvas *c2 = new TCanvas("2","muon eta",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();
  xaxis->SetTitle("#eta");
  yaxis->SetTitle("a.u.");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0,0.04);
  c2->cd();
  hist3->SetLineColor(4);
  hist6->SetLineColor(4);
  hist9->SetLineColor(4);
  hist3->SetLineStyle(1);
  hist6->SetLineStyle(2);
  hist9->SetLineStyle(5);
  hist3->SetFillStyle(0);
  hist6->SetFillStyle(0);
  hist9->SetFillStyle(0);
  hist3->Draw();
  hist6->Draw("SAME");
  hist9->Draw("SAME");
  // legend = new TLegend(0.15,0.75,0.4,0.85);
  //legend->AddEntry(hist3,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist6,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist9,"Tau -> 3 Mu MC","f");
  //legend->Draw();
  c2->Print("eta.png");

  TCanvas *c3 = new TCanvas("3","muon phi",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();
  xaxis->SetTitle("#varphi");
  yaxis->SetTitle("a.u. / 0.28 rad");
  yaxis->SetTitleOffset(1.3);
  yaxis->SetRangeUser(0,0.06);
  c3->cd();
  hist4->SetLineColor(4);
  hist7->SetLineColor(4);
  hist10->SetLineColor(4);
  hist4->SetLineStyle(1);
  hist7->SetLineStyle(2);
  hist10->SetLineStyle(5);
  hist4->SetFillStyle(0);
  hist7->SetFillStyle(0);
  hist10->SetFillStyle(0);
  hist4->Draw();
  hist7->Draw("SAME");
  hist10->Draw("SAME");
  //legend = new TLegend(0.15,0.75,0.4,0.85);
  //legend->AddEntry(hist4,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist7,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist10,"Tau -> 3 Mu MC","f");
  //legend->Draw();

  c3->Print("phi.png");

  gBenchmark->Show("select3Mu");

}
