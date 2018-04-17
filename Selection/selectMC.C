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

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TH1F* hist0 = new TH1F("tau -> 3 muon 0","GEN m_{#mu^{#pm}#mu^{pm}#mu^{#pm}}",50,0,5);
  TH1F* hist2 = new TH1F("tau -> 3 muon 2","GEN p_{T}",75,0,15);
  TH1F* hist3 = new TH1F("tau -> 3 muon 3","GEN #eta",50,-3,3);
  TH1F* hist4 = new TH1F("tau -> 3 muon 4","GEN #varphi",25,-3.5,3.5);
  TH1F* hist5 = new TH1F("tau -> 3 muon 5","GEN p_{T}",75,0,15);
  TH1F* hist6 = new TH1F("tau -> 3 muon 6","GEN #eta",50,-3,3);
  TH1F* hist7 = new TH1F("tau -> 3 muon 7","GEN #varphi",25,-3.5,3.5);
  TH1F* hist8 = new TH1F("tau -> 3 muon 8","GEN p_{T}",75,0,15);
  TH1F* hist9 = new TH1F("tau -> 3 muon 9","GEN #eta",50,-3,3);
  TH1F* hist10 = new TH1F("tau -> 3 muon 10","GEN #varphi",25,-3.5,3.5);
  UInt_t count1=0, count2=0, count3=0;
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
  //TClonesArray *muon_pt  = new TClonesArray(int);
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
  std::vector<int>   *muon_hltMatchBits = new std::vector<int>();
  std::vector<float> *vf_tC = new std::vector<float>();
  std::vector<float> *vf_dOF = new std::vector<float>();
  std::vector<float> *vf_nC = new std::vector<float>();
  std::vector<float> *vf_Prob = new std::vector<float>();
  std::vector<int>   *tri_category = new std::vector<int>();
  TFile *infile=0;
  TTree *eventTree=0;

  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    Bool_t isData=kFALSE;
    if(isam==0 && !hasData) continue;
    else if (isam==0) isData=kTRUE;
    Bool_t isSignal = (snamev[isam].CompareTo("dstau",TString::kIgnoreCase)==0);
    CSample* samp = samplev[isam];
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
      Bool_t hasGen = eventTree->GetBranchStatus("GenEvtInfo");

      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
        if(ientry%20000==0) cout << "Processing event " << ientry << ". " << (double)ientry/(double)eventTree->GetEntries()*100 << " percent done with this file." << endl;

	infoBr->GetEntry(ientry);        geninfoBr->GetEntry(ientry);
	genPartArr->Clear();             genPartBr->GetEntry(ientry);
	vertexArr->Clear();              vertexBr->GetEntry(ientry);
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
	tri_category->clear();           triCategoryBr->GetEntry(ientry);

	cout<<muon_pt->size()<<endl;
	cout<<vf_Prob->size()<<endl;
	cout<<tri_category->size()<<endl;
	cout<<endl;
      }//end of event loop
      
      hist0->Scale(1/hist0->GetEntries());
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
  
  //THStack *all = new THStack("tau -> 3 muons","invariant mass");
  TCanvas *c0 = new TCanvas("c0","invariant mass",1200,900);
  TAxis *xaxis = hist0->GetXaxis();
  TAxis *yaxis = hist0->GetYaxis();

  xaxis->SetTitle("m_{#mu^{#pm}#mu^{pm}#mu^{#pm}} (GeV)");
  yaxis->SetTitle("a.u. / 0.1 GeV");
  c0->cd();
  

  hist0->SetFillColor(5);
  hist0->SetFillStyle(0);
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
  TF1 *f1 = new TF1("m1","gaus",1.5,2);
  TF1 *f2 = new TF1("m2","gaus",0,5);
  TF1 *total = new TF1("mstotal","gaus(0)+gaus(3)",0,5);
  Double_t par[6]={50,1.5,0.1,5,2.5,2};
  hist0->Fit(f1,"R0");
  //hist2->Fit(f2,"R0+");
  //f1->GetParameters(&par[0]);
  //f2->GetParameters(&par[3]);
  total->SetParameters(par);
  hist0->Fit(total,"R0+");
  total->SetLineColor(2);
  total->SetLineWidth(2);
  //total->Draw("SAME");
  //f1->SetLineColor(3);
  //f1->Draw("SAME");

  auto legend = new TLegend(0.5,0.7,0.7,0.8);
  //legend->AddEntry(hist1,"3 muons with opposite signs","f");
  legend->AddEntry(hist0,"Tau -> 3 Mu MC","f");
  legend->Draw();
  

  c0->Print("invariant mass.png");
  cout<<count1<<" "<<count2<<" "<<count3<<endl;
  
  TCanvas *c1 = new TCanvas("1","muon pT",1200,900);
  xaxis = hist2->GetXaxis();
  yaxis = hist2->GetYaxis();
  xaxis->SetTitle("p_{T} (GeV)");
  yaxis->SetTitle("a.u / 0.2 GeV");
  c1->cd();
  hist2->SetLineColor(2);
  hist5->SetLineColor(6);
  hist8->SetLineColor(4);
  hist2->SetFillStyle(0);
  hist5->SetFillStyle(0);
  hist8->SetFillStyle(0);
  hist2->Draw();
  //hist5->Draw("SAME");
  //hist8->Draw("SAME");
  legend = new TLegend(0.35,0.75,0.6,0.85);
  legend->AddEntry(hist2,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist5,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist8,"Tau -> 3 Mu MC","f");
  //legend->Draw();
  c1->Print("pt.png");

  TCanvas *c2 = new TCanvas("2","muon eta",1200,900);
  xaxis = hist3->GetXaxis();
  yaxis = hist3->GetYaxis();
  xaxis->SetTitle("#eta");
  yaxis->SetTitle("a.u");
  c2->cd();
  hist3->SetLineColor(2);
  hist6->SetLineColor(6);
  hist9->SetLineColor(4);
  hist3->SetFillStyle(0);
  hist6->SetFillStyle(0);
  hist9->SetFillStyle(0);
  hist3->Draw();
  //hist6->Draw("SAME");
  //hist9->Draw("SAME");
  legend = new TLegend(0.15,0.75,0.4,0.85);
  legend->AddEntry(hist3,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist6,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist9,"Tau -> 3 Mu MC","f");
  //legend->Draw();
  c2->Print("eta.png");

  TCanvas *c3 = new TCanvas("3","muon phi",1200,900);
  xaxis = hist4->GetXaxis();
  yaxis = hist4->GetYaxis();
  xaxis->SetTitle("#varphi");
  yaxis->SetTitle("a.u / 0.28 rad");
  c3->cd();
  hist4->SetLineColor(2);
  hist7->SetLineColor(6);
  hist10->SetLineColor(4);
  hist4->Draw();
  //hist7->Draw("SAME");
  //hist10->Draw("SAME");
  legend = new TLegend(0.15,0.75,0.4,0.85);
  legend->AddEntry(hist4,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist7,"Tau -> 3 Mu MC","f");
  //legend->AddEntry(hist10,"Tau -> 3 Mu MC","f");
  //legend->Draw();

  c3->Print("phi.png");

  gBenchmark->Show("select3Mu");

}
