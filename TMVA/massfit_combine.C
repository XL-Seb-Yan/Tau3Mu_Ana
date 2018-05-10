#include <TMath.h>
#include <TLegend.h>
void massfit_combine(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37);

  // --- Create Variable --- 
  RooRealVar mass("mass","mass",1.65,1.89,""); //the name "BDT_mass" will be used by RooDataSet to import data
  mass.setBins(48);

  //--- signal PDF ---
  RooRealVar* mu1 = new RooRealVar("DG_mu1","#mu1",1.77682,1.7,1.8,"");
  RooRealVar* mu2 = new RooRealVar("DG_mu2","#mu2",1.77682,1.7,1.81,"");
  RooRealVar* sigma1 = new RooRealVar("DG_sigma1","#sigma_{1}",0.02,0.015,0.05,"");
  RooRealVar* sigma2 = new RooRealVar("DG_sigma2","#sigma_{2}",0.03,0.02,0.05,"");
  mu1->setConstant(kFALSE);
  mu2->setConstant(kFALSE);
  sigma1->setConstant(kFALSE);
  sigma2->setConstant(kFALSE);
  RooRealVar* frac = new RooRealVar("DG_frac","frac",0.5,0,1); 
  RooRealVar* Ns = new RooRealVar("DG_Ns","N_{s}",3000,0,10000,"events");
  Ns->setConstant(kFALSE);
  RooGaussian* gauss1 = new RooGaussian("G1","",mass,*mu1,*sigma1);
  RooGaussian* gauss2 = new RooGaussian("G2","",mass,*mu2,*sigma2);
  RooAddPdf* doublegauss = new RooAddPdf("SummedG1G2","",RooArgList(*gauss1,*gauss2),*frac);
  //RooAddPdf* ex_doublegauss = new RooAddPdf("DG_ex","extDgauss",RooArgList(*doublegauss),RooArgList(*Ns));

  //--- background PDF ---
  RooRealVar *pC = new RooRealVar("pol2_pC","C",5,0,10,"");
  pC->setConstant(kFALSE);
  RooRealVar *p0 = new RooRealVar("pol2_p0","p_0",5,0,10,"");
  p0->setConstant(kFALSE);
  RooRealVar *p1 = new RooRealVar("pol2_p1","p_1",5,0,10,"");
  p1->setConstant(kFALSE);
  RooRealVar *p2 = new RooRealVar("pol2_p2","p_2",0.5,0,2,"");
  p2->setConstant(kFALSE);
  RooRealVar *p3 = new RooRealVar("pol2_p3","p_3",0.5,0,1,"");
  p3->setConstant(kFALSE);
  RooRealVar *p4 = new RooRealVar("pol2_p4","p_4",0.5,0,1,"");
  p4->setConstant(kFALSE);
  RooRealVar *Nbkg   = new RooRealVar("pol2_Nbkg","N_{bkg}",2000,0,10000,"events");
  Nbkg->setConstant(kFALSE);
  RooBernstein* bern = new RooBernstein("Bern","",mass,RooArgList(*pC,*p0));
  //RooAddPdf* ex_bern = new RooAddPdf("Bern_ex","",RooArgList(*bern),RooArgList(*Nbkg));

  //--- combine PDF ---
  //RooAddPdf* ex_sum = new RooAddPdf("SUM_ex","",RooArgList(*bern,*doublegauss),RooArgList(*Nbkg,*Ns));
  
  // --- Import unbinned dataset ---
  TFile fileS("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/Tau3Mu_Ana/TMVA/InvMassS.root");
  TFile fileB("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/Tau3Mu_Ana/TMVA/InvMassB.root");
  TTree* treeS = (TTree*)fileS.Get("Events");
  TTree* treeB = (TTree*)fileB.Get("Events");
  RooArgList imarglist(mass);
  RooArgSet imargset(imarglist);
  RooDataSet dataS("MC","MC",imargset,Import(*treeS));//import branches with names match the "variable name" (not variable) listed in imargset
  RooDataSet dataB("DATA Sidebands","DATA Sidebands",imargset,Import(*treeB));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  //w.pdf("model_s")->fitTo(data);
  mass.setRange("SB1",1.65,1.72);
  mass.setRange("SB2",1.82,1.89);
  mass.setRange("Signal",1.72,1.82);
  RooFitResult *rS = doublegauss->fitTo(dataS,Range("Signal"),Save());
  RooFitResult *rB = bern->fitTo(dataB,Range("SB2,SB1"),Save());

  // --- Plot ---
  RooPlot *frame = mass.frame();
  bern->plotOn(frame,LineColor(kBlue),RooFit::Name("sidebandfun"));
  doublegauss->plotOn(frame,LineColor(kRed),RooFit::Name("signal"));

  //axis,log scale and range setting functions must be called after all plotOn functions being called
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->cd();
  TAxis* xaxis = frame->GetXaxis();
  TAxis* yaxis = frame->GetYaxis();
  xaxis->SetTitle("Mass (GeV)");
  xaxis->SetTitleOffset(1.2);
  yaxis->SetTitle("Events / 10 MeV");
  yaxis->SetTitleOffset(1.2);
  //frame->SetMaximum(600);
  frame->SetMinimum(0);
  //c1->SetLogy();
  frame->Draw();
  TLegend *l =  new TLegend(0.6,0.7,0.8,0.78);
  l->AddEntry(frame->findObject("sidebandfun"),"Signal Fit","l");
  //l->AddEntry(frame->findObject("bkgfun"),"Background Fit","l");
  l->AddEntry(frame->findObject("err"),"Fit Error 1 #sigma","f");
  //l->Draw("same");
  c1->Print("data.png");

  // --- Output root file ---
  RooWorkspace *wUP = new RooWorkspace("wUP","wUP");
  wUP->var("mass[1.65,1.89]");
  //wUP->import(data,Rename("data_obs"));
  wUP->import(*doublegauss);
  wUP->import(*bern);
  wUP->writeToFile("data-shapes-UnbinnedParam.root");
}
