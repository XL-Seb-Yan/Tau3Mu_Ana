#include <TMath.h>
#include <TLegend.h>
void massfit_side(int seed=37)
{
  using namespace RooFit;
  RooRandom::randomGenerator()->SetSeed(37); 

  // --- Create Variable --- 
  RooRealVar mass("mass","mass",1.65,1.89,""); //the name "BDT_mass" will be used by RooDataSet to import data
  mass.setBins(48);

  //--- signal PDF ---
  /*
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
  RooAddPdf* ex_doublegauss = new RooAddPdf("DG_ex","extDgauss",RooArgList(*doublegauss),RooArgList(*Ns));
  */

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
  TFile file("/afs/cern.ch/work/x/xuyan/work4/CMSSW_8_0_27/src/Tau3Mu_Ana/TMVA/InvMassB3.root");
  TTree* tree = (TTree*)file.Get("Events");
  RooArgList imarglist(mass);
  RooArgSet imargset(imarglist);
  RooDataSet data("DATA Sidebands","DATA Sidebands",imargset,Import(*tree));//import branches with names match the "variable name" (not variable) listed in imargset

  // --- Perform extended ML fit of composite PDF to toy data ---
  //w.pdf("model_s")->fitTo(data);
  mass.setRange("SB1",1.65,1.72);
  mass.setRange("SB2",1.82,1.89);
  RooFitResult *r = bern->fitTo(data,Range("SB1,SB2"),Save());
  //RooFitResult *r = bern->fitTo(data,Save());
  //RooFitResult *r = ex_doublegauss->fitTo(data,Save());
  //RooFitResult *r = gauss1->fitTo(data,Save());
  //RooFitResult *r = doublegauss->fitTo(data,Save());

  // --- Plot ---
  RooPlot *frame = mass.frame();
  data.plotOn(frame,RooFit::Name("data"));
  frame->SetTitle("DATA Sidebands");
  //ex_sum->plotOn(frame,LineColor(kRed));
  bern->plotOn(frame,LineColor(kBlue),RooFit::Name("sidebandfun"));
  //gauss1->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //ex_sum->plotOn(frame,Components("Bern"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  bern->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  //gauss1->plotOn(frame,VisualizeError(*r,1),FillColor(kOrange),RooFit::Name("err"));
  bern->plotOn(frame,LineColor(kBlue),RooFit::Name("sidebandfun"));
  //ex_bern->plotOn(frame,LineColor(kRed),RooFit::Name("sigfun"));
  //doublegauss->plotOn(frame,Components("G1"),LineStyle(kDashed),RooFit::Name("bkgfun"));
  //doublegauss->plotOn(frame,Components("G2"),LineStyle(kDashed),LineColor(kGreen),RooFit::Name("bkgfun"));
  data.plotOn(frame);

  // residuals
  RooHist* hresid = frame->residHist("data","sidebandfun");
  RooPlot* frame2 = mass.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;
  //w.pdf("model_s")->plotOn(frame,Components("Bern"),LineStyle(kDashed));

  //axis,log scale and range setting functions must be called after all plotOn functions being called
  gStyle->SetOptStat(111111);
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
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
  c2->cd();
  TAxis* xaxis2 = frame2->GetXaxis();
  TAxis* yaxis2 = frame2->GetYaxis();
  xaxis2->SetTitle("Mass (GeV)");
  xaxis2->SetTitleOffset(1.2);
  yaxis2->SetTitle("Residuals / 10 MeV");
  yaxis2->SetTitleOffset(1.2);
  frame2->SetMaximum(50);
  frame2->SetMinimum(-50);
  frame2->Draw();
  c2->Print("residual.png");
}
