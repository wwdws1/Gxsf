#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TDirectory.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooKeysPdf.h"
#include "RooGaussian.h"

#include "TMath.h"
#include "RooMinuit.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
using namespace std;
gSystem->Load("libRooFit");
using namespace RooFit;
int data_mcpdf_xi()
{
	gStyle->SetLabelSize(0.04, "x");
	gStyle->SetLabelSize(0.04, "y");
	gStyle->SetTitleSize(0.08, "xyz");
	gStyle->SetTitleOffset(0.8, "xyz");
	gStyle->SetTitleOffset(1.0, "y");
	gStyle->SetLabelOffset(0.01, "xyz");
	gStyle->SetPadLeftMargin(0.18);
	gStyle->SetPadBottomMargin(0.2);
	gStyle->SetPadRightMargin(0.10);
	gStyle->SetPadTopMargin(0.10);

	Double_t nbin(140), xlow(1.29), xup(1.36);

	TString BranName("xi_M"), DaHist(">>DaHknew");
	TString DaName = BranName + DaHist;

	TFile *Daf, *Inf, *Daf;
	TH1F *DaHknew, *InHknew, *DaHknew, *FitHist;
	TTree *DaOutTree, *DaOutTree, *InOutTree, *FitTree;

	TCut cut_per = "lam_M<=1.120683&&lam_M>=1.110683&&xi_M<=1.4&&xi_M>=1.2";
	Daf = new TFile("/besfs/groups/jpsi/jpsigroup/user/chendy/akmplxi/root_kmplxi/onlycut_kmp_delerr_kvxy_dedx_data.root");
	DaHknew = new TH1F("DaHknew", "DaHknew", nbin, xlow, xup);
	DaOutTree = (TTree *)Daf->Get("main");
	DaOutTree->Draw(DaName, cut_per);
	cout << " Da gaussal " << DaHknew->GetEntries() << endl;

	//data histogram
	FitHist = DaHknew;
	RooRealVar xi_M(BranName, "xi_M", 1.321, xlow, xup);
	RooRealVar lam_M("lam_M", "lam_M", -100, 100);
	RooArgSet ntupleVarSetData(xi_M, lam_M);

	oridata = new RooDataSet("oridata", "KnewMass", DaOutTree, ntupleVarSetData);
	RooDataSet data = (RooDataSet)oridata->reduce(cut_per);

	//double gauss
	TFile *File = new TFile("/besfs/groups/jpsi/jpsigroup/user/chendy/akmplxi/root_kmplxi/xi_MCPDF.root");
	TTree *Tree = File->Get("tree");
	RooDataSet MCXi("MCXI", "MCXI", Tree, xi_M);
	RooKeysPdf pdf("pdf", "pdf", xi_M, MCXi, RooKeysPdf::MirrorBoth, 1);

	RooRealVar mean("mean", "mean", 0, -1, 1);
	RooRealVar sigma("sigma", "sigma", 0.0001, 0.00, 0.3);
	RooGaussModel gaussmcore("gaussmcore", "coregauss", xi_M, mean, sigma);
	RooFFTConvPdf sig("sig", "sig", xi_M, pdf, gaussmcore);

	//Chebychev
	RooRealVar c0("c0", "coefficient #0", -100000, 100000);
	RooRealVar c1("c1", "coefficient #1", -100000, 100000);
	RooRealVar c2("c2", "coefficient #2", -100000, 100000);
	RooRealVar c3("c3", "coefficient #3", -100000, 100000);
	RooChebychev bgfunc("bgfunc", "background pdf", xi_M, RooArgList(c0, c1));

	//4.Events
	RooRealVar nsig("nsig", "signal fraction", 5000, 4000, 8000);
	RooRealVar nbkg("nbkg", "bg fraction", 2000, 1500, 4000);

	//PDF
	RooAddPdf *sum = 0;
	sum = new RooAddPdf("sum", "sum", RooArgList(sig, bgfunc), RooArgList(nsig, nbkg));

	//Fit
	RooFitResult *result = sum->fitTo(data, Extended(kTRUE), Save());

	// 7. Frame
	RooPlot *frame = xi_M.frame();
	frame->SetXTitle("M-#Xi [GeV/c^{2}]");
	frame->GetXaxis()->CenterTitle(kTRUE);
	frame->GetXaxis()->SetLabelFont(42);
	frame->SetTitleFont(42, "x");
	frame->SetTitleFont(42, "Y");
	frame->SetTitleOffset(1.0, "x");
	frame->GetYaxis()->CenterTitle(kTRUE);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetTitleOffset(1.0);
	frame->GetYaxis()->SetLabelSize(0.05);
	frame->GetXaxis()->SetTitleOffset(0.74);
	frame->GetYaxis()->SetTitleOffset(0.68);
	frame->GetXaxis()->SetLabelSize(0.05);
	frame->GetXaxis()->SetTitleSize(0.09);
	frame->GetYaxis()->SetTitle(Form("Events/%.1f MeV/c^{2}", (xup - xlow) / nbin * 1000));

	data.plotOn(frame, MarkerStyle(8), MarkerSize(0.8), Binning(52));
	sum->plotOn(frame, RooFit::LineColor(kRed));
	sum->plotOn(frame, Components("sig"), LineColor(kViolet), LineStyle(kDashed));
	sum->plotOn(frame, Components("bgfunc"), LineColor(kBlue), LineStyle(kDashed));
	TCanvas *Canvas = new TCanvas("Canvas", "");
	//.plotOn( frame, LineColor(3));
	frame->Draw();

	Int_t nParsToFit = (result->floatParsFinal()).getSize();
	Int_t nBinX = frame->GetNbinsX();
	Int_t ndof = nBinX - nParsToFit;
	RooCurve *curve = (RooCurve *)frame->getObject(1);
	RooHist *histo = (RooHist *)frame->getObject(0);
	Double_t chisq_test = curve->chiSquare(*histo, nParsToFit); //get the value of (chi2/ndf)
	Double_t chi2 = chisq_test * ndof;
	cout << "#chi^{2}2/ndf: " << chisq_test << endl;
	cout << "chi=   " << chi2 << endl;

	TPaveText *pt = new TPaveText(0.65, 0.45, 0.79, 0.88, "BRNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(10);
	pt->SetTextAlign(12);
	pt->SetTextSize(0.05);
	pt->SetLineColor(0);
	TString Par1V = Form("%5.1f", nbkg.getVal());
	TString Par1E = Form("%3.1f", nbkg.getError());
	TString Par1 = "N_{bkg} = " + Par1V + " #pm " + Par1E;
	TString Par2V = Form("%5.1f", nsig.getVal());
	TString Par2E = Form("%3.1f", nsig.getError());
	TString Par2 = "N_{sig} = " + Par2V + " #pm " + Par2E;

	TString Par10C = Form("%5.2f", curve->chiSquare(*histo, nParsToFit));
	TString Par10 = "#chi^{2}/ndf = " + Par10C;

	TText *text;
	text = pt->AddText(Par1);
	text = pt->AddText(Par2);
	text = pt->AddText(Par10);

	cout << "a = " << chisq_test << endl;
	cout << "b = " << chi2 << endl;
	pt->Draw();
	Canvas->Update();
	Canvas->Print("Data_mcpdf_xi.eps");
}
