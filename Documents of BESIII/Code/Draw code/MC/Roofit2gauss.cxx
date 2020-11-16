#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooAddModel.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

void fit_doubleguass_lambda()
{

	using namespace RooFit;
	//------Observvable--------
	RooRealVar mpi0("mpi0", "m() (Gev)", 1.111, 1.12);

	//------Build Gaussian signal PDF --------
	RooRealVar sigmean("sigmean", "B^{#pm} mass", 1.1157, 1.111, 1.12);
	//RooRealVar sigwidth("sigwidth","B^{#pm} width",0.00096,0.00001,0.01);
	//RooRealVar sigwidth("sigwidth","B^{#pm} width",0.0066,0.,0.01);
	RooRealVar sigwidth("sigwidth", "B^{#pm} width", 0.063, 0., 0.1);
	RooGaussian gauss("gauss", "title name", mpi0, sigmean, sigwidth);

	RooRealVar sigmean1("sigmean1", "B^{#pm} mass", 1.1157, 1.111, 1.12);
	RooRealVar sigwidth1("sigwidth1", "B^{#pm} width", 0.006, 0., 0.01);
	RooGaussModel gaussm1("gaussm1", "title name", mpi0, sigmean1, sigwidth1);

	// Build a composite resolution model f*gaussm+(1-f)*gaussm1
	//	RooRealVar gaussmfrac("gaussmfrac","fraction of gaussm",0.5) ;
	RooRealVar gaussmfrac("gaussmfrac", "fraction of gaussm", 0.3, 0, 1);
	RooAddModel gmsum("gmsum", "sum of gaussm and gaussm1", RooArgList(gauss, gaussm1), gaussmfrac);

	RooRealVar c0("c0", "c0", -0., -0.5, 0.5);
	RooRealVar c1("c1", "c1", -0., -0.5, 0.5.);
	RooChebychev bkg("bkg", "bkg", mpi0, RooArgList(c0, c1));

	//-------CONstruct signal+background PDF---------
	RooRealVar nsig("nsig", "#signal events", 1, 0., 1000000);
	RooRealVar nbkg("nbkg", "#background events", 1, 0., 1000000);
	RooAddPdf sum("sum", "g+b", RooArgList(gmsum, bkg), RooArgList(nsig, nbkg));
	//	RooAddPdf sum("sum","g+b",RooArgList(gauss),RooArgList(nsig));
	//--------Genrate a toyMC sample from composite PDF -------
	Double_t lo(1.111), up(1.12);
	int Nbin(80);

	//	TH1F *mass = new TH1F("mass","mass",Nbin, lo, up);

	//TFile *file = new TFile("/besfs/users/fengjh/7.0.3/workarea/Analysis/AAnalysis/XiXiAlg/XiXiAlg-00-00-01/ana_mgg/MXiXi_09_mc.root");  //第一次
	TFile *file = new TFile("/besfs/users/fengjh/7.0.3/workarea/Analysis/AAnalysis/XiXiAlg/XiXiAlg-00-00-01/ana_mgg/MMXiXi_09_mc.root"); //第二次
	TTree *r_data = (TTree *)file->Get("m3");
	TH1F *h1 = new TH1F("h1", "h1", Nbin, lo, up);
	r_data->Draw("mpi0>>h1");
	RooDataHist *data = new RooDataHist("data", "data", mpi0, h1);

	//-------Perform extended ML fit of composite PDF to toy data ------
	sum.fitTo(*data, Extended());

	//-------Plot toy data and composite PDF overlaid--------
	RooPlot *mpi0frame = mpi0.frame();
	//  mpi0frame->GetYaxis()->SetTitle("Events GeV/c^{2}");
	mpi0frame->GetYaxis()->SetTitle(Form("EVENTS/%.3f (GeV/c^{2})", (up - lo) / Nbin));
	mpi0frame->SetXTitle("M_{#Lambda} (GeV/c^{2})");
	mpi0frame->SetTitle("");
	mpi0frame->GetYaxis()->CenterTitle();
	mpi0frame->GetXaxis()->CenterTitle();
	data->plotOn(mpi0frame);
	sum.paramOn(mpi0frame);
	sum.plotOn(mpi0frame);
	sum.plotOn(mpi0frame, Components(bkg), LineStyle(kDashed));
	sum.plotOn(mpi0frame, Components(gmsum), LineStyle(kDashed));
	mpi0frame->Draw();

	double a = 0, b = 0, c = 0, frac_sig = 0, frac_bkg = 0, frac_sigmean = 0, frac_sigwidth = 0, frac_sigmean1 = 0, frac_sigwidth1 = 0, sig_err = 0, bkg_err = 0, sigmean_err = 0, sigwidth_err = 0, sigmean1_err = 0, sigwidth_err = 0;
	frac_sigmean = sigmean.getVal();
	sigmean_err = sigmean.getError();
	frac_sigwidth = sigwidth.getVal();
	sigwidth_err = sigwidth.getError();
	frac_sigmean1 = sigmean1.getVal();
	sigmean1_err = sigmean1.getError();
	frac_sigwidth1 = sigwidth1.getVal();
	sigwidth1_err = sigwidth1.getError();

	TLatex lt;
	lt.SetNDC();
	lt.SetTextAngle(0);
	lt.SetTextSize(0.04);
	lt.DrawLatex(0.612, 0.85, Form("N_{sig} = %1.0f #pm %1.0f", nsig.getVal(), nsig.getError()));
	lt.DrawLatex(0.612, 0.80, Form("N_{bkg} = %1.0f #pm %1.0f", nbkg.getVal(), nbkg.getError()));
	lt.DrawLatex(0.612, 0.75, Form("mean = %5.2f #pm %3.3f MeV", frac_sigmean * 1000, sigmean_err * 1000));
	lt.DrawLatex(0.612, 0.70, Form("sigma = %5.2f #pm %3.3f MeV", frac_sigwidth * 1000, sigwidth_err * 1000));
	lt.DrawLatex(0.612, 0.65, Form("mean_{1} = %5.2f #pm %3.3f MeV", frac_sigmean1 * 1000., sigmean1_err * 1000.));
	lt.DrawLatex(0.612, 0.60, Form("sigma_{1} = %5.2f #pm %3.3f MeV", frac_sigwidth1 * 1000., sigwidth1_err * 1000.));
	lt.DrawLatex(0.612, 0.55, Form("frac = %5.3f #pm %3.3f ", gaussmfrac.getVal(), gaussmfrac.getError()));
}
