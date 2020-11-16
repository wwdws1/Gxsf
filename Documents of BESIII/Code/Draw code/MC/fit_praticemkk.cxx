{

	using namespace RooFit;
	//------Observvable--------
	RooRealVar mkk("mkk", "m_{kk} (Gev)", 1, 1.04);

	//------Build Gaussian signal PDF --------
	RooRealVar sigmean("sigmean", "B^{#pm} mass", 1.019, 1, 1.04);
	RooRealVar sigwidth("sigwidth", "B^{#pm} width", 0.0027, 0.00001, 1);
	RooGaussian gauss("gauss", "title name", mkk, sigmean, sigwidth);

	//-------Build Argus background PDF ---------
	//	RooRealVar argpar("argpar","argus shape parameter",-20.0,-100.,-1);
	//	RooArgusBG argus("argus","Argus PDF",mkk,RooConst(1.020),argpar);

	//--build background-Chebychev background PDF---
	RooRealVar c0("c0", "c0", 1, -10., 10);
	//RooRealVar c1("c1","c1", 1., -10.,10);
	RooChebychev bkg("bkg", "bkg", mkk, RooArgList(c0));

	//-------CONstruct signal+background PDF---------
	RooRealVar nsig("nsig", "#signal events", 10, 0., 10000);
	RooRealVar nbkg("nbkg", "#background events", 5, 0., 10000);
	//	RooAddPdf sum("sum","g+a",RooArgList(gauss,argus),RooArgList(nsig,nbkg));
	RooAddPdf sum("sum", "g+b", RooArgList(gauss, bkg), RooArgList(nsig, nbkg));
	//RooAddPdf sum("sum","g+b",RooArgList(bkg),RooArgList(nbkg));
	//	RooAddPdf sum("sum","g+b",RooArgList(gauss),RooArgList(nsig));

	//--------Genrate a toyMC sample from composite PDF -------
	Double_t lo(1), up(1.04);
	int Nbin(30);

	TH1F *mass = new TH1F("mass", "mass", Nbin, lo, up);
	//	TString fname = "mkk";

	TFile *file = new TFile("/besfs/users/fengjh/7.0.3/workarea/Analysis/AAnalysis/kkggAlg/kkggAlg-00-00-01/ana2/4230_EtaPhi_data_wuMkkCut_ana_10.root");
	TTree *r_data = (TTree *)file->Get("m3");
	TH1F *h1 = new TH1F("h1", "h1", Nbin, lo, up);
	r_data->Draw("mkk>>h1");
	RooDataHist *data = new RooDataHist("data", "data", mkk, h1);

	//-------Perform extended ML fit of composite PDF to toy data ------
	sum.fitTo(*data, Extended());

	//-------Plot toy data and composite PDF overlaid--------
	RooPlot *mkkframe = mkk.frame();
	data->plotOn(mkkframe);
	sum.plotOn(mkkframe);
	sum.plotOn(mkkframe, Components(bkg), LineStyle(kDashed));
	mkkframe->Draw();

	//      h1->Draw();
}
