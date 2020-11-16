{
	using namespace RooFit;
	//------Observvable--------
	RooRealVar meta("meta", "m(#eta) (Gev/c^{2})", 0.1, 0.17);

	//------Build Gaussian signal PDF --------
	RooRealVar sigmean("sigmean", "B^{#pm} mass", 0.135, 0.01, 0.3);
	RooRealVar sigwidth("sigwidth", "B^{#pm} width", 0.027, 0.00001, 1);
	RooGaussian gauss("gauss", "title name", meta, sigmean, sigwidth);

	//-------Build Argus background PDF ---------
	//	RooRealVar argpar("argpar","argus shape parameter",-20.0,-100.,-1);
	//	RooArgusBG argus("argus","Argus PDF",meta,RooConst(1.020),argpar);

	//--build background-Chebychev background PDF---
	RooRealVar c0("c0", "c0", 1, -10., 10);
	//RooRealVar c1("c1","c1", 1., -10.,10);
	RooChebychev bkg("bkg", "bkg", meta, RooArgList(c0));

	//-------CONstruct signal+background PDF---------
	RooRealVar nsig("nsig", "#signal events", 10, 0., 10000);
	RooRealVar nbkg("nbkg", "#background events", 1, 0., 10000);
	//	RooAddPdf sum("sum","g+a",RooArgList(gauss,argus),RooArgList(nsig,nbkg));
	RooAddPdf sum("sum", "g+b", RooArgList(gauss, bkg), RooArgList(nsig, nbkg));
	//RooAddPdf sum("sum","g+b",RooArgList(bkg),RooArgList(nbkg));
	//	RooAddPdf sum("sum","g+b",RooArgList(gauss),RooArgList(nsig));

	//--------Genrate a toyMC sample from composite PDF -------
	Double_t lo(0.1), up(0.17);
	int Nbin(70);

	TH1F *mass = new TH1F("mass", "mass", Nbin, lo, up);
	//	TString fname = "meta";

	TFile *file = new TFile("../../../0912MCmix/root/etaMC4_NoGam.root");
	TTree *r_data = (TTree *)file->Get("trkRes");
	TH1F *h1 = new TH1F("h1", "h1", Nbin, lo, up);
	r_data->Draw("orig_m2gamma>>h1");
	RooDataHist *data = new RooDataHist("data", "data", meta, h1);

	//-------Perform extended ML fit of composite PDF to toy data ------
	sum.fitTo(*data, Extended());

	//-------Plot toy data and composite PDF overlaid--------
	RooPlot *metaframe = meta.frame();
	data->plotOn(metaframe);
	sum.plotOn(metaframe);
	sum.plotOn(metaframe, Components(bkg), LineStyle(kDashed));
	sum.plotOn(metaframe, Components(gauss), LineStyle(kDashed));
	metaframe->Draw();

	//      h1->Draw();
}
