{
    using namespace RooFit;

    TFile *f1 = new TFile("../../../0912MCmix/root/etaMC4_NoGam.root");
    TTree *t1 = (TTree *)f1->Get("trkRes");

    TFile *f2 = new TFile("../../root/data4_eta_NoGam.root");
    TTree *t2 = (TTree *)f2->Get("trkRes");

    Double_t xlow1 = 0.46;
    Double_t xup1 = 0.63;
    Int_t nbins1 = 68;

    //------Data Histogram--------

    RooRealVar orig_m2gamma("orig_m2gamma", "orig_m2gamma", xlow1, xup1);
    RooDataSet data("data", "data", t2, orig_m2gamma);

    //------Build signal PDF--------

    RooDataSet exmc("exmc", "exmc", t1, orig_m2gamma);
    RooKeysPdf exmcpdf("exmcpdf", "exmcpdf", orig_m2gamma, exmc, RooKeysPdf::MirrorBoth, 1);

    RooRealVar mean_gauss("mean_gauss", "mean_gauss", -100.0, 100.0);
    RooRealVar width_gauss("width_gauss", "#width_gauss", 0.0, 10.0);
    RooGaussModel gausspdf("gausspdf", "gausspdf", orig_m2gamma, mean_gauss, width_gauss);

    RooFFTConvPdf sigpdf("sigpdf", "sigpdf", orig_m2gamma, exmcpdf, gausspdf);

    //------Build background pdf------

    RooRealVar c0("c0", "c0", -10000.0, 10000.0);
    RooRealVar c1("c1", "c1", -10000.0, 10000.0);
    RooChebychev bkgpdf("bkgpdf", "bkgpdf", orig_m2gamma, RooArgList(c0, c1));

    //-----Construct signal+background PDF------

    RooRealVar nsig("nsig", "nsig", 200, 0, 500);
    RooRealVar nbkg("nbkg", "nbkg", 200, 0, 500);
    RooAddPdf sum("sum", "sum", RooArgList(sigpdf, bkgpdf), RooArgList(nsig, nbkg));

    //------Fit------

    RooFitResult *result = sum.fitTo(data, Extended(kTRUE), Save());

    //------Frame------

    TCanvas *Canvas = new TCanvas("Canvas", "Fit_M_eta.eps", 250, 50, 800, 600);

    RooPlot *metaframe = orig_m2gamma.frame();
    metaframe->GetXaxis()->SetTitle("M(#gamma#gamma) (GeV/c^{2})");
    metaframe->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
    metaframe->GetXaxis()->CenterTitle();
    metaframe->GetYaxis()->CenterTitle();

    data.plotOn(metaframe, MarkerStyle(8), MarkerSize(0.8), Binning(nbins1));
    sum.plotOn(metaframe, RooFit::LineColor(kRed));
    sum.plotOn(metaframe, Components("sigpdf"), LineColor(kViolet), LineStyle(kDashed));
    sum.plotOn(metaframe, Components("bkgpdf"), LineColor(kBlue), LineStyle(kDashed));

    metaframe->Draw();

    Int_t nParsToFit = (result->floatParsFinal()).getSize();
    Int_t nBinX = metaframe->GetNbinsX();
    Int_t ndof = nBinX - nParsToFit;
    RooCurve *curve = (RooCurve *)metaframe->getObject(1);
    RooHist *histo = (RooHist *)metaframe->getObject(0);
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
}