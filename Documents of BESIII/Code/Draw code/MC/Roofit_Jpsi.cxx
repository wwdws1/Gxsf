{
    using namespace RooFit;

    TFile *f1 = new TFile("/mnt/e/Work/IHEPBOX/root/psip/mc/jpsiMC9_pi0.root");
    TTree *t1 = (TTree *)f1->Get("trkRes");

    Double_t xlow1 = 3.08;
    Double_t xup1 = 3.12;
    Int_t nbins1 = 50;

    //------Data Histogram--------

    RooRealVar orig_recoilmpippim("orig_recoilmpippim", "orig_recoilmpippim", xlow1, xup1);
    RooDataSet jpsi("jpsi", "jpsi", t1, orig_recoilmpippim);

    //------Build signal PDF--------

    RooRealVar mean_gauss("mean_gauss", "mean_gauss", 3.097, 3.08, 3.12);
    RooRealVar width_gauss("width_gauss", "#width_gauss", 0.0, 0.1);
    RooGaussian gausspdf("gausspdf", "gausspdf", orig_recoilmpippim, mean_gauss, width_gauss);

    //------Build background pdf------

    RooRealVar c0("c0", "c0", -10000.0, 10000.0);
    RooRealVar c1("c1", "c1", -10000.0, 10000.0);
    RooChebychev bkgpdf("bkgpdf", "bkgpdf", orig_recoilmpippim, RooArgList(c0, c1));

    //-----Construct signal+background PDF------

    RooRealVar nsig("nsig", "nsig", 0, 500);
    RooRealVar nbkg("nbkg", "nbkg", 200, 0, 500);
    RooAddPdf sum("sum", "sum", RooArgList(gausspdf), RooArgList(nsig));

    //------Fit------

    RooFitResult *result = sum.fitTo(jpsi, Extended(kTRUE), Save());

    //------Frame------

    TCanvas *Canvas = new TCanvas("Canvas", "Fit_M_jpsi.eps", 250, 50, 800, 600);

    RooPlot *metaframe = orig_recoilmpippim.frame();
    metaframe->GetXaxis()->SetTitle("M(J/#psi) (GeV/c^{2})");
    metaframe->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
    metaframe->GetXaxis()->CenterTitle();
    metaframe->GetYaxis()->CenterTitle();

    jpsi.plotOn(metaframe, MarkerStyle(8), MarkerColor(kBlack), MarkerSize(0.8), Binning(nbins1), LineColor(kBlack));
    sum.plotOn(metaframe, RooFit::LineColor(kBlack));
    sum.plotOn(metaframe, Components("sigpdf"), LineColor(kViolet)); //, LineStyle(kDashed));
    //sum.plotOn(metaframe, Components("bkgpdf"), LineColor(kBlue));   //, LineStyle(kDashed));

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

    TString Par1V = Form("%5.1f", nsig.getVal());
    TString Par1E = Form("%3.1f", nsig.getError());
    TString Par1 = "N_{sig} = " + Par1V + " #pm " + Par1E;

    //TString Par2V = Form("%5.1f", nbkg.getVal());
    //TString Par2E = Form("%3.1f", nbkg.getError());
    //TString Par2 = "N_{bkg} = " + Par2V + " #pm " + Par2E;

    TString Par10C = Form("%5.2f", curve->chiSquare(*histo, nParsToFit));
    TString Par10 = "#chi^{2}/ndf = " + Par10C;

    TText *text;
    text = pt->AddText(Par1);
    ((TText *)pt->GetListOfLines()->Last())->SetTextColor(kViolet);

    //text = pt->AddText(Par2);
    //((TText *)pt->GetListOfLines()->Last())->SetTextColor(kBlue);

    text = pt->AddText(Par10);

    cout << "a = " << chisq_test << endl;
    cout << "b = " << chi2 << endl;
    pt->Draw();
    Canvas->Update();
}