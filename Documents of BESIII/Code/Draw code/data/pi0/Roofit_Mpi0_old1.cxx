{
    using namespace RooFit;

    TFile *f1 = new TFile("/mnt/e/Work/IHEPBOX/root/0912MCmix/pi0MC5_NoGam.root");
    TTree *t1 = (TTree *)f1->Get("trkRes");

    TFile *f2 = new TFile("/mnt/e/Work/IHEPBOX/root/0912data/data5_pi0_NoGam.root");
    TTree *t2 = (TTree *)f2->Get("trkRes");

    TFile *f4 = new TFile("/mnt/e/Work/IHEPBOX/root/0912MCmix/old/2020.09.28/LambdaSigmapi0MC5_pi0_NoGam.root");
    TTree *t4 = (TTree *)f4->Get("trkRes");

    TFile *f5 = new TFile("/mnt/e/Work/IHEPBOX/root/0912MCmix/old/2020.09.28/LambdaSigmapimMC5_pi0_NoGam.root");
    TTree *t5 = (TTree *)f5->Get("trkRes");

    TFile *f6 = new TFile("/mnt/e/Work/IHEPBOX/root/0912MCmix/old/2020.09.28/jpsiMC5_pi0_NoGam.root");
    TTree *t6 = (TTree *)f5->Get("trkRes");

    Double_t xlow1 = 0.08;
    Double_t xup1 = 0.19;
    Int_t nbins1 = 44;

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

    RooDataSet lspmc("lspmc", "lspmc", t4, orig_m2gamma);
    RooKeysPdf lspmcpdf1("lspmcpdf1", "lspmcpdf1", orig_m2gamma, lspmc, RooKeysPdf::MirrorBoth, 1);
    //RooRealVar mean_lspmc2("mean_lspmc2", "mean_lspmc2", -100.0, 100.0);
    //RooRealVar width_lspmc2("width_lspmc2", "#width_lspmc2", 0.0, 10.0);
    //RooGaussModel lspmcpdf2("lspmcpdf2", "lspmcpdf2", orig_m2gamma, mean_lspmc2, width_lspmc2);
    RooFFTConvPdf lspmcpdf("lspmcpdf", "lspmcpdf", orig_m2gamma, lspmcpdf1, gausspdf);

    RooDataSet lsp1mc1("lsp1mc1", "lsp1mc1", t5, orig_m2gamma);
    RooKeysPdf lsp1mcpdf1("lsp1mcpdf1", "lsp1mcpdf1", orig_m2gamma, lsp1mc1, RooKeysPdf::MirrorBoth, 1);
    //RooRealVar mean_lsp1mc2("mean_lsp1mc2", "mean_lsp1mc2", -100.0, 100.0);
    //RooRealVar width_lsp1mc2("width_lsp1mc2", "#width_lsp1mc2", 0.0, 10.0);
    //RooGaussModel lsp1mcpdf2("lsp1mcpdf2", "lsp1mcpdf2", orig_m2gamma, mean_lsp1mc2, width_lsp1mc2);
    RooFFTConvPdf lsp1mcpdf("lsp1mcpdf", "lsp1mcpdf", orig_m2gamma, lsp1mcpdf1, gausspdf);

    RooDataSet jpsimc1("jpsimc1", "jpsimc1", t6, orig_m2gamma);
    RooKeysPdf jpsimcpdf1("jpsimcpdf1", "jpsimcpdf1", orig_m2gamma, jpsimc1, RooKeysPdf::MirrorBoth, 1);
    //RooRealVar mean_jpsimc2("mean_jpsimc2", "mean_jpsimc2", -100.0, 100.0);
    //RooRealVar width_jpsimc2("width_jpsimc2", "#width_jpsimc2", 0.0, 10.0);
    //RooGaussModel jpsimcpdf2("jpsimcpdf2", "jpsimcpdf2", orig_m2gamma, mean_jpsimc2, width_jpsimc2);
    RooFFTConvPdf jpsimcpdf("jpsimcpdf", "jpsimcpdf", orig_m2gamma, jpsimcpdf1, gausspdf);

    //-----Construct signal+background PDF------

    RooRealVar nsig("nsig", "nsig", 30, 0, 121);
    RooRealVar nbkg("nbkg", "nbkg", 30, 0, 121);
    RooRealVar nlsp("nlsp", "nlsp", 43, 43, 43);
    RooRealVar nlsp1("nlsp1", "nlsp1", 5, 5, 5);
    RooRealVar njpsi("njpsi", "njpsi", 13, 13, 13);

    RooArgList allpdf(sigpdf, bkgpdf, lspmcpdf1, lsp1mcpdf1, jpsimcpdf1);
    RooArgList alln(nsig, nbkg, nlsp, nlsp1, njpsi);

    RooAddPdf sum("sum", "sum", allpdf, alln);

    //------Fit------

    RooFitResult *result = sum.fitTo(data, Extended(kTRUE), Save());

    //------Frame------

    TCanvas *Canvas = new TCanvas("Canvas", "Fit_M_pi0.eps", 250, 50, 800, 600);

    RooPlot *mpi0frame = orig_m2gamma.frame();
    mpi0frame->GetXaxis()->SetTitle("M(#gamma#gamma) (GeV/c^{2})");
    mpi0frame->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
    mpi0frame->GetXaxis()->CenterTitle();
    mpi0frame->GetYaxis()->CenterTitle();

    data.plotOn(mpi0frame, MarkerStyle(8), MarkerColor(kBlack), MarkerSize(0.8), Binning(nbins1), LineColor(kBlack));
    sum.plotOn(mpi0frame, RooFit::LineColor(kBlack));
    sum.plotOn(mpi0frame, Components("sigpdf"), LineColor(kRed));        //, LineStyle(kDashed));
    sum.plotOn(mpi0frame, Components("bkgpdf"), LineColor(kGreen));      //, LineStyle(kDashed));
    sum.plotOn(mpi0frame, Components("lspmcpdf1"), LineColor(kBlue));    //, LineStyle(kDashed));
    sum.plotOn(mpi0frame, Components("lsp1mcpdf1"), LineColor(kCyan));   //, LineStyle(kDashed));
    sum.plotOn(mpi0frame, Components("jpsimcpdf1"), LineColor(kOrange)); //, LineStyle(kDashed));

    mpi0frame->Draw();

    Int_t nParsToFit = (result->floatParsFinal()).getSize();
    Int_t nBinX = mpi0frame->GetNbinsX();
    Int_t ndof = nBinX - nParsToFit;
    RooCurve *curve = (RooCurve *)mpi0frame->getObject(1);
    RooHist *histo = (RooHist *)mpi0frame->getObject(0);
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

    TString Par2V = Form("%5.1f", nbkg.getVal());
    TString Par2E = Form("%3.1f", nbkg.getError());
    TString Par2 = "N_{bkg} = " + Par2V + " #pm " + Par2E;

    TString Par4V = Form("%5.1f", nlsp.getVal());
    TString Par4E = Form("%3.1f", nlsp.getError());
    TString Par4 = "N_{#Lambda#Sigma#pi^{0}} = " + Par4V;

    TString Par5V = Form("%5.1f", nlsp1.getVal());
    TString Par5E = Form("%3.1f", nlsp1.getError());
    TString Par5 = "N_{#Lambda#Sigma#pi^{-}} = " + Par5V;

    TString Par6V = Form("%5.1f", njpsi.getVal());
    TString Par6E = Form("%3.1f", njpsi.getError());
    TString Par6 = "N_{#pi^{+}#pi^{-}J/#psi} = " + Par6V;

    TString Par10C = Form("%5.2f", curve->chiSquare(*histo, nParsToFit));
    TString Par10 = "#chi^{2}/ndf = " + Par10C;

    TText *text;
    text = pt->AddText(Par1);
    ((TText *)pt->GetListOfLines()->Last())->SetTextColor(kRed);

    text = pt->AddText(Par2);
    ((TText *)pt->GetListOfLines()->Last())->SetTextColor(kGreen);

    text = pt->AddText(Par4);
    ((TText *)pt->GetListOfLines()->Last())->SetTextColor(kBlue);

    text = pt->AddText(Par5);
    ((TText *)pt->GetListOfLines()->Last())->SetTextColor(kCyan);

    text = pt->AddText(Par6);
    ((TText *)pt->GetListOfLines()->Last())->SetTextColor(kOrange);

    text = pt->AddText(Par10);

    cout << "a = " << chisq_test << endl;
    cout << "b = " << chi2 << endl;
    pt->Draw();
    Canvas->Update();
}