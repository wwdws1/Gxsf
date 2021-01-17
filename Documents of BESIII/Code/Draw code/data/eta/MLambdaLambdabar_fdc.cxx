#include <math.h>
void MLambdaLambdabar_fdc()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_LambdaLambdabar.eps", 250, 50, 800, 600);

	TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/dplot_eta.root");

	TFile *f2 = new TFile("E:/Work/IHEPBOX/root/psip/mc/dplot_eta_bkg.root");

	TFile *f3 = new TFile("E:/Work/IHEPBOX/root/psip/mc/mplot_eta_1670.root");

	Double_t xlow1 = 2.0;
	Double_t xup1 = 3.8;
	Int_t nbins1 = 90;

	TH1F *h1a = (TH1F *)f1->Get("h3");
	TH1F *h1b = (TH1F *)f2->Get("h3");
	TH1F *h1c = (TH1F *)f3->Get("h3");
	TH1F *h1d = (TH1F *)f3->Get("h27");
	TH1F *h1e = (TH1F *)f3->Get("h39");

	h1a->GetXaxis()->SetTitle("M(p#pi^{-}#bar{p}#pi^{+}) (GeV/c^{2})");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	//h1a->GetYaxis()->SetRangeUser(0, 70);

	h1a->Sumw2();
	h1b->Scale(0.14 * h1a->Integral() / h1b->Integral(), "nosw2");
	h1c->Scale(0.86 * h1a->Integral() / h1c->Integral(), "nosw2");
	h1d->Scale(0.86 * 1.24 * h1a->Integral() / h1d->Integral(), "nosw2");
	h1e->Scale(0.86 * 0.56 * h1a->Integral() / h1e->Integral(), "nosw2");

	//h1d->Scale(h1c->GetMaximum() / h1d->GetMaximum(), "nosw2");
	//h1a->Scale(h1b->Integral()/h1a->Integral());

	h1c->Add(h1b, 1);

	h1a->SetMarkerStyle(8);
	h1a->SetMarkerColor(kBlack);
	h1a->SetMarkerSize(1);
	h1a->SetLineColor(kBlack);
	h1a->SetLineWidth(1);

	h1b->SetLineColor(kGreen);

	h1c->SetLineColor(kBlue);

	h1d->SetLineColor(kRed);

	h1e->SetLineColor(kOrange);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum

	h1a->Draw("");
	h1b->Draw("SAME");
	h1c->Draw("SAME");
	h1d->Draw("SAME");
	h1e->Draw("SAME");

	TLegend *lg1 = new TLegend(0.68, 0.58, 0.86, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "DATA", "pl");
	lg1->AddEntry(h1b, "BKG(Sideband)", "l");
	lg1->AddEntry(h1c, "Model(1670)+BKG", "l");
	lg1->AddEntry(h1d, "PHSP Only", "l");
	lg1->AddEntry(h1e, "#Lambda(1670) Only", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(3.4, 0, 3.4, 300, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x1 = new TArrow(3.077, 0, 3.077, 590, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(3.117, 0, 3.117, 590, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
	c1->Update();

	//c1->SaveAs("./M_Lambda_2.eps");
}