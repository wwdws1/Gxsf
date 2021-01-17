#include <math.h>
void Angel_fdc()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_Lambdaeta_fdc.eps", 250, 50, 800, 600);

	TFile *f1 = new TFile("E:/Work/IHEPBOX/root/0912MCmix/dplot_eta.root");

	TFile *f2 = new TFile("E:/Work/IHEPBOX/root/0912MCmix/dplot_eta_bkg.root");

	TFile *f3 = new TFile("E:/Work/IHEPBOX/root/0912MCmix/mplot_eta_1670.root");

	double xlow1 = -1.0;
	double xup1 = 1.0;
	int nbins1 = 50;

	TH1F *h1a = (TH1F *)f1->Get("h7");
	TH1F *h1b = (TH1F *)f2->Get("h7");
	TH1F *h1c = (TH1F *)f3->Get("h7");
	TH1F *h1d = (TH1F *)f4->Get("h7");

	TH1F *h1e = (TH1F *)f3->Get("h31");
	TH1F *h1f = (TH1F *)f4->Get("h43");

	h1a->GetXaxis()->SetTitle("cos(#theta)");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->GetYaxis()->SetRangeUser(0, 30);

	h1a->Sumw2();
	h1b->Scale(0.1 * h1a->Integral() / h1b->Integral(), "nosw2");
	h1c->Scale(0.8 * h1a->Integral() / h1c->Integral(), "nosw2");
	h1d->Scale(0.8 * h1a->Integral() / h1d->Integral()); //, "nosw2");
	h1e->Scale(0.8 * 0.31 * h1a->Integral() / h1e->Integral(), "nosw2");
	h1f->Scale(0.8 * 0.23467186 * h1a->Integral() / h1f->Integral(), "nosw2");

	//h1d->Scale(h1c->GetMaximum() / h1d->GetMaximum(), "nosw2");
	//h1a->Scale(h1b->Integral()/h1a->Integral());

	h1c->Add(h1b, 1);
	h1d->Add(h1b, 1);

	h1a->SetMarkerStyle(8);
	h1a->SetMarkerColor(kBlack);
	h1a->SetMarkerSize(1);
	h1a->SetLineColor(1);
	h1a->SetLineWidth(1);

	h1b->SetLineColor(kGreen);

	//h1c->SetMarkerStyle(8);
	//h1c->SetMarkerColor(kBlue);
	h1c->SetLineColor(kRed);

	h1d->SetMarkerStyle(8);
	h1d->SetMarkerColor(kRed);
	h1d->SetLineColor(kRed);

	h1e->SetLineColor(kBlue);

	h1f->SetLineColor(kRed);

	h1a->Draw("");
	h1b->Draw("SAME");
	h1c->Draw("SAME");
	//h1d->Draw("SAME");
	h1e->Draw("SAME");
	//h1f->Draw("SAME");

	TLegend *lg1 = new TLegend(0.68, 0.68, 0.86, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "DATA", "pl");
	lg1->AddEntry(h1b, "BKG(Sideband)", "l");
	lg1->AddEntry(h1c, "Model(1670)+BKG", "pl");
	//lg1->AddEntry(h1d, "Model(1690)+BKG", "pl");
	lg1->AddEntry(h1e, "#Lambda(1670)Only", "l");
	//lg1->AddEntry(h1f, "#Lambda(1690)Only", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(1.7, 0, 1.7, 100, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(1.121, 0, 1.121, 200, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
	c1->Update();

	//c1->SaveAs("./M_Lambda_2.eps");
}