#include <math.h>
void MgammaLambda()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_gammaLambda.eps", 250, 50, 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912MCmix/pi0MC6.root");
	//ch1->Add("../root/MC_eta.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/INMC6_pi0.root");

	TChain *ch3 = new TChain("trkRes");
	ch3->Add("E:/Work/IHEPBOX/root/0912data/data6_pi0.root");

	TChain *ch4 = new TChain("trkRes");
	ch4->Add("E:/Work/IHEPBOX/root/0912MCmix/LambdaSigmapimMC6_pi0.root");

	TChain *ch5 = new TChain("trkRes");
	ch5->Add("E:/Work/IHEPBOX/root/0912MCmix/LambdaSigmapi0MC6_pi0.root");

	TChain *ch6 = new TChain("trkRes");
	ch6->Add("E:/Work/IHEPBOX/root/0912MCmix/jpsiMC6_pi0.root");

	TChain *ch7 = new TChain("trkRes");
	ch7->Add("E:/Work/IHEPBOX/root/0912MCmix/gammaMC6_pi0.root");

	int clr;

	clr = 2;

	Double_t xlow1;
	Double_t xup1;
	Int_t nbins1;

	if (clr == 1)
	{
		xlow1 = 1.0;
		xup1 = 2.6;
		nbins1 = 80;
	}
	else if (clr == 2)
	{
		xlow1 = 1.0;
		xup1 = 2.0;
		nbins1 = 100;
	}

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
	TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);
	TH1F *h1d = new TH1F("h1d", "", nbins1, xlow1, xup1);
	TH1F *h1e = new TH1F("h1e", "", nbins1, xlow1, xup1);
	TH1F *h1f = new TH1F("h1f", "", nbins1, xlow1, xup1);
	TH1F *h1g = new TH1F("h1g", "", nbins1, xlow1, xup1);

	if (clr == 1)
	{
		ch1->Draw("orig_mgamma1lambda>>h1a");
		ch2->Draw("orig_mgamma1lambda>>h1b");
		ch3->Draw("orig_mgamma1lambda>>h1c");
		ch4->Draw("orig_mgamma1lambda>>h1d");
		ch5->Draw("orig_mgamma1lambda>>h1e");
		ch6->Draw("orig_mgamma1lambda>>h1f");
		ch7->Draw("orig_mgamma1lambda>>h1g");
	}
	else if (clr == 2)
	{
		ch1->Draw("orig_mgamma2lambda>>h1a");
		ch2->Draw("orig_mgamma2lambda>>h1b");
		ch3->Draw("orig_mgamma2lambda>>h1c");
		ch4->Draw("orig_mgamma2lambda>>h1d");
		ch5->Draw("orig_mgamma2lambda>>h1e");
		ch6->Draw("orig_mgamma2lambda>>h1f");
		ch7->Draw("orig_mgamma2lambda>>h1g");
	}

	h1c->GetXaxis()->SetTitle("M(#gammap#pi^{-}) (GeV/c^{2})");
	h1c->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1c->GetXaxis()->CenterTitle();
	h1c->GetYaxis()->CenterTitle();
	h1a->SetLineColor(kBlue);
	h1a->SetLineWidth(1);

	h1a->Scale(0.35 * h1c->Integral() / h1a->Integral(), "nosw2");
	h1b->Scale(h1c->Integral() / h1b->Integral());
	h1c->Sumw2();
	h1d->Scale(0.05 * h1c->Integral() / h1d->Integral(), "nosw2");
	h1e->Scale(0.35 * h1c->Integral() / h1e->Integral(), "nosw2");
	h1f->Scale(0.1 * h1c->Integral() / h1f->Integral(), "nosw2");
	h1g->Scale(0.15 * h1c->Integral() / h1g->Integral(), "nosw2");
	//h1d->Scale(h1c->GetMaximum() / h1d->GetMaximum(), "nosw2");
	//h1a->Scale(h1b->Integral()/h1a->Integral());

	h1a->Add(h1d, 1);
	h1a->Add(h1e, 1);
	h1a->Add(h1f, 1);
	h1a->Add(h1g, 1);

	h1b->SetMarkerStyle(22);
	h1b->SetMarkerColor(2);
	h1b->SetMarkerSize(1);
	h1b->SetLineColor(2);
	h1b->SetLineWidth(1);

	h1c->SetMarkerStyle(8);
	h1c->SetMarkerColor(1);
	h1c->SetMarkerSize(1);
	h1c->SetLineColor(1);
	h1c->SetLineWidth(1);

	h1d->SetLineColor(kCyan);
	h1e->SetLineColor(kPink);
	h1f->SetLineColor(kOrange);
	h1g->SetLineColor(kGreen);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1c->Draw("");
	h1b->Draw("SAME");
	h1a->Draw("SAME");
	h1d->Draw("SAME");
	h1e->Draw("SAME");
	h1f->Draw("SAME");
	h1g->Draw("SAME");

	TLegend *lg1 = new TLegend(0.68, 0.68, 0.86, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "EX MC(S+B)", "l");
	lg1->AddEntry(h1b, "IN MC", "pl");
	lg1->AddEntry(h1c, "DATA", "pl");
	lg1->AddEntry(h1d, "BKG(#Sigma^{+}#bar{#Lambda}#pi^{-})", "l");
	lg1->AddEntry(h1e, "BKG(#Sigma^{0}#bar{#Lambda}#pi^{0})", "l");
	lg1->AddEntry(h1f, "BKG(#pi^{+}#pi^{-}J/#psi)", "l");
	lg1->AddEntry(h1g, "BKG(#Lambda#bar{#Lambda}#gamma)", "l");
	//lg1->AddEntry(h1d, "BKG(#gamma#Lambda#bar{#Lambda})", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(1.111, 0, 1.111, 300, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(1.121, 0, 1.121, 300, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
	c1->Update();

	//c1->SaveAs("./M_Lambda_2.eps");
}