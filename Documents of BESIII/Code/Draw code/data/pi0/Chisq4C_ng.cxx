#include <math.h>
void Chisq4C_ng()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "Chisq_4C_nGamma.eps", 250, 50, 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912MCmix/pi0MC6.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/gammaMC6_pi0.root");

	TChain *ch3 = new TChain("trkRes");
	ch3->Add("E:/Work/IHEPBOX/root/0912data/data6_pi0.root");

	Double_t xlow1 = -300.0;
	Double_t xup1 = 300.0;
	Int_t nbins1 = 400;

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
	TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);

	ch1->Draw("2m3_4C_chisq>>h1a");
	ch2->Draw("2m3_4C_chisq>>h1b");
	ch3->Draw("2m3_4C_chisq>>h1c");

	h1a->GetXaxis()->SetTitle("#diff(chi_{4C}^{2})");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->SetLineColor(kBlue);
	h1a->SetLineWidth(1);

	h1a->Scale(0.6 * h1c->Integral() / h1a->Integral(), "nosw2");
	h1b->Scale(0.2 * h1c->Integral() / h1b->Integral(), "nosw2");
	h1c->Sumw2();
	//h1b->Scale(h1c->GetMaximum() / h1b->GetMaximum());
	//h1b->Scale(h1a->Integral()/h1b->Integral());

	h1c->SetMarkerStyle(8);
	h1c->SetMarkerColor(kBlack);
	h1c->SetMarkerSize(0.75);
	h1c->SetLineColor(kBlack);
	h1c->SetLineWidth(1);

	//h1b->SetMarkerStyle(22);
	//h1b->SetMarkerColor(2); // 2
	//h1b->SetMarkerSize(0.75);
	h1b->SetLineColor(kGreen); // 2
	h1b->SetLineWidth(1);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1c->Draw("");
	h1b->Draw("SAME");
	h1a->Draw("SAME");

	TLegend *lg1 = new TLegend(0.68, 0.68, 0.86, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "Signal", "l");
	lg1->AddEntry(h1b, "BKG(#gamma)", "l");
	lg1->AddEntry(h1c, "DATA", "pl");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();

	//TArrow *x1 = new TArrow(15, 0, 15, 120, 0.03, "<");
	//x1->SetLineColor(2);
	//x1->SetLineWidth(1);
	//x1->Draw();

	c1->Update();

	//c1->SaveAs("./Chisq_4C_2.eps");
}