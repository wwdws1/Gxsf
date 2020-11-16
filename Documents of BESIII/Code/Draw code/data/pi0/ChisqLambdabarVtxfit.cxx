#include <math.h>
void ChisqLambdabarVtxfit()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "Chisq_LambdabarVtxfit.eps", 100, 100, 800, 500);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("../../0912MCmix/root/pi0MC4.root");
	//ch1->Add("../root/MC_eta.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("../../0912MCmix/root/INMC4_pi0.root");

	TChain *ch3 = new TChain("trkRes");
	ch3->Add("../root/data4_pi0.root");

	Double_t xlow1 = 0.0;
	Double_t xup1 = 200.0;
	Int_t nbins1 = 100;

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
	TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);

	ch1->Draw("lambdabar_chisq>>h1a");
	ch2->Draw("lambdabar_chisq>>h1b");
	ch3->Draw("lambdabar_chisq>>h1c");
	//ch1->Draw("lambdabar_chi>>h1a", "orig_m2gamma>0.121&&orig_m2gamma<0.149");
	//ch2->Draw("lambdabar_chi>>h1b", "orig_m2gamma>0.121&&orig_m2gamma<0.149");
	h1a->GetXaxis()->SetTitle("#chi_{#bar{#Lambda} Vtxfit}^{2}");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	h1a->Scale(h1c->GetMaximum() / h1a->GetMaximum(), "nosw2");
	h1b->Scale(h1c->GetMaximum() / h1b->GetMaximum());
	h1c->Sumw2();
	//h1b->Scale(h1a->Integral()/h1b->Integral());

	h1b->SetMarkerStyle(8);
	h1b->SetMarkerColor(1);
	h1b->SetMarkerSize(0.75);
	h1b->SetLineColor(1);
	h1b->SetLineWidth(1);

	h1c->SetMarkerStyle(22); // 22
	h1c->SetMarkerColor(2);	 // 2
	h1c->SetMarkerSize(0.75);
	h1c->SetLineColor(2); // 2
	h1c->SetLineWidth(1);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1a->Draw("");
	h1b->Draw("SAME");
	h1c->Draw("SAME");

	TLegend *lg1 = new TLegend(0.68, 0.68, 0.86, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "ExMC", "l");
	lg1->AddEntry(h1b, "InMC", "pl");
	lg1->AddEntry(h1c, "DATA", "pl");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();

	/*

	TArrow *x1 = new TArrow(1.00918, 0, 1.00918, 4.45, 0.03, "<");
	x1->SetLineColor(4);
	x1->SetLineWidth(3);
	x1->Draw();
	//  TArrow *x2= new TArrow(1.0265,2,1.0265,20,0.03,"<");
	TArrow *x2 = new TArrow(1.02882, 0, 1.02882, 4.45, 0.03, "<");
	x2->SetLineColor(4);
	x2->SetLineWidth(3);
	x2->Draw();

	*/
	c1->Update();

	//c1->SaveAs("./Chisq_4C_2.eps");
}