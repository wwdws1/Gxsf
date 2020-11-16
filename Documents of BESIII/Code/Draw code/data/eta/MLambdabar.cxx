#include <math.h>
void MLambdabar()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_Lambdabar.eps", 250, 50, 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912MCmix/etaMC4_NoLamb.root");
	//ch1->Add("../root/MC_eta.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/INMC5_eta.root");

	//TChain *ch2 = new TChain("trkRes");
	//ch2->Add("../root/0912data.root");

	TChain *ch3 = new TChain("trkRes");
	ch3->Add("E:/Work/IHEPBOX/root/0912data/data4_eta_NoLamb.root");

	Double_t xlow1 = 1.1;
	Double_t xup1 = 1.13;
	Int_t nbins1 = 60;
	/*
	Double_t xlow1 = 1.08;
	Double_t xup1 = 1.15;
	Int_t nbins1 = 140;
	*/
	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
	TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);

	ch1->Draw("orig_mlambdabar>>h1a");
	ch2->Draw("orig_mlambdabar>>h1b");
	ch3->Draw("orig_mlambdabar>>h1c");
	//ch1->Draw("orig_mlambdabar>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
	//ch2->Draw("orig_mlambdabar>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
	//ch1->Draw("orig_mlambdabar>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_m2gamma > 0.121 && orig_m2gamma < 0.149");
	//ch2->Draw("orig_mlambdabar>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_m2gamma > 0.121 && orig_m2gamma < 0.149");
	//ch1->Draw("orig_mlambdabar>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_m2gamma > 0.121 && orig_m2gamma < 0.149 && 4C_chisq < 41");
	//ch2->Draw("orig_mlambdabar>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_m2gamma > 0.121 && orig_m2gamma < 0.149 && 4C_chisq < 41");
	h1c->GetXaxis()->SetTitle("M(#bar{p}#pi^{+}) (GeV/c^{2})");
	h1c->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1c->GetXaxis()->CenterTitle();
	h1c->GetYaxis()->CenterTitle();
	//h1c->GetYaxis()->SetRangeUser(0, 65);
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	h1a->Scale(0.8 * h1c->Integral() / h1a->Integral(), "nosw2");
	h1b->Scale(h1c->Integral() / h1b->Integral());
	h1c->Sumw2();
	//h1b->Scale(h1c->GetMaximum() / h1b->GetMaximum());
	//h1b->Scale(h1a->Integral()/h1b->Integral());

	h1b->SetMarkerStyle(22); // 22
	h1b->SetMarkerColor(2);	 // 2
	h1b->SetMarkerSize(1);
	h1b->SetLineColor(2); // 2
	h1b->SetLineWidth(1);

	h1c->SetMarkerStyle(8);
	h1c->SetMarkerColor(1);
	h1c->SetMarkerSize(1);
	h1c->SetLineColor(1);
	h1c->SetLineWidth(1);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1c->Draw("");
	//h1b->Draw("SAME");
	h1a->Draw("SAME");

	TLegend *lg1 = new TLegend(0.68, 0.68, 0.86, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "EX MC", "l");
	//lg1->AddEntry(h1b, "IN MC", "pl");
	lg1->AddEntry(h1c, "DATA", "pl");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();

	TArrow *x1 = new TArrow(1.111, 0, 1.111, 50, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(1.121, 0, 1.121, 50, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();

	c1->Update();

	//c1->SaveAs("./M_Lambdabar_2.eps");
}