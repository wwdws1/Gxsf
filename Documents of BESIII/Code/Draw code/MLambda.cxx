#include <math.h>
void MLambda()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_Lambda.eps", 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("../root/pi0MC.root");
	//ch1->Add("../root/etaMC.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("../root/0912INMC1_pi0.root");

	Double_t xlow1 = 1.1;
	Double_t xup1 = 1.13;
	Int_t nbins1 = 100;

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);

	ch1->Draw("orig_mlambda>>h1a");
	ch2->Draw("orig_mlambda>>h1b");
	//ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	//ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	//ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569 && 4C_chisq < 37");
	//ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569 && 4C_chisq < 37");
	h1a->GetXaxis()->SetTitle("M(#Lambda) (GeV/c^{2})");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	//h1a->Scale(h1b->GetMaximum() / h1a->GetMaximum(), "nosw2");
	h1b->Sumw2();
	//h1a->Scale(h1b->Integral()/h1a->Integral());
	h1b->SetMarkerStyle(8);
	h1b->SetMarkerColor(1);
	h1b->SetMarkerSize(1);
	h1b->SetLineColor(1);
	h1b->SetLineWidth(1);
	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1a->Draw("");
	//h1b->Draw("SAME");

	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "Exclusive MC", "l");
	lg1->AddEntry(h1b, "Inclusive MC", "pl");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(1.1118, 0, 1.1118, 600, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(1.1194, 0, 1.1194, 600, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
	c1->Update();

	//c1->SaveAs("./M_Lambda_2.eps");
}