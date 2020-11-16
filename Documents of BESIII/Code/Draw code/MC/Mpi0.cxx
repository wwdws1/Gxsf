#include <math.h>
void Mpi0()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_pi0.eps", 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("../root/pi0MC.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("../root/0912INMC.root");

	Double_t xlow1 = 0.08;
	Double_t xup1 = 0.190;
	Int_t nbins1 = 44;

	TH1F *h0a = new TH1F("h0a", "", nbins1, xlow1, xup1);
	TH1F *h0b = new TH1F("h0b", "", nbins1, xlow1, xup1);

	ch1->Draw("orig_m2gamma>>h0a");
	ch2->Draw("orig_m2gamma>>h0b");
	//ch1->Draw("orig_m2gamma>>h0a", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
	//ch2->Draw("orig_m2gamma>>h0b", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
	//ch1->Draw("orig_m2gamma>>h0a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch2->Draw("orig_m2gamma>>h0b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch1->Draw("orig_m2gamma>>h0a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 41");
	//ch2->Draw("orig_m2gamma>>h0b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 41");
	h0a->GetXaxis()->SetTitle("M(#pi^{0}) (GeV/c^{2})");
	h0a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h0a->GetXaxis()->CenterTitle();
	h0a->GetYaxis()->CenterTitle();
	h0a->SetLineColor(4);
	h0a->SetLineWidth(1);

	h0a->Scale(h0b->GetMaximum() / h0a->GetMaximum(), "nosw2");
	h0b->Sumw2();
	//h0b->Scale(h0a->Integral()/h0b->Integral());
	h0b->SetMarkerStyle(8);
	h0b->SetMarkerColor(1);
	h0b->SetMarkerSize(1);
	h0b->SetLineColor(1);
	h0b->SetLineWidth(1);
	//h0b->SetLineStyle(1);
	//h0a->SetMinimum(0);//Outer shaft border minimum
	//h0b->SetMinimum(0);//Outer shaft border minimum
	h0a->Draw("");
	h0b->Draw("same");

	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h0a, "Exclusive MC", "l");
	lg1->AddEntry(h0b, "Inclusive MC", "pl");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(0.120, 0, 0.120, 9000, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(0.150, 0, 0.150, 9000, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
	c1->Update();
	
	//c1->SaveAs("./M_pi0_1.eps");
}