#include <math.h>
void Meta()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "M_eta.eps", 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("../root/MC_eta.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("../root/INMC.root");

	Double_t xlow1 = 0.46;
	Double_t xup1 = 0.63;
	Int_t nbins1 = 68;

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);

	//ch1->Draw("orig_m2gamma>>h1a");
	//ch2->Draw("orig_m2gamma>>h1b");
	//ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119");
	//ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119");
	//ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 37");
	ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 37");
	h1a->GetXaxis()->SetTitle("M_{#eta} (GeV/c^{2})");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	h1a->Scale(h1b->GetMaximum() / h1a->GetMaximum(), "nosw2");
	h1b->Sumw2();
	//h1b->Scale(h1a->Integral()/h1b->Integral());
	h1b->SetMarkerStyle(8);
	h1b->SetMarkerColor(1);
	h1b->SetMarkerSize(1);
	h1b->SetLineColor(1);
	h1b->SetLineWidth(1);
	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1a->Draw("");
	h1b->Draw("same");

	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "MC", "pl");
	lg1->AddEntry(h1b, "INMC", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(0.527, 0, 0.527, 1000, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(0.569, 0, 0.569, 1000, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
	c1->Update();
	
	//c1->SaveAs("./M_pi0_1.eps");
}