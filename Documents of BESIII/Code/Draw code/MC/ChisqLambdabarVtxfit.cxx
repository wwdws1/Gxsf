#include <math.h>
void ChisqLambdabarVtxfit()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "Chisq_LambdabarVtxfit.eps", 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("./MC2.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("./INMC2.root");

	Double_t xlow1 = 0.0;
	Double_t xup1 = 200.0;
	Int_t nbins1 = 100;

	TH1F *h0a = new TH1F("h0a", "", nbins1, xlow1, xup1);
	TH1F *h0b = new TH1F("h0b", "", nbins1, xlow1, xup1);

	ch1->Draw("lambdab_chi>>h0a");
	ch2->Draw("lambdab_chi>>h0b");
	//ch1->Draw("lambdab_chi>>h0a", "orig_m2gamma>0.121&&orig_m2gamma<0.149");
	//ch2->Draw("lambdab_chi>>h0b", "orig_m2gamma>0.121&&orig_m2gamma<0.149");
	h0a->GetXaxis()->SetTitle("#chi_{#bar{#Lambda} Vtxfit}^{2}(GeV/c^{2})");
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
	h0b->SetMarkerSize(0.75);
	h0b->SetLineColor(1);
	h0b->SetLineWidth(1);
	//h0b->SetLineStyle(1);
	//h0a->SetMinimum(0);//Outer shaft border minimum
	//h0b->SetMinimum(0);//Outer shaft border minimum
	h0a->Draw("");
	h0b->Draw("same");

	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h0a, "MC", "pl");
	lg1->AddEntry(h0b, "INMC", "l");
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
	c1->SaveAs("./Chisq_LambdabarVtxfit_2.eps");
}