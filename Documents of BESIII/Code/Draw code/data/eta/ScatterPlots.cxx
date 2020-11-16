#include <math.h>
void ScatterPlots()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "Scatter Plots.eps", 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912data/data4_eta_FitGam.root");
	//ch1->Add("../root/MC_eta.root");
	//ch1->Add("../root/INMC.root");

	Double_t xlow1 = 0.46;
	Double_t xup1 = 0.63;
	Int_t nbinsx1 = 68;
	Double_t ylow1 = 1.1;
	Double_t yup1 = 2.0;
	Int_t nbinsy1 = 45;

	TH2F *h2a = new TH2F("h2a", "", nbinsx1, xlow1, xup1, nbinsy1, ylow1, yup1);

	ch1->Draw("orig_mgamma2lambda:orig_m2gamma>>h2a");
	//ch1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h2a", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
	//ch1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h2a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h2a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.121 && orig_m2gamma < 0.149");
	//ch1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h2a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.121 && orig_m2gamma < 0.149 && 4C_chisq < 41");
	h2a->GetXaxis()->SetTitle("M(#gamma#gamma) (GeV/c^{2})");
	h2a->GetYaxis()->SetTitle("M(#gamma_{2}p#pi^{-}) (GeV/c^{2})");
	//h2a->GetYaxis()->SetTitle("#chi_{4C}^{2}(GeV/c^{2})");
	h2a->GetXaxis()->CenterTitle();
	h2a->GetYaxis()->CenterTitle();
	//h2a->SetLineColor(4);
	//h2a->SetLineWidth(1);

	//h2a->Scale(h0b->GetMaximum() / h2a->GetMaximum(), "nosw2");
	//h0b->Sumw2();
	//h0b->Scale(h2a->Integral()/h0b->Integral());
	//h0b->SetMarkerStyle(8);
	//h0b->SetMarkerColor(1);
	//h0b->SetMarkerSize(1);
	//h0b->SetLineColor(1);
	//h0b->SetLineWidth(1);
	//h0b->SetLineStyle(1);
	//h2a->SetMinimum(0);//Outer shaft border minimum
	//h0b->SetMinimum(0);//Outer shaft border minimum
	h2a->Draw("SCAT");
	//h2a->Draw("COLZ");
	//h0b->Draw("same");

	/*
	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h2a, "MC", "pl");
	lg1->AddEntry(h0b, "INMC", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	*/

	/*
	TArrow *x1 = new TArrow(92, 0, 92, 2000, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();
	*/
	c1->Update();
	//c1->SaveAs("./Chisq_4C_2.eps");
}