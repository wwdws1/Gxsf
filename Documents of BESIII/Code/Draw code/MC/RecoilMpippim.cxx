#include <math.h>
void RecoilMpippim()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "Recoil_M_pippim.eps", 250, 50, 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912MCmix/pi0MC6.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/jpsiMC6_pi0.root");

	Double_t xlow1 = 2.9;
	Double_t xup1 = 3.4;
	Int_t nbins1 = 100;

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);

	ch1->Draw("orig_recoilmpippim>>h1a", "2m1_4C_chisq < 0 && 2m3_4C_chisq < 0 && (orig_mgamma2lambda < 1.19 || orig_mgamma2lambda > 1.195) && orig_mgamma2lambdabar < 1.19 || orig_mgamma2lambdabar > 1.195");
	ch2->Draw("orig_recoilmpippim>>h1b", "2m1_4C_chisq < 0 && 2m3_4C_chisq < 0 && (orig_mgamma2lambda < 1.19 || orig_mgamma2lambda > 1.195) && orig_mgamma2lambdabar < 1.19 || orig_mgamma2lambdabar > 1.195");
	//ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	//ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	//ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569 && 4C_chisq < 37");
	//ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569 && 4C_chisq < 37");
	h1a->GetXaxis()->SetTitle("M_{Recoil}(#pi^{+}#pi^{-}) GeV/c^{2}");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->GetYaxis()->SetRangeUser(0, 16000);
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	//h1a->Scale(h1c->GetMaximum() / h1a->GetMaximum(), "nosw2");
	//h1b->Scale(h1a->GetMaximum() / h1b->GetMaximum(), "nosw2");
	h1b->Scale(0.325 * h1a->Integral() / h1b->Integral(), "nosw2");

	//h1b->SetMarkerStyle(8);
	//h1b->SetMarkerColor(1);
	//h1b->SetMarkerSize(1);
	h1b->SetLineColor(1);
	h1b->SetLineWidth(1);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1a->Draw("");
	h1b->Draw("SAME");

	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "Signal", "l");
	lg1->AddEntry(h1b, "BKG(J/#psi)", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(1.111, 0, 1.111, 200, 0.03, "<");
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