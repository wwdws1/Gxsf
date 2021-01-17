#include <math.h>
void Egamma()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "E_gamma.eps", 250, 50, 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912MCmix/pi0MC6.root");

	TChain *ch2 = new TChain("trkRes");
	//ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/gammaMC6_pi0.root");
	ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/LambdaSigmapi0MC6_pi0.root");

	Double_t xlow1;
	Double_t xup1;
	Int_t nbins1;

	int n;
	n = 1;

	if (!n)
	{
		xlow1 = 0.0;
		xup1 = 1.2;
		nbins1 = 120;
	}
	else if (n)
	{
		xlow1 = 0.0;
		xup1 = 0.8;
		nbins1 = 80;
	}

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);

	if (!n)
	{
		ch1->Draw("emc_en[0]>>h1a", "2m1_4C_chisq < 0 && 2m3_4C_chisq < 0 && (orig_mgamma2lambda < 1.19 || orig_mgamma2lambda > 1.195) && orig_mgamma2lambdabar < 1.19 || orig_mgamma2lambdabar > 1.195");
		ch2->Draw("emc_en[0]>>h1b", "2m1_4C_chisq < 0 && 2m3_4C_chisq < 0 && (orig_mgamma2lambda < 1.19 || orig_mgamma2lambda > 1.195) && orig_mgamma2lambdabar < 1.19 || orig_mgamma2lambdabar > 1.195");
		h1a->GetYaxis()->SetRangeUser(0, 2000);
		h1a->GetXaxis()->SetTitle("E(#gamma_{1}) (GeV)");
	}
	else if (n)
	{
		ch1->Draw("emc_en[1]>>h1a", "2m1_4C_chisq < 0 && 2m3_4C_chisq < 0 && (orig_mgamma2lambda < 1.19 || orig_mgamma2lambda > 1.195) && orig_mgamma2lambdabar < 1.19 || orig_mgamma2lambdabar > 1.195");
		ch2->Draw("emc_en[1]>>h1b", "2m1_4C_chisq < 0 && 2m3_4C_chisq < 0 && (orig_mgamma2lambda < 1.19 || orig_mgamma2lambda > 1.195) && orig_mgamma2lambdabar < 1.19 || orig_mgamma2lambdabar > 1.195");
		h1a->GetYaxis()->SetRangeUser(0, 4000);
		h1a->GetXaxis()->SetTitle("E(#gamma_{2}) (GeV)");
	}
	//ch1->Draw("4C_chisq>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119");
	//ch2->Draw("4C_chisq>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119");
	//ch1->Draw("4C_chisq>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch2->Draw("4C_chisq>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch1->Draw("4C_chisq>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	//ch2->Draw("4C_chisq>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");

	h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	//h1b->Scale(h1a->GetMaximum() / h1b->GetMaximum(), "nosw2");
	//h1b->Sumw2();
	h1b->Scale(0.65 * h1a->Integral() / h1b->Integral(), "nosw2");

	//h1b->SetMarkerStyle(8);
	//h1b->SetMarkerColor(1);
	//h1b->SetMarkerSize(0.75);
	h1b->SetLineColor(1);
	h1b->SetLineWidth(1);
	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1a->Draw("");
	h1b->Draw("same");

	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h1a, "Signal", "l");
	lg1->AddEntry(h1b, "BKG(#Sigma^{0})", "l");
	lg1->SetFillColor(0);
	lg1->SetTextFont(42);
	lg1->SetTextSize(0.05);
	lg1->SetBorderSize(0);
	lg1->Draw();
	/*
	TArrow *x1 = new TArrow(37, 0, 37, 2000, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();
	*/
	c1->Update();

	//c1->SaveAs("./Chisq_4C_2.eps");
}