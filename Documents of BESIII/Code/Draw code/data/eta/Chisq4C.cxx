#include <math.h>
void Chisq4C()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "Chisq_4C.eps", 250, 50, 800, 600);

	TChain *ch1 = new TChain("trkRes");
	ch1->Add("E:/Work/IHEPBOX/root/0912MCmix/etaMC4_NoChisq.root");
	//ch1->Add("../root/MC_eta.root");

	TChain *ch2 = new TChain("trkRes");
	ch2->Add("E:/Work/IHEPBOX/root/0912MCmix/INMC5_eta.root");

	TChain *ch3 = new TChain("trkRes");
	ch3->Add("E:/Work/IHEPBOX/root/0912data/data4_eta_NoChisq.root");

	Double_t xlow1 = 0.0;
	Double_t xup1 = 200.0;
	Int_t nbins1 = 100;

	TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
	TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
	TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);

	ch1->Draw("4C_chisq>>h1a");
	ch2->Draw("4C_chisq>>h1b");
	ch3->Draw("4C_chisq>>h1c");
	//ch1->Draw("4C_chisq>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119");
	//ch2->Draw("4C_chisq>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119");
	//ch1->Draw("4C_chisq>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch2->Draw("4C_chisq>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
	//ch1->Draw("4C_chisq>>h1a", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	//ch2->Draw("4C_chisq>>h1b", "orig_mlambda > 1.113 && orig_mlambda < 1.119 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
	h1a->GetXaxis()->SetTitle("#chi_{4C}^{2}");
	h1a->GetYaxis()->SetTitle(Form("Events/%.5f", (xup1 - xlow1) / nbins1));
	h1a->GetXaxis()->CenterTitle();
	h1a->GetYaxis()->CenterTitle();
	h1a->SetLineColor(4);
	h1a->SetLineWidth(1);

	h1a->Scale(0.5 * h1c->Integral() / h1a->Integral(), "nosw2");
	h1b->Scale(h1c->Integral() / h1b->Integral());
	h1c->Sumw2();
	//h1b->Scale(h1c->GetMaximum() / h1b->GetMaximum());
	//h1b->Scale(h1a->Integral()/h1b->Integral());

	h1b->SetMarkerStyle(22);
	h1b->SetMarkerColor(2); // 2
	h1b->SetMarkerSize(0.75);
	h1b->SetLineColor(2); // 2
	h1b->SetLineWidth(1);

	h1c->SetMarkerStyle(8);
	h1c->SetMarkerColor(1);
	h1c->SetMarkerSize(0.75);
	h1c->SetLineColor(1);
	h1c->SetLineWidth(1);

	//h1b->SetLineStyle(1);
	//h1a->SetMinimum(0);//Outer shaft border minimum
	//h1b->SetMinimum(0);//Outer shaft border minimum
	h1a->Draw("");
	//h1b->Draw("SAME");
	h1c->Draw("SAME");

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

	TArrow *x1 = new TArrow(40, 0, 40, 40, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	c1->Update();

	//c1->SaveAs("./Chisq_4C_2.eps");
}