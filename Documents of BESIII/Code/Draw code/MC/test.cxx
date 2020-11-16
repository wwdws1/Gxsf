#include <math.h>
void test()
{
	gStyle->SetOptStat(kFALSE);
	TCanvas *c1 = new TCanvas("c1", "test.eps", 600, 600);

	/*
	TTree *tr1 = new TTree();
	TFile *fi1 = TFile::Open("../root/MC_pi0.root");
	fi1->GetObject("trkRes", tr1);
	*/
	TChain *ch1 = new TChain("trkRes");
	TChain *ch2 = new TChain("trkRes");
	ch1->Add("../root/MC_pi0.root");
	//ch1->Add("../root/MC_eta.root");
	//ch1->Add("../root/INMC.root");
	ch2->Add("../root/INMC.root");
	
	Double_t xlow1 = 1.00;
	Double_t xup1 = 2.60;
	Int_t nbinsx1 = 100;
	Double_t ylow1 = 1.00;
	Double_t yup1 = 2.60;
	Int_t nbinsy1 = 100;

	TH1F *h1a = new TH1F("h1a", "", )

	TH2F *h2a = new TH2F("h0a", "", nbinsx1, xlow1, xup1, nbinsy1, ylow1, yup1);

	//tr1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h0a");
	//tr1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h0a", "orig_m2gamma > 0.121 && orig_m2gamma < 0.149");
	ch1->Draw("orig_m2gammalambda:orig_m2gammalambdabar>>h0a", "orig_m2gamma > 0.121 && orig_m2gamma < 0.149");
	h0a->GetXaxis()->SetTitle("M_{#Lambda#pi^{0}}(GeV/c^{2})");
	h0a->GetYaxis()->SetTitle("M_{#bar{#Lambda}#pi^{0}}(GeV/c^{2})");
	h0a->GetXaxis()->CenterTitle();
	h0a->GetYaxis()->CenterTitle();
	//h0a->SetLineColor(4);
	//h0a->SetLineWidth(1);

	//h0a->Scale(h0b->GetMaximum() / h0a->GetMaximum(), "nosw2");
	//h0b->Sumw2();
	//h0b->Scale(h0a->Integral()/h0b->Integral());
	//h0b->SetMarkerStyle(8);
	//h0b->SetMarkerColor(1);
	//h0b->SetMarkerSize(1);
	//h0b->SetLineColor(1);
	//h0b->SetLineWidth(1);
	//h0b->SetLineStyle(1);
	//h0a->SetMinimum(0);//Outer shaft border minimum
	//h0b->SetMinimum(0);//Outer shaft border minimum
	h0a->Draw("SCAT");
	//h0a->Draw("COLZ");
	//h0b->Draw("same");

	/*
	TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
	lg1->SetHeader("#sqrt{S}=3.686GeV");
	lg1->AddEntry(h0a, "MC", "pl");
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