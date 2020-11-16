#include <math.h>
void test1()
{
    gStyle->SetOptStat(kFALSE);
    TCanvas *c1 = new TCanvas("c1", "M_Lambda.eps", 800, 600);

    TChain *ch1 = new TChain("trkRes");
    ch1->Add("../root/0912data.root");

    TChain *ch2 = new TChain("trkRes");
    ch2->Add("../root/0912data1_pi0.root");

    TChain *ch3 = new TChain("trkRes");
    ch3->Add("../root/0912data2_pi0.root");

    TChain *ch4 = new TChain("trkRes");
    ch4->Add("../root/0912data3_pi0.root");

    TChain *ch5 = new TChain("trkRes");
    ch5->Add("../root/0912data4_pi0.root");

    TChain *ch6 = new TChain("trkRes");
    ch6->Add("../../0912MCmix/root/pi0MC4.root");

    Double_t xlow1 = 0.2;
    Double_t xup1 = 1.8;
    Int_t nbins1 = 80;

    TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
    TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
    TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);
    TH1F *h1d = new TH1F("h1d", "", nbins1, xlow1, xup1);
    TH1F *h1e = new TH1F("h1e", "", nbins1, xlow1, xup1);
    TH1F *h1f = new TH1F("h1f", "", nbins1, xlow1, xup1);

    ch1->Draw("orig_mpippim2gamma>>h1a");
    ch2->Draw("orig_mpippim2gamma>>h1b");
    ch3->Draw("orig_mpippim2gamma>>h1c");
    ch4->Draw("orig_mpippim2gamma>>h1d");
    ch5->Draw("orig_mpippim2gamma>>h1e");
    ch6->Draw("orig_mpippim2gamma>>h1f");
    //ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
    //ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
    //ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
    //ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569");
    //ch1->Draw("orig_mlambda>>h1a", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569 && 4C_chisq < 37");
    //ch2->Draw("orig_mlambda>>h1b", "orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && orig_m2gamma > 0.527 && orig_m2gamma < 0.569 && 4C_chisq < 37");
    h1a->GetXaxis()->SetTitle("M(#pi^{+}#pi^{-}#pi^{0}) GeV/c^{2}");
    h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
    h1a->GetXaxis()->CenterTitle();
    h1a->GetYaxis()->CenterTitle();
    h1a->SetLineColor(4);
    h1a->SetLineWidth(1);
    /*
    h1a->Scale(h1e->GetMaximum() / h1a->GetMaximum(), "nosw2");
    h1b->Scale(h1e->GetMaximum() / h1b->GetMaximum(), "nosw2");
    h1c->Scale(h1e->GetMaximum() / h1c->GetMaximum(), "nosw2");
    h1d->Scale(h1e->GetMaximum() / h1d->GetMaximum(), "nosw2");
    h1f->Scale(h1e->GetMaximum() / h1f->GetMaximum(), "nosw2");
    */
    //h1b->Sumw2();
    //h1a->Scale(h1b->Integral()/h1a->Integral());
    //h1b->SetMarkerStyle(8);
    //h1b->SetMarkerColor(1);
    //h1b->SetMarkerSize(1);
    h1a->SetLineColor(1);
    h1b->SetLineColor(2);
    h1c->SetLineColor(3);
    h1d->SetLineColor(4);
    h1e->SetLineColor(6);
    h1f->SetLineColor(7);
    //h1b->SetLineWidth(1);
    //h1b->SetLineStyle(1);
    //h1a->SetMinimum(0);//Outer shaft border minimum
    //h1b->SetMinimum(0);//Outer shaft border minimum
    h1a->Draw("");
    h1b->Draw("SAME");
    h1c->Draw("SAME");
    h1d->Draw("SAME");
    h1e->Draw("SAME");
    h1f->Draw("SAME");

    TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
    lg1->SetHeader("#sqrt{S}=3.686GeV");
    lg1->AddEntry(h1a, "No cut", "l");
    lg1->AddEntry(h1b, "#pi^{0}", "l");
    lg1->AddEntry(h1c, "#Lambda", "l");
    lg1->AddEntry(h1d, "#bar{#Lambda}", "l");
    lg1->AddEntry(h1e, "#chi^{2}_{4C}", "l");
    lg1->AddEntry(h1f, "ExMC", "l");
    lg1->SetFillColor(0);
    lg1->SetTextFont(42);
    lg1->SetTextSize(0.03);
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