#include <math.h>
void test3()
{
    gStyle->SetOptStat(kFALSE);
    TCanvas *c1 = new TCanvas("c1", "M_pi0.eps", 250, 50, 800, 600);

    TChain *ch1 = new TChain("trkRes");
    ch1->Add("../../../0912MCmix/root/pi0MC5_NoGam_test.root");

    TChain *ch2 = new TChain("trkRes");
    ch2->Add("../../../0912MCmix/root/INMC5_pi0_NoGam_test.root");

    TChain *ch3 = new TChain("trkRes");
    ch3->Add("../../root/data5_pi0_NoGam_test.root");

    TChain *ch4 = new TChain("trkRes");
    ch4->Add("../../../0912MCmix/root/LambdaSigmapi0MC5_pi0_NoGam_test.root");
    //ch4->Add("../../../0912MCmix/root/jpsiMC5_pi0_NoGam_test.root");

    Double_t xlow1 = 0.08;
    Double_t xup1 = 0.190;
    Int_t nbins1 = 44;

    TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
    TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
    TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);
    TH1F *h1d = new TH1F("h1d", "", nbins1, xlow1, xup1);

    ch1->Draw("orig_m2gamma>>h1a");
    ch2->Draw("orig_m2gamma>>h1b");
    ch3->Draw("orig_m2gamma>>h1c");
    ch4->Draw("orig_m2gamma>>h1d");
    //ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
    //ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
    //ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
    //ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
    //ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 41");
    //ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 41");
    h1a->GetXaxis()->SetTitle("M(#gamma#gamma) (GeV/c^{2})");
    h1a->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
    h1a->GetXaxis()->CenterTitle();
    h1a->GetYaxis()->CenterTitle();
    h1a->SetLineColor(4);
    h1a->SetLineWidth(1);

    h1a->Scale(h1c->Integral() / h1a->Integral(), "nosw2");
    h1b->Scale(h1c->Integral() / h1b->Integral());
    h1c->Sumw2();
    h1d->Scale(0.5 * (h1c->Integral() / h1d->Integral()), "nosw2");
    //h1b->Scale(h1c->GetMaximum() / h1b->GetMaximum());
    //h1b->Scale(h1a->Integral()/h1b->Integral());

    h1b->SetMarkerStyle(8);
    h1b->SetMarkerColor(1);
    h1b->SetMarkerSize(1);
    h1b->SetLineColor(1);
    h1b->SetLineWidth(1);

    h1c->SetMarkerStyle(22); // 22
    h1c->SetMarkerColor(2);  // 2
    h1c->SetMarkerSize(1);
    h1c->SetLineColor(2); // 2
    h1c->SetLineWidth(1);

    h1d->SetLineColor(8);
    h1d->SetLineWidth(1);
    //h1d->SetLineStyle(7);

    //h1b->SetLineStyle(1);
    //h1a->SetMinimum(0);//Outer shaft border minimum
    //h1b->SetMinimum(0);//Outer shaft border minimum
    h1a->Draw("");
    h1b->Draw("SAME");
    h1c->Draw("SAME");
    h1d->Draw("SAME");

    TLegend *lg1 = new TLegend(0.68, 0.68, 0.86, 0.88);
    lg1->SetHeader("#sqrt{S}=3.686GeV");
    lg1->AddEntry(h1a, "EX MC", "l");
    lg1->AddEntry(h1b, "IN MC", "pl");
    lg1->AddEntry(h1c, "DATA", "pl");
    lg1->AddEntry(h1d, "BKG(#Sigma^{0}#bar{#Sigma}^{0})", "l");
    //lg1->AddEntry(h1d, "BKG(#gamma#Lambda#bar{#Lambda})", "l");
    lg1->SetFillColor(0);
    lg1->SetTextFont(42);
    lg1->SetTextSize(0.05);
    lg1->SetBorderSize(0);
    lg1->Draw();
    /*
	TArrow *x1 = new TArrow(0.115, 0, 0.115, 100000, 0.03, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(1);
	x1->Draw();

	TArrow *x2 = new TArrow(0.155, 0, 0.155, 100000, 0.03, "<");
	x2->SetLineColor(2);
	x2->SetLineWidth(1);
	x2->Draw();
	*/
    c1->Update();

    //c1->SaveAs("./M_pi0_1.eps");
}