{
    gStyle->SetOptStat(kFALSE);
    TCanvas *c1 = new TCanvas("c1", "M_pi0.eps", 250, 50, 800, 600);

    TChain *ch1 = new TChain("trkRes");
    ch1->Add("../root/pi0MC3_NoGam_test.root");

    TChain *ch2 = new TChain("trkRes");
    ch2->Add("../root/Sigma0MC3_pi0_NoGam_test.root");

    int clr = 2;

    Double_t xlow1;
    Double_t xup1;
    Int_t nbins1;

    if (clr == 1)
    {
        xlow1 = 1.0;
        xup1 = 2.6;
        nbins1 = 80;
    }
    else if (clr == 2)
    {
        xlow1 = 1.0;
        xup1 = 2.0;
        nbins1 = 100;
    }

    TH1F *h1a = new TH1F("h1a", "", nbins1, xlow1, xup1);
    TH1F *h1b = new TH1F("h1b", "", nbins1, xlow1, xup1);
    TH1F *h1c = new TH1F("h1c", "", nbins1, xlow1, xup1);
    TH1F *h1d = new TH1F("h1d", "", nbins1, xlow1, xup1);

    if (clr == 1)
    {
        ch1->Draw("orig_mgamma1lambda>>h1a");
        ch2->Draw("orig_mgamma1lambda>>h1b");
        //ch3->Draw("orig_mgamma1lambda>>h1c");
        //ch4->Draw("orig_mgamma1lambda>>h1d");
    }
    else if (clr == 2)
    {
        ch1->Draw("orig_mgamma2lambda>>h1a");
        ch2->Draw("orig_mgamma2lambda>>h1b");
        //ch3->Draw("orig_mgamma2lambda>>h1c");
        //ch4->Draw("orig_mgamma2lambda>>h1d");
    }
    //ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
    //ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120");
    //ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
    //ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119");
    //ch1->Draw("orig_m2gamma>>h1a", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 41");
    //ch2->Draw("orig_m2gamma>>h1b", "orig_mlambda > 1.112 && orig_mlambda < 1.120 && orig_mlambdabar > 1.113 && orig_mlambdabar < 1.119 && 4C_chisq < 41");
    h1b->GetXaxis()->SetTitle("M(#pi^{0}) (GeV/c^{2})");
    h1b->GetYaxis()->SetTitle(Form("Events/%.5f GeV/c^{2}", (xup1 - xlow1) / nbins1));
    h1b->GetXaxis()->CenterTitle();
    h1b->GetYaxis()->CenterTitle();
    h1a->SetLineColor(4);
    h1a->SetLineWidth(1);

    h1a->Scale(h1b->Integral() / h1a->Integral(), "nosw2");
    //h1b->Sumw2();
    //h1a->Scale(h1b->GetMaximum() / h1a->GetMaximum(), "nosw2");
    //h1b->Scale(h1a->Integral()/h1b->Integral());
    h1b->SetMarkerStyle(8);
    h1b->SetMarkerColor(1);
    h1b->SetMarkerSize(1);
    h1b->SetLineColor(1);
    h1b->SetLineWidth(1);
    //h1b->SetLineStyle(1);
    //h1a->SetMinimum(0);//Outer shaft border minimum
    //h1b->SetMinimum(0);//Outer shaft border minimum
    h1b->Draw("");
    h1a->Draw("same");

    TLegend *lg1 = new TLegend(0.60, 0.68, 0.89, 0.88);
    lg1->SetHeader("#sqrt{S}=3.686GeV");
    lg1->AddEntry(h1a, "Exclusive MC", "l");
    lg1->AddEntry(h1b, "Inclusive MC", "pl");
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