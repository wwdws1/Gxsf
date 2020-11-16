#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "math.h"
using namespace std;
void OptMLambdabar()
{
	gStyle->SetPaperSize(20, 30);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.08);
	gStyle->SetPadBottomMargin(0.18);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetOptStat(0);
	TFile *f1 = new TFile("../../0912MCmix/root/pi0MC2.root");
	//TFile *f1 = new TFile("../root/MC_eta.root");
	TFile *f2 = new TFile("../root/0912data2_pi0.root");

	TTree *chain1 = (TTree *)f1->Get("trkRes");
	TTree *chain2 = (TTree *)f2->Get("trkRes");
	const int n = 200;
	const double h = 0.0002;
	double chisq1, chisq2;
	double NEntry, nentry1, nentry2;
	double s[n] = {0}, sb[n] = {0};
	double snr[n] = {0}, max = 0;
	chain1->SetBranchAddress("orig_mlambdabar", &chisq1);
	chain2->SetBranchAddress("orig_mlambdabar", &chisq2);
	nentry1 = chain1->GetEntries();
	nentry2 = chain2->GetEntries();
	for (int i = 0; i < nentry1; i++)
	{
		chain1->GetEntry(i);
		for (int j = 0; j < n; j++)
		{
			if (fabs(chisq1 - 1.115683) < j * h + h)
				s[j]++;
		}
	}
	for (int i = 0; i < nentry2; i++)
	{
		chain2->GetEntry(i);
		for (int j = 0; j < n; j++)
		{
			if (fabs(chisq2 - 1.115683) < j * h + h)
				sb[j]++;
		}
	}

	for (int j = 0; j < n; j++)
	{
		cout << "j = " << j << "\t signal = " << s[j] << "\t signal + bkg = " << sb[j] << endl;
		if (sb[j] == 0)
			snr[j] = 0;
		else
			snr[j] = s[j] / sqrt(sb[j]);
		if (snr[j] > max)
		{
			max = snr[j];
			NEntry = j;
		}
	}
	TH2F *n1 = new TH2F("", "", n, 0, n * h, 110, 0, max * 1.011);
	n1->SetXTitle(Form("M(#bar{#Lambda}) (GeV/c^{2})"));
	n1->SetYTitle(Form("#frac{s}{#sqrt{s+b}}"));
	n1->SetMarkerStyle(8);
	n1->SetMarkerSize(0.65);

	for (int i = 0; i < n; i++)
	{
		double z = snr[i];
		n1->Fill(i * h, z);
	}
	const double m = NEntry * h + h;
	cout << "NEntry= " << NEntry * h + h << endl;
	TCanvas *c1 = new TCanvas("c1", "", 100, 100, 800, 600); //Length:Width=7:5, else 3:2
	c1->Divide(1, 1);
	c1->cd(1);
	n1->Draw();

	TArrow *x1 = new TArrow(m, 0.0, m, max, 0.01, "<");
	x1->SetLineColor(2);
	x1->SetLineWidth(2);
	x1->Draw();
	TPaveText *pt = new TPaveText(0.6, 0.8, 0.85, 0.85, "BRNDC");
	pt->SetFillStyle(0);
	pt->SetBorderSize(0);
	pt->SetTextAlign(10);
	pt->SetTextSize(0.04);
	pt->Draw();
	char str[200];
	sprintf(str, "BEST = %5.5f", m);
	pt->AddText(str);
	c1->Update();
	//c1->Print("./Opt_M_Lambdabar_2.eps");
}