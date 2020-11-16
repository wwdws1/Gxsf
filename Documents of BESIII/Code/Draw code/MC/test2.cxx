#include <math.h>

void test2()
{
	TFile *f1 = new TFile("../root/test.root");
	TFile *newfile = new TFile("../root/test1.root", "recreate");
	//======================= for mu infor ==================================
	TTree *chain = (TTree *)f1->Get("mctruth");
	if (chain == 0)
	{
		cout << "Problem of Opening the root - tree of pid!" << endl;
	}

	/*
	TCanvas *myCanvas = new TCanvas("myCanvas");
	TPad *c1_1 = new TPad("c1_1", "c1_1", 0.01, 0.01, 0.95, 0.95);
	c1_1->Divide(1, 1);
	c1_1->Draw();
	c1_1->cd(1);
	c1_1->Range(0.425458, -114.842, 0.674993, 802.951);
	c1_1->SetBorderSize(2);
	c1_1->SetBottomMargin(0.125129);
	c1_1->SetFrameFillColor(0);
	*/

	/*
	double mdc_px[4], mdc_py[4], mdc_pz[4], mdc_p3[4], mdc_en[4];
	double emc_px[2], emc_py[2], emc_pz[2], emc_p3[2], emc_en[2];
	*/

	//double mc_pp_px, mc_pp_py, mc_pp_pz;
	//double mc_pm_px, mc_pm_py, mc_pm_pz;
	double mc_pip1_px, mc_pip1_py, mc_pip1_pz;
	double mc_pim2_px, mc_pim2_py, mc_pim2_pz;
	double mc_pip1_p3, mc_pim2_p3;

	/*
	chain->SetBranchAddress("mdc_px", mdc_px);
	chain->SetBranchAddress("mdc_py", mdc_py);
	chain->SetBranchAddress("mdc_pz", mdc_pz);
	chain->SetBranchAddress("mdc_p3", mdc_p3);
	chain->SetBranchAddress("mdc_en", mdc_en);
	chain->SetBranchAddress("emc_px", emc_px);
	chain->SetBranchAddress("emc_py", emc_py);
	chain->SetBranchAddress("emc_pz", emc_pz);
	chain->SetBranchAddress("emc_p3", emc_p3);
	chain->SetBranchAddress("emc_en", emc_en);
	*/

	//chain->SetBranchAddress("mc_pp_px", &mc_pp_px);
	//chain->SetBranchAddress("mc_pp_py", &mc_pp_py);
	//chain->SetBranchAddress("mc_pp_pz", &mc_pp_pz);
	//chain->SetBranchAddress("mc_pm_px", &mc_pm_px);
	//chain->SetBranchAddress("mc_pm_py", &mc_pm_py);
	//chain->SetBranchAddress("mc_pm_pz", &mc_pm_pz);
	chain->SetBranchAddress("mc_pip1_px", &mc_pip1_px);
	chain->SetBranchAddress("mc_pip1_py", &mc_pip1_py);
	chain->SetBranchAddress("mc_pip1_pz", &mc_pip1_pz);
	chain->SetBranchAddress("mc_pim2_px", &mc_pim2_px);
	chain->SetBranchAddress("mc_pim2_py", &mc_pim2_py);
	chain->SetBranchAddress("mc_pim2_pz", &mc_pim2_pz);

	TTree *t1 = chain->CloneTree(0);

	//t1->Branch("mc_pp_p3", &mc_pp_p3, "mc_pp_p3/D");
	//t1->Branch("mc_pm_p3", &mc_pm_p3, "mc_pm_p3/D");
	t1->Branch("mc_pip1_p3", &mc_pip1_p3, "mc_pip1_p3/D");
	t1->Branch("mc_pim2_p3", &mc_pim2_p3, "mc_pim2_p3/D");

	/*
	t1->Branch("orig_m2gamma", &orig_m2gamma, "orig_m2gamma/D");
	t1->Branch("orig_mlambda", &orig_mlambda, "orig_mlambda/D");
	t1->Branch("orig_mlambdabar", &orig_mlambdabar, "orig_mlambdabar/D");
	t1->Branch("orig_mtotal", &orig_mtotal, "orig_mtotal/D");
	t1->Branch("lambda_chisq", &lambda_chisq, "lambda_chisq/D");
	t1->Branch("lambdab_chisq", &lambdab_chisq, "lambdab_chisq/D");
	t1->Branch("ki_chi2", &ki_chi2, "ki_chi2/D");
	t1->Branch("lambda_lth", &lambda_lth, "lambda_lth/D");
	t1->Branch("lambda_lth_e", &lambda_lth_e, "lambda_lth_e/D");
	t1->Branch("lambdab_lth", &lambdab_lth, "lambdab_lth/D");
	t1->Branch("lambdab_lth_e", &lambdab_lth_e, "lambdab_lth_e/D");
	t1->Branch("lambdab_lth_re", &lambdab_lth_re, "lambdab_lth_re/D");
	t1->Branch("lambda_lth_re", &lambdab_lth_re, "lambda_lth_re/D");
	*/

	int Pass0 = 0, Pass1 = 0, Pass2 = 0, Pass3 = 0, Pass4 = 0, Pass5 = 0, Pass6 = 0, Pass7 = 0;
	int nevent = (int)chain->GetEntries();
	//cout << "Events" << nevent << endl;
	for (int i = 0; i < nevent; i++)
	{
		chain->GetEntry(i);
		//_pp_p3 = sqrt(pow(mc_pp_px, 2) + pow(mc_pp_py, 2) + pow(mc_pp_pz, 2));
		//mc_pm_p3 = sqrt(pow(mc_pm_px, 2) + pow(mc_pm_py, 2) + pow(mc_pm_pz, 2));
		mc_pip1_p3 = sqrt(pow(mc_pip1_px, 2) + pow(mc_pip1_py, 2) + pow(mc_pip1_pz, 2));
		mc_pim2_p3 = sqrt(pow(mc_pim2_px, 2) + pow(mc_pim2_py, 2) + pow(mc_pim2_pz, 2));

		//mc_pp_p3 = mc_pp_px;
		//mc_pm_p3 = mc_pm_px;
		//mc_pip_p3 = mc_pip_px;
		//mc_pim_p3 = mc_pim_px;

		t1->Fill();
	}
	//cout << "Pass0==" << Pass0 << endl;
	//cout << "Pass1==" << Pass1 << endl;
	//cout << "Pass2==" << Pass2 << endl;
	//cout << "Pass3==" << Pass3 << endl;
	//cout << "Pass4==" << Pass4 << endl;
	//cout << "Pass5==" << Pass5 << endl;
	//cout << "Pass6==" << Pass6 << endl;

	/*
	c1_1->Update();
	//   myCanvas->cd(1);
	M_sigma_pi0->Draw("colz");
	M_sigma_pi0->GetXaxis()->SetTitle("m_{#Sigma^{*+}} (GeV)");
	M_sigma_pi0->GetYaxis()->SetTitle("m_{#pi^{0}} (GeV)");
	M_sigma_pi0->SetTitle("");
	*/

	t1->Write();
	f1->Close();
}