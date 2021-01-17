{
	TFile *f1 = new TFile("E:/Work/IHEPBOX/root/0912data/data5_eta_fdc.root");
	//======================= for mu infor ==================================
	TTree *chain = (TTree *)f1->Get("trkRes");
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

	double orig_m2gamma, orig_mlambda, orig_mlambdabar, orig_mtotal;
	double lambda_chisq, lambdabar_chisq, chisq_4C;
	double orig_m2gammalambda, orig_m2gammalambdabar, orig_mlambdalambdabar, orig_mgamma2lambdalambdabar;
	//double mdc_px[4], mdc_py[4], mdc_pz[4], mdc_p3[4], mdc_en[4];
	//double emc_px[2], emc_py[2], emc_pz[2], emc_p3[2], emc_en[2];
	double trk_px[2], trk_py[2], trk_pz[2], trk_p3[2], trk_en[2];

	chain->SetBranchAddress("trk_px", trk_px);
	chain->SetBranchAddress("trk_py", trk_py);
	chain->SetBranchAddress("trk_pz", trk_pz);
	chain->SetBranchAddress("trk_p3", trk_p3);
	chain->SetBranchAddress("trk_en", trk_en);
	/*	
	chain->SetBranchAddress("emc_px", emc_px);
	chain->SetBranchAddress("emc_py", emc_py);
	chain->SetBranchAddress("emc_pz", emc_pz);
	chain->SetBranchAddress("emc_p3", emc_p3);
	chain->SetBranchAddress("emc_en", emc_en);
	*/

	chain->SetBranchAddress("orig_m2gamma", &orig_m2gamma);
	chain->SetBranchAddress("orig_mlambda", &orig_mlambda);
	chain->SetBranchAddress("orig_mlambdabar", &orig_mlambdabar);
	chain->SetBranchAddress("lambda_chisq", &lambda_chisq);
	chain->SetBranchAddress("lambdabar_chisq", &lambdabar_chisq);
	chain->SetBranchAddress("4C_chisq", &chisq_4C);
	chain->SetBranchAddress("orig_m2gammalambda", &orig_m2gammalambda);
	chain->SetBranchAddress("orig_mlambdalambdabar", &orig_mlambdalambdabar);
	chain->SetBranchAddress("orig_mgamma2lambdalambdabar", &orig_mgamma2lambdalambdabar);

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
	int Ncut0;
	int nevent = (int)chain->GetEntries();
	for (int i = 0; i < nevent; i++)
	{
		chain->GetEntry(i);
		if (orig_mlambdabar == sqrt(trk_en[1] * trk_en[1] - (trk_px[1] * trk_px[1] + trk_py[1] * trk_py[1] + trk_pz[1] * trk_pz[1])))
		{
			continue;
		}
		cout << orig_mlambdabar << "      " << sqrt(trk_en[1] * trk_en[1] - (trk_px[1] * trk_px[1] + trk_py[1] * trk_py[1] + trk_pz[1] * trk_pz[1])) << endl;
		Ncut0++;
	}
	cout << Ncut0 << endl;
	f1->Close();
}