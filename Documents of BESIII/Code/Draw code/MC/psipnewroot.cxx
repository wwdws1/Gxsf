{
	//---------pi0----------
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/pi0MC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/pi0MC10_FitGam.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/jpsiMC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/jpsiMC10_pi0.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/LambdaSigmapi0MC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/LambdaSigmapi0MC10_pi0_FitGam.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/LambdaSigmapimMC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/LambdaSigmapimMC11_pi0.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/gammaMC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/gammaMC9_pi0.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/INMC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/INMC10_pi0_FitGam.root", "recreate");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/INMC5_pi0_NoGam_test.root", "recreate");
	//---------eta----------
	TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/etaMC.root");
	TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/etaMC8_fdc.root", "recreate");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/etafdcMC5.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/Lambda1670MC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/Lambda1670MC5_eta.root", "recreate");
	//TFile *f1 = new TFile("E:/Work/IHEPBOX/root/psip/mc/INMC.root");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/INMC4_eta_FitGam.root", "recreate");
	//TFile *newfile = new TFile("E:/Work/IHEPBOX/root/psip/mc/INMC4_eta_NoGam_test.root", "recreate");
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

	/*
	double mdc_px[4], mdc_py[4], mdc_pz[4], mdc_p3[4], mdc_en[4];
	double emc_px[2], emc_py[2], emc_pz[2], emc_p3[2], emc_en[2];
	*/

	double orig_m2gamma, orig_mlambda, orig_mlambdabar, orig_mtotal;
	double lambda_chisq, lambdabar_chisq, chisq_4C;
	double orig_m2gammalambda, orig_m2gammalambdabar, orig_mlambdalambdabar;
	double orig_mgamma2lambda, orig_mgamma2lambdabar;
	double chisq_4C_2m1, chisq_4C_2m3;
	double orig_recoilmpippim;
	double orig_m2gammapp, orig_m2gammapm;

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

	chain->SetBranchAddress("orig_m2gamma", &orig_m2gamma);
	chain->SetBranchAddress("orig_mlambda", &orig_mlambda);
	chain->SetBranchAddress("orig_mlambdabar", &orig_mlambdabar);
	chain->SetBranchAddress("lambda_chisq", &lambda_chisq);
	chain->SetBranchAddress("lambdabar_chisq", &lambdabar_chisq);
	chain->SetBranchAddress("4C_chisq", &chisq_4C);
	chain->SetBranchAddress("orig_m2gammalambda", &orig_m2gammalambda);
	chain->SetBranchAddress("orig_mlambdalambdabar", &orig_mlambdalambdabar);
	chain->SetBranchAddress("orig_mgamma2lambda", &orig_mgamma2lambda);
	chain->SetBranchAddress("orig_mgamma2lambdabar", &orig_mgamma2lambdabar);
	chain->SetBranchAddress("2m1_4C_chisq", &chisq_4C_2m1);
	chain->SetBranchAddress("2m3_4C_chisq", &chisq_4C_2m3);
	chain->SetBranchAddress("orig_recoilmpippim", &orig_recoilmpippim);
	chain->SetBranchAddress("orig_m2gammapp", &orig_m2gammapp);
	chain->SetBranchAddress("orig_m2gammapm", &orig_m2gammapm);

	//mgamma2lambdalambdabar

	TTree *t1 = chain->CloneTree(0);

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

	int Pass0 = 0, Pass1 = 0, Pass2 = 0, Pass3 = 0, Pass4 = 0, Pass5 = 0, Pass6 = 0, Pass7 = 0, Pass8 = 0, Pass9 = 0, Pass10 = 0, Pass11 = 0;
	int nevent = (int)chain->GetEntries();
	//cout << "Events" << nevent << endl;
	for (int i = 0; i < nevent; i++)
	{
		chain->GetEntry(i);

		Pass0++;

		if (orig_m2gamma < 0.115 || orig_m2gamma > 0.155) // pi0
		{
			//continue;
		}
		if (orig_m2gamma < 0.515 || orig_m2gamma > 0.570) // eta
		{
			//continue;
		}
		if (orig_m2gamma < 0.08 || orig_m2gamma > 0.19) // fit pi0
		{
			//continue;
		}
		if (orig_m2gamma < 0.46 || orig_m2gamma > 0.63) // fit eta
		{
			//continue;
		}
		if (orig_m2gamma < 0.525 || orig_m2gamma > 0.56) // fdc eta
		{
			continue;
		}

		Pass1++;

		if (orig_mlambda < 1.111 || orig_mlambda > 1.121)
		{
			continue;
		}

		Pass2++;

		if (orig_mlambdabar < 1.111 || orig_mlambdabar > 1.121)
		{
			continue;
		}

		Pass3++;

		if (chisq_4C > 15) // pi0
		{
			//continue;
		}
		if (chisq_4C > 40) // eta
		{
			continue;
		}

		Pass4++;

		if (orig_mlambdalambdabar > 3.4) // pi0
		{
			//continue;
		}

		Pass5++;

		if (orig_mlambdalambdabar > 3.077 && orig_mlambdalambdabar < 3.117)
		{
			continue;
		}

		Pass6++;

		if (chisq_4C_2m1 >= 0)
		{
			continue;
		}

		Pass7++;

		if (chisq_4C_2m3 >= 0)
		{
			continue;
		}

		Pass8++;

		if ((orig_mgamma2lambda > 1.183 && orig_mgamma2lambda < 1.203) || (orig_mgamma2lambdabar > 1.183 && orig_mgamma2lambdabar < 1.203))
		{
			//continue; // pi0
		}

		Pass9++;

		if (orig_recoilmpippim > 3.087 && orig_recoilmpippim < 3.107)
		{
			continue;
		}

		Pass10++;

		if ((orig_m2gammapp > 1.17 && orig_m2gammapp < 1.2) || (orig_m2gammapm > 1.17 && orig_m2gammapm < 1.2))
		{
			//continue; // pi0
		}

		Pass11++;

		t1->Fill();
	}
	cout << "Pass0==" << Pass0 << endl;
	cout << "Pass1==" << Pass1 << endl;
	cout << "Pass2==" << Pass2 << endl;
	cout << "Pass3==" << Pass3 << endl;
	cout << "Pass4==" << Pass4 << endl;
	cout << "Pass5==" << Pass5 << endl;
	cout << "Pass6==" << Pass6 << endl;
	cout << "Pass7==" << Pass7 << endl;
	cout << "Pass8==" << Pass8 << endl;
	cout << "Pass9==" << Pass9 << endl;
	cout << "Pass10==" << Pass10 << endl;
	cout << "Pass11==" << Pass11 << endl;

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
	newfile->Close();
}