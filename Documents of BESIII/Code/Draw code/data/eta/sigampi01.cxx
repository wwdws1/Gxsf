{
	TFile *f1 = new TFile("/scratchfs/bes/lufx/fengjh/sigmapi0_tag/4700/ana/root/MC.root");
	TFile* newfile = new TFile("mc2.root","recreate");
	//======================= for mu infor ================================== 
	TTree *chain = (TTree *)f1 -> Get("trkRes");
	if (chain == 0){
		cout << "Problem of Opening the root - tree of pid!" << endl;
	}

	TCanvas *myCanvas = new TCanvas("myCanvas");
	TPad *c1_1 = new TPad("c1_1", "c1_1",0.01,0.01,0.95,0.95);
	c1_1->Divide(1,1);
	c1_1->Draw();
	c1_1->cd(1);
	c1_1->Range(0.425458,-114.842,0.674993,802.951);
	c1_1->SetBorderSize(2);
	c1_1->SetBottomMargin(0.125129);
	c1_1->SetFrameFillColor(0);

	double ki4c_sig_px[5],ki4c_sig_py[5],ki4c_sig_pz[5],ki4c_sig_en[5];
	double ki1c_pi0_chi1;
	double chisq_lambda, DL_DLerr;
	double M_BC, Delta_E;

	chain->SetBranchAddress("ki4c_sig_px", ki4c_sig_px);
	chain->SetBranchAddress("ki4c_sig_py", ki4c_sig_py);
	chain->SetBranchAddress("ki4c_sig_pz", ki4c_sig_pz);
	chain->SetBranchAddress("ki4c_sig_en", ki4c_sig_en);
	chain->SetBranchAddress("ki1c_pi0_chi1", &ki1c_pi0_chi1);
	chain->SetBranchAddress("chisq_lambda", &chisq_lambda);
	chain->SetBranchAddress("DL_DLerr", &DL_DLerr);
	chain->SetBranchAddress("M_BC", &M_BC);
	chain->SetBranchAddress("Delta_E", &Delta_E);

        TBranch        *b_indexmc;
        TBranch        *b_pdgid;
        TBranch        *b_motheridx;
        int  mypdgid[100],mymotheridx[100];
        Int_t pdgid[100],motheridx[100],indexmc;  //[indexmc]
        Int_t kkindexmc;
         
	chain->SetBranchAddress("indexmc", &kkindexmc);
	chain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
	chain->SetBranchAddress("motheridx", motheridx, &b_motheridx);


	TTree* t1 = new TTree("mytree","clone tree");
	double M_lambda, M_sigma, M_pi0, M_lambdac;
	t1->Branch("ki1c_pi0_chi1",&ki1c_pi0_chi1,"ki1c_pi0_chi1/D");
	t1->Branch("chisq_lambda",&chisq_lambda,"chisq_lambda/D");
	t1->Branch("DL_DLerr",&DL_DLerr,"DL_DLerr/D");
	t1->Branch("M_BC",&M_BC,"M_BC/D");
	t1->Branch("Delta_E",&Delta_E,"Delta_E/D");
	t1->Branch("M_lambda",&M_lambda,"M_lambda/D");
	t1->Branch("M_sigma",&M_sigma,"M_sigma/D");
	t1->Branch("M_pi0",&M_pi0,"M_pi0/D");
	t1->Branch("M_lambdac",&M_lambdac,"M_lambdac/D");

        t1->Branch("indexmc"  ,&indexmc,"indexmc/I");
        t1->Branch("pdgid"    ,mypdgid,"mypdgid[indexmc]/I");
        t1->Branch("motheridx",mymotheridx,"mymotheridx[indexmc]/I");

	TH2F *M_sigma_pi0  = new TH2F("M_sigma_pi0","",100,1.2837,1.4837,  100,0.115,0.150);
	TLorentzVector son3,son4,son5,son6,son7;
        int Pass0 =0,Pass1 =0,Pass2 =0,Pass3 =0,Pass4 =0, Pass5 =0, Pass6 =0, Pass7 =0;
	int nevent =(int) chain->GetEntries(); cout<<"Events"<<nevent<<endl;
	for(int i=0; i<nevent; i++ ){
		chain->GetEntry(i);
		//      ppp.SetPxPyPzE(m_px_ki4c_sig[0],m_py_ki4c_sig[0],m_pz_ki4c_sig[0],m_en_ki4c_sig[0]);
		//      ppim.SetPxPyPzE(m_px_ki4c_sig[1],m_py_ki4c_sig[1],m_pz_ki4c_sig[1],m_en_ki4c_sig[1]);
		//      ppm.SetPxPyPzE(m_px_ki4c_sig[2],m_py_ki4c_sig[2],m_pz_ki4c_sig[2],m_en_ki4c_sig[2]);
		//      ppip.SetPxPyPzE(m_px_ki4c_sig[3],m_py_ki4c_sig[3],m_pz_ki4c_sig[3],m_en_ki4c_sig[3]);
		//      pgam.SetPxPyPzE(m_px_ki4c_sig[4],m_py_ki4c_sig[4],m_pz_ki4c_sig[4],m_en_ki4c_sig[4]);

		son3.SetPx(ki4c_sig_px[0]); son3.SetPy(ki4c_sig_py[0]); son3.SetPz(ki4c_sig_pz[0]); son3.SetE(ki4c_sig_en[0]);//p
		son4.SetPx(ki4c_sig_px[1]); son4.SetPy(ki4c_sig_py[1]); son4.SetPz(ki4c_sig_pz[1]); son4.SetE(ki4c_sig_en[1]);//pi-
		son5.SetPx(ki4c_sig_px[2]); son5.SetPy(ki4c_sig_py[2]); son5.SetPz(ki4c_sig_pz[2]); son5.SetE(ki4c_sig_en[2]);//pi+
		son6.SetPx(ki4c_sig_px[3]); son6.SetPy(ki4c_sig_py[3]); son6.SetPz(ki4c_sig_pz[3]); son6.SetE(ki4c_sig_en[3]);//pi0
		son7.SetPx(ki4c_sig_px[4]); son7.SetPy(ki4c_sig_py[4]); son7.SetPz(ki4c_sig_pz[4]); son7.SetE(ki4c_sig_en[4]);//lambda
                Pass0++;
		// cout<<" son3 "<<" px = "<<son3.Px()<<"\t  py = "<<son3.Py()<<"\t  pz = "<<son3.Pz()<<"\t  en = "<<son3.E()<<endl;
		// cout<<" son4 "<<" px = "<<son4.Px()<<"\t  py = "<<son4.Py()<<"\t  pz = "<<son4.Pz()<<"\t  en = "<<son4.E()<<endl;
		// cout<<" son5 "<<" px = "<<son5.Px()<<"\t  py = "<<son5.Py()<<"\t  pz = "<<son5.Pz()<<"\t  en = "<<son5.E()<<endl;
		// cout<<" son6 "<<" px = "<<son6.Px()<<"\t  py = "<<son6.Py()<<"\t  pz = "<<son6.Pz()<<"\t  en = "<<son6.E()<<endl;
		// cout<<" son7 "<<" px = "<<son7.Px()<<"\t  py = "<<son7.Py()<<"\t  pz = "<<son7.Pz()<<"\t  en = "<<son7.E()<<endl;
	
         	if(ki1c_pi0_chi1>200) continue; Pass1++;
                if(chisq_lambda>100) continue;  Pass2++;
                if(DL_DLerr<2) continue;        Pass3++;
                M_lambda = (son3 + son4).M();   
                if(M_lambda<1.111 || M_lambda>1.121) continue; Pass4++;
		M_sigma = (son3 + son4 + son5).M() - (son3 + son4).M() + 1.115683;
                if(M_sigma<1.2828 || M_sigma>1.4828) continue;   Pass5++;
		M_pi0 = son6.M();
                TLorentzVector p4sigLc = son3 + son4 + son5 + son6;
                p4sigLc.Boost(-0.011,0,0); double Ebeam = 2.35022;
                M_lambdac = sqrt( Ebeam*Ebeam - p4sigLc.P()*p4sigLc.P() );
		//M_lambdac = (son3 + son4 + son5 + son6).M() - (son3 + son4 + son5).M() + 1.3828;
                if(M_lambdac<2.25 || M_lambdac>2.36) continue;  Pass6++;
		M_sigma_pi0->Fill(M_sigma,M_pi0);

                for(int mcnum=0; mcnum<kkindexmc; mcnum++){
                        mypdgid[mcnum]     = pdgid[mcnum];
                        mymotheridx[mcnum] = motheridx[mcnum];
                } //topology block
                indexmc=kkindexmc;
          
               t1->Fill();
	}
        cout<<"Pass0=="<<Pass0<<endl;
        cout<<"Pass1=="<<Pass1<<endl;
        cout<<"Pass2=="<<Pass2<<endl;
        cout<<"Pass3=="<<Pass3<<endl;
        cout<<"Pass4=="<<Pass4<<endl;
        cout<<"Pass5=="<<Pass5<<endl;
        cout<<"Pass6=="<<Pass6<<endl;

	c1_1->Update();
	//   myCanvas->cd(1);
	M_sigma_pi0->Draw("colz");
        M_sigma_pi0->GetXaxis()->SetTitle("m_{#Sigma^{*+}} (GeV)");
        M_sigma_pi0->GetYaxis()->SetTitle("m_{#pi^{0}} (GeV)");
	M_sigma_pi0->SetTitle("");

	t1->Write();
	f1->Close(); 
}

