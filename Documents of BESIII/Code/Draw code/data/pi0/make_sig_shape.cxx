{
	TH1F *hcom = new TH1F("hcom", "",290,1.89,2.035);
	TChain *chain = new TChain("ksa");
	chain->Add("../all_DsSTDs_pwa.root"); //signal

	Int_t  mode_tag,match_sig;
	Double_t  mass_sig(0),mrec_sig(0),mrec_tag(0),kschi2(0),mass_tag(0);
	chain->SetBranchAddress("match_sig",   &match_sig);
	chain->SetBranchAddress("mass_sig",    &mass_sig);
	chain->SetBranchAddress("mass_tag",   &mass_tag);
	chain->SetBranchAddress("mrec_sig",   &mrec_sig);
	chain->SetBranchAddress("mrec_tag",   &mrec_tag);

	Int_t count=0;
	Long64_t nevent = chain->GetEntries();
	for(int j=0; j<nevent; j++ ) 
	{
		chain->GetEntry(j);
		//if (mass_sig>2.035||mass_sig<1.89)continue;
		if( abs(mode_tag)==442 )continue;
		if( abs(mode_tag)==400 && (mass_tag>1.991 || mass_tag<1.948 ) )continue;
		if( abs(mode_tag)==401 && (mass_tag>1.986 || mass_tag<1.950 ) )continue;
		if( abs(mode_tag)==402 && (mass_tag>1.987 || mass_tag<1.946 ) )continue;
		if( abs(mode_tag)==502 && (mass_tag>1.983 || mass_tag<1.953 ) )continue;
		if( abs(mode_tag)==460 && (mass_tag>1.996 || mass_tag<1.940 ) )continue;
		if( abs(mode_tag)==406 && (mass_tag>1.983 || mass_tag<1.953 ) )continue;
		if( abs(mode_tag)==404 && (mass_tag>1.982 || mass_tag<1.947) )continue;
		if( abs(mode_tag)==421 && (mass_tag>1.982 || mass_tag<1.952 ) )continue;
		if( abs(mode_tag)==440 && (mass_tag>2.000 || mass_tag<1.930 ) )continue;
		if(match_sig!=1)continue;
		count++;
		hcom->Fill(mass_sig);
	}

	cout<<"count= "<<count<<endl;
	TFile f("ds_Shape.root","RECREATE");
	hcom->Write();
	f.Close();
}
