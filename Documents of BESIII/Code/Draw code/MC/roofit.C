#include <fstream> 
#include <sstream>
#include <iostream>
#include <math.h>
#include "TGraph.h" 
#include "TLatex.h" 
#include "TLegend.h" 
#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "myfit.h"  
void print_4momentums(TLorentzVector son);
int main(int argc, char *argv[]) 
{ 
  Int_t dataflag = atoi(argv[1]); //dataflag = 1 (data), 0 (MC) 
  // ---- define trees  
  Double_t ki4c_sig_chi2,ki4c_2pi0_pi02cos; 
  Double_t mphi1,mphi2,mphi3,metac, gammarec,gammaE,mchicj; 
  Double_t cos1,cos2,cos3,phicos1,phicos2,phicos3; 
  Double_t RunNo, EvtRec; 
  Double_t PxMCtruth[100], PyMCtruth[100], PzMCtruth[100], EnMCtruth[100];
  Int_t xtopo,sigidx,xtag; 

  Double_t ECMS; // energy in center of mass system

  stringstream InputROOT,  OutputROOT, TopoFile;  
  InputROOT.str("");
  InputROOT.clear();
  OutputROOT.str("");
  OutputROOT.clear();
  TopoFile.str("");
  TopoFile.clear();
  // --- chain the file  
  TChain *chain0 = new TChain("trkRes");  
  TChain *chain1 = new TChain("topo"); 
  chain0->AddFriend("topo"); 

  if(dataflag ==4009 ){              //  MC@4180MeV  
    ECMS = 4.009;
    chain0->Add("/scratchfs/bes/peijh/OmegaEta/4009/mc/root/*.root"); 
    std::cout<<"End of data reading for @etac4k "<<dataflag<<std::endl; 
    OutputROOT << "omega_exmc.root";
  }else if(dataflag ==40091 ){      //  Data@4180MeV  
    ECMS = 4.009;
    chain0->Add("/scratchfs/bes/peijh/4C_data/4009/root/*.root");
    std::cout<<"End of data reading for @etac4k "<<dataflag<<std::endl;
    OutputROOT << "omega_data.root";
  } 
  //Making several pictures in the same Postscript file  
  TFile file(OutputROOT.str().c_str(),"recreate");  
  //================ style sets don't touch it 
#include "sty.C" 
  /////////// end of style set     
  Double_t ki4c_sig_px[14], ki4c_sig_py[14], ki4c_sig_pz[14], ki4c_sig_en[14];  
  chain0->SetBranchAddress("ki4c_sig_px",ki4c_sig_px);  
  chain0->SetBranchAddress("ki4c_sig_py",ki4c_sig_py);  
  chain0->SetBranchAddress("ki4c_sig_pz",ki4c_sig_pz);  
  chain0->SetBranchAddress("ki4c_sig_en",ki4c_sig_en);  
  chain0->SetBranchAddress("ki4c_sig_chi2",&ki4c_sig_chi2);  
  chain0->SetBranchAddress("ki4c_2pi0_pi02cos",&ki4c_2pi0_pi02cos); 
  // To read out the flexible dimesion array  
  TBranch        *b_indexmc;  
  TBranch        *b_pdgid;    
  TBranch        *b_motheridx;
  //--- tree to writ out
  int  mypdgid[100],mymotheridx[100];
  double mydang[100], dang[100], mygam_TDC1, mygam_TDC2; //mydang: the minimum of the angles between one neutralTrack and all the chargedTracks
  Int_t pdgid[100],motheridx[100],indexmc;  //[indexmc]
  Int_t kkindexmc;
  //   ecMdcTrack
  //-- 
  //chain0->SetBranchAddress("m_best_indexmc", &kkindexmc);
  //chain0->SetBranchAddress("pdgid", pdgid, &b_pdgid);
  //chain0->SetBranchAddress("motherIdx", motheridx, &b_motheridx);

  chain0->SetBranchAddress("m_best_indexmc", &kkindexmc);
  chain0->SetBranchAddress("m_best_pdgid", pdgid, &b_pdgid);
  chain0->SetBranchAddress("m_best_motherIdx", motheridx, &b_motheridx);
  

  double Meta,Mpi0,Momega,coseta_ecms,eta_ecmsAng,eta_cms;
  TLorentzVector son1,son2,son3,son4,son5,son6,son7,son8,son9,son10,phi1,phi2,etac;  
  TTree m3("m3","invariant mass of two body");  
  m3.Branch("ki4c_sig_chi2", &ki4c_sig_chi2, "ki4c_sig_chi2/D");
  m3.Branch("Momega", &Momega, "Momega/D"); // before 4C
  m3.Branch("Meta",&Meta,"Meta/D");  
  m3.Branch("Mpi0",&Mpi0,"Mpi0/D");  // after 4C
  m3.Branch("indexmc",&indexmc,"indexmc/I");
  m3.Branch("pdgid",mypdgid,"mypdgid[100]/I");  
  m3.Branch("motheridx",mymotheridx,"mymotheridx[100]/I");
  //m3.Branch("coseta_ecms", &coseta_ecms, "coseta_ecms/D"); 
  // for data  
  TH1F *hchisq = new TH1F("ki4c_sig_chi2","",200,0.0,200);  
  TH1F *hcoseta_ecms = new TH1F("coseta_ecms","",0.0,0.0,1.0);  
  TH1F *hgg34 = new TH1F("Meta","",50,0.0,0.3);  
  TH1F *hgg56 = new TH1F("Mpi0","",50,0.4,0.7);  
  TH1F *homega = new TH1F("Momega","",80,0.6,0.9);  
  TH2F *hb = new TH2F("MetaMomega",  "",100,0.4,0.7,100,0.6,0.9);  

  //   for output data to a text file 
  TTree *newtree; 
  //---------- 
  std::vector<int> Vsig;  
  Int_t nentries = (Int_t)chain0->GetEntries();     //Data  
  ofstream TopoRec;
  TopoRec.open(TopoFile.str().c_str(), ios_base::out);
  int Pass =0,Pass1 =0,Pass2 =0,Pass3 =0,Pass4 =0;
  //execute
  for(Int_t j=0;j<nentries;j++)  
  {  
    chain0->GetEntry(j);  
    Vsig.clear(); 
    // mometum for track list: K+ K- photon photon 
    son1.SetPx(ki4c_sig_px[0]); son1.SetPy(ki4c_sig_py[0]); son1.SetPz(ki4c_sig_pz[0]); son1.SetE(ki4c_sig_en[0]);//pi+
    son2.SetPx(ki4c_sig_px[1]); son2.SetPy(ki4c_sig_py[1]); son2.SetPz(ki4c_sig_pz[1]); son2.SetE(ki4c_sig_en[1]);//pi-
    son3.SetPx(ki4c_sig_px[2]); son3.SetPy(ki4c_sig_py[2]); son3.SetPz(ki4c_sig_pz[2]); son3.SetE(ki4c_sig_en[2]);//gamma1
    son4.SetPx(ki4c_sig_px[3]); son4.SetPy(ki4c_sig_py[3]); son4.SetPz(ki4c_sig_pz[3]); son4.SetE(ki4c_sig_en[3]);//gamma2
    son5.SetPx(ki4c_sig_px[4]); son5.SetPy(ki4c_sig_py[4]); son5.SetPz(ki4c_sig_pz[4]); son5.SetE(ki4c_sig_en[4]);//gamma3
    son6.SetPx(ki4c_sig_px[5]); son6.SetPy(ki4c_sig_py[5]); son6.SetPz(ki4c_sig_pz[5]); son6.SetE(ki4c_sig_en[5]);//gamma4
    //good charged track selection
    //
    bool cflag = TrackSelection(son1,son2);
    if(!cflag) continue;
    Mpi0 = (son3+son4).M();
    Meta = (son5+son6).M();
    Momega = (son1+son2+son3+son4).M();

      if(ki4c_sig_chi2>50) continue;
      Pass1++;
      if(TMath::Abs(Meta-0.5479)>0.03) continue;
      if(TMath::Abs(Mpi0-0.135)>0.02) continue;
      Pass2++;
      //if(Momega>0.9) continue;
      //if(Momega<0.65) continue;
      Pass3++;

    hchisq->Fill(ki4c_sig_chi2);
    hgg34->Fill(Mpi0);         
    hgg56->Fill(Meta);         
    homega->Fill(Momega);              
    hb->Fill(Meta,Momega);
    hcoseta_ecms->Fill(coseta_ecms);
    Pass++;
    m3.Fill();
    //=========== fill MC truth
    if(dataflag<0){
      //TopoRec << "RunNo.: " << Run <<"; EventNo.: " << rec << endl;
       for(int i=0;i<kkindexmc;i++){
        mypdgid[i]=pdgid[i];mymotheridx[i]=motheridx[i];
        TopoRec << "indexmc= " << kkindexmc
          <<  "; pdgid[" <<i<< "]=" << pdgid[i]
          << "; motheridx[" << i << "]=" << motheridx[i]
          << "; (" << PxMCtruth[i]
          << ", " << PyMCtruth[i]
          << ", " << PzMCtruth[i]
          << ", " << EnMCtruth[i]
          << ")" << endl;
      }
      indexmc = kkindexmc;
    }
    for(int mcnum=0; mcnum<kkindexmc; mcnum++){
      mypdgid[mcnum]     = pdgid[mcnum];
      mymotheridx[mcnum] = motheridx[mcnum];
    } //topology block
    indexmc=kkindexmc;

  } //nentries  (师妹，你的事例挑选从这里就循环完了，所以，你的所以的cut条件都要在这里之前完成)
  //end of looping over event
  TopoRec << "Total number of Events which passed Cuts: " << Pass <<endl; 
  TopoRec.close();
  m3.Write();

  double a =Pass1*1.0/nentries; double b =Pass2*1.0/Pass1; double c =Pass3*1.0/Pass2; double d =Pass4*1.0/Pass3;
  double efficiency = Pass/95000.;
  std::cout<<"total events: "<<nentries<<std::endl; 
  std::cout<<"chi  cut event :   "<<Pass1<< "   mkk/src  :"<<a<< std::endl; 
  std::cout<<"pi0eta cut event :   "<<Pass2<< "    chi/mkk  :"<<b<<std::endl; 
  std::cout<<"Momega cut event :   "<<Pass3<< "    Meta/chi   :" <<c<<std::endl; 
//  std::cout<<"gammaAng cut event :   "<<Pass4<< "     gammaAng/Meta :"<<d<<std::endl; 
  std::cout<<"pass cut event :   "<<Pass<<std::endl; 
  std::cout<<"efficiency  :   "<<efficiency<<std::endl; 

  return 1; 
}  

void print_4momentums(TLorentzVector son){

  cout<<"px = "<<son.Px();
  cout<<"\t  py = "<<son.Py();
  cout<<"\t  pz = "<<son.Pz();
  cout<<"\t  en = "<<son.E();
  cout<<endl;

}

