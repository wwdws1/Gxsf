#ifndef Physics_Analysis_psip_H
#define Physics_Analysis_psip_H

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "VertexFit/WTrackParameter.h"
#include <vector>

class psip : public Algorithm
{

public:
    psip(const std::string &name, ISvcLocator *pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

private:
    void InitVar();
    void corgen(HepMatrix &, HepVector &, int);
    void corset(HepSymMatrix &, HepMatrix &, int);
    void calibration(RecMdcKalTrack *, HepVector &, int);

    double m_vr0cut;
    double m_vz0cut;

    //Declare energy, dphi, dthe cuts for fake gamma's
    double m_energyThresholdB;
    double m_energyThresholdE;
    //
    int m_test4C;
    int m_test5C;
    bool m_datacont;

    //ReadBeamParFromDb m_reader;
    double m_ecms;
    // Declare r0, z0 cut for charged tracks
    //Declare energy, dphi, dthe cuts for fake gamma's
    double m_energyThreshold;
    double m_gammaPhiCut;
    double m_gammaThetaCut;
    double m_gammaAngleCut;
    //
    int m_checkDedx;
    int m_checkTof;
    //
    int m_test4k4C;
    int m_test4k5C;
    int m_test4k5C_4C;
    int m_testmctruth;

    NTuple::Tuple *m_tuple11;
    NTuple::Item<double> parti_px_mc1;
    NTuple::Item<double> parti_py_mc1;
    NTuple::Item<double> parti_pz_mc1;
    NTuple::Item<double> parti_e_mc1;
    NTuple::Item<double> parti_m_mc1;

    NTuple::Item<double> parti_px_mc2;
    NTuple::Item<double> parti_py_mc2;
    NTuple::Item<double> parti_pz_mc2;
    NTuple::Item<double> parti_e_mc2;
    NTuple::Item<double> parti_m_mc2;

    NTuple::Item<double> parti_px_mc3;
    NTuple::Item<double> parti_py_mc3;
    NTuple::Item<double> parti_pz_mc3;
    NTuple::Item<double> parti_e_mc3;
    NTuple::Item<double> parti_m_mc3;

    NTuple::Item<double> parti_px_mc4;
    NTuple::Item<double> parti_py_mc4;
    NTuple::Item<double> parti_pz_mc4;
    NTuple::Item<double> parti_e_mc4;
    NTuple::Item<double> parti_m_mc4;

    NTuple::Item<double> parti_px_mc11;
    NTuple::Item<double> parti_py_mc11;
    NTuple::Item<double> parti_pz_mc11;
    NTuple::Item<double> parti_e_mc11;
    NTuple::Item<double> parti_m_mc11;

    NTuple::Item<double> parti_px_mc12;
    NTuple::Item<double> parti_py_mc12;
    NTuple::Item<double> parti_pz_mc12;
    NTuple::Item<double> parti_e_mc12;
    NTuple::Item<double> parti_m_mc12;

    NTuple::Item<double> parti_px_mc21;
    NTuple::Item<double> parti_py_mc21;
    NTuple::Item<double> parti_pz_mc21;
    NTuple::Item<double> parti_e_mc21;
    NTuple::Item<double> parti_m_mc21;

    NTuple::Item<double> parti_px_mc22;
    NTuple::Item<double> parti_py_mc22;
    NTuple::Item<double> parti_pz_mc22;
    NTuple::Item<double> parti_e_mc22;
    NTuple::Item<double> parti_m_mc22;

    NTuple::Item<double> parti_px_mcgamma;
    NTuple::Item<double> parti_py_mcgamma;
    NTuple::Item<double> parti_pz_mcgamma;
    NTuple::Item<double> parti_e_mcgamma;
    NTuple::Item<double> parti_m_mcgamma;

    NTuple::Item<double> m_tstlambar;
    NTuple::Item<double> m_tstlam;

    NTuple::Item<int> m_mc_nlam;
    NTuple::Item<int> m_mc_nlambar;
    NTuple::Item<int> m_mc_npi0;
    NTuple::Item<int> m_mc_neta;
    NTuple::Item<int> m_mc_npp;
    NTuple::Item<int> m_mc_npim;
    NTuple::Item<int> m_mc_npm;
    NTuple::Item<int> m_mc_npip;
    NTuple::Item<int> m_mc_ngamma;

    NTuple::Tuple *m_tuple6; // 4c

    NTuple::Item<long> m_run; // run number
    NTuple::Item<long> m_rec;

    NTuple::Item<long> m_nps;
    NTuple::Array<long> m_pid;
    NTuple::Array<long> m_midx;
    NTuple::Array<long> m_motherid;

    NTuple::Item<long> m_idxmc1;
    NTuple::Array<long> m_pdgid1;
    NTuple::Array<long> m_motherid1;

    NTuple::Array<long> mc_tag;
    NTuple::Array<long> vert_mcpar1;

    NTuple::Item<long> m_ichannel;

    NTuple::Item<double> m_mass_mc;
    NTuple::Item<double> m_phi;

    NTuple::Array<double> m_vtx_Rxy;
    NTuple::Array<double> m_vtx_z0;
    NTuple::Array<double> m_vtx_cos;
    NTuple::Array<double> m_vtx_p;
    NTuple::Item<double> m_totE_ch;
    NTuple::Item<int> m_nKp;
    NTuple::Item<int> m_nKm;
    NTuple::Item<int> m_npp;
    NTuple::Item<int> m_npm;
    NTuple::Item<int> m_npip;
    NTuple::Item<int> m_npim;

    NTuple::Array<double> m_R;
    NTuple::Item<long> chargTrk_index;

    NTuple::Item<double> chi2_lamb;
    NTuple::Item<double> chi2_lambar;
    NTuple::Array<double> m_px_mdc;
    NTuple::Array<double> m_py_mdc;
    NTuple::Array<double> m_pz_mdc;
    NTuple::Array<double> m_p3_mdc;
    NTuple::Array<double> m_en_mdc;
    NTuple::Array<double> m_px_emc;
    NTuple::Array<double> m_py_emc;
    NTuple::Array<double> m_pz_emc;
    NTuple::Array<double> m_p3_emc;
    NTuple::Array<double> m_en_emc;
    NTuple::Array<double> m_px_mdc_lamfit;
    NTuple::Array<double> m_py_mdc_lamfit;
    NTuple::Array<double> m_pz_mdc_lamfit;
    NTuple::Array<double> m_p3_mdc_lamfit;
    NTuple::Array<double> m_en_mdc_lamfit;
    NTuple::Array<double> m_px_mdc_ori;
    NTuple::Array<double> m_py_mdc_ori;
    NTuple::Array<double> m_pz_mdc_ori;
    NTuple::Array<double> m_p3_mdc_ori;
    NTuple::Array<double> m_en_mdc_ori;
    NTuple::Array<double> m_pid_e;
    NTuple::Array<double> m_pid_u;
    NTuple::Array<double> m_pid_pi;
    NTuple::Array<double> m_pid_k;
    NTuple::Array<double> m_pid_p;

    NTuple::Array<double> m_TDC;
    NTuple::Array<double> m_dthe;
    NTuple::Array<double> m_dphi;
    NTuple::Array<double> m_dang;
    NTuple::Array<double> m_eraw;
    NTuple::Array<double> m_tagg;

    NTuple::Item<long> gamma_index;

    NTuple::Item<double> dthe_m_4c;
    NTuple::Item<double> dang_m_4c;
    NTuple::Item<double> dphi_m_4c;
    NTuple::Item<double> eraw_m_4c;
    NTuple::Item<double> TDC_m_4c;
    NTuple::Item<double> tagg_m_4c;

    NTuple::Array<double> m_pull_1;
    NTuple::Array<double> m_pull_2;
    NTuple::Array<double> m_pull_3;
    NTuple::Array<double> m_pull_4;
    NTuple::Array<double> m_pull_5;

    NTuple::Item<int> m_ngLam1;
    NTuple::Item<int> m_ngLam2;

    NTuple::Item<double> m_mlam_vert;
    NTuple::Item<double> m_mlambar_vert;
    NTuple::Item<double> m_mlamlambar_vert;
    NTuple::Item<double> m_reco_g;

    NTuple::Item<double> m_secchis1;
    NTuple::Item<double> m_mlengthlam;
    NTuple::Item<double> m_mlengthlam_err;
    NTuple::Item<double> m_mlengthlam_rerr;

    NTuple::Item<double> m_secchis2;
    NTuple::Item<double> m_mlengthlamb;
    NTuple::Item<double> m_mlengthlamb_err;
    NTuple::Item<double> m_mlengthlamb_rerr;

    NTuple::Item<double> m_secchis12diff1;
    NTuple::Item<double> m_secchis12diff2;
    NTuple::Item<double> m_secchis12diff3;
    NTuple::Item<double> m_secchis12diff4;
    NTuple::Item<double> m_secchis12diff5;

    NTuple::Item<double> m_chi1;
    NTuple::Item<double> m_chi1_1g;
    NTuple::Item<double> m_chi1_0g;
    NTuple::Item<double> m_chi1_2g;

    NTuple::Array<double> m_px_trk;
    NTuple::Array<double> m_py_trk;
    NTuple::Array<double> m_pz_trk;
    NTuple::Array<double> m_p3_trk;
    NTuple::Array<double> m_en_trk;

    NTuple::Array<double> m_px_ptc;
    NTuple::Array<double> m_py_ptc;
    NTuple::Array<double> m_pz_ptc;
    NTuple::Array<double> m_p3_ptc;
    NTuple::Array<double> m_en_ptc;

    NTuple::Array<double> m_px_trk_1g;
    NTuple::Array<double> m_py_trk_1g;
    NTuple::Array<double> m_pz_trk_1g;
    NTuple::Array<double> m_p3_trk_1g;
    NTuple::Array<double> m_en_trk_1g;

    NTuple::Item<double> m_pppm;
    NTuple::Item<double> m_pippim;
    NTuple::Item<double> m_m2p_ki;
    NTuple::Item<double> m_m2pi_ki;
    NTuple::Item<double> m_phi1_true;
    NTuple::Item<double> m_phi2_true;
    NTuple::Item<double> m_etac_true;
    NTuple::Item<double> m_gphi1;
    NTuple::Item<double> m_gphi2;

    NTuple::Item<double> m_pipm_fit;
    NTuple::Item<double> m_pppm_fit;
    NTuple::Item<double> m_2lam_fit;
    NTuple::Item<double> m_lam_fit1;
    NTuple::Item<double> m_lam_fit2;
    NTuple::Item<double> m_2phi_bk;
    NTuple::Item<double> m_phi1_bk1;
    NTuple::Item<double> m_phi2_bk1;
    NTuple::Item<double> m_totm_bk;
    NTuple::Item<double> m_2gam_bk1;
    NTuple::Item<double> m_lam0_bk1;
    NTuple::Item<double> m_lam1_bk1;
    NTuple::Item<double> m_gphi1_bk1;
    NTuple::Item<double> m_gphi2_bk1;
    NTuple::Item<double> m_phi1_bk;
    NTuple::Item<double> m_phi2_bk;

    NTuple::Item<double> m_mix0_bk1;
    NTuple::Item<double> m_mix1_bk1;
    NTuple::Item<double> m_mix2_bk1;
    NTuple::Item<double> m_mix3_bk1;
    NTuple::Item<double> m_mix4_bk1;
    NTuple::Item<double> m_mix5_bk1;
    NTuple::Item<double> m_mix6_bk1;
    NTuple::Item<double> m_mix7_bk1;
    NTuple::Item<double> m_mix8_bk1;
    NTuple::Item<double> m_mix9_bk1;
    NTuple::Item<double> m_mix10_bk1;

    NTuple::Array<double> m_gphi1_2g;
    NTuple::Array<double> m_gphi2_2g;

    NTuple::Item<double> m_cos_ll1;
    NTuple::Item<double> m_cos_ll2;
    NTuple::Item<double> m_cos_ll3;
    NTuple::Item<double> m_etac;
    NTuple::Item<double> m_etac_sigma;
    NTuple::Item<double> m_etac_sigma_lamb;
    NTuple::Item<double> m_etac_cos_glamb;
    NTuple::Item<double> m_etac_sigma_bk, m_etac_lamb_bk;
    NTuple::Item<int> m_etac_tag;

    NTuple::Item<double> hel_cos_pp;
    NTuple::Item<double> hel_cos_pim;
    NTuple::Item<double> hel_cos_pm;
    NTuple::Item<double> hel_cos_pip;
    NTuple::Item<double> hel_cos_pppm;
    NTuple::Item<double> hel_cos_lam;
    NTuple::Item<double> hel_cos_lambar;
    NTuple::Item<double> hel_cos_g;

    NTuple::Item<int> r4c_nGoodcharge;
    NTuple::Item<int> r4c_ngam;
    NTuple::Item<int> r4c_charge;
    NTuple::Item<int> m_idxmc;
    NTuple::Array<int> m_pdgid;
    NTuple::Array<int> m_motheridx;
    NTuple::Array<int> m_trkidx;

    NTuple::Item<double> m_chisq;
    //phik+-
    NTuple::Array<long> r4c_phi_kid;
    NTuple::Array<double> r4c_mphi_k;
    NTuple::Array<double> r4c_phi_kpx;
    NTuple::Array<double> r4c_phi_kpy;
    NTuple::Array<double> r4c_phi_kpz;
    NTuple::Array<double> r4c_phi_ke;
    //K+K- form Jpsi
    NTuple::Array<long> r4c_jpsi_kid;
    NTuple::Array<double> r4c_mjpsi_k;
    NTuple::Array<double> r4c_jpsi_kpx;
    NTuple::Array<double> r4c_jpsi_kpy;
    NTuple::Array<double> r4c_jpsi_kpz;
    NTuple::Array<double> r4c_jpsi_ke;
    //Kmiss and K
    NTuple::Array<long> r4c_KKmpid;
    NTuple::Array<double> r4c_mKKmp;
    NTuple::Array<double> r4c_KKmppx;
    NTuple::Array<double> r4c_KKmppy;
    NTuple::Array<double> r4c_KKmppz;
    NTuple::Array<double> r4c_KKmpe;
    //gamma
    NTuple::Array<long> r4c_gid;
    //    NTuple::Array<double> r4c_mg;
    NTuple::Array<double> r4c_gpx;
    NTuple::Array<double> r4c_gpy;
    NTuple::Array<double> r4c_gpz;
    NTuple::Array<double> r4c_ge;
    //2K
    NTuple::Array<long> m_kid;
    NTuple::Array<double> m_mk;
    NTuple::Array<double> m_kpx;
    NTuple::Array<double> m_kpy;
    NTuple::Array<double> m_kpz;
    NTuple::Array<double> m_ke;
    //2gamma--eta
    NTuple::Item<double> m_etaid;
    NTuple::Item<double> m_meta;
    NTuple::Item<double> m_etapx;
    NTuple::Item<double> m_etapy;
    NTuple::Item<double> m_etapz;
    NTuple::Item<double> m_etae;

    NTuple::Item<double> m_m2gam;
    NTuple::Item<double> m_m4k;
    NTuple::Item<double> m_mKKKmp;
    NTuple::Item<double> m_mKKphikp;
    NTuple::Item<double> m_mKKphikm;
    NTuple::Item<double> m_mphikp;
    NTuple::Item<double> m_mphikm;
    NTuple::Item<double> m_mkp0km0;
    NTuple::Item<double> m_mkp1km1;
    NTuple::Item<double> m_mkp0km1;
    NTuple::Item<double> m_mkp1km0;

    //--------------------------------------------------------------------------
    NTuple::Tuple *m_tuple; //mctruth

    ///K+ from J/psi
    NTuple::Item<long> mt_kpid;
    NTuple::Item<double> mt_mkp;
    NTuple::Item<double> mt_kppx;
    NTuple::Item<double> mt_kppy;
    NTuple::Item<double> mt_kppz;
    NTuple::Item<double> mt_kpe;
    ///K- from J/psi
    NTuple::Item<long> mt_kmid;
    NTuple::Item<double> mt_mkm;
    NTuple::Item<double> mt_kmpx;
    NTuple::Item<double> mt_kmpy;
    NTuple::Item<double> mt_kmpz;
    NTuple::Item<double> mt_kme;
    ///K+ from phi
    NTuple::Item<long> mt_phi_kpid;
    NTuple::Item<double> mt_phi_mkp;
    NTuple::Item<double> mt_phi_kppx;
    NTuple::Item<double> mt_phi_kppy;
    NTuple::Item<double> mt_phi_kppz;
    NTuple::Item<double> mt_phi_kpe;
    ///K- from phi
    NTuple::Item<long> mt_phi_kmid;
    NTuple::Item<double> mt_phi_mkm;
    NTuple::Item<double> mt_phi_kmpx;
    NTuple::Item<double> mt_phi_kmpy;
    NTuple::Item<double> mt_phi_kmpz;
    NTuple::Item<double> mt_phi_kme;
    //gamma
    NTuple::Array<double> mt_gid;
    NTuple::Array<double> mt_gpx;
    NTuple::Array<double> mt_gpy;
    NTuple::Array<double> mt_gpz;
    NTuple::Array<double> mt_ge;
    NTuple::Item<double> mt_mphi;
    NTuple::Item<double> mt_m2g;
    NTuple::Array<double> mt_m2k;
    NTuple::Item<double> mt_m4k;
};

#endif