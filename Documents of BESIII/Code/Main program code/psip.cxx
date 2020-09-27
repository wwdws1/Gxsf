#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"

#include "DstEvent/TofHitStatus.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "TMath.h"
using CLHEP::Hep2Vector;
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "ParticleID/ParticleID.h"
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/VertexFit.h"
#include "psipAlg/psip.h"

#include "McTruth/DecayMode.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/McEvent.h"
#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/MucMcHit.h"
#include "McTruth/TofMcHit.h"

#include "CLHEP/Random/RandGauss.h"
#include "RooFit.h"
#include "TRandom.h"

#include <vector>
const double PI = 3.1415927;
const double me = 0.000511;
const double mmu = 0.105658;
const double mpi = 0.139570;
const double mk = 0.493677;
const double mp = 0.938272;

const double velc = 299.792458;            // tof path unit in mm
typedef std::vector<int> Vint;             // 定义int形矢量Vint
typedef std::vector<double> Vdouble;       // 定义double形矢量Vdouble
typedef std::vector<WTrackParameter> WTk;  // 定义径迹参数矢量WTk
typedef std::vector<HepLorentzVector> Vp4; // 定义洛伦兹矢量Vp4

int Ncut0, Ncut1, Ncut2, Ncut3, Ncut4, Ncut5, Ncut6, Ncut7, Ncut8, Ncut9, Ncut10, Ncut11, Ncut12;
int counter[6] = {0, 0, 0, 0, 0, 0};
int evt_mc = 0;
int couter_ph = 0;
int pdgid58 = 0;
int pdgid59 = 0;
int pdgid441 = 0;
int pdgid9080221 = 0;
int pdgid333 = 0;

//////////////////////////////////////////////////////////////////////////////////////////

psip::psip(const std::string &name, ISvcLocator *pSvcLocator) : Algorithm(name, pSvcLocator)
{
    declareProperty("Vr0cut", m_vr0cut = 1.0);                    // 在r方向上截取1cm
    declareProperty("Vz0cut", m_vz0cut = 10.0);                   // 在z方向上截取10cm
    declareProperty("EnergyThreshold", m_energyThreshold = 0.04); // 能量阈值0.04GeV
    declareProperty("GammaPhiCut", m_gammaPhiCut = 20.0);         // 光子phi角截取20度
    declareProperty("GammaThetaCut", m_gammaThetaCut = 20.0);     // 光子theta角截取20度
    declareProperty("GammaAngleCut", m_gammaAngleCut = 10.0);     // 光子angle截取10度
    declareProperty("CheckDedx", m_checkDedx = 1);
    declareProperty("CheckTof", m_checkTof = 1); // Tof
    declareProperty("Test4k4C", m_test4k4C = 1); // 4C
    declareProperty("Test4k5C", m_test4k5C = 1); // 5C
    declareProperty("Test4k5C_4C", m_test4k5C_4C = 1);
    declareProperty("Testmctruth", m_testmctruth = 1);
}

//////////////////////////////////////////////////////////////////////////////////////////

StatusCode psip::initialize() // 初始化
{
    MsgStream log(msgSvc(), name());

    log << MSG::INFO << "in initialize()" << endmsg;
    StatusCode status;
    if (m_testmctruth == 1)
    {
        NTuplePtr nt6(ntupleSvc(), "FILE1/trkRes");
        if (nt6)
            m_tuple6 = nt6;
        else
        {
            m_tuple6 = ntupleSvc()->book("FILE1/trkRes", CLID_ColumnWiseTuple, "ks N-Tuple example");
            if (m_tuple6)
            {
                status = m_tuple6->addItem("indexmc", m_idxmc, 0, 500);
                status = m_tuple6->addIndexedItem("pdgid", m_idxmc, m_pdgid);
                status = m_tuple6->addIndexedItem("motheridx", m_idxmc, m_motheridx);
                status = m_tuple6->addIndexedItem("trkidx", m_idxmc, m_trkidx);

                status = m_tuple6->addIndexedItem("mc_vertex", m_idxmc, vert_mcpar1);

                status = m_tuple6->addItem("trk_px", 2, m_px_trk);
                status = m_tuple6->addItem("trk_py", 2, m_py_trk);
                status = m_tuple6->addItem("trk_pz", 2, m_pz_trk);
                status = m_tuple6->addItem("trk_p3", 2, m_p3_trk);
                status = m_tuple6->addItem("trk_en", 2, m_en_trk);
                status = m_tuple6->addItem("chi2_lamb", chi2_lamb);
                status = m_tuple6->addItem("chi2_lambar", chi2_lambar);
                status = m_tuple6->addItem("mdc_px", 4, m_px_mdc);
                status = m_tuple6->addItem("mdc_py", 4, m_py_mdc);
                status = m_tuple6->addItem("mdc_pz", 4, m_pz_mdc);
                status = m_tuple6->addItem("mdc_p3", 4, m_p3_mdc);
                status = m_tuple6->addItem("mdc_en", 4, m_en_mdc);
                status = m_tuple6->addItem("emc_px", 2, m_px_emc);
                status = m_tuple6->addItem("emc_py", 2, m_py_emc);
                status = m_tuple6->addItem("emc_pz", 2, m_pz_emc);
                status = m_tuple6->addItem("emc_p3", 2, m_p3_emc);
                status = m_tuple6->addItem("emc_en", 2, m_en_emc);
                status = m_tuple6->addItem("pull_drho", 4, m_pull_1);
                status = m_tuple6->addItem("pull_phi", 4, m_pull_2);
                status = m_tuple6->addItem("pull_kappa", 4, m_pull_3);
                status = m_tuple6->addItem("pull_dz", 4, m_pull_4);
                status = m_tuple6->addItem("pull_tanlamda", 4, m_pull_5);

                status = m_tuple6->addItem("trk_nCh", chargTrk_index, 0, 100);
                status = m_tuple6->addIndexedItem("trk_Rxy", chargTrk_index, m_vtx_Rxy);
                status = m_tuple6->addIndexedItem("trk_cos", chargTrk_index, m_vtx_cos);
                status = m_tuple6->addIndexedItem("trk_z0", chargTrk_index, m_vtx_z0);
                status = m_tuple6->addIndexedItem("trk_p0", chargTrk_index, m_vtx_p);
                status = m_tuple6->addItem("trk_ach", m_totE_ch);
                status = m_tuple6->addItem("trk_npp", m_npp);
                status = m_tuple6->addItem("trk_npm", m_npm);
                status = m_tuple6->addItem("trk_npip", m_npip);
                status = m_tuple6->addItem("trk_npim", m_npim);

                status = m_tuple6->addItem("gam_nGam", gamma_index, 0, 100);
                status = m_tuple6->addIndexedItem("gam_dthe", gamma_index, m_dthe);
                status = m_tuple6->addIndexedItem("gam_dphi", gamma_index, m_dphi);
                status = m_tuple6->addIndexedItem("gam_dang", gamma_index, m_dang);
                status = m_tuple6->addIndexedItem("gam_eraw", gamma_index, m_eraw);
                status = m_tuple6->addIndexedItem("gam_TDC", gamma_index, m_TDC);
                status = m_tuple6->addIndexedItem("gam_tagg", gamma_index, m_tagg);

                status = m_tuple6->addItem("orig_m2pi", m_phi1_bk);
                status = m_tuple6->addItem("orig_m2p", m_phi2_bk);
                status = m_tuple6->addItem("orig_mpppim", m_phi1_bk1);
                status = m_tuple6->addItem("orig_mpmpip", m_phi2_bk1);
                status = m_tuple6->addItem("orig_m2gamma", m_2gam_bk1);
                status = m_tuple6->addItem("orig_mlambda", m_lam0_bk1);
                status = m_tuple6->addItem("orig_mlambdabar", m_lam1_bk1);

                status = m_tuple6->addItem("orig_m2gammalambda", m_mix0_bk1);
                status = m_tuple6->addItem("orig_m2gammalambdabar", m_mix1_bk1);
                status = m_tuple6->addItem("orig_mlambdalambdabar", m_mix2_bk1);
                status = m_tuple6->addItem("orig_mgamma1lambda", m_mix3_bk1);
                status = m_tuple6->addItem("orig_mgamma1lambdabar", m_mix4_bk1);
                status = m_tuple6->addItem("orig_mgamma2lambda", m_mix5_bk1);
                status = m_tuple6->addItem("orig_mgamma2lambdabar", m_mix6_bk1);
                status = m_tuple6->addItem("orig_mpippim2gamma", m_mix7_bk1);
                status = m_tuple6->addItem("orig_recoilmpppm", m_mix8_bk1);
                status = m_tuple6->addItem("orig_mgamma1lambdalambdabar", m_mix9_bk1);
                status = m_tuple6->addItem("orig_mgamma2lambdalambdabar", m_mix10_bk1);

                status = m_tuple6->addItem("orig_mtotal", m_totm_bk);
                status = m_tuple6->addItem("orig_m4trk", m_2phi_bk);
                status = m_tuple6->addItem("ki_mpppim", m_phi1_true);
                status = m_tuple6->addItem("ki_mpmpip", m_phi2_true);
                status = m_tuple6->addItem("ki_metac", m_etac_true);
                status = m_tuple6->addItem("ki_dthe", dthe_m_4c);
                status = m_tuple6->addItem("ki_dang", dang_m_4c);
                status = m_tuple6->addItem("ki_dphi", dphi_m_4c);
                status = m_tuple6->addItem("ki_eraw", eraw_m_4c);
                status = m_tuple6->addItem("ki_tdc", TDC_m_4c);
                status = m_tuple6->addItem("ki_tagg", tagg_m_4c);

                status = m_tuple6->addItem("lambda_n", m_ngLam1);
                status = m_tuple6->addItem("lambdab_n", m_ngLam2);
                status = m_tuple6->addItem("lambda_chisq", m_secchis1);
                status = m_tuple6->addItem("lambda_lth", m_mlengthlam);
                status = m_tuple6->addItem("lambda_lth_e", m_mlengthlam_err);
                status = m_tuple6->addItem("lambda_lth_re", m_mlengthlam_rerr);
                status = m_tuple6->addItem("lambdabar_chisq", m_secchis2);
                status = m_tuple6->addItem("lambdabar_lth", m_mlengthlamb);
                status = m_tuple6->addItem("lambdabar_lth_e", m_mlengthlamb_err);
                status = m_tuple6->addItem("lambdabar_lth_re", m_mlengthlamb_rerr);

                status = m_tuple6->addItem("2vtxfit_minchisq_mlambda", m_secchis12diff1);
                status = m_tuple6->addItem("2vtxfit_minchisq_mlambdabar", m_secchis12diff2);
                status = m_tuple6->addItem("2vtxfit_pass4c_elambdalambdabar", m_secchis12diff3);
                status = m_tuple6->addItem("2vtxfit_nopass4c_elambdalambdabar", m_secchis12diff4);
                status = m_tuple6->addItem("2vtxfit_elambdalambdabar_diff", m_secchis12diff5);

                status = m_tuple6->addItem("4C_chisq", m_chi1);
                status = m_tuple6->addItem("ki0g_chi2", m_chi1_0g);

            } //if m_tuple6
            else
            {
                log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple6) << endmsg;
                return StatusCode::FAILURE;
            }
        }
        // for mctruth

        NTuplePtr nt11(ntupleSvc(), "FILE1/mctruth");
        if (nt11)
            m_tuple11 = nt11;
        else
        {
            m_tuple11 = ntupleSvc()->book("FILE1/mctruth", CLID_ColumnWiseTuple, "NTuple");
            if (m_tuple11)
            {
                status = m_tuple11->addItem("mc_nlambda", m_mc_nlam);
                status = m_tuple11->addItem("mc_lambda_px", parti_px_mc1);
                status = m_tuple11->addItem("mc_lambda_py", parti_py_mc1);
                status = m_tuple11->addItem("mc_lambda_pz", parti_pz_mc1);
                status = m_tuple11->addItem("mc_lambda_e", parti_e_mc1);
                status = m_tuple11->addItem("mc_lambda_m", parti_m_mc1);

                status = m_tuple11->addItem("mc_nlambdabar", m_mc_nlambar);
                status = m_tuple11->addItem("mc_lambdabar_px", parti_px_mc2);
                status = m_tuple11->addItem("mc_lambdabar_py", parti_py_mc2);
                status = m_tuple11->addItem("mc_lambdabar_pz", parti_pz_mc2);
                status = m_tuple11->addItem("mc_lambdabar_e", parti_e_mc2);
                status = m_tuple11->addItem("mc_lambdabar_m", parti_m_mc2);

                status = m_tuple11->addItem("mc_npi0", m_mc_npi0);
                status = m_tuple11->addItem("mc_pi0_px", parti_px_mc3);
                status = m_tuple11->addItem("mc_pi0_py", parti_py_mc3);
                status = m_tuple11->addItem("mc_pi0_pz", parti_pz_mc3);
                status = m_tuple11->addItem("mc_pi0_e", parti_e_mc3);
                status = m_tuple11->addItem("mc_pi0_m", parti_m_mc3);

                status = m_tuple11->addItem("mc_neta", m_mc_neta);
                status = m_tuple11->addItem("mc_eta_px", parti_px_mc4);
                status = m_tuple11->addItem("mc_eta_py", parti_py_mc4);
                status = m_tuple11->addItem("mc_eta_pz", parti_pz_mc4);
                status = m_tuple11->addItem("mc_eta_e", parti_e_mc4);
                status = m_tuple11->addItem("mc_eta_m", parti_m_mc4);

                status = m_tuple11->addItem("mc_npp", m_mc_npp);
                status = m_tuple11->addItem("mc_pp_px", parti_px_mc11);
                status = m_tuple11->addItem("mc_pp_py", parti_py_mc11);
                status = m_tuple11->addItem("mc_pp_pz", parti_pz_mc11);
                status = m_tuple11->addItem("mc_pp_e", parti_e_mc11);
                status = m_tuple11->addItem("mc_pp_m", parti_m_mc11);

                status = m_tuple11->addItem("mc_npim", m_mc_npim);
                status = m_tuple11->addItem("mc_pim_px", parti_px_mc12);
                status = m_tuple11->addItem("mc_pim_py", parti_py_mc12);
                status = m_tuple11->addItem("mc_pim_pz", parti_pz_mc12);
                status = m_tuple11->addItem("mc_pim_e", parti_e_mc12);
                status = m_tuple11->addItem("mc_pim_m", parti_m_mc12);

                status = m_tuple11->addItem("mc_npm", m_mc_npm);
                status = m_tuple11->addItem("mc_pm_px", parti_px_mc21);
                status = m_tuple11->addItem("mc_pm_py", parti_py_mc21);
                status = m_tuple11->addItem("mc_pm_pz", parti_pz_mc21);
                status = m_tuple11->addItem("mc_pm_e", parti_e_mc21);
                status = m_tuple11->addItem("mc_pm_m", parti_m_mc21);

                status = m_tuple11->addItem("mc_npip", m_mc_npip);
                status = m_tuple11->addItem("mc_pip_px", parti_px_mc22);
                status = m_tuple11->addItem("mc_pip_py", parti_py_mc22);
                status = m_tuple11->addItem("mc_pip_pz", parti_pz_mc22);
                status = m_tuple11->addItem("mc_pip_e", parti_e_mc22);
                status = m_tuple11->addItem("mc_pip_m", parti_m_mc22);

                status = m_tuple11->addItem("mc_ngamma", m_mc_ngamma);
                status = m_tuple11->addItem("mc_gamma_px", parti_px_mcgamma);
                status = m_tuple11->addItem("mc_gamma_py", parti_py_mcgamma);
                status = m_tuple11->addItem("mc_gamma_pz", parti_pz_mcgamma);
                status = m_tuple11->addItem("mc_gamma_e", parti_e_mcgamma);
                status = m_tuple11->addItem("mc_gamma_m", parti_m_mcgamma);
            }
            else
            {
                log << MSG::ERROR << "Cannot book N-tuple:"
                    << long(m_tuple11) << endmsg;
                return StatusCode::FAILURE;
            }
        }
    }
    return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////////////

StatusCode psip::execute() // 执行
{
    psip::InitVar();

    //cout << __LINE__<<endl;

    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in execute()" << endreq;

    HepLorentzVector ecms(3.686 * sin(0.011), 0, 0, 3.686); // 质心系能量3.686

    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    int runNum = eventHeader->runNumber();
    int eventNum = eventHeader->eventNumber();

    Ncut0++; // 如果runNum > 0 ，事例累加（总的事例数）

    //cout << "******************************" << Ncut0 << "******************************" << endl;

    //cout << __LINE__ << endl;

    SmartDataPtr<EvtRecEvent>
        evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    log << MSG::DEBUG << "ncharg, nneu, tottks = " // 带电径迹数，中性径迹数，总径迹数
        << evtRecEvent->totalCharged() << " , "
        << evtRecEvent->totalNeutral() << " , "
        << evtRecEvent->totalTracks() << endreq;
    if (evtRecEvent->totalTracks() > 99)
    {
        return StatusCode::SUCCESS; // 如果带电径迹数目大于99，返回SUCCESS。(总径迹数过多，即本底过多，舍弃)
    }

    //////////////////////////////////////////////////////////////////////////////////

    // 下面一段代码用于使用Inclusive MC的数据进行验证

    if (runNum < 0)
    {
        int numParticle = 0;
        int mcnlam = 0;
        int mcnlambar = 0;
        int mcnpi0 = 0;
        int mcneta = 0;
        int mcnpp = 0;
        int mcnpm = 0;
        int mcnpip = 0;
        int mcnpim = 0;
        int mcng = 0;

        bool Decay = false;
        int rootIndex = -1;
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        if (!mcParticleCol)
        {
            std::cout << "Could not retrieve McParticelCol" << std::endl;
            return StatusCode::FAILURE;
        }
        else
        {
            Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            bool strange = false;
            for (; iter_mc != mcParticleCol->end(); iter_mc++)
            {
                int temp_id = (*iter_mc)->particleProperty();
                int motherpdg = ((*iter_mc)->mother()).particleProperty();

                //HepLorentzVector mctrue_track = (*iter_mc)->initialFourMomentum();    // 原始四动量
                //HepLorentzVector mctrack_iniposition = (*iter_mc)->initialPosition(); // 原始位置
                //HepLorentzVector mctrack_finposition = (*iter_mc)->finalPosition();   // 最后位置

                if ((*iter_mc)->primaryParticle() && temp_id == 11 && motherpdg == 11)
                {
                    strange = true;
                }
                if ((*iter_mc)->primaryParticle())
                    continue;
                if (!(*iter_mc)->decayFromGenerator())
                    continue;

                if (temp_id == 100443)
                {
                    rootIndex = (*iter_mc)->trackIndex();
                    Decay = true;
                }
                if (!Decay)
                    continue;
                int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
                int trkidx = (*iter_mc)->trackIndex() - rootIndex;
                if (strange && motherpdg != 100443)
                {
                    mcidx--;
                }
                else if (temp_id != 100443)
                {
                    trkidx--;
                }

                m_pdgid[numParticle] = temp_id;
                m_motheridx[numParticle] = mcidx;
                m_trkidx[numParticle] = trkidx;
                vert_mcpar1[numParticle] = (*iter_mc)->vertexIndex0();
                numParticle++;

                //cout << "trkidx,temp_id,mcidx = " << setw(8) << trkidx << setw(8) << temp_id << setw(8) << mcidx << endl;

                //cout << " mother " << motherpdg << " itself " << temp_id << endl;

                double mcpx = (*iter_mc)->initialFourMomentum().px(); //最初px
                double mcpy = (*iter_mc)->initialFourMomentum().py(); //最初py
                double mcpz = (*iter_mc)->initialFourMomentum().pz(); //最初pz
                double mcen = (*iter_mc)->initialFourMomentum().e();  //最初能量
                HepLorentzVector p4mc;
                p4mc.setPx(mcpx);
                p4mc.setPy(mcpy);
                p4mc.setPz(mcpz);
                p4mc.setE(mcen);
                double mmcpar = p4mc.mag();

                //cout << " mass of the particle " << p4mc.mag() << endl;

                if (motherpdg == 100443 && temp_id == 3122)
                {
                    parti_px_mc1 = mcpx;
                    parti_py_mc1 = mcpy;
                    parti_pz_mc1 = mcpz;
                    parti_e_mc1 = mcen;
                    parti_m_mc1 = mmcpar;
                    mcnlam++;
                } // lambda数量、px、py、pz、能量
                else if (motherpdg == 100443 && temp_id == -3122)
                {
                    parti_px_mc2 = mcpx;
                    parti_py_mc2 = mcpy;
                    parti_pz_mc2 = mcpz;
                    parti_e_mc2 = mcen;
                    parti_m_mc2 = mmcpar;
                    mcnlambar++;
                } // lambdabar数量、px、py、pz、能量
                else if (motherpdg == 100443 && temp_id == 111)
                {
                    parti_px_mc3 = mcpx;
                    parti_py_mc3 = mcpy;
                    parti_pz_mc3 = mcpz;
                    parti_e_mc3 = mcen;
                    parti_m_mc3 = mmcpar;
                    mcnpi0++;
                } // pi0数量、px、py、pz、能量
                else if (motherpdg == 100443 && temp_id == 221)
                {
                    parti_px_mc4 = mcpx;
                    parti_py_mc4 = mcpy;
                    parti_pz_mc4 = mcpz;
                    parti_e_mc4 = mcen;
                    parti_m_mc4 = mmcpar;
                    mcneta++;
                } // eta数量、px、py、pz、能量
                else if (motherpdg == 3122 && temp_id == 2212)
                {
                    parti_px_mc11 = mcpx;
                    parti_py_mc11 = mcpy;
                    parti_pz_mc11 = mcpz;
                    parti_e_mc11 = mcen;
                    parti_m_mc11 = mmcpar;
                    mcnpp++;
                } // pp数量、px、py、pz、能量
                else if (motherpdg == 3122 && temp_id == -211)
                {
                    parti_px_mc12 = mcpx;
                    parti_py_mc12 = mcpy;
                    parti_pz_mc12 = mcpz;
                    parti_e_mc12 = mcen;
                    parti_m_mc12 = mmcpar;
                    mcnpim++;
                } // pim数量、px、py、pz、能量
                else if (motherpdg == -3122 && temp_id == -2212)
                {
                    parti_px_mc21 = mcpx;
                    parti_py_mc21 = mcpy;
                    parti_pz_mc21 = mcpz;
                    parti_e_mc21 = mcen;
                    parti_m_mc21 = mmcpar;
                    mcnpm++;
                } // pm数量、px、py、pz、能量
                else if (motherpdg == -3122 && temp_id == 211)
                {
                    parti_px_mc22 = mcpx;
                    parti_py_mc22 = mcpy;
                    parti_pz_mc22 = mcpz;
                    parti_e_mc22 = mcen;
                    parti_m_mc22 = mmcpar;
                    mcnpip++;
                } // pip数量、px、py、pz、能量
                else if (motherpdg == 111 && temp_id == 22)
                {
                    parti_px_mcgamma = mcpx;
                    parti_py_mcgamma = mcpy;
                    parti_pz_mcgamma = mcpz;
                    parti_e_mcgamma = mcen;
                    parti_m_mcgamma = mmcpar;
                    mcng++;
                } // 来自pi0的gamma数量、px、py、pz、能量(不区分来源)
                else if (motherpdg == 221 && temp_id == 22)
                {
                    parti_px_mcgamma = mcpx;
                    parti_py_mcgamma = mcpy;
                    parti_pz_mcgamma = mcpz;
                    parti_e_mcgamma = mcen;
                    parti_m_mcgamma = mmcpar;
                    mcng++;
                } // 来自eta的gamma数量、px、py、pz、能量(不区分来源)
                else
                {
                    continue;
                }
            }
        }

        m_idxmc = numParticle;    // 总粒子数
        m_mc_nlam = mcnlam;       // lamda数
        m_mc_nlambar = mcnlambar; // lamdabar数
        m_mc_npi0 = mcnpi0;       // pi0数
        m_mc_neta = mcneta;       // eta数
        m_mc_npp = mcnpp;         // pp数
        m_mc_npm = mcnpm;         // pm数
        m_mc_npip = mcnpip;       // pip数
        m_mc_npim = mcnpim;       // pim数
        m_mc_ngamma = mcng;       // gamma数
        m_tuple11->write();
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);

    // for Total charged
    // 定义iGood和iCharge(好的带电事例数和总电荷数)并清空数据
    Vint iGood, iCharge, iGoodP, iGoodM;
    iGood.clear();
    iCharge.clear();
    iGoodP.clear();
    iGoodM.clear();

    int nCharge = 0;

    Hep3Vector xorigin(0, 0, 0);
    Hep3Vector eorigin(0, 0, 0); // 定义原点
    IVertexDbSvc *vtxsvc;
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if (vtxsvc->isVertexValid())
    {
        double *dbv = vtxsvc->PrimaryVertex();
        double *vv = vtxsvc->SigmaPrimaryVertex();
        xorigin.setX(dbv[0]);
        xorigin.setY(dbv[1]);
        xorigin.setZ(dbv[2]);
        eorigin.setX(vv[0]);
        eorigin.setY(vv[1]);
        eorigin.setZ(vv[2]);
    } // 将X,Y,Z保存在数组中

    //cout << " charged track reconstruction " << endl;

    int iv = 0;
    double etot = 0.0; // E Total
    for (int i = 0; i < evtRecEvent->totalCharged(); i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!((*itTrk)->isMdcKalTrackValid()))
            continue;                                        // 不是MDC的轨迹不要
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack(); // 重建MDC轨迹

        HepVector a = mdcKalTrk->helix();
        HepSymMatrix Ea = mdcKalTrk->err();
        HepPoint3D point0(0., 0., 0.);
        HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]); // 定义IP值
        VFHelix helixip(point0, a, Ea);                    // 螺旋
        helixip.pivot(IP);
        HepVector vecipa = helixip.a();
        double Rvxy0 = fabs(vecipa[0]); // the distance to IP in xy plane
        double Rvz0 = vecipa[3];        // the distance to IP in z direction
        double Rvphi0 = vecipa[1];
        double cos_1 = cos(mdcKalTrk->theta());
        double pch = mdcKalTrk->p();
        if (fabs(cos_1) > 0.93)
            continue; // 如果cos theta绝对值的值大于0.93就不要（不得超过MDC测量范围）
        if (pch > 2.0)
            continue; // 如果粒子动量大于2.0也不要（总的动量是3.686，一半不可能大于2）

        m_vtx_cos[iv] = cos_1;
        m_vtx_Rxy[iv] = Rvxy0;
        m_vtx_z0[iv] = Rvz0;
        m_vtx_p[iv] = pch;
        iv++;

        iGood.push_back(i);
        iCharge.push_back(mdcKalTrk->charge());
        nCharge += mdcKalTrk->charge(); // 总电量
        if (mdcKalTrk->charge() > 0)
        {
            iGoodP.push_back(i);
        }
        else
        {
            iGoodM.push_back(i);
        }
    }

    //m_totE_ch = etot; // 总能量

    int nGood = iGood.size();
    chargTrk_index = nGood;
    log << MSG::DEBUG << "ngood = " << nGood << " , totcharge = " << nCharge << endreq;

    if (nGood != 4 || nCharge != 0)
    {
        return StatusCode::SUCCESS;
    }
    Ncut2++; // 如果总带电径迹少于4个，不要

    //
    // Finish Good Charged Track Selection
    // 完成好带电径迹筛选
    //

    Vint iGam;
    iGam.clear();
    int ii = 0;

    HepLorentzVector emctrk;

    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) // 挑选出EMC的事例
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) // 非EMC事例不要
            continue;

        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        int tsttag = 0;
        for (int j = 0; j < nGood; j++)
        {
            EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + iGood[j];
            if (!(*jtTrk)->isExtTrackValid())
                continue;
            RecExtTrack *extTrk = (*jtTrk)->extTrack();
            if (extTrk->emcVolumeNumber() == -1)
                continue;
            Hep3Vector extpos = extTrk->emcPosition();

            double angd = extpos.angle(emcpos);
            double thed = extpos.theta() - emcpos.theta();
            double phid = extpos.deltaPhi(emcpos); // 定义angd,thed和phid角

            if (angd < dang) // 如果angd小于200
            {
                dang = angd;
                dthe = thed;
                dphi = phid;
                tsttag = iGood[j];
            }
        }
        dthe = dthe * 180 / (CLHEP::pi);
        dphi = dphi * 180 / (CLHEP::pi);
        dang = dang * 180 / (CLHEP::pi); // 将弧度转化为角度

        //if (dang < 10) continue;

        double eraw = emcTrk->energy(); // 能量
        double tdc = emcTrk->time();    // 时间

        //if (eraw < m_energyThreshold) continue;

        // iGam位置筛选

        double costh = cos(emcTrk->theta());
        bool ok1 = 1;
        bool ok2 = 1;
        if ((fabs(costh) > 0.8 && fabs(costh) < 0.86) || fabs(costh) > 0.92)
            continue;          // 不得超出电磁量能器的测量范围
        if (fabs(costh) < 0.8) // 桶部测量范围
        {
            if (eraw < 0.025) // 桶部能量最小分辨率
                ok1 = 0;
        }
        if (fabs(costh) > 0.86 && fabs(costh) < 0.92) // 端盖测量范围
        {
            if (eraw < 0.050) // 端盖能量最小分辨率
                ok2 = 0;
        }

        if (!ok1)
            continue;
        if (!ok2)
            continue;
        // 超过测量范围的都不要

        if (fabs(dang) < m_gammaAngleCut) // dang角的绝对值小于10不要
            continue;
        if (tdc < 0 || tdc > 15) // 飞行时间要在0到15ns之间
            continue;
        m_dthe[ii] = dthe;
        m_dphi[ii] = dphi;
        m_dang[ii] = dang;
        m_eraw[ii] = eraw;
        m_TDC[ii] = emcTrk->time(); // 记下dthe、dphi、dang角、能量和时间

        ii++;
        iGam.push_back(i);
    }

    //
    // Finish Good Photon Selection
    // 完成好光子的筛选
    //

    int nGam = iGam.size();
    gamma_index = nGam;

    log << MSG::DEBUG << "num Good Photon = " << nGam << " , "
        << evtRecEvent->totalNeutral() << endreq;

    if (nGam < 2) // 如果光子数少于2就不要
        return StatusCode::SUCCESS;
    Ncut3++;

    //
    // Give each photon four momentum
    // 为每个光子赋予四动量
    //

    Vp4 pGam;
    pGam.clear();

    for (int i = 0; i < nGam; i++)
    {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
        RecEmcShower *emcTrk = (*itTrk)->emcShower(); // EMC束流
        double eraw = emcTrk->energy();               // 能量
        double phi = emcTrk->phi();                   // phi角
        double the = emcTrk->theta();                 // theta角
        HepLorentzVector ptrk;
        ptrk.setPx(eraw * sin(the) * cos(phi));
        ptrk.setPy(eraw * sin(the) * sin(phi));
        ptrk.setPz(eraw * cos(the));
        ptrk.setE(eraw);
        pGam.push_back(ptrk); // 根据能量和角度可以得到光子的四动量
    }

    HepLorentzVector p4c[4], ptk[4], p_pmpip[2], p_pppim[2]; //定义洛伦兹矢量
    // p4c顺序：1.Lambda 2.Lambdabar 3.gamma1 4.gamma2
    // ptk(未过4C)顺序：1.pp 2.pm 3.pip 4.pim

    //
    // begin 4C fit
    // 开始4C拟合
    //

    if (m_test4k4C == 1)
    {

        HepPoint3D vx(0., 0., 0.); // ini vertex
        HepSymMatrix Evx(3, 0);
        double bx = 1E+6;
        double by = 1E+6;
        double bz = 1E+6;
        Evx[0][0] = bx * bx;
        Evx[1][1] = by * by;
        Evx[2][2] = bz * bz;
        VertexParameter vxpar;
        vxpar.setVx(vx);
        vxpar.setEvx(Evx);

        HepPoint3D pvx = xorigin; // for every run
        HepSymMatrix pEvx(3, 0);
        pEvx[0][0] = eorigin[0] * eorigin[0];
        pEvx[1][1] = eorigin[1] * eorigin[1];
        pEvx[2][2] = eorigin[2] * eorigin[2];
        VertexParameter pvxpar;
        pvxpar.setVx(pvx);
        pvxpar.setEvx(pEvx);

        // Test vertex fit

        VertexFit *vtxfit = VertexFit::instance();              // 顶点拟合
        SecondVertexFit *svtxfit = SecondVertexFit::instance(); // 次级顶点拟合

        VertexFit *vtxfit1 = VertexFit::instance();
        SecondVertexFit *svtxfit1 = SecondVertexFit::instance();

        Vdouble temp_chi2VT1, temp_chi2VT2, temp_chi2VT12;
        temp_chi2VT1.clear();
        temp_chi2VT2.clear();
        temp_chi2VT12.clear();

        WTk temp_wtkLam1, temp_wtkLam2; // lambda & lambdabar 径迹矢量
        temp_wtkLam1.clear();
        temp_wtkLam2.clear();

        HepLorentzVector p4_pp, p4_pm, p4_pip, p4_pim; // pp pm pip pim 拟合前四动量

        Vp4 temp_p4_pp, temp_p4_pm, temp_p4_pip, temp_p4_pim; // pp pm pip pim 拟合后四动量
        temp_p4_pp.clear();
        temp_p4_pm.clear();
        temp_p4_pip.clear();
        temp_p4_pim.clear();

        double p3;
        int count = 0;

        WTrackParameter wvppTrk, wvpmTrk, wvpipTrk, wvpimTrk; // pp pm pip pim 径迹参数

        Vdouble temp_length_lam, temp_length_lam_err; // lambda 衰变长度&误差
        temp_length_lam.clear();
        temp_length_lam_err.clear();

        Vdouble temp_length_lam_bar, temp_length_lam_bar_err; // lambdabar 衰变长度&误差
        temp_length_lam_bar.clear();
        temp_length_lam_bar_err.clear();
        for (int n0 = 0; n0 < 2; n0++)
        {
            RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
            RecMdcKalTrack *ppTrk = (*(evtRecTrkCol->begin() + iGoodP[n0]))->mdcKalTrack();
            p4_pp.setPx(ppTrk->px());
            p4_pp.setPy(ppTrk->py());
            p4_pp.setPz(ppTrk->pz());
            p3 = ppTrk->p();
            p4_pp.setE(sqrt(p3 * p3 + mp * mp));
            wvppTrk = WTrackParameter(mp, ppTrk->getZHelixP(), ppTrk->getZErrorP());
            for (int n1 = 0; n1 < 2; n1++)
            {
                RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
                RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + iGoodM[n1]))->mdcKalTrack();
                p4_pim.setPx(pimTrk->px());
                p4_pim.setPy(pimTrk->py());
                p4_pim.setPz(pimTrk->pz());
                p3 = pimTrk->p();
                p4_pim.setE(sqrt(p3 * p3 + mpi * mpi));
                wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());
                for (int n2 = 0; n2 < 2; n2++)
                {
                    if (n2 == n0)
                    {
                        continue;
                    }
                    RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
                    RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + iGoodP[n2]))->mdcKalTrack();
                    p4_pip.setPx(pipTrk->px());
                    p4_pip.setPy(pipTrk->py());
                    p4_pip.setPz(pipTrk->pz());
                    p3 = pipTrk->p();
                    p4_pip.setE(sqrt(p3 * p3 + mpi * mpi));
                    wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
                    for (int n3 = 0; n3 < 2; n3++)
                    {
                        if (n3 == n1)
                        {
                            continue;
                        }
                        RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
                        RecMdcKalTrack *pmTrk = (*(evtRecTrkCol->begin() + iGoodM[n3]))->mdcKalTrack();
                        p4_pm.setPx(pmTrk->px());
                        p4_pm.setPy(pmTrk->py());
                        p4_pm.setPz(pmTrk->pz());
                        p3 = pmTrk->p();
                        p4_pm.setE(sqrt(p3 * p3 + mp * mp));
                        wvpmTrk = WTrackParameter(mp, pmTrk->getZHelixP(), pmTrk->getZErrorP());

                        //////////////////////////////////////////////////
                        // Lambda vertex fit and secondary vertex fit
                        //////////////////////////////////////////////////

                        vtxfit->init();
                        vtxfit->AddTrack(0, wvppTrk);      // 添加顶点拟合参数-径迹：pp径迹
                        vtxfit->AddTrack(1, wvpimTrk);     // 添加顶点拟合参数-径迹：pim径迹
                        vtxfit->AddVertex(0, vxpar, 0, 1); // 添加顶点拟合参数-顶点位置：初级顶点
                        if (!vtxfit->Fit(0))
                        {
                            continue;
                        }
                        if (vtxfit->chisq(0) < 0)
                        {
                            continue;
                        }

                        vtxfit->Swim(0);
                        vtxfit->BuildVirtualParticle(0);
                        WTrackParameter wALambda = vtxfit->wVirtualTrack(0); // 返回lambda径迹参数
                        VertexParameter vtxlambda = vtxfit->vpar(0);         // 返回lambda初级顶点参数

                        // secondary vertex fit
                        // 次级顶点拟合

                        svtxfit->init();
                        svtxfit->AddTrack(0, wALambda);    // 添加顶点拟合参数-径迹：lambda径迹
                        svtxfit->setVpar(vtxlambda);       // 设置顶点拟合参数-顶点参数：lambda初级顶点参数
                        svtxfit->setPrimaryVertex(pvxpar); // 设置顶点拟合参数-初级顶点：lambda初级顶点

                        if (!svtxfit->Fit())
                        {
                            continue;
                        }
                        if (svtxfit->chisq() < 0)
                        {
                            continue;
                        }

                        //////////////////////////////////////////////////
                        // Lambdabar vertex fit and secondary vertex fit
                        //////////////////////////////////////////////////

                        vtxfit1->init();
                        vtxfit1->AddTrack(0, wvpmTrk);      // 添加顶点拟合参数-径迹：pm径迹
                        vtxfit1->AddTrack(1, wvpipTrk);     // 添加顶点拟合参数-径迹：pip径迹
                        vtxfit1->AddVertex(0, vxpar, 0, 1); // 添加顶点拟合参数-顶点位置：初级顶点
                        if (!vtxfit1->Fit(0))
                        {
                            continue;
                        }
                        if (vtxfit1->chisq() < 0)
                        {
                            continue;
                        }

                        vtxfit1->Swim(0);
                        vtxfit1->BuildVirtualParticle(0);
                        WTrackParameter wALambda_bar = vtxfit1->wVirtualTrack(0); // 返回lambdabar径迹参数
                        VertexParameter vtxlambda_bar = vtxfit1->vpar(0);         // 返回lambdabar初级顶点参数

                        // secondary vertex fit
                        // 次级顶点拟合

                        svtxfit1->init();
                        svtxfit1->AddTrack(0, wALambda_bar); // 添加顶点拟合参数-径迹：lambdabar径迹
                        svtxfit1->setVpar(vtxlambda_bar);    // 设置顶点拟合参数-顶点参数：lambdabar初级顶点参数
                        svtxfit1->setPrimaryVertex(pvxpar);  // 设置顶点拟合参数-初级顶点：lambdabar初级顶点

                        if (!svtxfit1->Fit())
                        {
                            continue;
                        }
                        if (svtxfit1->chisq() < 0)
                        {
                            continue;
                        }

                        //////////////////////////////////////////////////
                        // Save vertex fit and secondary vertex fit data
                        //////////////////////////////////////////////////

                        temp_chi2VT1.push_back(svtxfit->chisq());
                        temp_chi2VT2.push_back(svtxfit1->chisq());
                        temp_chi2VT12.push_back(svtxfit->chisq() + svtxfit1->chisq());

                        temp_wtkLam1.push_back(wALambda);                           // lambda径迹
                        temp_p4_pp.push_back(p4_pp);                                // pp动量
                        temp_p4_pim.push_back(p4_pim);                              // pim动量
                        temp_length_lam.push_back(svtxfit->decayLength());          // lambda衰变长度
                        temp_length_lam_err.push_back(svtxfit->decayLengthError()); // lambda衰变长度误差

                        temp_wtkLam2.push_back(wALambda_bar);                            // lambdabar径迹
                        temp_p4_pm.push_back(p4_pm);                                     // pm四动量
                        temp_p4_pip.push_back(p4_pip);                                   // pip四动量
                        temp_length_lam_bar.push_back(svtxfit1->decayLength());          // lambdabar衰变长度
                        temp_length_lam_bar_err.push_back(svtxfit1->decayLengthError()); // lambdabar衰变长度误差

                        count++;
                    }
                }
            }
        }

        int ncddLam1 = temp_wtkLam1.size(); // lambda数量
        int ncddLam2 = temp_wtkLam2.size(); // lambdabar数量

        double chi2VT12 = 1999.9;
        double chi2VT12min = 1999.9;
        bool vtxfit12fail = true;

        if (count < 1)
        {
            return StatusCode::SUCCESS;
        }

        for (int i = 0; i < count; i++)
        {
            if (temp_chi2VT12[i] < chi2VT12)
            {
                vtxfit12fail = false;
            }
        }

        if (vtxfit12fail)
        {
            return StatusCode::SUCCESS;
        }

        Ncut4++;

        double t_chi2VT1, t_chi2VT2;

        for (int i = 0; i < count; i++)
        {
            if (temp_chi2VT12[i] < chi2VT12min)
            {
                chi2VT12min = temp_chi2VT12[i];
                t_chi2VT1 = i;
                t_chi2VT2 = i;
            }
        }

        HepLorentzVector p4_lam, p4_lamb, p4_gam1, p4_gam2, p4_emc1, p4_emc2;

        Vp4 temp_p4_lam, temp_p4_lamb, temp_p4_gam1, temp_p4_gam2, temp_p4_emc1, temp_p4_emc2;
        temp_p4_lam.clear();
        temp_p4_lamb.clear();
        temp_p4_gam1.clear();
        temp_p4_gam2.clear();

        int kmfit_gn1, kmfit_gn2;
        Vint temp_count;
        temp_count.clear();

        Vdouble temp_chi2VT3;
        temp_chi2VT3.clear();

        int count1 = 0;
        double chi2VT3 = 999.9;
        bool kmfit_fail = true;

        KalmanKinematicFit *kmfit_2g = KalmanKinematicFit::instance();
        for (int i = 0; i < count; i++)
        {
            chi2VT3 = 999.9;
            kmfit_fail = true;
            for (int gn1 = 0; gn1 < nGam - 1; gn1++)
            {
                RecEmcShower *g1trk = (*(evtRecTrkCol->begin() + iGam[gn1]))->emcShower();
                for (int gn2 = gn1 + 1; gn2 < nGam; gn2++)
                {
                    RecEmcShower *g2trk = (*(evtRecTrkCol->begin() + iGam[gn2]))->emcShower();

                    kmfit_2g->init();
                    kmfit_2g->AddTrack(0, temp_wtkLam1[i]); // 添加运动学拟合参数-径迹：Lambda径迹
                    kmfit_2g->AddTrack(1, temp_wtkLam2[i]); // 添加运动学拟合参数-径迹：Lambdabar径迹
                    kmfit_2g->AddTrack(2, 0.0, g1trk);      // 添加运动学拟合参数-径迹：gamma1径迹
                    kmfit_2g->AddTrack(3, 0.0, g2trk);      // 添加运动学拟合参数-径迹：gamma2径迹
                    kmfit_2g->AddFourMomentum(0, ecms);     // 添加运动学拟合参数-四动量：质心系四动量 ;
                    if (!kmfit_2g->Fit())
                    {
                        continue;
                    }
                    if (kmfit_2g->chisq() < 0)
                    {
                        continue;
                    }
                    if (kmfit_2g->chisq() < chi2VT3)
                    {
                        chi2VT3 = kmfit_2g->chisq();

                        p4_lam = kmfit_2g->pfit(0);  // 传递Lambda四动量
                        p4_lamb = kmfit_2g->pfit(1); // 传递Lambdabar四动量
                        p4_gam1 = kmfit_2g->pfit(2); // 传递gamma1四动量
                        p4_gam2 = kmfit_2g->pfit(3); // 传递gamma2四动量

                        p4_emc1 = pGam[gn1];
                        p4_emc2 = pGam[gn2];
                        kmfit_fail = false;
                    }
                }
            }
            if (kmfit_fail)
            {
                continue;
            }
            temp_chi2VT3.push_back(chi2VT3);
            temp_p4_lam.push_back(p4_lam);
            temp_p4_lamb.push_back(p4_lamb);
            temp_p4_gam1.push_back(p4_gam1);
            temp_p4_gam2.push_back(p4_gam2);
            temp_p4_emc1.push_back(p4_emc1);
            temp_p4_emc2.push_back(p4_emc2);
            temp_count.push_back(i);
            count1++;
        }

        if (count1 < 1)
        {
            return StatusCode::SUCCESS;
        }

        Vdouble length_lam, length_lam_err; // lambda 衰变长度&误差
        length_lam.clear();
        length_lam_err.clear();
        Vdouble length_lam_bar, length_lam_bar_err; // lambdabar 衰变长度&误差
        length_lam_bar.clear();
        length_lam_bar_err.clear();

        double chi2VT1 = 999.9;
        double chi2VT2 = 999.9;

        chi2VT12 = 1999.9;

        int champion = -1;
        int champion1 = -1;

        for (int i = 0; i < count1; i++)
        {
            if (temp_chi2VT12[temp_count[i]] < chi2VT12)
            {
                chi2VT12 = temp_chi2VT12[temp_count[i]];
                champion = temp_count[i];
                champion1 = i;
            }
        }

        if (champion < 0)
        {
            return StatusCode::SUCCESS;
        }

        Ncut10++;

        length_lam.push_back(temp_length_lam[champion]);
        length_lam_err.push_back(temp_length_lam_err[champion]);
        length_lam_bar.push_back(temp_length_lam_bar[champion]);
        length_lam_bar_err.push_back(temp_length_lam_bar_err[champion]);

        chi2VT1 = temp_chi2VT1[champion];
        chi2VT2 = temp_chi2VT2[champion];
        chi2VT3 = temp_chi2VT3[champion1];

        ptk[0] = temp_p4_pp[champion];
        ptk[1] = temp_p4_pm[champion];
        ptk[2] = temp_p4_pip[champion];
        ptk[3] = temp_p4_pim[champion];
        ptk[4] = temp_p4_emc1[champion1];
        ptk[5] = temp_p4_emc2[champion1];

        p4c[0] = temp_p4_lam[champion1];
        p4c[1] = temp_p4_lamb[champion1];
        p4c[2] = temp_p4_gam1[champion1];
        p4c[3] = temp_p4_gam2[champion1];

        m_ngLam1 = ncddLam1;
        m_secchis1 = chi2VT1;
        m_ngLam2 = ncddLam2;
        m_secchis2 = chi2VT2;
        m_chi1 = chi2VT3;

        m_secchis12diff1 = 9.9;
        m_secchis12diff2 = 9.9;
        m_secchis12diff3 = 9.9;
        m_secchis12diff4 = 9.9;
        m_secchis12diff5 = 9.9;

        m_secchis12diff3 = (p4c[0] + p4c[1]).e();
        m_secchis12diff4 = (temp_p4_pp[t_chi2VT1] + temp_p4_pim[t_chi2VT2] + temp_p4_pip[t_chi2VT1] + temp_p4_pm[t_chi2VT2]).e();
        if (!(chi2VT12 - chi2VT12min))
        {
            m_secchis12diff2 = (temp_p4_pp[t_chi2VT1] + temp_p4_pim[t_chi2VT2] - p4c[0]).mag(); // Delta M(Lambda)
            m_secchis12diff1 = (temp_p4_pip[t_chi2VT1] + temp_p4_pm[t_chi2VT2] - p4c[1]).mag(); // Delta M(anti-Lambda)
            m_secchis12diff5 = m_secchis12diff3 - m_secchis12diff4;
        }

        m_mlengthlam = length_lam[0];                                   // Lambda衰变长度
        m_mlengthlamb = length_lam_bar[0];                              // Lambdabar衰变长度
        m_mlengthlam_err = length_lam_err[0];                           // Lambda衰变长度误差
        m_mlengthlamb_err = length_lam_bar_err[0];                      // Lambdabar衰变长度误差
        m_mlengthlam_rerr = length_lam[0] / length_lam_err[0];          // Lambda相对衰变长度
        m_mlengthlamb_rerr = length_lam_bar[0] / length_lam_bar_err[0]; // Lambdabar相对衰变长度

        m_px_mdc[0] = ptk[0].px();
        m_py_mdc[0] = ptk[0].py();
        m_pz_mdc[0] = ptk[0].pz();
        m_en_mdc[0] = ptk[0].e();
        m_p3_mdc[0] = ptk[0].v().mag(); // pp

        m_px_mdc[1] = ptk[1].px();
        m_py_mdc[1] = ptk[1].py();
        m_pz_mdc[1] = ptk[1].pz();
        m_en_mdc[1] = ptk[1].e();
        m_p3_mdc[1] = ptk[1].v().mag(); // pm

        m_px_mdc[2] = ptk[2].px();
        m_py_mdc[2] = ptk[2].py();
        m_pz_mdc[2] = ptk[2].pz();
        m_en_mdc[2] = ptk[2].e();
        m_p3_mdc[2] = ptk[2].v().mag(); // pip

        m_px_mdc[3] = ptk[3].px();
        m_py_mdc[3] = ptk[3].py();
        m_pz_mdc[3] = ptk[3].pz();
        m_en_mdc[3] = ptk[3].e();
        m_p3_mdc[3] = ptk[3].v().mag(); // pim

        m_px_emc[0] = p4c[2].px();
        m_py_emc[0] = p4c[2].py();
        m_pz_emc[0] = p4c[2].pz();
        m_en_emc[0] = p4c[2].e();
        m_p3_emc[0] = p4c[2].v().mag(); // gamma1

        m_px_emc[1] = p4c[3].px();
        m_py_emc[1] = p4c[3].py();
        m_pz_emc[1] = p4c[3].pz();
        m_en_emc[1] = p4c[3].e();
        m_p3_emc[1] = p4c[3].v().mag(); // gamma2

        m_2gam_bk1 = (p4c[2] + p4c[3]).mag(); // 双光子不变质量

        // 保存lambda lambdabar的四动量

        m_px_trk[0] = p4c[0].px();
        m_py_trk[0] = p4c[0].py();
        m_pz_trk[0] = p4c[0].pz();
        m_en_trk[0] = p4c[0].e();
        m_p3_trk[0] = p4c[0].v().mag();
        m_lam0_bk1 = p4c[0].mag();

        m_px_trk[1] = p4c[1].px();
        m_py_trk[1] = p4c[1].py();
        m_pz_trk[1] = p4c[1].pz();
        m_en_trk[1] = p4c[1].e();
        m_p3_trk[1] = p4c[1].v().mag();
        m_lam1_bk1 = p4c[1].mag();

        // 保存粒子组合的不变质量

        m_mix0_bk1 = (p4c[0] + p4c[2] + p4c[3]).mag();          // 2gamma & Lambda
        m_mix1_bk1 = (p4c[1] + p4c[2] + p4c[3]).mag();          // 2gamma & Lambdabar
        m_mix2_bk1 = (p4c[0] + p4c[1]).mag();                   // Lambda & Lambdabar
        m_mix3_bk1 = (p4c[2] + p4c[0]).mag();                   // gamma1 & Lambda
        m_mix4_bk1 = (p4c[2] + p4c[1]).mag();                   // gamma1 & Lambdabar
        m_mix5_bk1 = (p4c[3] + p4c[0]).mag();                   // gamma2 & Lambda
        m_mix6_bk1 = (p4c[3] + p4c[1]).mag();                   // gamma2 & Lambdabar
        m_mix7_bk1 = (ptk[2] + ptk[3] + p4c[2] + p4c[3]).mag(); // pip & pim & pi0
        m_mix8_bk1 = (ecms - ptk[0] - ptk[1]).mag();            // pppm recoil
        m_mix9_bk1 = (p4c[0] + p4c[1] + p4c[2]).mag();          // gamma1 & Lambda & Lambdabar
        m_mix10_bk1 = (p4c[0] + p4c[1] + p4c[3]).mag();         // gamma2 & Lambda & Lambdabar

        m_totm_bk = (p4c[0] + p4c[1] + p4c[2] + p4c[3]).mag();

    } // 4c

    //counter[5]++;
    m_tuple6->write();
    return StatusCode::SUCCESS;
}

//////////////////////////////////////////////////////////////////////////////////////////

StatusCode psip::finalize()
{
    cout << "total number        : " << Ncut0 << endl;
    cout << "nGood >= 4          : " << Ncut2 << endl;
    //cout << "4 tracks for lamb   : " << Ncut2 << endl;
    cout << "nGam >= 2           : " << Ncut3 << endl;
    cout << "2 vertex fit ok     : " << Ncut4 << endl;
    //cout << "lambda fit          : " << Ncut5 << endl;
    //cout << "pm pip vertex fit   : " << Ncut6 << endl;
    //cout << "lambda_bar fit      : " << Ncut7 << endl;
    //cout << "lambda and lambda_bar fit     : " << Ncut8 << endl;
    //cout << "Nlam >= 1           : " << counter[0] << endl;
    //cout << "Nlambar >= 1        : " << counter[1] << endl;
    //cout << "pass 4c 0g          : " << counter[2] << endl;
    //cout << "pass 4c 1g          : " << Ncut9 << endl;
    cout << "pass 4c             : " << Ncut10 << endl;
    //cout << "p pi- to lamb       : " << counter[3] << endl;
    //cout << "pbar pi+ to lambar  : " << counter[4] << endl;
    //cout << "fill tree           : " << counter[5] << endl;
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in finalize()" << endmsg;
    return StatusCode::SUCCESS;
}

void psip::InitVar()
{
    m_chi1 = 9999.0;
    m_chi1_0g = 9999.0;
    m_secchis1 = 9999.0;
    m_secchis2 = 9999.0;
    for (int i2 = 0; i2 < 2; i2++)
    {
        m_px_trk[i2] = -10;
        m_py_trk[i2] = -10;
        m_pz_trk[i2] = -10;
        m_en_trk[i2] = -10;
        m_p3_trk[i2] = -10;
    }
    for (int i3 = 0; i3 < 4; i3++)
    {
        m_px_mdc[i3] = -10;
        m_py_mdc[i3] = -10;
        m_pz_mdc[i3] = -10;
        m_en_mdc[i3] = -10;
        m_p3_mdc[i3] = -10;
        m_pull_1[i3] = -999.0;
        m_pull_2[i3] = -999.0;
        m_pull_3[i3] = -999.0;
        m_pull_4[i3] = -999.0;
        m_pull_5[i3] = -999.0;
    }
}

// corset(): sets up the generation by calculating C from V.

void psip::corset(HepSymMatrix &V, HepMatrix &C, int n)
{

    //cout<<"v="<<V<<endl;
    //cout<<"c="<<C<<endl;

    double sum;

    // Compute square root of matrix sigma

    for (int j = 0; j < n; j++)
    {
        sum = 0;
        for (int k = 0; k < j; k++)
        {
            sum = sum + C[j][k] * C[j][k];
        }

        //cout<<"sum="<<sum<<endl;
        //cout<<"v("<<j<<","<<j<<")="<<V[j][j]<<endl;

        C[j][j] = sqrt(abs(V[j][j] - sum));

        //cout<<"c("<<j<<","<<j<<")="<<C[j][j]<<endl;

        // Off Diagonal terms

        for (int i = j + 1; i < n; i++)
        {
            sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum = sum + C[i][k] * C[j][k];
            }
            C[i][j] = (V[i][j] - sum) / C[j][j];
        }
    }
}

// Calibration for helix parameters
// corgen(): generates a set of n random numbers Gaussian-distributed with covariance
// matrix V (V = C*C') and mean values zero.
void psip::corgen(HepMatrix &C, HepVector &x, int n)
{
    int i, j;
    int nmax = 100;

    if (n > nmax)
    {
        printf("Error in corgen: array overflown");
    }

    double tmp[3];
    for (int p = 0; p < n; p++)
    {
        tmp[p] = gRandom->Gaus(0, 1);

        //cout<<"tmp["<<p<<"]="<<tmp[p]<<endl;
    }
    for (i = 0; i < n; i++)
    {
        x[i] = 0.0;
        for (j = 0; j <= i; j++)
        {
            x[i] = x[i] + C[i][j] * tmp[j];
        }
    }
}

void psip::calibration(RecMdcKalTrack *trk, HepVector &wtrk_zHel, int n)
{

    HepVector pip_calerr_d2(5, 0);
    HepVector pim_calerr_d2(5, 0);
    HepVector pp_calerr_d2(5, 0);
    HepVector pm_calerr_d2(5, 0);

    HepVector pip_calmean_d2(5, 0);
    HepVector pim_calmean_d2(5, 0);
    HepVector pp_calmean_d2(5, 0);
    HepVector pm_calmean_d2(5, 0);

    pip_calerr_d2[0] = 1.0;
    pip_calerr_d2[1] = 1.18;
    pip_calerr_d2[2] = 1.21;
    pip_calerr_d2[3] = 1.0;
    pip_calerr_d2[4] = 1.13;

    pim_calerr_d2[0] = 1.0;
    pim_calerr_d2[1] = 1.18;
    pim_calerr_d2[2] = 1.21;
    pim_calerr_d2[3] = 1.0;
    pim_calerr_d2[4] = 1.12;
    pm_calmean_d2[1] = 0.045;
    pm_calmean_d2[2] = -0.033;
    pm_calmean_d2[3] = 0;
    pm_calmean_d2[4] = 0.170;
    if (trk->charge() > 0 && n == 0)
    {
        // pip
        HepSymMatrix wpip_zerr(5, 0);
        wpip_zerr = trk->getZError();
        HepSymMatrix wpip_zcal(3, 0);

        wpip_zcal[0][0] = (pip_calerr_d2[1] * pip_calerr_d2[1] - 1) * wpip_zerr[1][1];
        wpip_zcal[1][1] = (pip_calerr_d2[2] * pip_calerr_d2[2] - 1) * wpip_zerr[2][2];
        wpip_zcal[2][2] = (pip_calerr_d2[4] * pip_calerr_d2[4] - 1) * wpip_zerr[4][4];

        HepMatrix wpip_zerrc(3, 3, 0);
        psip::corset(wpip_zcal, wpip_zerrc, 3);
        HepVector wpip_zgen(3, 0);
        psip::corgen(wpip_zerrc, wpip_zgen, 3);

        wtrk_zHel[0] = trk->getZHelix()[0];
        wtrk_zHel[1] = trk->getZHelix()[1] + pip_calmean_d2[1] * sqrt(wpip_zerr[1][1]) + wpip_zgen[0];
        wtrk_zHel[2] = trk->getZHelix()[2] + pip_calmean_d2[2] * sqrt(wpip_zerr[2][2]) + wpip_zgen[1];
        wtrk_zHel[3] = trk->getZHelix()[3];
        wtrk_zHel[4] = trk->getZHelix()[4] + pip_calmean_d2[4] * sqrt(wpip_zerr[4][4]) + wpip_zgen[2];
    }
    if (trk->charge() < 0 && n == 0)
    {
        // pim
        HepSymMatrix wpim_zerr(5, 0);
        wpim_zerr = trk->getZError();

        HepSymMatrix wpim_zcal(3, 0);

        wpim_zcal[0][0] = (pim_calerr_d2[1] * pim_calerr_d2[1] - 1) * wpim_zerr[1][1];
        wpim_zcal[1][1] = (pim_calerr_d2[2] * pim_calerr_d2[2] - 1) * wpim_zerr[2][2];
        wpim_zcal[2][2] = (pim_calerr_d2[4] * pim_calerr_d2[4] - 1) * wpim_zerr[4][4];

        HepMatrix wpim_zerrc(3, 3, 0);
        psip::corset(wpim_zcal, wpim_zerrc, 3);
        HepVector wpim_zgen(3, 0);
        psip::corgen(wpim_zerrc, wpim_zgen, 3);

        wtrk_zHel[0] = trk->getZHelix()[0];
        wtrk_zHel[1] = trk->getZHelix()[1] + pim_calmean_d2[1] * sqrt(wpim_zerr[1][1]) + wpim_zgen[0];
        wtrk_zHel[2] = trk->getZHelix()[2] + pim_calmean_d2[2] * sqrt(wpim_zerr[2][2]) + wpim_zgen[1];
        wtrk_zHel[3] = trk->getZHelix()[3];
        wtrk_zHel[4] = trk->getZHelix()[4] + pim_calmean_d2[4] * sqrt(wpim_zerr[4][4]) + wpim_zgen[2];
    }
    if (trk->charge() > 0 && n == 1)
    {
        // kp
        HepSymMatrix wpp_zerr(5, 0);
        wpp_zerr = trk->getZErrorK();

        HepSymMatrix wpp_zcal(3, 0);

        wpp_zcal[0][0] = (pp_calerr_d2[1] * pp_calerr_d2[1] - 1) * wpp_zerr[1][1];
        wpp_zcal[1][1] = (pp_calerr_d2[2] * pp_calerr_d2[2] - 1) * wpp_zerr[2][2];
        wpp_zcal[2][2] = (pp_calerr_d2[4] * pp_calerr_d2[4] - 1) * wpp_zerr[4][4];

        HepMatrix wpp_zerrc(3, 3, 0);
        psip::corset(wpp_zcal, wpp_zerrc, 3);
        HepVector wpp_zgen(3, 0);
        psip::corgen(wpp_zerrc, wpp_zgen, 3);

        wtrk_zHel[0] = trk->getZHelixK()[0];
        wtrk_zHel[1] = trk->getZHelixK()[1] + pp_calmean_d2[1] * sqrt(wpp_zerr[1][1]) + wpp_zgen[0];
        wtrk_zHel[2] = trk->getZHelixK()[2] + pp_calmean_d2[2] * sqrt(wpp_zerr[2][2]) + wpp_zgen[1];
        wtrk_zHel[3] = trk->getZHelixK()[3];
        wtrk_zHel[4] = trk->getZHelixK()[4] + pp_calmean_d2[4] * sqrt(wpp_zerr[4][4]) + wpp_zgen[2];
    }
    if (trk->charge() < 0 && n == 1)
    {
        // km
        HepSymMatrix wpm_zerr(5, 0);
        wpm_zerr = trk->getZErrorK();

        HepSymMatrix wpm_zcal(3, 0);

        wpm_zcal[0][0] = (pm_calerr_d2[1] * pm_calerr_d2[1] - 1) * wpm_zerr[1][1];
        wpm_zcal[1][1] = (pm_calerr_d2[2] * pm_calerr_d2[2] - 1) * wpm_zerr[2][2];
        wpm_zcal[2][2] = (pm_calerr_d2[4] * pm_calerr_d2[4] - 1) * wpm_zerr[4][4];

        HepMatrix wpm_zerrc(3, 3, 0);
        psip::corset(wpm_zcal, wpm_zerrc, 3);
        HepVector wpm_zgen(3, 0);
        psip::corgen(wpm_zerrc, wpm_zgen, 3);

        wtrk_zHel[0] = trk->getZHelixK()[0];
        wtrk_zHel[1] = trk->getZHelixK()[1] + pm_calmean_d2[1] * sqrt(wpm_zerr[1][1]) + wpm_zgen[0];
        wtrk_zHel[2] = trk->getZHelixK()[2] + pm_calmean_d2[2] * sqrt(wpm_zerr[2][2]) + wpm_zgen[1];
        wtrk_zHel[3] = trk->getZHelixK()[3];
        wtrk_zHel[4] = trk->getZHelixK()[4] + pm_calmean_d2[4] * sqrt(wpm_zerr[4][4]) + wpm_zgen[2];
    }
}
