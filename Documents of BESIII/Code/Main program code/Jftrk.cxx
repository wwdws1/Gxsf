#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "JftrkAlg/Jftrk.h"
#include "ParticleID/ParticleID.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/Helix.h"

#include "McTruth/McParticle.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/DecayMode.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/MucMcHit.h"
#include "McTruth/McEvent.h"

#include "RooFit.h"
#include "TRandom.h"
#include "CLHEP/Random/RandGauss.h"

#include <vector>
const double PI   = 3.1415927;
const double me = 0.000511;
const double mmu= 0.105658;
const double mpi= 0.139570;
const double mk = 0.493677;
const double mp = 0.938272;

const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;           //定义int形矢量Vint
typedef std::vector<double> Vdouble;     //定义double形矢量Vdouble
typedef std::vector<WTrackParameter> WTk;//定义径迹参数矢量WTk
typedef std::vector<HepLorentzVector> Vp4;//定义洛伦兹矢量Vp4

int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6, Ncut7, Ncut8, Ncut9, Ncut10, Ncut11, Ncut12;
int counter[6] = {0,0,0,0,0,0};  //数组存放变量
int evt_mc = 0;
int couter_ph = 0;
int pdgid58=0;
int pdgid59=0;
int pdgid441=0;int pdgid9080221=0;int pdgid333=0;


Jftrk::Jftrk(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {
		declareProperty("Vr0cut", m_vr0cut=1.0);                     //在r方向上截取1cm
		declareProperty("Vz0cut", m_vz0cut=10.0);                    //在z方向上截取10cm
		declareProperty("EnergyThreshold", m_energyThreshold=0.04);  //能量阈值0.04GeV
		declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);         //光子phi角截取20度
		declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);     //光子theta角截取20度
		declareProperty("GammaAngleCut", m_gammaAngleCut=10.0);    //光子angle截取10度
		declareProperty("CheckDedx", m_checkDedx = 1);   
		declareProperty("CheckTof",  m_checkTof = 1);        //Tof
		declareProperty("Test4k4C", m_test4k4C = 1);        //4C
		declareProperty("Test4k5C", m_test4k5C = 1);        //5C
		declareProperty("Test4k5C_4C", m_test4k5C_4C = 1);
		declareProperty("Testmctruth", m_testmctruth = 1);
	}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$(^*_*^)
StatusCode Jftrk::initialize(){           //初始化
	MsgStream log(msgSvc(), name());

	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;
	if( m_testmctruth == 1) 
	{
		NTuplePtr nt6(ntupleSvc(), "FILE1/trkRes");
		if ( nt6 ) m_tuple6 = nt6;
		else 
		{
			m_tuple6 = ntupleSvc()->book ("FILE1/trkRes", CLID_ColumnWiseTuple, "ks N-Tuple example");
			if ( m_tuple6 )    
			{
				status = m_tuple6->addItem("indexmc",          m_idxmc, 0, 500);
				status = m_tuple6->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
				status = m_tuple6->addIndexedItem("motheridx", m_idxmc, m_motheridx);
				status = m_tuple6->addIndexedItem("trkidx",     m_idxmc, m_trkidx);

				status = m_tuple6->addIndexedItem("mc_vertex",    m_idxmc, vert_mcpar1);

				status = m_tuple6->addItem ("trk_px", 3,m_px_trk);
				status = m_tuple6->addItem ("trk_py", 3,m_py_trk);
				status = m_tuple6->addItem ("trk_pz", 3,m_pz_trk);
				status = m_tuple6->addItem ("trk_p3", 3,m_p3_trk);
				status = m_tuple6->addItem ("trk_en", 3,m_en_trk);
				status = m_tuple6->addItem ("chi2_lamb", chi2_lamb);
				status = m_tuple6->addItem ("chi2_lambar", chi2_lambar);
				//////事例选择需要比较的变量mdc_p3
				status = m_tuple6->addItem ("mdc_px", 4,m_px_mdc); ///p pi- p-bar pi+
				status = m_tuple6->addItem ("mdc_py", 4,m_py_mdc); 
				status = m_tuple6->addItem ("mdc_pz", 4,m_pz_mdc);
				status = m_tuple6->addItem ("mdc_p3", 4,m_p3_mdc);
				status = m_tuple6->addItem ("mdc_en", 4,m_en_mdc);

				status = m_tuple6->addItem ("trk_nCh", chargTrk_index, 0,100);
				status = m_tuple6->addIndexedItem ("trk_Rxy", chargTrk_index, m_vtx_Rxy);
				status = m_tuple6->addIndexedItem ("trk_cos", chargTrk_index, m_vtx_cos);  
				status = m_tuple6->addIndexedItem ("trk_z0",  chargTrk_index,m_vtx_z0);
				status = m_tuple6->addIndexedItem ("trk_p0",  chargTrk_index,m_vtx_p);
				status = m_tuple6->addItem ("trk_ach",m_totE_ch);
				status = m_tuple6->addItem ("trk_npp",m_npp);
				status = m_tuple6->addItem ("trk_npm",m_npm);
				status = m_tuple6->addItem ("trk_npip",m_npip);
				status = m_tuple6->addItem ("trk_npim",m_npim);

				status = m_tuple6->addItem ("gam_nGam", gamma_index, 0,100);
				status = m_tuple6->addIndexedItem ("gam_dthe", gamma_index,m_dthe);
				status = m_tuple6->addIndexedItem ("gam_dphi", gamma_index,m_dphi);
				status = m_tuple6->addIndexedItem ("gam_dang", gamma_index,m_dang);
				status = m_tuple6->addIndexedItem ("gam_eraw", gamma_index,m_eraw);
				status = m_tuple6->addIndexedItem ("gam_TDC",  gamma_index,m_TDC);
				////需要通过gam_tag的值判断该光子与哪条径迹的夹角
				status = m_tuple6->addIndexedItem ("gam_tagg", gamma_index,m_tagg);
				///////gamVspbar////
				status = m_tuple6->addItem ("orig_m2pi",    m_phi1_bk);
				status = m_tuple6->addItem ("orig_m2p",     m_phi2_bk);
				status = m_tuple6->addItem ("orig_mpppim",    m_phi1_bk1);
				status = m_tuple6->addItem ("orig_mpmpip",    m_phi2_bk1);
				status = m_tuple6->addItem ("orig_m4trk",     m_2phi_bk);
				/////事例挑选需要比较的量
				status = m_tuple6->addItem ("ki_mpppim",  m_phi1_true);
				status = m_tuple6->addItem ("ki_mpmpip",  m_phi2_true);
				status = m_tuple6->addItem ("ki_metac",   m_etac_true);
				status = m_tuple6->addItem ("ki_chi2",      m_chi1);
				status = m_tuple6->addItem ("ki0g_chi2",    m_chi1_0g); 

				status = m_tuple6->addItem ("lambda_n",     m_ngLam1);
				status = m_tuple6->addItem ("lambdab_n",    m_ngLam2);
				status = m_tuple6->addItem ("lambda_chi",   m_secchis1);
				//////需要比较的量
				status = m_tuple6->addItem ("lambda_lth",   m_mlengthlam);
				status = m_tuple6->addItem ("lambdab_lth",  m_mlengthlamb);

				status = m_tuple6->addItem ("lambda_lth_e", m_mlengthlam_err);
				status = m_tuple6->addItem ("lambdab_lth_e",m_mlengthlamb_err);
				status = m_tuple6->addItem ("lambdab_chi",  m_secchis2);

			} ///// if m_tuple6
			else    
			{
				log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple6) << endmsg;
				return StatusCode::FAILURE;
			}
		}
		//////for mctruth
		NTuplePtr nt11(ntupleSvc(), "FILE1/mctruth");  
		if ( nt11 ) m_tuple11 = nt11;
		else { 
			m_tuple11 = ntupleSvc()->book ("FILE1/mctruth", CLID_ColumnWiseTuple, "NTuple");
			if ( m_tuple11 )    {  
				status = m_tuple11->addItem("mc_ng",        m_mc_ng);
				status = m_tuple11->addItem("mc_g_px", parti_px_mc1);
				status = m_tuple11->addItem("mc_g_py", parti_py_mc1);
				status = m_tuple11->addItem("mc_g_pz", parti_pz_mc1);
				status = m_tuple11->addItem("mc_g_e",  parti_e_mc1);
				status = m_tuple11->addItem("mc_g_m",  parti_m_mc1);

				status = m_tuple11->addItem("mc_nlam",   m_mc_nlam);
				status = m_tuple11->addItem("mc_lam_px", parti_px_mc2);
				status = m_tuple11->addItem("mc_lam_py", parti_py_mc2);
				status = m_tuple11->addItem("mc_lam_pz", parti_pz_mc2);
				status = m_tuple11->addItem("mc_lam_e",  parti_e_mc2);
				status = m_tuple11->addItem("mc_lam_m",  parti_m_mc2);

				status = m_tuple11->addItem("mc_nlambar",   m_mc_nlambar);
				status = m_tuple11->addItem("mc_lamb_px", parti_px_mc3);
				status = m_tuple11->addItem("mc_lamb_py", parti_py_mc3);
				status = m_tuple11->addItem("mc_lamb_pz", parti_pz_mc3);
				status = m_tuple11->addItem("mc_lamb_e",  parti_e_mc3);
				status = m_tuple11->addItem("mc_lamb_m",  parti_m_mc3);
				/*
				   status = m_tuple11->addItem("mc_npp",        m_mc_npp);
				   status = m_tuple11->addItem("mc_pp_px", parti_px_mc4);
				   status = m_tuple11->addItem("mc_pp_py", parti_py_mc4);
				   status = m_tuple11->addItem("mc_pp_pz", parti_pz_mc4);
				   status = m_tuple11->addItem("mc_pp_e",  parti_e_mc4);
				   status = m_tuple11->addItem("mc_pp_m",  parti_m_mc4);

				   status = m_tuple11->addItem("mc_npim",   m_mc_npim);
				   status = m_tuple11->addItem("mc_pim_px", parti_px_mc5);
				   status = m_tuple11->addItem("mc_pim_py", parti_py_mc5);
				   status = m_tuple11->addItem("mc_pim_pz", parti_pz_mc5);
				   status = m_tuple11->addItem("mc_pim_e",  parti_e_mc5);
				   status = m_tuple11->addItem("mc_pim_m",  parti_m_mc5);

				   status = m_tuple11->addItem("mc_npm",   m_mc_npm);
				   status = m_tuple11->addItem("mc_pm_px", parti_px_mc6);
				   status = m_tuple11->addItem("mc_pm_py", parti_py_mc6);
				   status = m_tuple11->addItem("mc_pm_pz", parti_pz_mc6);
				   status = m_tuple11->addItem("mc_pm_e",  parti_e_mc6);
				   status = m_tuple11->addItem("mc_pm_m",  parti_m_mc6);

				   status = m_tuple11->addItem("mc_npip",   m_mc_npip);
				   status = m_tuple11->addItem("mc_pip_px", parti_px_mc7);
				   status = m_tuple11->addItem("mc_pip_py", parti_py_mc7);
				   status = m_tuple11->addItem("mc_pip_pz", parti_pz_mc7);
				   status = m_tuple11->addItem("mc_pip_e",  parti_e_mc7);
				   status = m_tuple11->addItem("mc_pip_m",  parti_m_mc7);
				 */
			}
			else    {
				log << MSG::ERROR << "    Cannot book N-tuple:"
					<< long(m_tuple11) << endmsg;
				return StatusCode::FAILURE;
			}


		}

	}
	return StatusCode::SUCCESS;
}

//************************************************************************(*^_^*)
StatusCode Jftrk::execute() {    //主要的程序

	Jftrk::InitVar();

	//	cout << __LINE__<<endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	HepLorentzVector ecms(3.097*sin(0.011), 0, 0, 3.097);   //质心系能量3.097
	//	HepLorentzVector ecms(0.034,0,0,3.097);  洛伦兹矢量
	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNum = eventHeader->runNumber();
	int eventNum = eventHeader->eventNumber();


	if (runNum == -27878 || runNum == -27875) return StatusCode::SUCCESS;

	if (runNum < 0 || runNum == 0){
		ecms.setPx(0.03407);ecms.setPy(0.0);ecms.setPz(0.0);ecms.setE(3.09725);
	}             // 如果runNum小于或等于零，则px=0.03407，py=0，pz=0，E=3.09725  
	Ncut0++;      //如果runNum > 0 ，事例累加（总的事例数）
	//	cout << __LINE__<<endl;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	log << MSG::DEBUG <<"ncharg, nneu, tottks = "
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;
	if ( evtRecEvent->totalTracks() > 99 ) return SUCCESS; 
	//如果总的带电径迹数目大于99，返回SUCCESS。（说明我要求有足够多的径迹，否则不要）

	////////////////////////////////////////////////////////////////////////
	if (runNum < 0){
		int numParticle=0;
		int mcng=0; int mcnlam=0; int mcnlambar=0; 
		int mcnpp=0; int mcnpm=0; int mcnpip=0; int mcnpim=0; 
		//如果runNum小于0，清空粒子数、mcng、lam、lambar、p+、p-、pi+、pi-

		bool Decay=false;
		int rootIndex=-1;
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		if (!mcParticleCol){
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}
		else{
			Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
			bool strange = false;
			for(; iter_mc!=mcParticleCol->end(); iter_mc++){
				int temp_id = (*iter_mc)->particleProperty();
				int motherpdg = ((*iter_mc)->mother()).particleProperty();

				HepLorentzVector  mctrue_track = (*iter_mc)->initialFourMomentum();  //原始四动量
				HepLorentzVector mctrack_iniposition = (*iter_mc)->initialPosition();//原始位置
				HepLorentzVector mctrack_finposition = (*iter_mc)->finalPosition();  //最后位置

				if((*iter_mc)->primaryParticle() && temp_id == 11 && motherpdg == 11) {strange = true;}
				if((*iter_mc)->primaryParticle()) continue;
				if(!(*iter_mc)->decayFromGenerator()) continue;

				if(temp_id == 443){
					rootIndex = (*iter_mc)->trackIndex();
					Decay = true;
				}
				if(!Decay) continue;
				int mcidx = ((*iter_mc)->mother()).trackIndex()-rootIndex;
				int trkidx=(*iter_mc)->trackIndex() - rootIndex;
				if(strange && motherpdg != 443)mcidx--;
				m_pdgid[numParticle] = temp_id;
				m_motheridx[numParticle] = mcidx;
				m_trkidx[numParticle] = trkidx;
				vert_mcpar1[numParticle]  = (*iter_mc)->vertexIndex0();
				numParticle++;
				//cout << " mother " << motherpdg << " itself " << temp_id << endl;

				double mcpx = (*iter_mc)->initialFourMomentum().px(); //最初px
				double mcpy = (*iter_mc)->initialFourMomentum().py(); //最初py
				double mcpz = (*iter_mc)->initialFourMomentum().pz(); //最初pz
				double mcen = (*iter_mc)->initialFourMomentum().e();  //最初能量
				HepLorentzVector p4mc;  
				p4mc.setPx(mcpx); p4mc.setPy(mcpy);
				p4mc.setPz(mcpz); p4mc.setE(mcen);
				double mmcpar = p4mc.mag();

				//cout << " mass of the particle " << p4mc.mag() << endl;

				if ( motherpdg == 443 && temp_id == 22 ){
					parti_px_mc1 = mcpx; 
					parti_py_mc1 = mcpy;
					parti_pz_mc1 = mcpz;
					parti_e_mc1  = mcen;
					parti_m_mc1 = mmcpar;
					mcng++;                      //mcng数量、px、py、pz、能量
				}
				else if ( motherpdg == 443 && temp_id == 3122 ){
					parti_px_mc2 = mcpx;
					parti_py_mc2 = mcpy;
					parti_pz_mc2 = mcpz;
					parti_e_mc2  = mcen;
					parti_m_mc2 = mmcpar;
					mcnlam++;                  //lam数量、px、py、pz、能量
				}
				else if ( motherpdg == 443 && temp_id == -3122 ){
					parti_px_mc3 = mcpx;
					parti_py_mc3 = mcpy;
					parti_pz_mc3 = mcpz;
					parti_e_mc3  = mcen;
					parti_m_mc3 = mmcpar;
					mcnlambar++;               //lambar数量、px、py、pz、能量
				}
				/*
				   else if ( motherpdg == 3122 && temp_id == 2212 ){
				   parti_px_mc4 = mcpx;
				   parti_py_mc4 = mcpy;
				   parti_pz_mc4 = mcpz;
				   parti_e_mc4 = mcen;
				   parti_m_mc4 = mmcpar;
				   mcnpp++;                 //p+数量、px、py、pz、能量
				   }

				   else if ( motherpdg == 3122 && temp_id == -211 ){
				   parti_px_mc5 = mcpx; 
				   parti_py_mc5 = mcpy; 
				   parti_pz_mc5 = mcpz;
				   parti_e_mc5 = mcen;
				   parti_m_mc5 = mmcpar;
				   mcnpim++;              //pi-数量、px、py、pz、能量
				   }

				   else if ( motherpdg == -3122 && temp_id == -2212 ){
				   parti_px_mc6 = mcpx;
				   parti_py_mc6 = mcpy;
				   parti_pz_mc6 = mcpz;
				   parti_e_mc6  = mcen;
				   parti_m_mc6 = mmcpar;
				   mcnpm++;                //p-数量、px、py、pz、能量
				   }
				   else if ( motherpdg == -3122 && temp_id == 211 ){
				   parti_px_mc7 = mcpx;
				   parti_py_mc7 = mcpy;
				   parti_pz_mc7 = mcpz;
				   parti_e_mc7 =  mcen;
				   parti_m_mc7 = mmcpar;
				   mcnpip++;                //pi+数量、px、py、pz、能量
				   }
				 */
				else continue;
			}
		}
		m_idxmc = numParticle;                                        //粒子数
		m_mc_ng = mcng; m_mc_nlam = mcnlam; m_mc_nlambar = mcnlambar; //lambar数
		//		m_mc_npp = mcnpp; m_mc_npim = mcnpim; m_mc_npm = mcnpm;       //p-数
		//		m_mc_npip = mcnpip;                                           //pi+数
		m_tuple11->write(); 
	}

	/////////

	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

	Vint iGood, icharg;  //for charged plus and minus  
	iGood.clear();
	icharg.clear();      // 定义iGood和icharg类型为整形并清空数据（好的事例和好的电量）

	Vint ichpip, ichpim, ichpp, ichpm;  //for charged plus and minus  
	ichpip.clear();
	ichpim.clear();
	ichpp.clear();
	ichpm.clear();  //定义整形的pi+、pi—、p+、p—(电量)个数并清空里面的数据

	Vp4 pchpip, pchpim, pchpp, pchpm;
	pchpip.clear();
	pchpim.clear();
	pchpp.clear();
	pchpm.clear();  //定义pi+、pi—、p+、p—的动量并清零

	int nCharge = 0;

	Hep3Vector xorigin(0,0,0);
	Hep3Vector eorigin(0,0,0);   //定义原点
	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
		double* dbv = vtxsvc->PrimaryVertex();
		double*  vv = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
		eorigin.setX(vv[0]);
		eorigin.setY(vv[1]);
		eorigin.setZ(vv[2]);   // 将X Y Z 保存在数组中
	}

	int iv=0;
	double etot = 0.0;
	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!((*itTrk)->isMdcKalTrackValid())) continue;       //不是MDC的轨迹不要
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();  //重建MDC轨迹

		HepVector a = mdcKalTrk->helix();
		HepSymMatrix Ea = mdcKalTrk->err();
		HepPoint3D point0(0.,0.,0.);
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);    //定义IP值
		VFHelix helixip(point0,a,Ea);   //螺旋
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]); //the distance to IP in xy plane  
		double  Rvz0=vecipa[3];    //the distance to IP in z direction  
		double  Rvphi0=vecipa[1];
		double cos_1 = cos(mdcKalTrk->theta());
		double pch = mdcKalTrk->p();    
		if(fabs(cos_1) > 0.93)continue;  // 如果cos theta绝对值的值大于0.93就不要（不得超过MDC测量范围）
		if(pch > 2.0)continue;    // 如果粒子动量大于2.0也不要（总的动量是3.097，一半不可能大于2）

		m_vtx_cos[iv] = cos_1;
		m_vtx_Rxy[iv] = Rvxy0;
		m_vtx_z0[iv] = Rvz0;
		m_vtx_p[iv]  = pch;  
		iv++;

		double trk_px = mdcKalTrk->px();
		double trk_py = mdcKalTrk->py();
		double trk_pz = mdcKalTrk->pz();  //定义轨迹的三动量
		HepLorentzVector ptrk;  //ptrk为洛伦兹矢量
		ptrk.setPx(trk_px);
		ptrk.setPy(trk_py);
		ptrk.setPz(trk_pz);    
		double p3 = ptrk.mag();

		if ( pch > 0.4 ){
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::proton);
			//			ptrk.setE(sqrt(p3*p3+mp*mp));
			ptrk.setE(sqrt(pch*pch+mp*mp));   
			//如果总动量大于0.4，可以确定粒子是质子，则可以算出能量E = sqrt(p3*p3+mp*mp） 
			if ( mdcKalTrk->charge() > 0 ) {
				ichpp.push_back(i);
				pchpp.push_back(ptrk);
			}
			//当带电量大于零（带正电）时，返回p+的个数和p+的径迹 
			else {
				ichpm.push_back(i);
				pchpm.push_back(ptrk);	
			}
			//当带电量为负时，返回p—的个数和p—的径迹
		}
		else {
			RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);
			//			ptrk.setE(sqrt(p3*p3+mpi*mpi));
			ptrk.setE(sqrt(pch*pch+mpi*mpi));
			//当总动量小于0.4时，确定粒子为pi
			if ( mdcKalTrk->charge() > 0 ) {
				ichpip.push_back(i);
				pchpip.push_back(ptrk);
			}
			//如果带电量大于零，返回pi+的个数和pi+的径迹
			else {
				ichpim.push_back(i);
				pchpim.push_back(ptrk);
			}
			//如果带电量小于零，返回pi—的个数和pi—的径迹
		}
		iGood.push_back(i);
		icharg.push_back(mdcKalTrk->charge());
		nCharge += mdcKalTrk->charge();       //总电量
		if (!(*itTrk)->isEmcShowerValid())continue;
		etot += (*itTrk)->emcShower()->energy();
	}

	m_totE_ch = etot;    //总能量

	int nGood = iGood.size();
	chargTrk_index = nGood;
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;

	//	if( nGood != 4 ||  nCharge != 0  ) return SUCCESS;
	//	if( nGood < 4 ) return SUCCESS;
	Ncut1++;      //Good数目
	if( ichpip.size() != 1 || ichpim.size() != 1 ) return SUCCESS;
	if( ichpp.size() != 1 || ichpm.size() != 1 ) return SUCCESS;
	Ncut2++;     
	//如果pi+，pi—，p+，p—各有一个，累加到Ncut2中
	m_npp = ichpp.size();  m_npm = ichpm.size();
	m_npip = ichpip.size();  m_npim = ichpim.size();


	Vint iGam;
	iGam.clear();
	int  ii=0;
	//cout << " event no. is " << eventNum << endl;

	//cout << " total track " << evtRecEvent->totalTracks() << endl;
	//cout << " total charged track " << evtRecEvent->totalCharged() << endl;

	for(int i=evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
		//挑选出EMC的事例
	{
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());

		double dthe = 200.;
		double dphi = 200.;
		double dang = 200.;
		int tsttag = 0;
		for(int j = 0; j < nGood; j++) {
			EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + iGood[j];
			if(!(*jtTrk)->isExtTrackValid()) continue;
			RecExtTrack *extTrk = (*jtTrk)->extTrack();
			if(extTrk->emcVolumeNumber() == -1) continue;
			Hep3Vector extpos = extTrk->emcPosition();

			double angd = extpos.angle(emcpos);
			double thed = extpos.theta() - emcpos.theta();
			double phid = extpos.deltaPhi(emcpos);	  //定义angd,thed和phid角
			if(angd < dang){          //如果angd小于200
				dang = angd;
				dthe = thed;
				dphi = phid;
				tsttag = iGood[j];
			}
		}
		dthe = dthe * 180 / (CLHEP::pi);
		dphi = dphi * 180 / (CLHEP::pi);
		dang = dang * 180 / (CLHEP::pi);  //将弧度转化为角度
		//		if( dang < 10) continue;
		double eraw = emcTrk->energy();   //能量
		double tdc = emcTrk->time();      //时间
		//     if(eraw < m_energyThreshold) continue;  
		double costh = cos(emcTrk->theta());
		bool ok1 = 1;
		bool ok2 = 1;
		if( (fabs(costh)>0.8&&fabs(costh)<0.86) || fabs(costh)>0.92 )continue; 
		//不得超出电磁量能器的测量范围
		if(fabs(costh)<0.8){            //桶部测量范围（costh绝对值小于0.8）
			if (eraw < 0.025) ok1 = 0;  //桶部能量最小分辨率
		}
		if (fabs(costh) > 0.86 && fabs(costh) < 0.92){   //端盖测量范围 
			if (eraw < 0.050) ok2 = 0;                   //端盖能量最小分辨率
		}
		if(!ok1)continue;
		if(!ok2)continue;   //超过测量范围的都不要
		if(fabs(dang) < m_gammaAngleCut) continue;   // dang角的绝对值小于10不要 
		if(tdc < 0 || tdc > 15) continue;            //  飞行时间要在0到15ns之间
		m_dthe[ii] = dthe;
		m_dphi[ii] = dphi;
		m_dang[ii] = dang;
		m_eraw[ii] = eraw;
		m_TDC[ii] = emcTrk->time();                // 记下dthe、dphi、dang角、能量和时间
		if ( fabs( tsttag - ichpp[0] ) < 0.001 ) m_tagg[ii] = 4;
		else if ( fabs( tsttag - ichpm[0] ) < 0.001 ) m_tagg[ii] = 3;
		else if ( fabs( tsttag - ichpip[0] ) < 0.001 ) m_tagg[ii] = 2;
		else if ( fabs( tsttag - ichpim[0] ) < 0.001 ) m_tagg[ii] = 1;
		else m_tagg[ii] = 0;                      //将p+、p-、pi+、pi-比较并保存在数组中
		ii++;
		iGam.push_back(i);
	}

	// Finish Good Photon Selection  

	int nGam = iGam.size();
	gamma_index=nGam;

	log << MSG::DEBUG << "num Good Photon " << nGam  << " , "
		<<evtRecEvent->totalNeutral()<<endreq;
	if( nGam > 10 ) return StatusCode::SUCCESS;
	// 如果光子数大于10就不要
	Ncut3++;
	//////////////////////////////gamVSpbar//////////////////////////////////////////

	/*
	   m_gammVspm_dang= -9;
	   m_gammVspm_dthe  = -9;
	   m_gammVspm_phi = -9;


	   double pmangd = -9;
	   double pmthed = -9;
	   double pmphid = -9;
	   double pmdthe = 200.;
	   double pmdphi = 200.;
	   double pmdang = 200.;
	   for(int i=evtRecEvent->totalCharged(); i<evtRecEvent->totalTracks(); i++)
	   {

	   EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
	   if(!(*itTrk)->isEmcShowerValid()) continue;
	   RecEmcShower *emcTrk = (*itTrk)->emcShower();
	   Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
	   for ( int j5 = 0; j5 < ichpm.size(); j5++ ){
	   EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + ichpm[j5];
	   if(!(*jtTrk)->isExtTrackValid()) continue;
	   RecExtTrack *extTrk = (*jtTrk)->extTrack();
	   if(extTrk->emcVolumeNumber() == -1) continue;//break
	   Hep3Vector extpos = extTrk->emcPosition();

	   pmangd = extpos.angle(emcpos);
	   pmthed = extpos.theta() - emcpos.theta();
	   pmphid = extpos.deltaPhi(emcpos);


	   if(pmangd < pmdang){
	   pmdang = pmangd;
	   pmdthe = pmthed;
	   pmdphi = pmphid;
	   }
	   }
	   pmdthe = pmdthe * 180 / (CLHEP::pi);
	   pmdphi = pmdphi * 180 / (CLHEP::pi);
	   pmdang = pmdang * 180 / (CLHEP::pi);

	   }
	   m_gammVspm_dang = pmdang;
	   m_gammVspm_dthe = pmdthe;
	   m_gammVspm_phi  = pmdphi;
	 */    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Vp4 pGam;
	pGam.clear();
	for(int i = 0; i < nGam; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
		RecEmcShower* emcTrk = (*itTrk)->emcShower();  //EMC束流
		double eraw = emcTrk->energy();   //能量
		double phi = emcTrk->phi();       //phi角
		double the = emcTrk->theta();     //the角
		HepLorentzVector ptrk;
		ptrk.setPx(eraw*sin(the)*cos(phi));
		ptrk.setPy(eraw*sin(the)*sin(phi));
		ptrk.setPz(eraw*cos(the));
		ptrk.setE(eraw);
		pGam.push_back(ptrk);          //根据能量和角度可以得到光子的四动量

	}


	HepLorentzVector p4c[11], ptk[7], p_pmpip[2], p_pppim[2];  //定义洛伦兹矢量

	double chisq   = 999.; 
	double chisq_1 = 999.;

	if( m_test4k4C == 1 ) {

		HepPoint3D vx(0., 0., 0.);  //ini vertex
		HepSymMatrix Evx(3, 0);
		double bx = 1E+6;
		double by = 1E+6;
		double bz = 1E+6;
		Evx[0][0] = bx*bx;
		Evx[1][1] = by*by;
		Evx[2][2] = bz*bz;
		VertexParameter vxpar;
		vxpar.setVx(vx);
		vxpar.setEvx(Evx);

		HepPoint3D pvx = xorigin;    //for every run
		HepSymMatrix pEvx(3,0);
		pEvx[0][0]=eorigin[0]*eorigin[0];
		pEvx[1][1]=eorigin[1]*eorigin[1];
		pEvx[2][2]=eorigin[2]*eorigin[2];
		VertexParameter pvxpar;
		pvxpar.setVx(pvx);
		pvxpar.setEvx(pEvx);

		//Test vertex fit
		VertexFit *vtxfit = VertexFit::instance();               //顶点拟合
		SecondVertexFit *svtxfit = SecondVertexFit::instance();  //次级顶点拟合
		double chi2VT1 = 999.9;
		WTk wtkLam1, wtkLam2;    //径迹矢量
		wtkLam1.clear();  wtkLam2.clear();

		Vint itpp, itpm, itpip, itpim;  //整形矢量
		itpp.clear();  itpm.clear(); 
		itpip.clear(); itpim.clear(); 

		Vdouble length_lam, length_lam_err; //长整形矢量：衰变长度
		length_lam.clear(); 
		length_lam_err.clear();

		Vp4 p4_pp, p4_pm, p4_pip, p4_pim;
		p4_pp.clear(); p4_pm.clear();
		p4_pip.clear(); p4_pim.clear();  //四动量

		WTrackParameter wvppTrk, wvpimTrk, wvpmTrk, wvpipTrk;     //径迹参数

		// reconstruct Lambda and anti-Lambda 

		for ( int j1 = 0; j1 < ichpp.size(); j1++ ){
			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton); // j1:重建质子径迹
			RecMdcKalTrack *ppTrk = (*(evtRecTrkCol->begin()+ichpp[j1]))->mdcKalTrack();//p+径迹
			wvppTrk = WTrackParameter(mp, ppTrk->getZHelixP(), ppTrk->getZErrorP());   //p+径迹Z轴误差
			for ( int j2 = 0; j2 < ichpim.size(); j2++ ){     
				RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);  //j2:重建pi的径迹
				RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin()+ichpim[j2]))->mdcKalTrack(); //pi-径迹
				wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());     //pi-径迹Z轴误差

				vtxfit->init();
				vtxfit->AddTrack(0,  wvppTrk);    //p+拟合径迹
				vtxfit->AddTrack(1,  wvpimTrk);   //p-拟合径迹
				vtxfit->AddVertex(0, vxpar,0,1);  //顶点拟合范围
				if( !vtxfit->Fit(0) ) continue;
				if( vtxfit->chisq(0) < 0 ) continue; //

				vtxfit->Swim(0);
				vtxfit->BuildVirtualParticle(0);    //
				WTrackParameter wALambda = vtxfit->wVirtualTrack(0); //lambda径迹参数
				VertexParameter vtxlambda = vtxfit->vpar(0);        //lambda拟合参数

				////secondary vertex fit
				svtxfit->init();
				svtxfit->AddTrack(0, wALambda);
				svtxfit->setVpar(vtxlambda);
				svtxfit->setPrimaryVertex(pvxpar);  //建立初级顶点拟合?
				if(!svtxfit->Fit()) continue;
				if(svtxfit->chisq() < 0) continue;
				if(svtxfit->chisq() < chi2VT1){
					chi2VT1 = svtxfit->chisq();
					//		wvppTrk = vtxfit->wtrk(0);
					//		wvpimTrk = vtxfit->wtrk(1);

					wtkLam1.push_back(wALambda);  //lambda
					itpp.push_back(ichpp[j1]);    //p+
					itpim.push_back(ichpim[j2]);  //pi-
					p4_pp.push_back(pchpp[j1]);  //p+动量
					p4_pim.push_back(pchpim[j2]); //pi-动量
					length_lam.push_back(svtxfit->decayLength());  //lambda衰变长度
					length_lam_err.push_back(svtxfit->decayLengthError());  //lambda衰变长度误差
				}
			}
		}

		int ncddLam1 = wtkLam1.size();
		if( ncddLam1 < 1 ) return StatusCode::SUCCESS;
		counter[0]++;         //lambda数量
		m_ngLam1 = ncddLam1;
		m_secchis1 = chi2VT1;

		//		svtxfit->ctau();

		double chi2VT2 = 999.9;
		Vdouble length_lam_bar, length_lam_bar_err;
		length_lam_bar.clear();
		length_lam_bar_err.clear();

		VertexFit *vtxfit1 = VertexFit::instance();
		SecondVertexFit *svtxfit1 = SecondVertexFit::instance();   //第二顶点拟合

		for ( int j3 = 0; j3 < ichpm.size(); j3++ ){
			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);   //j3:重建MDC质子径迹
			RecMdcKalTrack *pmTrk = (*(evtRecTrkCol->begin()+ichpm[j3]))->mdcKalTrack();  //p-径迹
			wvpmTrk = WTrackParameter(mp, pmTrk->getZHelixP(), pmTrk->getZErrorP());     //p-径迹Z轴误差

			for ( int j4 = 0; j4 < ichpip.size(); j4++ ){
				RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);   //j4:重建pi径迹
				RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin()+ichpip[j4]))->mdcKalTrack(); //pi+径迹
				wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());    //pi+径迹Z轴误差

				////primary vertex fit for lambda_bar
				vtxfit1->init();
				vtxfit1->AddTrack(0,  wvpmTrk);  //p-径迹
				vtxfit1->AddTrack(1,  wvpipTrk); //pi+径迹
				vtxfit1->AddVertex(0, vxpar,0,1); //顶点范围
				if(!vtxfit1->Fit(0)) continue;
				if(vtxfit1->chisq() < 0) continue;

				//if(!vtxfit1->Fit(0)) continue;
				vtxfit1->Swim(0);
				vtxfit1->BuildVirtualParticle(0);   
				WTrackParameter wALambda_bar = vtxfit1->wVirtualTrack(0); //FIXME  //径迹参数
				VertexParameter vtxlambda_bar = vtxfit1->vpar(0);     //拟合参数
				//		HepLorentzVector  pppm = vtxfit->pfit(0);

				svtxfit1->init();
				svtxfit1->AddTrack(0, wALambda_bar); //lambdabar径迹
				svtxfit1->setVpar(vtxlambda_bar);
				svtxfit1->setPrimaryVertex(pvxpar);  //初级顶点拟合
				if(!svtxfit1->Fit()) continue;
				if(svtxfit1->chisq() < 0) continue;
				Ncut7++;   //lambdabar拟合

				if(svtxfit1->chisq() < chi2VT2){
					chi2VT2 = svtxfit1->chisq();
					//	wvpmTrk = vtxfit1->wtrk(0);
					//	wvpipTrk = vtxfit1->wtrk(1);

					wtkLam2.push_back(wALambda_bar);
					itpm.push_back(ichpm[j3]);    //p-
					itpip.push_back(ichpip[j4]);  //pi+
					p4_pm.push_back(pchpm[j3]);   //p-四动量
					p4_pip.push_back(pchpip[j4]); //pi+四动量
					length_lam_bar.push_back(svtxfit1->decayLength());  //lambdabar衰变长度
					length_lam_bar_err.push_back(svtxfit1->decayLengthError());   //lambdabar衰变长度误差
				}
			}
		}

		int ncddLam2 = wtkLam2.size();
		if( ncddLam2 < 1 ) return StatusCode::SUCCESS;
		counter[1]++;   //lam数量
		m_ngLam2 = ncddLam2;

		m_secchis2 = chi2VT2;

		Ncut8++;
		////for J/psi-->Lambda Lambda_bar
		KalmanKinematicFit * kmfit_0g = KalmanKinematicFit::instance();   //卡尔曼滤波拟合
		for ( int n5 = 0; n5 < wtkLam1.size(); n5++ ){
			for ( int n6 = 0; n6 < wtkLam2.size(); n6++ ){
				kmfit_0g->init();
				kmfit_0g->AddTrack(0, wtkLam1[n5]);//K  //lam1径迹
				kmfit_0g->AddTrack(1, wtkLam2[n6]);//K  //lam2径迹
				kmfit_0g->AddFourMomentum(0, ecms);     //四动量
				bool oksq_0g = kmfit_0g->Fit();        
				if (oksq_0g) {
					double chi22 = kmfit_0g->chisq();
					m_chi1_0g = chi22;
					for (int ipull = 0; ipull < 2; ipull++){
						p4c[ipull]= kmfit_0g->pfit(ipull);   //动量拟合
					}
				}
			}
		}

		if ( nGam > 0 ){
			KalmanKinematicFit *kmfit_1g = KalmanKinematicFit::instance();
			for(int ni = 0;ni < nGam-0;ni++) {
				RecEmcShower *g1trk = (*(evtRecTrkCol->begin()+iGam[ni]))->emcShower();
				for ( int n1 = 0; n1 < wtkLam1.size(); n1++ ){
					for ( int n2 = 0; n2 < wtkLam2.size(); n2++ ){
						kmfit_1g->init();
						kmfit_1g->AddTrack(0, wtkLam1[n1]); //lam1径迹拟合
						kmfit_1g->AddTrack(1, wtkLam2[n2]); //lam2径迹拟合
						kmfit_1g->AddTrack(2, 0.0, g1trk);  //g1trk径迹拟合
						kmfit_1g->AddFourMomentum(0, ecms); //四动量
						bool oksq = kmfit_1g->Fit();
						if(oksq) {
							double chi2 = kmfit_1g->chisq();
							if(chi2 < chisq) {
								chisq = chi2;
							}//for chi2  
						}///ok
					}//n2
				}//n1 
			}     //for ni nGam
			m_chi1 = chisq;
		}

		////for J/psi-->Lambda Lambda_bar

		if ( m_chi1_0g < 200 ) {
			counter[2]++;
			ptk[0] = p4_pp[0];     //p+
			ptk[0].setE(sqrt(p4_pp[0].v().mag2()+mp*mp));

			ptk[1] = p4_pim[0];    //pi-
			ptk[1].setE(sqrt(p4_pim[0].v().mag2()+mpi*mpi));

			ptk[2] = p4_pm[0];     //p-
			ptk[2].setE(sqrt(p4_pm[0].v().mag2()+mp*mp));

			ptk[3] = p4_pip[0];    //pi+
			ptk[3].setE(sqrt(p4_pip[0].v().mag2()+mpi*mpi));

			m_phi1_bk1 = ( ptk[0] + ptk[1] ).mag(); //p+ + pi-
			m_phi2_bk1 = ( ptk[2] + ptk[3] ).mag();//p- + pi+

			m_phi1_bk = ( ptk[3] + ptk[1] ).mag();  //pi+ + pi-
			m_phi2_bk = ( ptk[2] + ptk[0] ).mag();  //p+ + p-

			m_2phi_bk = ( ptk[0] + ptk[1] + ptk[2] + ptk[3] ).mag(); //p+ + pi- + p- + pi+

			m_mlengthlam = length_lam[0];               //lam衰变长度
			m_mlengthlamb = length_lam_bar[0];         //lamdarbar衰变长度
			m_mlengthlam_err = length_lam_err[0];      //lamdar衰变长度误差
			m_mlengthlamb_err = length_lam_bar_err[0]; //lamdarbar衰变长度误差

			////// correct the information of p pi- p-bar pi+

			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);  //重建MDC质子径迹
			RecMdcKalTrack *lamb_pp = (*(evtRecTrkCol->begin()+itpp[0]))->mdcKalTrack();  //p+径迹
			WTrackParameter wvlamb_pp = WTrackParameter(mp, lamb_pp->getZHelixP(), lamb_pp->getZErrorP());//p+径迹Z轴误差

			RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);  //pi径迹
			RecMdcKalTrack *lamb_pim = (*(evtRecTrkCol->begin()+itpim[0]))->mdcKalTrack();  //pi-径迹
			WTrackParameter wvlamb_pim = WTrackParameter(mpi, lamb_pim->getZHelix(), lamb_pim->getZError());//pi-径迹误差

			HepLorentzVector p4lamb;
			p4lamb.setPx(p4c[0].px()); 
			p4lamb.setPy(p4c[0].py());                
			p4lamb.setPz(p4c[0].pz());                
			p4lamb.setE(p4c[0].e());    //四个洛伦兹矢量

			KalmanKinematicFit *kmfit_lamb = KalmanKinematicFit::instance();  //卡尔曼滤波拟合
			kmfit_lamb->init();
			kmfit_lamb->AddTrack(0, wvlamb_pp);   //p+
			kmfit_lamb->AddTrack(1, wvlamb_pim);  //pi-
			kmfit_lamb->AddFourMomentum(0, p4lamb); //四动量
			chi2_lamb = kmfit_lamb->chisq();
			if ( chi2_lamb < 200 ) Ncut11++; 
			if ( !kmfit_lamb->Fit() ) return StatusCode::SUCCESS;
			counter[3]++;            

			bool oksq_lamb = kmfit_lamb->Fit();
			for (int it1 = 0; it1 < 2; it1++){
				p_pppim[it1]= kmfit_lamb->pfit(it1);
			}

			RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);  //质子
			RecMdcKalTrack *lamb_pm = (*(evtRecTrkCol->begin()+itpm[0]))->mdcKalTrack();  //p-
			WTrackParameter wvlamb_pm = WTrackParameter(mp, lamb_pm->getZHelixP(), lamb_pm->getZErrorP()); //误差

			RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);  //pi
			RecMdcKalTrack *lamb_pip = (*(evtRecTrkCol->begin()+itpip[0]))->mdcKalTrack();  //pi+
			WTrackParameter wvlamb_pip = WTrackParameter(mpi, lamb_pip->getZHelix(), lamb_pip->getZError()); //误差

			KalmanKinematicFit *kmfit_lambbar = KalmanKinematicFit::instance();  //卡尔曼滤波拟合
			HepLorentzVector p4lambbar;  
			p4lambbar.setPx(p4c[1].px());
			p4lambbar.setPy(p4c[1].py());
			p4lambbar.setPz(p4c[1].pz());
			p4lambbar.setE(p4c[1].e());  //lambbar四动量

			kmfit_lambbar->init();
			kmfit_lambbar->AddTrack(0, wvlamb_pm);  //p-
			kmfit_lambbar->AddTrack(1, wvlamb_pip); //pi+
			kmfit_lambbar->AddFourMomentum(0, p4lambbar);  //lambbar四动量
			chi2_lambar = kmfit_lambbar->chisq();
			if ( !kmfit_lambbar->Fit() ) return StatusCode::SUCCESS;
			counter[4]++;  //p-、pi+顶点拟合

			bool oksq_lambbar = kmfit_lambbar->Fit();
			for (int it2 = 0; it2 < 2; it2++){
				p_pmpip[it2]= kmfit_lambbar->pfit(it2);
			}

			m_px_mdc[0] = p_pppim[0].px();
			m_py_mdc[0] = p_pppim[0].py();
			m_pz_mdc[0] = p_pppim[0].pz();
			m_en_mdc[0] = p_pppim[0].e();
			m_p3_mdc[0] = p_pppim[0].v().mag();  //p

			m_px_mdc[1] = p_pppim[1].px();
			m_py_mdc[1] = p_pppim[1].py();
			m_pz_mdc[1] = p_pppim[1].pz();
			m_en_mdc[1] = p_pppim[1].e();
			m_p3_mdc[1] = p_pppim[1].v().mag();  //pi-

			m_px_mdc[2] = p_pmpip[0].px();
			m_py_mdc[2] = p_pmpip[0].py();
			m_pz_mdc[2] = p_pmpip[0].pz();
			m_en_mdc[2] = p_pmpip[0].e();
			m_p3_mdc[2] = p_pmpip[0].v().mag();  //p-

			m_px_mdc[3] = p_pmpip[1].px();
			m_py_mdc[3] = p_pmpip[1].py();
			m_pz_mdc[3] = p_pmpip[1].pz();
			m_en_mdc[3] = p_pmpip[1].e();
			m_p3_mdc[3] = p_pmpip[1].v().mag();  //pi+

			for (int itk=0; itk<2; itk++){
				m_px_trk[itk] = p4c[itk].px();
				m_py_trk[itk] = p4c[itk].py();
				m_pz_trk[itk] = p4c[itk].pz();
				m_en_trk[itk] = p4c[itk].e();
				m_p3_trk[itk] = p4c[itk].v().mag();  //四动量
			}

			HepLorentzVector p4phi1,p4phi2;
			p4phi1 = p4c[0]; p4phi2 = p4c[1];

			m_phi1_true = p4phi1.mag();
			m_phi2_true = p4phi2.mag();
			m_etac_true = (p4phi1+p4phi2).mag();
		}
	}// 4c

	counter[5]++;
	m_tuple6->write();
	return StatusCode::SUCCESS;
}
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&(0_0)
StatusCode Jftrk::finalize() {
	cout<<"total number        : "<<Ncut0<<endl;
	cout<<"nGood >= 4          : "<<Ncut1<<endl;
	cout<<"four track for lamb : "<<Ncut2<<endl;
	cout<<"nGam > 1 < 10       : "<<Ncut3<<endl;
	cout<<"p+ pi- vertex fit   : "<<Ncut4<<endl;
	cout<<"lambda fit          : "<<Ncut5<<endl;
	cout<<"p- pi+ vertex fit   : "<<Ncut6<<endl;
	cout<<"lambda_bar fit      : "<<Ncut7<<endl;
	cout<<"lambda and lambda_bar fit     : "<<Ncut8<<endl;
	cout<<"Nlam >= 1           : "<<counter[0]<<endl;
	cout<<"Nlambar >= 1        : "<<counter[1]<<endl;
	cout<<"pass 4c 0g          : "<<counter[2]<<endl;
	cout<<"pass 4c 1g          : "<<Ncut9<<endl;
	cout<<"pass 4c 2g          : "<<Ncut10<<endl;
	cout<<"p pi- to lamb       : "<<counter[3]<<endl;
	cout<<"pbar pi+ to lambar  : "<<counter[4]<<endl;
	cout<<"fill tree           : "<<counter[5]<<endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}

void Jftrk::InitVar(){
	m_chi1 = 9999.0; m_chi1_0g = 9999.0;
	m_secchis1 = 9999.0; m_secchis2 = 9999.0;
	for ( int i2 = 0 ; i2 < 2; i2++){
		m_px_trk[i2] = -10; m_py_trk[i2] = -10;
		m_pz_trk[i2] = -10; m_en_trk[i2] = -10; m_p3_trk[i2] = -10;
	}
	for ( int i3 = 0 ; i3 < 4; i3++){
		m_px_mdc[i3] = -10; m_py_mdc[i3] = -10;
		m_pz_mdc[i3] = -10; m_en_mdc[i3] = -10; m_p3_mdc[i3] = -10;
	}
}
