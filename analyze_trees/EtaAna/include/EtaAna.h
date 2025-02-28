#ifndef _ETAGG_ANA_
#define _ETAGG_ANA_

#define MAX_BEAM 400
#define MAX_FCAL 100
#define MAX_BCAL 200
#define MAX_TOF  100
#define MAX_SC   100
#define MAX_MC   5

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRandom3.h"
#include "particleType.h"

class EtaAna {
	private:
		
		int m_runNumber;
		TRandom3 *m_random;
		
		// Geometry:
		
		int m_phaseVal;
		Particle_t m_Target;
		double m_targetLength, m_targetDensity, m_targetAtten;
		
		TVector3 m_fcalFace;
		TVector3 m_fcalCorrection;
		TVector3 m_vertex;
		
		// Constants:
		
		static constexpr double m_c = TMath::C() * 1.e-7; // [cm/ns]
		static constexpr double m_e = 0.510998928e-3;     // [GeV]
		static constexpr double m_fcalBlockSize = 4.0157; // [cm]
		
		// Variables:
		
		string m_outputFileName;
		
		TFile *m_inputFile;
		TTree *m_tree;
		
		// decides whether to create and fill thrown histograms:
		bool m_FillThrown = false;
		
		int m_event;
		
		vector<vector<Particle_t>> m_reaction_types;
		
		// Common Functions:
		
		int SetGeometry();
		
		int LoadTree();
		int CheckEventMultiplicities();
		int AcceptRejectEvent();
		
		void ReadEvent();
		
		void PlotThrown(double, double);
		void FillAngularMatrix(double thrownEnergy, double thrownAngle, 
			double recAngle, double weight);
		void FillAngularMatrix_vetos(int vetoOption, double thrownEnergy, double thrownAngle, 
			double recAngle, double weight);
		void FillInvmassMatrix(double theta, double mgg, double beamEnergy, double weight);
		
		int GetFCALShowerList(vector<int> &goodShowers, int &nFCALShowers_EnergyCut, 
			double energyCut, double extraEnergyCut, double fiducialCut, double timingCut);
		int GetBCALShowerList(vector<int> &goodShowers, double energyCut, double timingCut);
		int GetSCHitList(vector<int> &goodHits);
		int GetBeamPhotonList(vector<pair<int,double>> &goodPhotons, double minEnergyCut, double maxEnergyCut);
		
		TVector3 GetFCALPosition(int index);
		TVector3 GetBCALPosition(int index);
		
		double CalcInvmass(double e1, double e2, TVector3 pos1, TVector3 pos2);
		double CalcInvmass(double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3);
		double CalcInvmassShifted(double deltaZ, double e1, double e2, TVector3 pos1, TVector3 pos2);
		double CalcInvmassShifted(double deltaZ, double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3);
		double CalcProdTheta(double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3);
		double CalcProdThetaShifted(double deltaZ, double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3);
		double CalculateCost(double deltaZ, double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3);
		
		int FCALFiducialCut(TVector3, double);
		
		void CheckTOFMatch(TVector3, double&, double&, double&, double);
		
		double GetEnergyAfterRecoil(double, double, double, double);
		double GetFCALEnergyResolution(double);
		void GetThrownEnergyAndAngle(double&, double&);
		void GetThrownEnergyAndAngleBGGEN(double&, double&);
		
		void SmearShowerEnergy(double&);
		
		int GetAcceptanceHistogram();
		
		bool IsElasticCut(double, double, double);
		bool IsEtaCut(double);
		bool IsCoplanarBCAL(double);
		bool IsCoplanarSC(double);
		
		//-------------------------------------------------------------------//
		
		// Different ways to do analysis (each function is defined in it's own .cc file):
		void EtaggAnalysis();
		void EtaggAnalysis_matrix();
		void EtaggAnalysis_FCAL();
		void EtaggAnalysis_TOF();
		void EtaggAnalysis_bggen();
		void Omega3gAnalysis();
		
		// functions for default analysis option:
		void InitializeDefaultHists();
		void      ResetDefaultHists();
		void      WriteDefaultHists();
		
		// functions for matrix analysis option:
		void InitializeMatrixHists();
		void      ResetMatrixHists();
		void      WriteMatrixHists();
		
		// functions for FCAL analysis option:
		void InitializeFCALHists();
		void      ResetFCALHists();
		void      WriteFCALHists();
		
		// functions for TOF analysis option:
		void InitializeTOFHists();
		void      ResetTOFHists();
		void      WriteTOFHists();
		
		// functions for analyzing bggen mc:
		void InitializeReactionTypes();
		void InitializeBGGENHists();
		void      ResetBGGENHists();
		void      WriteBGGENHists();
		int  GetFinalState_bggen(int debug=0);
		
		// functions for omega->pi0+gamma analysis option:
		void InitializeOmegaHists();
		void      ResetOmegaHists();
		void      WriteOmegaHists();
		
		//-------------------------------------------------------------------//
		// TTree Variables:
		
		int    m_eventNum;
		double m_rfTime;
		int    m_nbeam;
		int    m_tagSystem[MAX_BEAM];
		int    m_tagCounter[MAX_BEAM];
		double m_beamE[MAX_BEAM];
		double m_beamT[MAX_BEAM];
		double m_accScaleFactor[MAX_BEAM];
		
		int    m_nfcal;
		double m_fcalE[MAX_FCAL];
		double m_fcalX[MAX_FCAL];
		double m_fcalY[MAX_FCAL];
		double m_fcalZ[MAX_FCAL];
		double m_fcalT[MAX_FCAL];
		int    m_fcalNblocks[MAX_FCAL];
		double m_fcalE1E9[MAX_FCAL];
		double m_fcalE9E25[MAX_FCAL];
		
		int    m_nbcal;
		double m_bcalE[MAX_BCAL];
		double m_bcalX[MAX_BCAL];
		double m_bcalY[MAX_BCAL];
		double m_bcalZ[MAX_BCAL];
		double m_bcalT[MAX_BCAL];
		
		int    m_ntof;
		double m_tofX[MAX_TOF];
		double m_tofY[MAX_TOF];
		double m_tofZ[MAX_TOF];
		double m_tofT[MAX_TOF];
		
		int    m_nsc;
		double m_scSector[MAX_SC];
		double m_scPhi[MAX_SC];
		double m_scdE[MAX_SC];
		double m_scT[MAX_SC];
		double m_scPulseHeight[MAX_SC];
		
		int    m_nmc;
		int    m_mcPDGType[MAX_MC];
		double m_mcX[MAX_MC];
		double m_mcY[MAX_MC];
		double m_mcZ[MAX_MC];
		double m_mcT[MAX_MC];
		double m_mcE[MAX_MC];
		double m_mcP[MAX_MC];
		double m_mcTheta[MAX_MC];
		double m_mcPhi[MAX_MC];
		
		double m_thrownBeamEnergy;
		
		//-------------------------------------------------------------------//
		// Histograms:
		
		// Defaults for histogram bin sizes:
		
		double m_minInvmassBin      =  0.300;
		double m_maxInvmassBin      =  1.100;
		double m_invmassBinSize     =  0.001;
		
		double m_minRecAngleBin     =  0.000;
		double m_maxRecAngleBin     =  5.500;
		double m_recAngleBinSize    =  0.010;
		
		double m_minThrownAngleBin  =  0.000;
		double m_maxThrownAngleBin  =  5.000;
		double m_thrownAngleBinSize =  0.010;
		
		double m_minBeamEnergyBin   =  7.000;
		double m_maxBeamEnergyBin   = 12.000;
		double m_beamEnergyBinSize  =  0.050;
		
		TH2F *h_acceptance; // acceptance as 2-d grid of theta vs. beam energy
		
		/////////////////////////////////////////
		// Defulat Analysis:
		
		TH1F *h_mcVertex, *h_mcVertexAccepted;
		TH2F *h_thrown;
		
		TH1F *h_fcalRFdt, *h_bcalRFdt, *h_tofRFdt, *h_scRFdt;
		TH1F *h_beamRFdt, *h_beamRFdt_cut;
		
		TH2F *h_tofMatch;
		
		static const int m_nVetos = 9;
		TH2F *h_elasticity[m_nVetos];
		TH2F *h_elasticityConstr[m_nVetos];
		TH2F *h_mgg[m_nVetos];
		TH2F *h_mggConstr[m_nVetos];
		TH2F *h_mggConstr_coh[m_nVetos];
		TH2F *h_mm[m_nVetos],      *h_mm_coh[m_nVetos];
		TH2F *h_mm_elas[m_nVetos], *h_mm_elas_coh[m_nVetos];
		TH2F *h_pt[m_nVetos],      *h_ptCoh[m_nVetos];
		
		TH2F *h_t_vs_theta;
		
		TH2F *h_mgg_vs_vertexZ, *h_mgg_vs_vertexR;
		
		TH2F *h_scDeltaPhi,           *h_bcalDeltaPhi;
		TH2F *h_scDeltaPhi_singleHit, *h_bcalDeltaPhi_singleHit;
		
		TH1F *h_nSC, *h_nSC_nobcal, *h_nSC_onebcal, *h_nSC_bcalVeto, *h_nSC_extra;
		
		TH2F *h_xy1, *h_xy2;
		
		vector<TH3F*> h_AngularMatrix_vetos;
		
		TH2F *h_ecomp1,       *h_ecomp2,       *h_ecomp;
		TH2F *h_ecomp1_clean, *h_ecomp2_clean, *h_ecomp_clean;
		
		/////////////////////////////////////////
		// Default Analysis with no cut on beam energy:
		
		TH3F *h_invmassMatrix;
		TH3F *h_invmassMatrix_prompt, *h_invmassMatrix_acc;
		TH3F *h_AngularMatrix;
		
		/////////////////////////////////////////
		// Varying FCAL Cuts:
		
		TH2F *h_mgg_FCAL;
		TH2F *h_mgg_FCALECut;
		TH2F *h_mgg_FCALFidCut;
		TH2F *h_mgg_FCALCuts;
		TH2F *h_mgg_FCALGoodMult;
		TH2F *h_mgg_FCALMult;
		
		vector<double> m_fcalEnergyCuts;
		vector<TH2F*> h_mgg_FCALECutVec,           h_mgg_FCALExtraECutVec;
		vector<TH3F*> h_AngularMatrix_FCALECutVec, h_AngularMatrix_FCALExtraECutVec;
		
		vector<double> m_fcalFiducialCuts;
		vector<TH2F*> h_mgg_FCALFidCutVec;
		vector<TH3F*> h_AngularMatrix_FCALFidCutVec;
		
		/////////////////////////////////////////
		// Varying TOF Cuts:
		
		vector<double> m_TOFTimingCuts;
		vector<TH2F*> h_mgg_TOFTimingCutVec;
		vector<TH3F*> h_AngularMatrix_TOFTimingCutVec;
		
		vector<double> m_TOFDistanceCuts;
		vector<TH2F*> h_mgg_TOFDistanceCutVec;
		vector<TH3F*> h_AngularMatrix_TOFDistanceCutVec;
		
		TH2F *h_mgg_noTOF, *h_mgg_TOF, *h_mgg_singleTOF;
		TH3F *h_AngularMatrix_noTOF, *h_AngularMatrix_TOF, *h_AngularMatrix_singleTOF;
		
		/////////////////////////////////////////
		// Varying Beam Cuts:
		
		/////////////////////////////////////////
		// BGGEN Analysis:
		
		TH2F *h_mgg_const_bggen_signal,   *h_elas_bggen_signal,   *h_elas_bggen_signal_cut;
		TH2F *h_mgg_const_bggen_etapion,  *h_elas_bggen_etapion,  *h_elas_bggen_etapion_cut;
		TH2F *h_mgg_const_bggen_eta2pion, *h_elas_bggen_eta2pion, *h_elas_bggen_eta2pion_cut;
		TH2F *h_mgg_const_bggen_omega,    *h_elas_bggen_omega,    *h_elas_bggen_omega_cut;
		TH2F *h_mgg_const_bggen_rho,      *h_elas_bggen_rho,      *h_elas_bggen_rho_cut;
		TH2F *h_mgg_const_bggen_bkgd,     *h_elas_bggen_bkgd,     *h_elas_bggen_bkgd_cut;
		
		TH1F *h_thrown_reactions_bggen;
		TH1F *h_rec_reactions_bggen;
		
		vector<TH2F*> h_bcalDeltaPhi_bggen, h_scDeltaPhi_bggen, h_scDeltaPhi_singleHit_bggen;
		
		TH2F *h_thrown_proton,        *h_thrown_neutron;
		TH3F *h_AngularMatrix_proton, *h_AngularMatrix_neutron;
		
		TH1F *h_thrown_proton_1d,  *h_rec_proton_1d;
		TH1F *h_thrown_neutron_1d, *h_rec_neutron_1d;
		
		/////////////////////////////////////////
		// Omega->pi0+gamma:
		
		TH1F *h_3gamma_m12,      *h_3gamma_m13,      *h_3gamma_m23,      *h_3gamma_m3g;
		TH1F *h_3gamma_m12_elas, *h_3gamma_m13_elas, *h_3gamma_m23_elas, *h_3gamma_m3g_elas;
		TH1F *h_3gamma_vz,       *h_3gamma_vz_elas;
		TH1F *h_3gamma_theta_targ, *h_3gamma_theta_fdc1, *h_3gamma_theta_fdc2, *h_3gamma_theta_fdc3;
		TH2F *h_xy_targ, *h_xy_fdc1, *h_xy_fdc2, *h_xy_fdc3;
		//----------------------------------------------------------------------------------------------//
		
	public:
		//-------------------------------------------------------------------//
		// Cuts:
		
		double m_FCALRFCut;
		double m_BCALRFCut;
		double m_BeamRFCut;
		double m_TOFRFCut;
		double m_FCALEnergyCut, m_FCALExtraEnergyCut;
		double m_BCALEnergyCut;
		double m_minBeamEnergyCut, m_maxBeamEnergyCut;
		double m_FCALTOFCut;
		double m_SCDeltaPhiCut;
		double m_BCALDeltaPhiCut;
		double m_ElasMean_p0;
		double m_ElasMean_p1;
		double m_ElasWidth;
		double m_ElasSigmaCut;
		
		double beamBunchesMain = 1.0;
		double beamBunchesAcc  = 5.0; // how many bunches to use on each side for accidental subtraction
		
		//-------------------------------------------------------------------//
		// Constructor/Destructor:
		
		EtaAna();
		~EtaAna(){};
		
		//-------------------------------------------------------------------//
		// Function Declarations:
		
		int  GetPrimexPhase(int runNumber);
		
		int  SetRunNumber(int runNumber);
		void SetOutputFileName(string fileName) { m_outputFileName = fileName; return; };
		
		int  SetCuts(TString configFileName);
		void DumpCuts();
		
		void RunAnalysis(TString inputFileName, int analysisOption=0);
		void  InitHistograms(int analysisOptions=0);
		void ResetHistograms(int analysisOptions=0);
		void WriteHistograms(int analysisOptions=0);
		
		void SetFillThrown(bool doFill) { m_FillThrown = doFill; return; };
};

#endif
