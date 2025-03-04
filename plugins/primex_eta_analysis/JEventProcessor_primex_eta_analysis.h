// $Id$
//
//    File: JEventProcessor_primex_eta_analysis.h
// Created: Tue Mar  4 10:14:54 AM EST 2025
// Creator: andrsmit (on Linux ifarm2401.jlab.org 5.14.0-503.19.1.el9_5.x86_64 x86_64)
//

/// For more information on the syntax changes between JANA1 and JANA2, visit: https://jeffersonlab.github.io/JANA2/#/jana1to2/jana1-to-jana2

#ifndef _JEventProcessor_primex_eta_analysis_
#define _JEventProcessor_primex_eta_analysis_

#include <JANA/JEventProcessor.h>
#include <JANA/Services/JLockService.h> // Required for accessing services

// Hall-D headers:
#include "DANA/DEvent.h"
#include "HDGEOMETRY/DGeometry.h"
#include "TRACKING/DMCThrown.h"
#include "TRIGGER/DTrigger.h"
#include "TRIGGER/DL1Trigger.h"
#include "FCAL/DFCALShower.h"
#include "BCAL/DBCALShower.h"
#include "CCAL/DCCALShower.h"
#include "TOF/DTOFPoint.h"
#include "START_COUNTER/DSCHit.h"
#include "PID/DMCReaction.h"
#include "PID/DBeamPhoton.h"
#include "PID/DEventRFBunch.h"
#include "DVector3.h"
#include "DLorentzVector.h"
#include "particleType.h"

// ROOT headers:
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"

class JEventProcessor_primex_eta_analysis:public JEventProcessor{
    public:
        JEventProcessor_primex_eta_analysis() { 
			SetTypeName(NAME_OF_THIS);
		}
        ~JEventProcessor_primex_eta_analysis() = default;
		
        const char* className(void){return "JEventProcessor_primex_eta_analysis";}

    private:
        void Init() override;
        void BeginRun(const std::shared_ptr<const JEvent>& event) override;
        void Process(const std::shared_ptr<const JEvent>& event) override;
        void EndRun() override {};
        void Finish() override {};

    	std::shared_ptr<JLockService> lockService; //Used to access all the services, its value should be set inside Init()
		
		//---------------------------------------//
		// Functions
		
		void InitializeHistograms();
		
		Particle_t GetTargetType(int32_t);
		int GetPrimExPhase(int32_t);
		
		int FCALFiducialCut(DVector3, double);
		
		double GetEnergyAfterRecoil(double, double, double, double);
		
		double GetFCALEnergyRes(double);
		
		double GetAccScalingFactor(double);
		
		bool IsElasticCut(double, double, double);
		bool IsEtaCut(double);
		bool IsCoplanarBCAL(double);
		bool IsCoplanarSC(double);
		
		void CheckTOFMatch(DVector3, double, vector<const DTOFPoint*>, double&, double&, double&, double);
		
		void EtaGGAnalysis(
			vector<const DFCALShower*>, vector<const DBeamPhoton*>, 
			vector<const DBCALShower*>, vector<const DTOFPoint*>, vector<const DSCHit*>, 
			int, int, double, int, int, double, double, double, bool, double, double
		);
		
		//---------------------------------------//
		// Geometry
		
		DVector3 m_beamSpot;
		DVector3 m_fcalFace, m_ccalFace;
		DVector3 m_fcalCorrection, m_ccalCorrection;
		
		vector<vector<DVector3>> m_scPos, m_scNorm;
		
		int m_phaseVal = 0;
		
		double m_HodoscopeHiFactor    = 1.0;
		double m_HodoscopeHiFactorErr = 1.0;
		double m_HodoscopeLoFactor    = 1.0;
		double m_HodoscopeLoFactorErr = 1.0;
		double m_MicroscopeFactor     = 1.0;
		double m_MicroscopeFactorErr  = 1.0;
		double m_TAGMEnergyBoundHi    = 1.0;
		double m_TAGMEnergyBoundLo    = 1.0;
		
		//---------------------------------------//
		// Cuts (defaults set inside Init)
		
		double m_BeamRFCut, m_FCALRFCut, m_BCALRFCut, m_CCALRFCut, m_TOFRFCut;
		
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
		
		int m_UseLogWeight, m_BypassTrigger;
		
		//---------------------------------------//
		// Constants 
		
		Particle_t m_Target;
		
		const double m_c = 29.9792458;   // [cm/ns]
		
		double m_beamBunchesMain = 1.0;
		double m_beamBunchesAcc  = 5.0;
		
		//---------------------------------------//
		// Histograms
		
		double m_minInvmassBin      =  0.300;
		double m_maxInvmassBin      =  1.100;
		double m_invmassBinSize     =  0.001;
		
		double m_minRecAngleBin     =  0.000;
		double m_maxRecAngleBin     =  5.500;
		double m_recAngleBinSize    =  0.010;
		
		double m_minThrownAngleBin  =  0.000;
		double m_maxThrownAngleBin  =  5.000;
		double m_thrownAngleBinSize =  0.010;
		
		// Thrown angular distributions for different energy ranges:
		
		TH1F *h_ThrownAngle[13];
		
		// Timing and energy sum for different trigger types:
		
		static const int m_nTrigs = 4;
		vector<string> m_triggerNames = {"FCAL+CCAL Trigger", 
			"Low-energy FCAL Trigger", "PS Trigger", "CCAL Trigger"};
		
		TH1F *h_fcalRFdt[m_nTrigs];
		TH1F *h_bcalRFdt[m_nTrigs];
		TH1F *h_ccalRFdt[m_nTrigs];
		TH1F *h_taghRFdt[m_nTrigs];
		TH1F *h_tagmRFdt[m_nTrigs];
		TH1F  *h_tofRFdt[m_nTrigs];
		TH1F   *h_scRFdt[m_nTrigs];
		
		TH1F *h_fcalEnergySum[m_nTrigs];
		
		// FCAL-TOF Matching:
		
		TH1F *h_fcalTOFdx, *h_fcalTOFdy, *h_fcalTOFdr;
		TH1F *h_fcalTOFdt, *h_fcalTOFdt_cut;
		TH1F *h_fcalTOFMatches;
		
		// RF timing after cuts are applies:
		
		TH1F *h_beamRFdt_cut, *h_scRFdt_cut;
		
		// Elasticity:
		
		TH2F *h_elas;
		TH2F *h_elasCorr, *h_elasCorrMain, *h_elasCorrSide;
		TH2F *h_elasCorr_coh;
		
		// Invariant mass:
		
		TH2F *h_mgg,       *h_mggMain,       *h_mggSide;
		TH2F *h_mggConstr, *h_mggConstrMain, *h_mggConstrSide;
		TH2F *h_mggConstr_coh;
		
		// Missing mass:
		
		TH2F *h_mm,     *h_mm_etaCut,     *h_mm_elasEtaCut;
		TH2F *h_mm_coh, *h_mm_coh_etaCut, *h_mm_coh_elasEtaCut;
		
		// Misc:
		
		TH2F *h_elasVSmgg;
		TH2F *h_xy_1, *h_xy_2;
		
		// For MC:
		
		TH2F *h_recVSthrown;
		TH2F *h_mggThrown, *h_mggConstrThrown;
		
		// Plots for different veto options:
		 
		static const int m_nVetos = 8;
		
		TH2F *h_elas_veto[m_nVetos];
		TH2F *h_mgg_veto[m_nVetos], *h_mggConstr_veto[m_nVetos], *h_mggConstr_coh_veto[m_nVetos];
		TH2F *h_mm_veto[m_nVetos],  *h_mm_coh_veto[m_nVetos];
};

#endif // _JEventProcessor_primex_eta_analysis_

