// $Id$
//
//    File: JEventProcessor_eta_gg_tree.h
// Created: Tue Mar  4 10:14:42 AM EST 2025
// Creator: andrsmit (on Linux ifarm2401.jlab.org 5.14.0-503.19.1.el9_5.x86_64 x86_64)
//

/// For more information on the syntax changes between JANA1 and JANA2, visit: https://jeffersonlab.github.io/JANA2/#/jana1to2/jana1-to-jana2

#ifndef _JEventProcessor_eta_gg_tree_
#define _JEventProcessor_eta_gg_tree_

#include <JANA/JEventProcessor.h>
// #include <JANA/Services/JLockService.h> // Required for accessing services

// Hall-D headers:
#include "DANA/DEvent.h"
#include "HDGEOMETRY/DGeometry.h"
#include "TRACKING/DMCThrown.h"
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
#include "include/particleType.h"
#include "ANALYSIS/DTreeInterface.h"

// ROOT headers:
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TTree.h"

class JEventProcessor_eta_gg_tree : public JEventProcessor {
    public:
        JEventProcessor_eta_gg_tree() {
			SetTypeName(NAME_OF_THIS);
		}
        ~JEventProcessor_eta_gg_tree() = default;
		
        const char* className(void){return "JEventProcessor_eta_gg_tree";}

    private:
        void Init() override;
        void BeginRun(const std::shared_ptr<const JEvent>& event) override;
        void Process(const std::shared_ptr<const JEvent>& event) override;
        void EndRun() override;
        void Finish() override;

    	// std::shared_ptr<JLockService> lockService; //Used to access all the services, its value should be set inside Init()
		
		//---------------------------------------//
		// Functions
		
		Particle_t GetTargetType(int32_t);
		int GetPrimExPhase(int32_t);
		
		double GetEnergyAfterRecoil(double, double, double, double);
		double GetFCALEnergyRes(double);
		double GetAccScalingFactor(double);
		
		void WriteEvent(uint64_t eventnumber, double rfTime, vector<const DMCThrown*> mc_thrown, const DMCReaction* mc_reaction);
		
		void WriteEvent(uint64_t eventnumber, double rfTime,
			vector<const DBeamPhoton*> beam_photons, 
			vector<const DFCALShower*> fcal_showers,
			vector<const DBCALShower*> bcal_showers,
			vector<const DTOFPoint*> tof_points,
			vector<const DSCHit*> sc_hits,
			vector<const DMCThrown*> mc_thrown, const DMCReaction* mc_reaction);
		
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
		
		int m_USE_LOG_WEIGHT, m_SAVE_MC_NOHITS;
		
		double m_FCAL_RF_CUT, m_BCAL_RF_CUT, m_TOF_RF_CUT, m_BEAM_RF_CUT;
		
		double m_MIN_BEAM_ENERGY;
		double m_DELTA_E_CUT;
		
		//---------------------------------------//
		// Constants 
		
		Particle_t m_Target;
		
		const double m_c = 29.9792458;   // [cm/ns]
		
		double m_beamBunchesMain = 1.0;
		double m_beamBunchesAcc  = 5.0;
		
		//---------------------------------------//
		// Tree
		
		DTreeInterface *dTreeInterface;
		static thread_local DTreeFillData dTreeFillData;
};

#endif // _JEventProcessor_eta_gg_tree_

