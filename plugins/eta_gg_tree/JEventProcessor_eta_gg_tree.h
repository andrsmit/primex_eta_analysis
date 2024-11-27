// $Id$
//
//    File: JEventProcessor_eta_gg_tree.h
// Created: Fri Aug 11 14:26:44 EDT 2023
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.92.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_eta_gg_tree_
#define _JEventProcessor_eta_gg_tree_

// JANA headers:
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>

// Hall-D headers:
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
#include "ANALYSIS/DTreeInterface.h"

// ROOT headers:
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TTree.h"

using namespace jana;
using namespace std;

class JEventProcessor_eta_gg_tree:public jana::JEventProcessor{
	public:
		JEventProcessor_eta_gg_tree();
		~JEventProcessor_eta_gg_tree(){};
		const char* className(void){return "JEventProcessor_eta_gg_tree";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void){return NOERROR;};
		jerror_t fini(void);
		
		//---------------------------------------//
		// Functions
		
		double energy_after_recoil(double eb, double theta, double m0, double mp);
		double fcal_energy_res(double e);
		double get_acc_scaling_factor(double eb);
		
		void write_events(uint64_t eventnumber, double rfTime, vector<const DMCThrown*> mc_thrown);
		
		void write_events(uint64_t eventnumber, double rfTime,
			vector<const DBeamPhoton*> beam_photons, 
			vector<const DFCALShower*> fcal_showers,
			vector<const DBCALShower*> bcal_showers,
			vector<const DTOFPoint*> tof_points,
			vector<const DSCHit*> sc_hits,
			vector<const DMCThrown*> mc_thrown);
		
		//---------------------------------------//
		// Geometry
		
		DVector3 m_beamSpot;
		DVector3 m_fcalFace, m_ccalFace;
		DVector3 m_fcal_correction, m_ccal_correction;
		
		vector<vector<DVector3>> m_sc_pos, m_sc_norm;
		
		int m_phase_val = 0;
		
		double m_HodoscopeHiFactor    = 1.0;
		double m_HodoscopeHiFactorErr = 1.0;
		double m_HodoscopeLoFactor    = 1.0;
		double m_HodoscopeLoFactorErr = 1.0;
		double m_MicroscopeFactor     = 1.0;
		double m_MicroscopeFactorErr  = 1.0;
		double m_TAGMEnergyBoundHi    = 1.0;
		double m_TAGMEnergyBoundLo    = 1.0;
		
		//---------------------------------------//
		// Cuts (defaults set inside constructor)
		
		int m_USE_LOG_WEIGHT, m_SAVE_MC_NOHITS;
		
		double m_FCAL_RF_CUT, m_BCAL_RF_CUT, m_TOF_RF_CUT, m_BEAM_RF_CUT;
		
		double m_MIN_BEAM_ENERGY;
		double m_DELTA_E_CUT;
		
		//---------------------------------------//
		// Constants 
		
		Particle_t m_Target;
		
		const double m_pi0      =  0.1349770;   // [GeV]
		const double m_eta      =  0.547862;    // [GeV]
		const double m_etap     =  0.95778;     // [GeV]
		
		const double m_Proton   =  0.938272046; // [GeV]
		const double m_Deuteron =  1.875612859; // [GeV]
		const double m_He4      =  3.727379238; // [GeV]
		const double m_Be9      =  8.39479;     // [GeV]
		const double m_C12      = 11.17793;     // [GeV]
		
		const double m_c        = 29.9792458;   // [cm/ns]
		
		double m_beam_bunches_main = 1.0;
		double m_beam_bunches_acc  = 5.0;
		
		//---------------------------------------//
		// Tree
		
		DTreeInterface *dTreeInterface;
		static thread_local DTreeFillData dTreeFillData;
};

#endif // _JEventProcessor_eta_gg_tree_
