// $Id$
//
//    File: JEventProcessor_primex_eta_analysis.h
// Created: Fri Aug 11 14:26:44 EDT 2023
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.92.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_primex_eta_analysis_
#define _JEventProcessor_primex_eta_analysis_

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

// ROOT headers:
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"

using namespace jana;
using namespace std;

class JEventProcessor_primex_eta_analysis:public jana::JEventProcessor{
	public:
		JEventProcessor_primex_eta_analysis();
		~JEventProcessor_primex_eta_analysis(){};
		const char* className(void){return "JEventProcessor_primex_eta_analysis";}

	private:
		jerror_t init(void);
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
		jerror_t erun(void){return NOERROR;};
		jerror_t fini(void){return NOERROR;};
		
		//---------------------------------------//
		// Functions
		
		int fcal_fiducial_cut(DVector3 pos, double layer_cut);
		
		double energy_after_recoil(double eb, double theta, double m0, double mp);
		
		double fcal_energy_res(double e);
		
		double get_acc_scaling_factor(double eb);
		
		void check_TOF_match(DVector3 pos, double rfTime, vector<const DTOFPoint*> tof_points, 
			double &dx_min, double &dy_min, double &dt_min, double rf_time_cut);
		
		void eta_gg_analysis(
			vector<const DFCALShower*> fcal_showers, 
			vector<const DBeamPhoton*> beam_photons, 
			vector<const DBCALShower*> bcal_showers, 
			vector<const DTOFPoint*> tof_points, vector<const DSCHit*> sc_hits, 
			int n_fcal_showers, int n_good_fcal_showers, 
			double bcal_energy_sum, int n_bcal_showers, int n_bcal_showers_1ns, double bcal_phi, double bcal_rfdt,
			double rfTime, bool is_mc, double thrown_beam_energy, double thrown_eta_angle
		);
		
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
		
		int m_USE_LOG_WEIGHT, m_BYPASS_TRIGGER;
		
		double m_FCAL_RF_CUT, m_CCAL_RF_CUT, m_BCAL_RF_CUT, m_TOF_RF_CUT, m_BEAM_RF_CUT;
		
		double m_MIN_BCAL_ENERGY;
		double m_MIN_FCAL_ENERGY;
		double m_MIN_CCAL_ENERGY;
		double m_MIN_BEAM_ENERGY;
		
		double m_FCAL_TOF_CUT;
		
		double m_ELAS_CUT_SIGMA, m_ELAS_CUT_WIDTH;
		double m_ELAS_CUT_MU_P0, m_ELAS_CUT_MU_P1;
		
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
		
		//---------------------------------------//
		// Histograms
		
		static const int N_TRIGS = 4;
		vector<string> trigger_names = {"FCAL+CCAL Trigger", 
			"Low-energy FCAL Trigger", "PS Trigger", "CCAL Trigger"};
		
		TH1F *h_fcal_energy_sum[N_TRIGS];
		
		TH1F *h_theta_thrown[13];
		
		TH1F *h_fcal_rf_dt[N_TRIGS];
		TH1F *h_bcal_rf_dt[N_TRIGS];
		TH1F *h_ccal_rf_dt[N_TRIGS];
		TH1F *h_tagh_rf_dt[N_TRIGS];
		TH1F *h_tagm_rf_dt[N_TRIGS];
		TH1F  *h_tof_rf_dt[N_TRIGS];
		TH1F   *h_sc_rf_dt[N_TRIGS];
		
		TH1F *h_fcal_tof_dx, *h_fcal_tof_dy, *h_fcal_tof_dr;
		TH1F *h_fcal_tof_dt, *h_fcal_tof_dt_cut;
		TH1F *h_fcal_tof_matches;
		
		TH1F *h_beam_rf_dt_cut;
		
		TH2F *h_elas_vs_mgg;
		TH2F *h_elas;
		TH2F *h_elas_corr, *h_elas_corr_main, *h_elas_corr_side;
		
		TH2F *h_mgg,       *h_mgg_main,       *h_mgg_side;
		TH2F *h_mgg_const, *h_mgg_const_main, *h_mgg_const_side;
		TH2F *h_mgg_const_corr;
		
		// plots 
		static const int m_n_vetos = 6;
		
		TH2F *h_elas_veto[m_n_vetos];
		TH2F *h_mgg_veto[m_n_vetos];
		TH2F *h_mgg_const_veto[m_n_vetos];
		
		TH2F *h_rec_vs_thrown;
		TH2F *h_mgg_thrown;
		TH2F *h_mgg_const_thrown;
		
		TH2F *h_hmass,              *h_mm_vs_theta;
		TH2F *h_hmass_eta_cut,      *h_mm_vs_theta_eta_cut;
		TH2F *h_hmass_eta_elas_cut, *h_mm_vs_theta_eta_elas_cut;
		
		TH2F *h_xy_1, *h_xy_2;
};

#endif // _JEventProcessor_primex_eta_analysis_
