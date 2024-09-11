#ifndef _EtaAna_
#define _EtaAna_

#define MAX_BEAM 1000
#define MAX_FCAL 100
#define MAX_BCAL 100
#define MAX_TOF  200
#define MAX_MC   20

using namespace std;

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
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
		
		int    m_phase_val;
		double m_target_length, m_target_density, m_target_atten;
		
		TVector3 m_fcal_face;
		TVector3 m_fcal_correction;
		TVector3 m_vertex;
		
		// Constants:
		
		static constexpr double m_c = TMath::C() * 1.e-7;   // [cm/ns]
		static constexpr double m_e = 0.510998928e-3;       // [GeV]
		static constexpr double m_fcal_block_size = 4.0157; // [cm]
		
		double m_Target;
		
		static constexpr double m_pi0      =  0.1349770;   // [GeV]
		static constexpr double m_eta      =  0.547862;    // [GeV]
		static constexpr double m_etap     =  0.95778;     // [GeV]
		
		static constexpr double m_Proton   =  0.938272046; // [GeV]
		static constexpr double m_Deuteron =  1.875612859; // [GeV]
		static constexpr double m_He4      =  3.727379238; // [GeV]
		static constexpr double m_Be9      =  8.39479;     // [GeV]
		static constexpr double m_C12      = 11.17793;     // [GeV]
		
		// Variables:
		
		string m_output_fname;
		
		TFile *m_infile;
		TTree *m_tree;
		
		int m_event;
		
		int loadTree();
		void setGeometry();
		void readEvent();
		void etaggAnalysis();
		
		void initializeReactionTypes();
		int getFinalState_bggen(int debug=0);
		TVector3 getFCALPosition(int index);
		
		int fcal_fiducial_cut(TVector3 pos, double cut_layer);
		
		void check_TOF_match(TVector3 pos, double &dx_min, double &dy_min, double &dz_min, double rf_time_cut);
		
		double energy_after_recoil(double eb, double theta, double m0, double mp);
		double fcal_energy_res(double e);
		
		// TTree Variables:
		int    m_eventNum;
		double m_rfTime;
		int    m_nbeam;
		int    m_tag_sys[MAX_BEAM];
		int    m_tag_counter[MAX_BEAM];
		double m_beam_e[MAX_BEAM];
		double m_beam_t[MAX_BEAM];
		double m_acc_scale_factor[MAX_BEAM];
		int    m_nfcal;
		double m_fcal_e[MAX_FCAL];
		double m_fcal_x[MAX_FCAL];
		double m_fcal_y[MAX_FCAL];
		double m_fcal_z[MAX_FCAL];
		double m_fcal_t[MAX_FCAL];
		int    m_fcal_nblocks[MAX_FCAL];
		int    m_nbcal;
		double m_bcal_e[MAX_BCAL];
		double m_bcal_x[MAX_BCAL];
		double m_bcal_y[MAX_BCAL];
		double m_bcal_z[MAX_BCAL];
		double m_bcal_t[MAX_BCAL];
		int    m_ntof;
		double m_tof_x[MAX_TOF];
		double m_tof_y[MAX_TOF];
		double m_tof_z[MAX_TOF];
		double m_tof_t[MAX_TOF];
		int    m_nmc;
		double m_mc_pdgtype[MAX_MC];
		double m_mc_x[MAX_MC];
		double m_mc_y[MAX_MC];
		double m_mc_z[MAX_MC];
		double m_mc_t[MAX_MC];
		double m_mc_e[MAX_MC];
		double m_mc_p[MAX_MC];
		double m_mc_theta[MAX_MC];
		double m_mc_phi[MAX_MC];
		
		// Histograms:
		
		vector<vector<Particle_t>> m_reaction_types;
		
		TH1F *h_thrown;
		TH1F *h_accepted;
		
		vector<TH1F*> h_theta;
		int m_n_bcal_vetos = 4;
		vector<vector<TH1F*>> h_theta_veto;
		vector<TH1F*> h_nbcal, h_bcal_energy, h_bcal_energy_single;
		vector<TH2F*> h_bcal_dt_vs_eta_angle;
		vector<TH1F*> h_bcal_deltaPhi;
		
		TH2F *h_elas_vs_mgg;
		TH2F *h_elas;
		TH2F *h_elas_corr;
		TH2F *h_mgg;
		TH2F *h_mgg_const;
		TH2F *h_mgg_const_corr;
		TH2F *h_mm_vs_theta;
		TH2F *h_mm_vs_theta_eta_cut;
		TH2F *h_mm_vs_theta_eta_elas_cut;
		TH2F *h_xy_1;
		TH2F *h_xy_2;
		
		//----------------------------------------------------------------------------------------------//
		
	public:
		vector<int> m_pdg_types;
		// Cuts:
		
		double m_FCAL_RF_CUT, m_BCAL_RF_CUT, m_TOF_RF_CUT, m_BEAM_RF_CUT;
		double m_MIN_BCAL_ENERGY;
		double m_MIN_FCAL_ENERGY;
		double m_MIN_BEAM_ENERGY;
		
		double m_FCAL_TOF_CUT;
		
		double m_ELAS_CUT_SIGMA, m_ELAS_CUT_WIDTH;
		double m_ELAS_CUT_MU_P0, m_ELAS_CUT_MU_P1;
		
		double m_beam_bunches_main = 1.0;
		double m_beam_bunches_acc  = 5.0; // how many bunches to use on each side for accidental subtraction
		
		EtaAna();
		~EtaAna(){};
		
		int getPrimexPhase(int run_number);
		
		// default analysis:
		void initHistograms();
		void runAnalysis(TString infname);
		void resetHistograms();
		void writeHistograms();
		
		void setRunNumber(int runNum);
		void setOutputFileName(string name);
};

#endif
