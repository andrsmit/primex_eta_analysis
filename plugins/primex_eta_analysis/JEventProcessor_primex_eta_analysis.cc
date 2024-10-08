// $Id$
//
//    File: JEventProcessor_primex_eta_analysis.cc
// Created: Fri Aug 11 14:26:44 EDT 2023
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.92.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_primex_eta_analysis.h"

extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_primex_eta_analysis());
	}
} // "C"

//------------------
// JEventProcessor_primex_eta_analysis (Constructor)
//------------------
JEventProcessor_primex_eta_analysis::JEventProcessor_primex_eta_analysis()
{
	// default values for the RF timing cuts for each sub-detector:
	m_BEAM_RF_CUT =  2.004;
	m_FCAL_RF_CUT =  2.0;
	m_BCAL_RF_CUT = 12.0;
	m_CCAL_RF_CUT =  2.0;
	m_TOF_RF_CUT  =  1.0;
	
	// default values for the minimum energy cuts:
	m_MIN_FCAL_ENERGY = 0.5; // energy of each FCAL shower used in the analysis
	m_MIN_BEAM_ENERGY = 8.0; // energy of the tagged photon energy
	m_MIN_BCAL_ENERGY = 0.;  // energy sum of BCAL showers within timing cut
	m_MIN_CCAL_ENERGY = 0.5; // energy sum of CCAL showers within timing cut
	
	// default cut value for selecting a match between the FCAL and TOF:
	m_FCAL_TOF_CUT = 8.0; // distance between fcal shower and closest DTOFPoint
	
	// default value for elasticity cut:
	m_ELAS_CUT_SIGMA = 0.031;
	m_ELAS_CUT_WIDTH = 3.0;
	m_ELAS_CUT_MU_P0 = 1.0; // mu = p0 + p1*E_gamma
	m_ELAS_CUT_MU_P1 = 0.0;
	
	// miscellaneous:
	m_USE_LOG_WEIGHT = 0; // force log-weighted FCAL position (should be automatically set in DFCALShower_factory)
	m_BYPASS_TRIGGER = 0; // determines whether or not to check the trigger bits set for each event
	
	//-------------------------------------------------------------------------------------//
	// allow for command-line overriding of the default values:
	
	gPARMS->SetDefaultParameter("primex_eta_analysis:FCAL_RF_CUT",     m_FCAL_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:BEAM_RF_CUT",     m_BEAM_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:BCAL_RF_CUT",     m_BCAL_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:CCAL_RF_CUT",     m_CCAL_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:TOF_RF_CUT",      m_TOF_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:MIN_FCAL_ENERGY", m_MIN_FCAL_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis:MIN_BEAM_ENERGY", m_MIN_BEAM_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis:MIN_BCAL_ENERGY", m_MIN_BCAL_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis:MIN_CCAL_ENERGY", m_MIN_CCAL_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis:FCAL_TOF_CUT",    m_FCAL_TOF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_SIGMA",  m_ELAS_CUT_SIGMA);
	gPARMS->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_WIDTH",  m_ELAS_CUT_WIDTH);
	gPARMS->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_MU_P0",  m_ELAS_CUT_MU_P0);
	gPARMS->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_MU_P1",  m_ELAS_CUT_MU_P1);
	gPARMS->SetDefaultParameter("primex_eta_analysis:USE_LOG_WEIGHT",  m_USE_LOG_WEIGHT);
	gPARMS->SetDefaultParameter("primex_eta_analysis:BYPASS_TRIGGER",  m_BYPASS_TRIGGER);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_primex_eta_analysis::init(void)
{
	TDirectory *dir_primex_eta = new TDirectoryFile("primex_eta_analysis", "primex_eta_analysis");
	dir_primex_eta->cd();
	
	// Distribution of the FCAL energy sum for each different trigger type:
	for(int itrig=0; itrig<N_TRIGS; itrig++) {
		h_fcal_energy_sum[itrig] = new TH1F(Form("fcal_energy_sum_%d",itrig), 
			Form("FCAL Shower Energy Sum (%s); E_{FCAL} [GeV]", trigger_names[itrig].c_str()), 
			1200, 0., 12.);
	}
	
	// Thrown angle distribution with different cuts on the beam photon energy:
	TDirectory *dir_thrown = new TDirectoryFile("thrown", "thrown");
	dir_thrown->cd();
	for(int icut=0; icut<13; icut++) {
		double eb_cut = 7.6 + 0.2*(double)(icut);
		h_theta_thrown[icut] = new TH1F(Form("theta_thrown_%02d", icut), 
			Form("Thrown Angle of #eta (E_{#gamma} > %.1f GeV)", eb_cut), 650, 0., 6.5);
	}
	dir_thrown->cd("../");
	
	//====================================================================================//
	
	TDirectory *dir_timing = new TDirectoryFile("rf_timing", "rf_timing");
	dir_timing->cd();
	// Timing distributions for each different trigger type:
	for(int itrig=0; itrig<N_TRIGS; itrig++) {
		h_fcal_rf_dt[itrig] = new TH1F(Form("fcal_rf_dt_%d",itrig), 
			Form("FCAL - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100., 100.);
		h_bcal_rf_dt[itrig] = new TH1F(Form("bcal_rf_dt_%d",itrig), 
			Form("BCAL - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100., 100.);
		h_ccal_rf_dt[itrig] = new TH1F(Form("ccal_rf_dt_%d",itrig), 
			Form("CCAL - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100., 100.);
		h_tof_rf_dt[itrig] = new TH1F(Form("tof_rf_dt_%d",itrig), 
			Form("TOF - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100., 100.);
		h_sc_rf_dt[itrig] = new TH1F(Form("sc_rf_dt_%d",itrig), 
			Form("SC - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100., 100.);
		h_tagh_rf_dt[itrig] = new TH1F(Form("tagh_rf_dt_%d",itrig), 
			Form("TAGH - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100.0, 100.0);
		h_tagm_rf_dt[itrig] = new TH1F(Form("tagm_rf_dt_%d",itrig), 
			Form("TAGM - RF Time (%s); [ns]", trigger_names[itrig].c_str()), 
			2000, -100.0, 100.0);
	}
	dir_timing->cd("../");
	
	//====================================================================================//
	
	TDirectory *dir_gg = new TDirectoryFile("eta_gg", "eta_gg");
	dir_gg->cd();
	
	// FCAL-TOF matching distributions:
	h_fcal_tof_dx = new TH1F("fcal_tof_dx", "x_{FCAL} - x_{TOF} (closest DTOFPoint); [cm]", 2000, -100., 100.);
	h_fcal_tof_dy = new TH1F("fcal_tof_dy", "y_{FCAL} - y_{TOF} (closest DTOFPoint); [cm]", 2000, -100., 100.);
	h_fcal_tof_dr = new TH1F("fcal_tof_dr", "Distance between FCAL Shower and closest DTOFPoint; [cm]", 1000, 0., 100.);
	
	h_fcal_tof_dt      = new TH1F("fcal_tof_dt",      "t_{FCAL} - t_{TOF}; [ns]", 2000, -100., 100.);
	h_fcal_tof_dt_cut  = new TH1F("fcal_tof_dt_cut",  "t_{FCAL} - t_{TOF}; [ns]", 2000, -100., 100.);
	h_fcal_tof_matches = new TH1F("fcal_tof_matches", "Number of TOF Matches per 2-#gamma Event", 3, -0.5, 2.5);
	
	// Plot the level of accidentals after elasticity cut:
	h_beam_rf_dt_cut = new TH1F("beam_rf_dt_cut", "; t_{#gamma} - t_{RF} (ns); counts / 0.1 ns", 2000, -100., 100.);
	
	// Plot the level of accidentals in the start counter after elasticity cut:
	h_sc_rf_dt_cut = new TH1F("sc_rf_dt_cut", "; t_{SC} - t_{RF} (ns); counts / 0.1 ns", 2000, -100., 100.);
	
	// Elasticity vs. mass ratio:
	h_elas_vs_mgg = new TH2F("elas_vs_mgg", 
		"; M_{#gamma#gamma}/M_{#eta}(PDG); #left(E_{1}+E_{2}#right)/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)", 
		1000, 0., 2., 1000, 0., 2.);
	
	//------------------------------------//
	
	// Elasticity with tagged photon:
	h_elas = new TH2F("elas", "Elasticity", 650, 0., 6.5, 1000, 0., 2.);
	h_elas->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elas->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#gamma}");
	
	// Elasticity with coherently-produced eta:
	h_elas_corr = new TH2F("elas_corr", "Elasticity (corrected for recoil)", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_corr->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elas_corr->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	h_elas_corr_main = new TH2F("elas_corr_main", "Elasticity (corrected for recoil)", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_corr_main->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elas_corr_main->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	h_elas_corr_side = new TH2F("elas_corr_side", "Elasticity (corrected for recoil)", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_corr_side->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elas_corr_side->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	//------------------------------------//
	
	// 2-photon invariant mass vs. angle:
	h_mgg = new TH2F("mgg", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg->Sumw2();
	
	// main RF bunch:
	h_mgg_main = new TH2F("mgg_main", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_main->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_main->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg_main->Sumw2();
	
	// accidental sidebands:
	h_mgg_side = new TH2F("mgg_side", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_side->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_side->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg_side->Sumw2();
	
	//------------------------------------//
	
	// Energy-constrained mass vs. angle:
	h_mgg_const = new TH2F("mgg_const", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_const->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const->Sumw2();
	
	// main RF bunch:
	h_mgg_const_main = new TH2F("mgg_const_main", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_main->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_const_main->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const_main->Sumw2();
	
	// Accidental side-bands:
	h_mgg_const_side = new TH2F("mgg_const_side", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_side->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_const_side->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const_side->Sumw2();
	
	// Energy-constrained mass vs. energy-constrained angle:
	h_mgg_const_corr = new TH2F("mgg_const_corr", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_corr->GetXaxis()->SetTitle("#theta_{#gamma#gamma}^{constr} [#circ]");
	h_mgg_const_corr->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const_corr->Sumw2();
	
	//------------------------------------//
	
	// x-y position of FCAL showers that survive all cuts:
	h_xy_1 = new TH2F("xy_1", "Postion of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -100., 100., 500, -100., 100.);
	h_xy_2 = new TH2F("xy_2", "Postion of Shower 2; x_{2} [cm]; y_{2} [cm]", 500, -100., 100., 500, -100., 100.);
	
	//------------------------------------//
	// Veto Plots:
	
	for(int ihist=0; ihist<m_n_vetos; ihist++) {
		
		h_elas_veto[ihist] = new TH2F(Form("elas_veto_%d",ihist), 
			Form("Elasticity (Veto Option %d)", ihist+1), 650, 0., 6.5, 1000, 0., 2.0);
		h_elas_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_elas_veto[ihist]->GetYaxis()->SetTitle(
			"E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
		
		h_mgg_veto[ihist] = new TH2F(Form("mgg_veto_%d",ihist), 
			Form("Two-Photon Invariant Mass (Veto Option %d)", ihist+1), 650, 0., 6.5, 600, 0., 1.2);
		h_mgg_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg_veto[ihist]->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
		h_mgg_veto[ihist]->Sumw2();
		
		h_mgg_const_veto[ihist] = new TH2F(Form("mgg_const_veto_%d",ihist), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)", ihist+1), 650, 0., 6.5, 600, 0., 1.2);
		h_mgg_const_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg_const_veto[ihist]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
		h_mgg_const_veto[ihist]->Sumw2();
		
		h_mgg_const_coh_veto[ihist] = new TH2F(Form("mgg_const_coh_veto_%d",ihist), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)", ihist+1), 650, 0., 6.5, 600, 0., 1.2);
		h_mgg_const_coh_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg_const_coh_veto[ihist]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
		h_mgg_const_coh_veto[ihist]->Sumw2();
		
		h_mm_veto[ihist] = new TH2F(Form("mm_veto_%d",ihist),
			Form("Missing Mass Squared (Veto Option %d)", ihist+1), 650, 0., 6.5, 1500, 0., 30.0);
		h_mm_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_veto[ihist]->GetYaxis()->SetTitle("M_{miss}^{2} [GeV^{2}/c^{4}]");
		h_mm_veto[ihist]->Sumw2();
		
		h_mm_const_veto[ihist] = new TH2F(Form("mm_const_veto_%d",ihist),
			Form("Missing Mass Squared (Energy-constrained) (Veto Option %d)", ihist+1), 650, 0., 6.5, 1500, 0., 30.0);
		h_mm_const_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_const_veto[ihist]->GetYaxis()->SetTitle("M_{miss}^{2} [GeV^{2}/c^{4}]");
		h_mm_const_veto[ihist]->Sumw2();
		
		h_mm_const_coh_veto[ihist] = new TH2F(Form("mm_const_coh_veto_%d",ihist),
			Form("Missing Mass Squared (Energy-constrained) (Veto Option %d)", ihist+1), 650, 0., 6.5, 1500, 0., 30.0);
		h_mm_const_coh_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_const_coh_veto[ihist]->GetYaxis()->SetTitle("M_{miss}^{2} [GeV^{2}/c^{4}]");
		h_mm_const_coh_veto[ihist]->Sumw2();
	}
	
	//------------------------------------//
	
	// Missing-mass vs. angle:
	
	h_mm_vs_theta = new TH2F("mm_vs_theta", 
		"Squared Missing Mass", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_vs_theta->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_vs_theta->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_vs_theta->Sumw2();
	
	h_mm_const_vs_theta = new TH2F("mm_const_vs_theta", 
		"Squared Missing Mass", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_const_vs_theta->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_const_vs_theta->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_const_vs_theta->Sumw2();
	
	h_mm_const_coh_vs_theta = new TH2F("mm_const_coh_vs_theta", 
		"Squared Missing Mass", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_const_coh_vs_theta->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_const_coh_vs_theta->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_const_coh_vs_theta->Sumw2();
	
	// With cut on invariant mass of 2-photon pair:
	
	h_mm_vs_theta_eta_cut = new TH2F("mm_vs_theta_eta_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} cut)", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_vs_theta_eta_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_vs_theta_eta_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_vs_theta_eta_cut->Sumw2();
	
	h_mm_const_vs_theta_eta_cut = new TH2F("mm_const_vs_theta_eta_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} cut)", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_const_vs_theta_eta_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_const_vs_theta_eta_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_const_vs_theta_eta_cut->Sumw2();
	
	h_mm_const_coh_vs_theta_eta_cut = new TH2F("mm_const_coh_vs_theta_eta_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} cut)", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_const_coh_vs_theta_eta_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_const_coh_vs_theta_eta_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_const_coh_vs_theta_eta_cut->Sumw2();
	
	// apply elasticity cut to both:
	
	h_mm_vs_theta_eta_elas_cut = new TH2F("mm_vs_theta_eta_elas_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} + elasticity cut)", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_vs_theta_eta_elas_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_vs_theta_eta_elas_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_vs_theta_eta_elas_cut->Sumw2();
	
	h_mm_const_vs_theta_eta_elas_cut = new TH2F("mm_const_vs_theta_eta_elas_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} + elasticity cut)", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_const_vs_theta_eta_elas_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_const_vs_theta_eta_elas_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_const_vs_theta_eta_elas_cut->Sumw2();
	
	h_mm_const_coh_vs_theta_eta_elas_cut = new TH2F("mm_const_coh_vs_theta_eta_elas_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} + elasticity cut)", 650, 0., 6.5, 1500, 0., 30.);
	h_mm_const_coh_vs_theta_eta_elas_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_const_coh_vs_theta_eta_elas_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_const_coh_vs_theta_eta_elas_cut->Sumw2();
	
	/*
	// Hybrid (rotated) mass:
	h_hmass = new TH2F("hmass", "Hybrid Mass", 
		650, 0., 6.5, 1000, -1.0, 1.0);
	h_hmass->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_hmass->GetYaxis()->SetTitle("Hybrid Mass");
	h_hmass->Sumw2();
	
	// apply eta invariant mass cut:
	
	h_hmass_eta_cut = new TH2F("hmass_eta_cut", "Hybrid Mass (m_{#gamma#gamma} cut)", 
		650, 0., 6.5, 1000, -1.0, 1.0);
	h_hmass_eta_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_hmass_eta_cut->GetYaxis()->SetTitle("Hybrid Mass");
	h_hmass_eta_cut->Sumw2();
	
	// apply elasticity cut:
	
	h_hmass_eta_elas_cut = new TH2F("hmass_eta_elas_cut", "Hybrid Mass (m_{#gamma#gamma} + elasticity cut)", 
		650, 0., 6.5, 1000, -1.0, 1.0);
	h_hmass_eta_elas_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_hmass_eta_elas_cut->GetYaxis()->SetTitle("Hybrid Mass");
	h_hmass_eta_elas_cut->Sumw2();
	*/
	//------------------------------------//
	
	// reconstructed angle of two-photon pair vs. thrown angle (only filled for MC):
	h_rec_vs_thrown = new TH2F("rec_vs_thrown", 
		"Reconstructed Angle vs. Thrown Angle; #theta_{thrown} [#circ]; #theta_{rec} [#circ]", 
		650, 0., 6.5, 650, 0., 6.5);
	h_rec_vs_thrown->Sumw2();
	
	// invariant mass vs. thrown angle of eta (only filled for MC):
	h_mgg_thrown = new TH2F("mgg_thrown", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_thrown->GetXaxis()->SetTitle("#theta_{thrown} [#circ]");
	h_mgg_thrown->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg_thrown->Sumw2();
	
	// energy-constrained invariant mass vs. thrown angle of eta (only filled for MC):
	h_mgg_const_thrown = new TH2F("mgg_const_thrown", "Energy-Constrained Inv Mass", 
		650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_thrown->GetXaxis()->SetTitle("#theta_{thrown} [#circ]");
	h_mgg_const_thrown->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const_thrown->Sumw2();
	
	dir_gg->cd("../");
	
	//====================================================================================//
	
	dir_primex_eta->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_primex_eta_analysis::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	//--------------------------------------------------------------//
	// Get geometry information for each run from CCDB:
	
	DGeometry*   dgeom = NULL;
	DApplication* dapp = dynamic_cast<DApplication*>(eventLoop->GetJApplication());
	if(dapp)     dgeom = dapp->GetDGeometry(runnumber);
	
	double loc_targetZ = 65.0;
	double loc_x, loc_y, loc_z;
	
	if(dgeom==NULL) {
		cerr << "No geometry accessbile to plugin." << endl;
		return RESOURCE_UNAVAILABLE;
	}
	
	// Get target position:
	dgeom->GetTargetZ(loc_targetZ);
	
	// Get position of FCAL face:
	dgeom->GetFCALPosition(loc_x, loc_y, loc_z);
	m_fcalFace.SetXYZ(loc_x, loc_y, loc_z);
	
	// Get position of CCAL face:
	dgeom->GetCCALPosition(loc_x, loc_y, loc_z);
	m_ccalFace.SetXYZ(loc_x, loc_y, loc_z);
	
	// Get beam spot on center of target:
	jana::JCalibration *jcalib = japp->GetJCalibration(runnumber);
	std::map<string, float> beam_spot;
	jcalib->Get("PHOTON_BEAM/beam_spot", beam_spot);
	m_beamSpot.SetXYZ(beam_spot.at("x"), beam_spot.at("y"), loc_targetZ);
	
	// Get start counter geomety:
	dgeom->GetStartCounterGeom(m_sc_pos, m_sc_norm);
	
	//--------------------------------------------------------------//
	// Set the target according to the run number:
	
	if(runnumber < 60000 || 
		( 70000 <= runnumber && runnumber <=  79999) || 
		(120000 <= runnumber && runnumber <= 129999)
	) {
		m_Target = Proton;
	} else if(
		( 60000 <= runnumber && runnumber <=  61354) || 
		( 80000 <= runnumber && runnumber <=  81381) || 
		(110000 <= runnumber && runnumber <= 110621)
	) { 
		m_Target = Be9;
	} else if(
		( 60000 <= runnumber && runnumber <=  69999) || 
		( 80000 <= runnumber && runnumber <=  81716) || 
		(110000 <= runnumber && runnumber <= 112001) || 
		( 90034 <= runnumber && runnumber <=  90200) || 
		( 90607 <= runnumber && runnumber <=  90660)) {
		m_Target = Helium;
	} else if(
		( 90207 <= runnumber && runnumber <=  90249) || 
		( 90558 <= runnumber && runnumber <=  90601)
	) {
		m_Target = Deuteron;
	} else if(90263 <= runnumber && runnumber <= 90536) {
		m_Target = C12;
	} else {
		m_Target = Proton;
	}
	
	//--------------------------------------------------------------//
	
	if(runnumber>60000 && runnumber<69999) 
	{
		m_phase_val = 1;
		
		//--------------------------------------------------------------------//
		// For phase 1, update the geometry from Compton alignment studies:
		
		// (2/4/2024): Correction to alignment after Simon updated beam spot with new CDC alignment:
		
		m_fcal_correction.SetXYZ(0.455 - m_fcalFace.X(), -0.032 - m_fcalFace.Y(), 0.0);
		m_fcalFace += m_fcal_correction;
		
		m_ccal_correction.SetXYZ(-0.082 - m_ccalFace.X(), 0.061 - m_ccalFace.Y(), 0.0);
		if(runnumber>=61483) m_ccal_correction.SetY(0.051 - m_ccalFace.Y());
		m_ccalFace += m_ccal_correction;
		
		if(runnumber<61483) 
		{
			m_beamSpot.SetX( 0.027);
			m_beamSpot.SetY(-0.128);
		} else if(runnumber<61774) 
		{
			m_beamSpot.SetX( 0.001);
			m_beamSpot.SetY(-0.077);
		} else 
		{
			m_beamSpot.SetX( 0.038);
			m_beamSpot.SetY(-0.095);
		}
	}
	else if(runnumber>80000 && runnumber < 89999) 
	{
		m_phase_val = 2;
		m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
		m_ccal_correction.SetXYZ(0.0, 0.0, 0.0);
	} 
	else if(runnumber>110000 && runnumber < 119999) 
	{
		m_phase_val = 3;
		m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
		m_ccal_correction.SetXYZ(0.0, 0.0, 0.0);
	}
	else 
	{
		m_phase_val = 0;
		m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
		m_ccal_correction.SetXYZ(0.0, 0.0, 0.0);
	}
	
	/*------------------------------------------------------------------------------------------------------*/
	// Code to obtain the scaling factors for accidental beam bunches 
	//     (copied from DAnalysisUtilities.cc in gluex_root_analysis)
	
	ostringstream locCommandStream;
	locCommandStream << "ccdb dump ANALYSIS/accidental_scaling_factor -r " << runnumber;
	FILE* locInputFile = gSystem->OpenPipe(locCommandStream.str().c_str(), "r");
	if(locInputFile == NULL) {
		
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		return NOERROR;
		/*
		cerr << "Could not load ANALYSIS/accidental_scaling_factor from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	//get the first line
	char buff[1024]; // I HATE char buffers
	if(fgets(buff, sizeof(buff), locInputFile) == NULL)
	{
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		gSystem->ClosePipe(locInputFile);
		return NOERROR;
		/*
		//vector<double> locCachedValues = { -1., -1., -1., -1., -1., -1., -1., -1. };
		//dAccidentalScalingFactor_Cache[runnumber] = locCachedValues;   // give up for this run
		gSystem->ClosePipe(locInputFile);
		cerr << "Could not parse ANALYSIS/accidental_scaling_factor from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	//get the second line (where the # is)
	if(fgets(buff, sizeof(buff), locInputFile) == NULL)
	{
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		gSystem->ClosePipe(locInputFile);
		return NOERROR;
		/*
		//vector<double> locCachedValues = { -1., -1., -1., -1., -1., -1., -1., -1. };
		//dAccidentalScalingFactor_Cache[runnumber] = locCachedValues;   // give up for this run
		gSystem->ClosePipe(locInputFile);
		cerr << "Could not parse ANALYSIS/accidental_scaling_factor from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	// catch some CCDB error conditions
	if(strncmp(buff, "Cannot", 6) == 0) 
	{
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_TAGMEnergyBoundHi    = 9.00;
		m_TAGMEnergyBoundLo    = 8.00;
		
		gSystem->ClosePipe(locInputFile);
		return NOERROR;
		/*
		// no assignment for this run
		//vector<double> locCachedValues = { -1., -1., -1., -1., -1., -1., -1., -1. };
		//dAccidentalScalingFactor_Cache[runnumber] = locCachedValues;   // give up for this run
		gSystem->ClosePipe(locInputFile);
		cerr << "No data available for ANALYSIS/accidental_scaling_factor, run " << runnumber << " from CCDB !" << endl;
		gSystem->Exit(1);        // make sure we don't fail silently
		return RESOURCE_UNAVAILABLE;    // sanity check, this shouldn't be executed!
		*/
	}
	
	istringstream locStringStream(buff);
	
	double locHodoscopeHiFactor    = -1.0;
	double locHodoscopeHiFactorErr = -1.0;
	double locHodoscopeLoFactor    = -1.0;
	double locHodoscopeLoFactorErr = -1.0;
	double locMicroscopeFactor     = -1.0;
	double locMicroscopeFactorErr  = -1.0;
	double locTAGMEnergyBoundHi    = -1.0;
	double locTAGMEnergyBoundLo    = -1.0;
	
	//extract it
	locStringStream >> locHodoscopeHiFactor >> locHodoscopeHiFactorErr >> locHodoscopeLoFactor
		>> locHodoscopeLoFactorErr >> locMicroscopeFactor >> locMicroscopeFactorErr
		>> locTAGMEnergyBoundHi >> locTAGMEnergyBoundLo;
	
	//Close the pipe
	gSystem->ClosePipe(locInputFile);
	
	m_HodoscopeHiFactor    = locHodoscopeHiFactor;
	m_HodoscopeHiFactorErr = locHodoscopeHiFactorErr;
	m_HodoscopeLoFactor    = locHodoscopeLoFactor;
	m_HodoscopeLoFactorErr = locHodoscopeLoFactorErr;
	m_MicroscopeFactor     = locMicroscopeFactor;
	m_MicroscopeFactorErr  = locMicroscopeFactorErr;
	m_TAGMEnergyBoundHi    = locTAGMEnergyBoundHi;
	m_TAGMEnergyBoundLo    = locTAGMEnergyBoundLo;
	
	/*------------------------------------------------------------------------------------------------------*/
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_primex_eta_analysis::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	//-----------------------------------------------------//
	// Get RF Time
	
	const DEventRFBunch *locRFBunch = NULL;
	try {
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch(...) {return NOERROR;}
	double locRFTime = locRFBunch->dTime;
	
	//-----------------------------------------------------//
	// Data objects
	
	vector<const DBeamPhoton*> locDBeamPhotons;
	eventLoop->Get(locDBeamPhotons);
	
	vector<const DFCALShower*> locDFCALShowers;
	eventLoop->Get(locDFCALShowers);
	
	vector<const DCCALShower*> locDCCALShowers;
	eventLoop->Get(locDCCALShowers);
	
	vector<const DBCALShower*> locDBCALShowers;
	eventLoop->Get(locDBCALShowers);
	
	vector<const DTOFPoint*> locDTOFPoints;
	eventLoop->Get(locDTOFPoints);
	
	vector<const DSCHit*> locDSCHits;
	eventLoop->Get(locDSCHits);
	
	vector<const DMCThrown*> locDMCThrown;	
	eventLoop->Get(locDMCThrown);
	
	//-----------------------------------------------------//
	// Trigger information
	
	bool trig_conditions[N_TRIGS];
	for(int itrig=0; itrig<N_TRIGS; itrig++) { trig_conditions[itrig] = false; }
	
	if(locDMCThrown.size() > 0) 
	{
		trig_conditions[0] = true;
		trig_conditions[1] = true;
	} else if(m_BYPASS_TRIGGER) {
		trig_conditions[0] = true;
		trig_conditions[1] = true;
	} else {
		const DL1Trigger *trig = NULL;
		try {
			eventLoop->GetSingle(trig);
		} catch (...) {}
		if (trig == NULL) { return NOERROR; }
		
		uint32_t trigmask    = trig->trig_mask;	
		uint32_t fp_trigmask = trig->fp_trig_mask;
		
		if(!trigmask)   return NOERROR;
		if(fp_trigmask) return NOERROR;
		
		if(trigmask & (1 <<  1)) trig_conditions[0] = true; // FCAL Energy Sum
		if(trigmask & (1 <<  2)) trig_conditions[1] = true; // FCAL Energy Sum (low-threshold, phase 3 only)
		if(trigmask & (1 <<  3)) trig_conditions[2] = true; // PS
		if(trigmask & (1 << 10)) trig_conditions[3] = true; // CCAL Energy Sum
	}
	
	//-----------------------------------------------------//
	// Apply fill lock for multi-threaded running:
	
	japp->RootFillLock(this);
	
	//-----------------------------------------------------//
	// Plot thrown distributions (if MC)
	
	bool   locIsMC         = false;
	double locThrownEnergy = 0.;
	double locThrownAngle  = 0.;
	
	if(locDMCThrown.size() > 0) {
		const DMCReaction* locDMCReactions = NULL;
		eventLoop->GetSingle(locDMCReactions);
		
		TLorentzVector locEtaMCP4(0, 0, 0, 0);
		vector<TLorentzVector> locPhotonsMCList; locPhotonsMCList.clear();
		vector<TLorentzVector> locPipsMCList; locPipsMCList.clear();
		vector<TLorentzVector> locPimsMCList; locPimsMCList.clear();
		//vector<TLorentzVector> locPsMCList; locPsMCList.clear();
		//vector<TLorentzVector> locNsMCList; locNsMCList.clear();
		
		for(unsigned int i = 0; i < locDMCThrown.size(); i++) {
			const DMCThrown *mcthrown = locDMCThrown[i];
			double p     = mcthrown->momentum().Mag();
			double theta = mcthrown->momentum().Theta();
			double phi   = mcthrown->momentum().Phi();
			double px    = p * sin(theta) * cos(phi);
			double py    = p * sin(theta) * sin(phi);
			double pz    = p * cos(theta);
			TLorentzVector thrownP4(px, py, pz, p);
			if(mcthrown->type ==  1) locPhotonsMCList.push_back(thrownP4); // photon
			if(mcthrown->type ==  8) locPipsMCList.push_back(thrownP4);    // pi+
			if(mcthrown->type ==  9) locPimsMCList.push_back(thrownP4);    // pi-
			//if(mcthrown->type == 13) locNsMCList.push_back(thrownP4);    // neutron
			//if(mcthrown->type == 14) locPsMCList.push_back(thrownP4);    // proton
		}
		if(locPhotonsMCList.size() == 2 && locPipsMCList.size() == 0 && locPimsMCList.size() == 0) {
			locEtaMCP4 = locPhotonsMCList[0] + locPhotonsMCList[1];
		}
		if(locPhotonsMCList.size() == 2 && locPipsMCList.size() == 1 && locPimsMCList.size() == 1) {
			locEtaMCP4 = locPhotonsMCList[0] + locPhotonsMCList[1] + locPipsMCList[0] + locPimsMCList[0];
		}
		if(locPhotonsMCList.size() == 6 && locPipsMCList.size() == 0 && locPimsMCList.size() == 0) {
			locEtaMCP4 = locPhotonsMCList[0] + locPhotonsMCList[1] + locPhotonsMCList[2] + 
				locPhotonsMCList[3] + locPhotonsMCList[4] + locPhotonsMCList[5];
		}
		
		locIsMC         = true;
		locThrownEnergy = locDMCReactions->beam.energy();
		locThrownAngle  = locEtaMCP4.Theta() * TMath::RadToDeg();
		
		for(int icut=0; icut<13; icut++) {
			double eb_cut = 7.6 + 0.2*(double)(icut);
			if(locThrownEnergy>=eb_cut) {
				h_theta_thrown[icut]->Fill(locThrownAngle);
			}
		}
		if(locThrownEnergy<m_MIN_BEAM_ENERGY) {
			japp->RootFillUnLock(this);
			return NOERROR;
		}
	}
	
	//-----------------------------------------------------//
	// RF Timing Histograms:
	
	for(vector<const DBeamPhoton*>::const_iterator gam = locDBeamPhotons.begin(); 
		gam != locDBeamPhotons.end(); gam++) {
		double loc_t = (*gam)->time() - locRFTime;
		if((*gam)->dSystem==SYS_TAGH) {
			for(int itrig=0; itrig<N_TRIGS; itrig++) {
				if(trig_conditions[itrig]) h_tagh_rf_dt[itrig]->Fill(loc_t);
			}
		} else {
			for(int itrig=0; itrig<N_TRIGS; itrig++) {
				if(trig_conditions[itrig]) h_tagm_rf_dt[itrig]->Fill(loc_t);
			}
		}
	}
	
	int    locNFCALShowers  = 0, locNGoodFCALShowers = 0;
	double locFCALEnergySum = 0.;
	for(vector<const DFCALShower*>::const_iterator show = locDFCALShowers.begin(); 
		show != locDFCALShowers.end(); show++) {
		
		DVector3 loc_pos;
		if(m_USE_LOG_WEIGHT) {
			loc_pos = (*show)->getPosition_log();
		} else {
			loc_pos = (*show)->getPosition();
		}
		loc_pos = loc_pos - m_beamSpot + m_fcal_correction;
		double loc_t = (*show)->getTime() - (loc_pos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<N_TRIGS; itrig++) {
			if(trig_conditions[itrig]) h_fcal_rf_dt[itrig]->Fill(loc_t);
		}
		if(fabs(loc_t) < m_FCAL_RF_CUT) {
			locFCALEnergySum += (*show)->getEnergy();
			locNFCALShowers++;
			if(!fcal_fiducial_cut(loc_pos, 2.0) && (*show)->getEnergy() > m_MIN_FCAL_ENERGY) {
				locNGoodFCALShowers++;
			}
		}
	}
	for(int itrig=0; itrig<N_TRIGS; itrig++) {
		if(trig_conditions[itrig]) h_fcal_energy_sum[itrig]->Fill(locFCALEnergySum);
	}
	
	int    locNBCALShowers  = 0, locNBCALShowers_1ns = 0;
	double locBCALEnergySum = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	for(vector<const DBCALShower*>::const_iterator show = locDBCALShowers.begin(); 
		show != locDBCALShowers.end(); show++) {
		DVector3 loc_pos((*show)->x, (*show)->y, (*show)->z);
		loc_pos -= m_beamSpot;
		double loc_t = (*show)->t - (loc_pos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<N_TRIGS; itrig++) {
			if(trig_conditions[itrig]) h_bcal_rf_dt[itrig]->Fill(loc_t);
		}
		if(fabs(loc_t) < m_BCAL_RF_CUT) {
			locBCALEnergySum += (*show)->E;
			locNBCALShowers++;
			locBCALRFDT = loc_t;
			locBCALPhi  = loc_pos.Phi() * (180./TMath::Pi());
			if(fabs(loc_t) < 1.0) {
				locNBCALShowers_1ns++;
			}
		}
	}
	
	int    locNCCALShowers  = 0;
	double locCCALEnergySum = 0.;
	for(vector<const DCCALShower*>::const_iterator show = locDCCALShowers.begin(); 
		show != locDCCALShowers.end(); show++) {
		DVector3 loc_pos((*show)->x1, (*show)->y1, (*show)->z);
		loc_pos = loc_pos - m_beamSpot + m_ccal_correction;
		double loc_t = (*show)->time - (loc_pos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<N_TRIGS; itrig++) {
			if(trig_conditions[itrig]) h_ccal_rf_dt[itrig]->Fill(loc_t);
		}
		if(fabs(loc_t) < m_CCAL_RF_CUT) {
			locCCALEnergySum += (*show)->E;
			locNCCALShowers++;
		}
	}
	
	for(vector<const DTOFPoint*>::const_iterator tof = locDTOFPoints.begin(); 
		tof != locDTOFPoints.end(); tof++) {
		DVector3 loc_pos = (*tof)->pos - m_beamSpot;
		double loc_t = (*tof)->t - (loc_pos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<N_TRIGS; itrig++) {
			if(trig_conditions[itrig]) h_tof_rf_dt[itrig]->Fill(loc_t);
		}
	}
	
	for(vector<const DSCHit*>::const_iterator sc = locDSCHits.begin(); 
		sc != locDSCHits.end(); sc++) {
		double loc_t = (*sc)->t - locRFTime;
		for(int itrig=0; itrig<N_TRIGS; itrig++) {
			if(trig_conditions[itrig]) h_sc_rf_dt[itrig]->Fill(loc_t);
		}
	}
	
	//-----------------------------------------------------//
	// eta->2gamma analysis:
	
	if(locNFCALShowers>1) {
		eta_gg_analysis(
			locDFCALShowers, locDBeamPhotons,     locDBCALShowers,  locDTOFPoints, locDSCHits,
			locNFCALShowers, locNGoodFCALShowers, 
			locBCALEnergySum, locNBCALShowers, locNBCALShowers_1ns, locBCALPhi, locBCALRFDT,
			locRFTime, locIsMC, locThrownEnergy, locThrownAngle
		);
	}
	
	//-----------------------------------------------------//
	
	japp->RootFillUnLock(this);
	
	return NOERROR;
}

void JEventProcessor_primex_eta_analysis::eta_gg_analysis(
	vector<const DFCALShower*> fcal_showers, 
	vector<const DBeamPhoton*> beam_photons, 
	vector<const DBCALShower*> bcal_showers, 
	vector<const DTOFPoint*> tof_points, 
	vector<const DSCHit*> sc_hits, 
	int n_fcal_showers, int n_good_fcal_showers, 
	double bcal_energy_sum, int n_bcal_showers, int n_bcal_showers_1ns, double bcal_phi, double bcal_rfdt,
	double rfTime, bool is_mc, double thrown_beam_energy, double thrown_eta_angle
)
{
	// Reject events with more than 1 shower in the BCAL:
	//if(n_bcal_showers>1) return;
	
	// Apply multiplicity cut on the number of FCAL showers:
	if(!(n_fcal_showers==2 && n_good_fcal_showers==2)) return;
	
	int n_fcal_showers_total = static_cast<int>(fcal_showers.size());
	
	//----------------------------------------------------------------------------------//
	// Loop over all possible combinations of pairs of FCAL showers that pass the cuts:
	
	for(int ishow=0; ishow<n_fcal_showers_total; ishow++) {
		
		const DFCALShower *show1 = fcal_showers[ishow];
		DVector3 pos1;
		if(m_USE_LOG_WEIGHT) {
			pos1 = show1->getPosition_log();
		} else {
			pos1 = show1->getPosition();
		}
		pos1 = pos1 - m_beamSpot + m_fcal_correction;
		
		double t1  = show1->getTime() - (pos1.Mag()/m_c) - rfTime;
		double e1  = show1->getEnergy();
		
		// apply minimum energy and RF timing cuts:
		if(e1 < m_MIN_FCAL_ENERGY || fabs(t1) >= m_FCAL_RF_CUT) continue;
		
		// apply fiducial cut to remove the innermost two FCAL layers:
		if(fcal_fiducial_cut(pos1, 2.0)) continue;
		
		double px1 = e1*pos1.X()/pos1.Mag();
		double py1 = e1*pos1.Y()/pos1.Mag();
		double pz1 = e1*pos1.Z()/pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		check_TOF_match(pos1, rfTime, tof_points, tof_dx1, tof_dy1, tof_dt1, m_TOF_RF_CUT);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		// plot FCAL-TOF matching distributions for monitoring:
		if(n_fcal_showers==2) {
			h_fcal_tof_dx->Fill(tof_dx1);
			h_fcal_tof_dy->Fill(tof_dy1);
			h_fcal_tof_dr->Fill(tof_dr1);
			h_fcal_tof_dt->Fill(t1-tof_dt1);
			if(tof_dr1 < m_FCAL_TOF_CUT) {
				h_fcal_tof_dt_cut->Fill(t1-tof_dt1);
			}
		}
		
		//-----------------------------------------------------//
		
		for(int jshow=ishow+1; jshow<n_fcal_showers_total; jshow++) {
			
			const DFCALShower *show2 = fcal_showers[jshow];
			DVector3 pos2;
			if(m_USE_LOG_WEIGHT) {
				pos2 = show2->getPosition_log();
			} else {
				pos2 = show2->getPosition();
			}
			pos2 = pos2 - m_beamSpot + m_fcal_correction;
			
			double t2  = show2->getTime() - (pos2.Mag()/m_c) - rfTime;
			double e2  = show2->getEnergy();
			
			// apply minimum energy and RF timing cuts:
			if(e2 < m_MIN_FCAL_ENERGY || fabs(t2) >= m_FCAL_RF_CUT) continue;
			
			// apply fiducial cut to remove the innermost two FCAL layers:
			if(fcal_fiducial_cut(pos2, 2.0)) continue;
			
			double px2 = e2*pos2.X()/pos2.Mag();
			double py2 = e2*pos2.Y()/pos2.Mag();
			double pz2 = e2*pos2.Z()/pos2.Mag();
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			check_TOF_match(pos2, rfTime, tof_points, tof_dx2, tof_dy2, tof_dt2, m_TOF_RF_CUT);
			double tof_dr2 = sqrt(pow(tof_dx2,2.0)+pow(tof_dy2,2.0));
			
			//-----------------------------------------------------//
			// TOF Veto
			
			// reject combinations of FCAL showers where both showers are near a TOF hit:
			bool tof_veto = false;
			if((tof_dr1 < m_FCAL_TOF_CUT) && (tof_dr2 < m_FCAL_TOF_CUT)) tof_veto = true;
			
			// count the number of FCAL-TOF matches for monitoring:
			int n_tof_matches = 0;
			if(tof_dr1 < m_FCAL_TOF_CUT) n_tof_matches++;
			if(tof_dr2 < m_FCAL_TOF_CUT) n_tof_matches++;
			h_fcal_tof_matches->Fill(n_tof_matches);
			
			if(tof_veto) continue;
			
			//-----------------------------------------------------//
			// Two-Photon kinematics:
			
			double Egg  =  e1 +  e2; // energy of 2-photon pair
			double pggx = px1 + px2; // momentum along x-axis
			double pggy = py1 + py2; // momentum along y-axis
			double pggz = pz1 + pz2; // momentum along z-axis
			
			// transverse momentum:
			double pggt = sqrt(pow(pggx,2.0) + pow(pggy,2.0));
			
			// polar angle:
			double prod_th = (180./TMath::Pi()) * atan2(pggt,pggz);
			
			// azimuthal angle:
			double prod_phi = (180./TMath::Pi()) * atan2(pggy,pggx);
			
			// opening angle:
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
			
			// invariant mass:
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			//-----------------------------------------------------//
			// check different veto options:
			
			//----------------------------------------------------//
			// Check SC Matches:
			
			int loc_n_sc_hits = 0;
			int loc_n_sc_hits_coplanar = 0;
			for(vector<const DSCHit*>::const_iterator sc = sc_hits.begin(); sc != sc_hits.end(); sc++) {
				
				// only check hits between 1ns < (t_sc - t_RF) < 7ns 
				//    and with dE > 0.0002 (from DNeutralShower_factory)
				
				double loc_t  = (*sc)->t - rfTime;
				double loc_dE = (*sc)->dE;
				
				if((0.5<invmass) && (invmass<0.60)) h_sc_rf_dt_cut->Fill(loc_t);
				
				if((1.0 < loc_t) && (loc_t < 7.0) && (loc_dE > 0.0002)) {
					int sector = (*sc)->sector;
					double phi = m_sc_pos[sector-1][0].Phi() * (180./TMath::Pi());
					loc_n_sc_hits++;
					if(fabs(fabs(phi-prod_phi)-180.0) < 36.0) loc_n_sc_hits_coplanar++;
				}
			}
			//----------------------------------------------------//
			
			vector<int> loc_veto_options;
			for(int iveto=0; iveto<m_n_vetos; iveto++) loc_veto_options.push_back(0);
			
			// Option 0 (no veto): No Veto is applied:
			loc_veto_options[0] = 1;
			
			// Option 1 (strict): Remove events with any BCAL shower within +/-12ns:
			if(n_bcal_showers==0) loc_veto_options[1] = 1;
			
			// Option 2 (loose): Remove events with any BCAL shower within +/-1ns:
			if(n_bcal_showers_1ns==0) loc_veto_options[2] = 1;
			
			// Option 3 (looser): Keep events where there is EITHER (i) no BCAL shower within +/-12ns, 
			//        OR (ii) 1 BCAL shower that has opposite phi angle to two-photon pair in FCAL:
			if((n_bcal_showers==0) ||
				(n_bcal_showers==1 && fabs(fabs(bcal_phi-prod_phi)-180.0) < 30.0)) loc_veto_options[3] = 1;
			
			// Option 4 (improvement on option 3): Keep events where there is EITHER (i) no BCAL shower within +/-12ns, 
			//        OR (ii) 1 BCAL shower that has opposite phi angle to two-photon pair in FCAL 
			//        and is more than 1ns removed from RF time:
			if((n_bcal_showers==0) || 
				(n_bcal_showers==1 && fabs(fabs(bcal_phi-prod_phi)-180.0) < 30.0 && bcal_rfdt>1.0)) {
				loc_veto_options[4] = 1;
				
				// Option 5 (add in SC Veto): Remove events where there is a hit in the SC outside of the range:
				//          150 < |phi_SC - phi_FCAL| < 210:
				if(loc_n_sc_hits_coplanar==loc_n_sc_hits) loc_veto_options[5] = 1;
			}
			
			// Option 6: Use tight veto on SC only:
			if(loc_n_sc_hits==0) loc_veto_options[6] = 1;
			
			// Option 7: tight SC + tight BCAL vetos:
			if(loc_n_sc_hits==0 && n_bcal_showers==0) loc_veto_options[7] = 1;
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin(); 
				gam != beam_photons.end(); gam++) {
				
				double eb    = (*gam)->lorentzMomentum().E(); // energy of beam photon
				double brfdt = (*gam)->time() - rfTime;
				
				// remove beam photons below the minimum energy cut:
				if(eb < m_MIN_BEAM_ENERGY) continue;
				
				// Accidental subtraction procedure:
				//   - Fill histograms with a weight of 1.0 for beam photons within main RF bunch
				//   - Select two side-bands to the left and two-sidebands to the right (4 in total)
				//   - Fill histograms with a weight of -1/4 for beam photons in these sidebands.
				//      - An extra scaling factor is applied to out-of-time beam photons to account
				//        for the non-uniformity of beam bunches.
				//      - This scaling factor comes from the /ANALYSIS/accidental_scaling_factor tabls in the CCDB.
				//      - reference: GlueX-doc-4122 (B. Zihlmann)
				
				double fill_weight = 0.0;
				if(fabs(brfdt) < m_BEAM_RF_CUT) 
				{
					fill_weight = 1.0;
				}
				else if((-30.060<=brfdt && brfdt<=-22.044) || (22.044<=brfdt && brfdt<= 30.060)) 
				{
					fill_weight = -0.25*get_acc_scaling_factor(eb);
				}
				//else { continue; }
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				double eeta_coh = energy_after_recoil(eb, prod_th, m_eta, ParticleMass(m_Target));
				
				// Calculate the energy of the eta meson, assuming production on a free nucleon:
				double eeta = energy_after_recoil(eb, prod_th, m_eta, m_Proton);
				
				// Apply a cut on the elasticity
				//  (ratio of measured energy of 2-photons, to the calculated energy above):
				
				bool elas_cut = false;
				double loc_elas_mean  = m_ELAS_CUT_MU_P0 + m_ELAS_CUT_MU_P1*prod_th;
				double loc_elas_width = m_ELAS_CUT_WIDTH * m_ELAS_CUT_SIGMA;
				if(fabs((Egg/eeta)-loc_elas_mean)<loc_elas_width) elas_cut = true;
				
				// set a variable to indicate if the two-photon mass is consistent with an eta meson:
				bool eta_cut = false;
				if(0.497862<invmass && invmass<0.597862) eta_cut = true;
				
				// Plot timing distribution of beam photons after elasticity cut to see the level of accidentals:
				if(elas_cut && loc_veto_options[5]) {
					h_beam_rf_dt_cut->Fill(brfdt);
				}
				
				// If the beam photon wasn't in the main RF bunch or selected sidebands, skip it:
				if(fill_weight==0.0) continue;
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a coherently-produced eta meson:
				
				double sig1 = fcal_energy_res(e1);
				double sig2 = fcal_energy_res(e2);
				double sigr = pow(sig1/sig2,2.0);
				
				// Energy-constrained invariant mass assuming production on free nucleon:
				double e1c = e1/(1.+sigr) + (eeta-e2)/(1.+(1./sigr));
				double e2c = eeta - e1c;
				double invmass_const = sqrt(2.*e1c*e2c*(1.-cos12));
				
				// Energy-constrained invariant mass assuming coherent production on nucleus:
				double e1c_coh = e1/(1.+sigr) + (eeta_coh-e2)/(1.+(1./sigr));
				double e2c_coh = eeta_coh - e1c;
				double invmass_const_coh = sqrt(2.*e1c_coh*e2c_coh*(1.-cos12)); // energy-constrained invariant mass
				
				// re-compute the polar angle of the two-photon pair using these adjusted energies:
				double px1c  = e1c*pos1.X()/pos1.Mag();
				double py1c  = e1c*pos1.Y()/pos1.Mag();
				double pz1c  = e1c*pos1.Z()/pos1.Mag();
				double px2c  = e2c*pos2.X()/pos2.Mag();
				double py2c  = e2c*pos2.Y()/pos2.Mag();
				double pz2c  = e2c*pos2.Z()/pos2.Mag();
				double pggxc = px1c + px2c;
				double pggyc = py1c + py2c;
				double pggzc = pz1c + pz2c;
				double pggtc = sqrt(pow(pggxc,2.0) + pow(pggyc,2.0));
				double prod_th_const = (180./TMath::Pi()) * atan2(pggtc,pggzc);
				
				//-----------------------------------------------------//
				// Hybrid Mass
				
				//double hmass = (invmass/m_eta)*cos(TMath::Pi()/4.0) - (Egg/eeta)*sin(TMath::Pi()/4.0);
				
				//-----------------------------------------------------//
				// Missing Mass 
				//   - This is calculated assuming the reaction took place on a quasifree nucleon
				//     However, we neglect the fermi-momentum of said nucleon
				
				//double mmsq = 2.0*m_Proton*eb - 2.0*eb*Egg + m_Proton*m_Proton + m_eta*m_eta 
				//	- 2.0*m_Proton*Egg + 2.0*eb*cos(prod_th*TMath::Pi()/180.)*sqrt(Egg*Egg - m_eta*m_eta);
				double sinth = pggt / sqrt(pow(pggt,2.0) + pow(pggz,2.0));
				double costh = sqrt(1.0 - pow(sinth,2.0));
				
				double  mmsq = pow(ParticleMass(m_Target),2.0) + pow(m_eta,2.0) 
					- 2.0*(Egg*(eb+ParticleMass(m_Target)) 
					- eb*(sqrt(pow(Egg,2.0)-pow(m_eta,2.0))*costh + ParticleMass(m_Target)));
				
				double mmsq_const = pow(ParticleMass(m_Target),2.0) + pow(m_eta,2.0) 
					- 2.0*(eeta*(eb+ParticleMass(m_Target)) 
					- eb*(sqrt(pow(eeta,2.0)-pow(m_eta,2.0))*costh + ParticleMass(m_Target)));
				
				double mmsq_const_coh = pow(ParticleMass(m_Target),2.0) + pow(m_eta,2.0) 
					- 2.0*(eeta_coh*(eb+ParticleMass(m_Target)) 
					- eb*(sqrt(pow(eeta_coh,2.0)-pow(m_eta,2.0))*costh + ParticleMass(m_Target)));
				
				//-----------------------------------------------------//
				// Default Cuts
				
				// plot elasticity vs. mass ratio:
				if(loc_veto_options[5]) {
					
					h_elas_vs_mgg->Fill(invmass/m_eta, Egg/eeta, fill_weight);
					
					// plot the elasticity distribution for events where the mass is close to the eta:
					if(eta_cut) {
						h_elas->Fill(prod_th, Egg/eb, fill_weight);
						h_elas_corr->Fill(prod_th, Egg/eeta, fill_weight);
						
						if(fill_weight==1.0) h_elas_corr_main->Fill(prod_th, Egg/eeta);
						else                 h_elas_corr_side->Fill(prod_th, Egg/eeta, -1.0*fill_weight);
					}
					
					// apply elasticity cut and plot the invariant mass distriubtion:
					if(elas_cut) {
						// invariant mass vs. polar angle:
						h_mgg->Fill(prod_th, invmass, fill_weight);
						
						// energy-constrained invariant mass vs. polar angle:
						h_mgg_const->Fill(prod_th, invmass_const, fill_weight);
						
						// for monitoring purposes, plot the invariant mass separately for 
						//   beam photons in the main RF bunch and for beam photons in the accidental sidebands:
						
						if(fill_weight==1.0) {
							h_mgg_main->Fill(prod_th, invmass);
							h_mgg_const_main->Fill(prod_th, invmass_const);
						} else {
							h_mgg_side->Fill(prod_th, invmass, -1.0*fill_weight);
							h_mgg_const_side->Fill(prod_th, invmass_const, -1.0*fill_weight);
						}
						
						// energy-constrained invariant mass vs. energy-constrained polar angle:
						h_mgg_const_corr->Fill(prod_th_const, invmass_const, fill_weight);
						
						// for MC, plot invariant mass vs thrown information:
						if(is_mc) {
							h_mgg_thrown->Fill(thrown_eta_angle, invmass, fill_weight);
							h_mgg_const_thrown->Fill(thrown_eta_angle, invmass_const, fill_weight);
							
							// plot reconstructed vs. thrown angle:
							if(0.5078 < invmass_const && invmass_const < 0.5978) {
								h_rec_vs_thrown->Fill(thrown_eta_angle, prod_th, fill_weight);
							}
						}
						
						if(eta_cut) {
							// plot x-y distribution of showers:
							h_xy_1->Fill(pos1.X(), pos1.Y(), fill_weight);
							h_xy_2->Fill(pos2.X(), pos2.Y(), fill_weight);
						}
					}
				}
				
				h_mm_vs_theta->Fill(prod_th, mmsq, fill_weight);
				h_mm_const_vs_theta->Fill(prod_th, mmsq_const, fill_weight);
				h_mm_const_coh_vs_theta->Fill(prod_th, mmsq_const_coh, fill_weight);
				if(eta_cut) {
					h_mm_vs_theta_eta_cut->Fill(prod_th, mmsq, fill_weight);
					h_mm_const_vs_theta_eta_cut->Fill(prod_th, mmsq_const, fill_weight);
					if(elas_cut) {
						h_mm_vs_theta_eta_elas_cut->Fill(prod_th, mmsq, fill_weight);
						h_mm_const_vs_theta_eta_elas_cut->Fill(prod_th, mmsq_const, fill_weight);
					}
				}
				
				for(int iveto=0; iveto<m_n_vetos; iveto++) {
					if(loc_veto_options[iveto]) {
						if(eta_cut) h_elas_veto[iveto]->Fill(prod_th, Egg/eeta, fill_weight);
						if(elas_cut) {
							h_mgg_veto[iveto]->Fill(prod_th, invmass, fill_weight);
							h_mgg_const_veto[iveto]->Fill(prod_th, invmass_const, fill_weight);
							h_mgg_const_coh_veto[iveto]->Fill(prod_th, invmass_const_coh, fill_weight);
							if(eta_cut) {
								h_mm_veto[iveto]->Fill(prod_th, mmsq, fill_weight);
								h_mm_const_veto[iveto]->Fill(prod_th, mmsq_const, fill_weight);
								h_mm_const_coh_veto[iveto]->Fill(prod_th, mmsq_const_coh, fill_weight);
							}
						}
					}
				}
				
				/*
				// hybrid mass:
				h_hmass->Fill(prod_th, hmass, fill_weight);
				
				// missing mass:
				h_mm_vs_theta->Fill(prod_th, mmsq, fill_weight);
				
				// plot the above two distributions with a cut around the eta mass:
				
				if(eta_cut) {
					h_hmass_eta_cut->Fill(prod_th, hmass, fill_weight);
					h_mm_vs_theta_eta_cut->Fill(prod_th, mmsq, fill_weight);
					
					// next, apply elasticity cut:
					if(elas_cut) {
						h_hmass_eta_elas_cut->Fill(prod_th, hmass, fill_weight);
						h_mm_vs_theta_eta_elas_cut->Fill(prod_th, mmsq, fill_weight);
					}
				}
				*/
				
			} // loop over DBeamPhotons
		} // inner loop over DFCALShowers
	} // outer loop over DFCALShowers
	
	
	return;
}

int JEventProcessor_primex_eta_analysis::fcal_fiducial_cut(DVector3 pos, double layer_cut) 
{
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = (1.5 + layer_cut) * 4.0157; // 4.0157 cm is the size of one FCAL block
	
	double fcal_face_x = m_beamSpot.X() + (pos.X() * (m_fcalFace.Z() - m_beamSpot.Z())/pos.Z());
	double fcal_face_y = m_beamSpot.Y() + (pos.Y() * (m_fcalFace.Z() - m_beamSpot.Z())/pos.Z());
	
	fcal_face_x -= m_fcalFace.X();
	fcal_face_y -= m_fcalFace.Y();
	
	if((-1.*fcal_inner_layer_cut < fcal_face_x && fcal_face_x < fcal_inner_layer_cut)
		&& (-1.*fcal_inner_layer_cut < fcal_face_y 
		&& fcal_face_y < fcal_inner_layer_cut)) fid_cut = 1;
	
	// apply fiducial cut for bad channels in phase 1:
	/*
	if(phase_val==1) {
		if(-32. < fcal_face_y && fcal_face_y < -21. 
			&& -18. < fcal_face_x && fcal_face_x < 6.) fid_cut = 1;
		
		if(15. < fcal_face_y && fcal_face_y < 25. 
			&& 48. < fcal_face_x && fcal_face_x < 58.) fid_cut = 1;
		
		if(-29. < fcal_face_y && fcal_face_y < -19. 
			&& 61. < fcal_face_x && fcal_face_x < 71.) fid_cut = 1;
		
		if(65. < fcal_face_y && fcal_face_y < 75. 
			&& 24. < fcal_face_x && fcal_face_x < 34.) fid_cut = 1;
	}
	*/
	return fid_cut;
}

void JEventProcessor_primex_eta_analysis::check_TOF_match(DVector3 pos1, double rfTime, 
	vector<const DTOFPoint*> tof_points, double &dx_min, double &dy_min, 
	double &dt_min, double rf_time_cut) {
	
	dx_min = 1000.;
	dy_min = 1000.;
	dt_min = 1000.;
	
	for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
		tof != tof_points.end(); tof++) {
		
		double xt = (*tof)->pos.X() - m_beamSpot.X();
		double yt = (*tof)->pos.Y() - m_beamSpot.Y();
		double zt = (*tof)->pos.Z() - m_beamSpot.Z();
		double rt = sqrt(xt*xt + yt*yt + zt*zt);
		double tt = (*tof)->t - (rt/m_c);
		double dt = tt - rfTime;
		xt *= pos1.Z() / zt;
		yt *= pos1.Z() / zt;
		double dx = pos1.X() - xt;
		double dy = pos1.Y() - yt;
		
		if(fabs(dt) < rf_time_cut) {
			if((dx*dx + dy*dy) < (dx_min*dx_min + dy_min*dy_min)) {
				dx_min = dx;
				dy_min = dy;
				dt_min = dt;
			}
		}
	}
	
	return;
}

double JEventProcessor_primex_eta_analysis::energy_after_recoil(double eb, double theta, 
	double m0, double mp) 
{
	theta *= (TMath::Pi()/180.);
  
	double t1 = eb*cos(theta);
	double t2 = mp+eb;
	double t3 = mp*eb + m0*m0*0.5;
	
	double a = t1*t1-t2*t2;
	double b = 2.*t2*t3;
	double c = -m0*m0*t1*t1-t3*t3;
	double d = b*b - 4.*a*c;
	
	if(d < 0. || a == 0.) {
		cout << "IMAGINARY ETA ENERGY!!!" << endl;
		return 0.;
	}
	
	double energy = (-b-sqrt(d))/2./a;
	return energy;
}

double JEventProcessor_primex_eta_analysis::fcal_energy_res(double e)
{
	// hard-coded values for the FCAL energy resolution (taken from GlueX NIM paper)
	
	double a = 0.062, b = 0.047;
	double sig = (a*a)/e + (b*b);
	sig = sqrt(sig) * e;
	return sig;
}

double JEventProcessor_primex_eta_analysis::get_acc_scaling_factor(double eb)
{
	if(eb > m_TAGMEnergyBoundHi)
		return m_HodoscopeHiFactor;
	else if(eb > m_TAGMEnergyBoundLo)
		return m_MicroscopeFactor;
	else
		return m_HodoscopeLoFactor;
}

