// $Id$
//
//    File: JEventProcessor_primex_eta_analysis_FCAL.cc
// Created: Fri Aug 11 14:26:44 EDT 2023
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.92.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_primex_eta_analysis_FCAL.h"

extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_primex_eta_analysis_FCAL());
	}
} // "C"

//------------------
// JEventProcessor_primex_eta_analysis_FCAL (Constructor)
//------------------
JEventProcessor_primex_eta_analysis_FCAL::JEventProcessor_primex_eta_analysis_FCAL()
{
	//-------------------------------------------------------------------------------------//
	// initialize cut vectors:
	
	// FCAL fiducial cuts:
	m_fiducial_cuts.clear();
	for(int icut=0; icut<13; icut++) {
		double loc_cut     = 0.0 + 0.25*(double)(icut);
		m_fiducial_cuts.push_back(loc_cut);
	}
	
	//-------------------------------------------------------------------------------------//
	// initialize default cut values:
	
	// default values for the RF timing cuts for each sub-detector:
	m_BEAM_RF_CUT =  2.004;
	m_FCAL_RF_CUT =  2.0;
	m_BCAL_RF_CUT = 10.0;
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
	m_USE_LOG_WEIGHT = 1; // use log-weighted FCAL position
	m_BYPASS_TRIGGER = 0; // determines whether or not to check the trigger bits set for each event
	
	//-------------------------------------------------------------------------------------//
	// allow for command-line overriding of the default values:
	
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:FCAL_RF_CUT",     m_FCAL_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:BEAM_RF_CUT",     m_BEAM_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:BCAL_RF_CUT",     m_BCAL_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:CCAL_RF_CUT",     m_CCAL_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:TOF_RF_CUT",      m_TOF_RF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:MIN_FCAL_ENERGY", m_MIN_FCAL_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:MIN_BEAM_ENERGY", m_MIN_BEAM_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:MIN_BCAL_ENERGY", m_MIN_BCAL_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:MIN_CCAL_ENERGY", m_MIN_CCAL_ENERGY);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:FCAL_TOF_CUT",    m_FCAL_TOF_CUT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:ELAS_CUT_SIGMA",  m_ELAS_CUT_SIGMA);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:ELAS_CUT_WIDTH",  m_ELAS_CUT_WIDTH);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:ELAS_CUT_MU_P0",  m_ELAS_CUT_MU_P0);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:ELAS_CUT_MU_P1",  m_ELAS_CUT_MU_P1);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:USE_LOG_WEIGHT",  m_USE_LOG_WEIGHT);
	gPARMS->SetDefaultParameter("primex_eta_analysis_FCAL:BYPASS_TRIGGER",  m_BYPASS_TRIGGER);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_primex_eta_analysis_FCAL::init(void)
{
	TDirectory *dir_primex_eta = new TDirectoryFile("primex_eta_analysis_FCAL", "primex_eta_analysis_FCAL");
	dir_primex_eta->cd();
	
	// invariant mass (with and without constraint) without FCAL cuts:
	
	h_mgg = new TH2F("mgg", "No FCAL multiplicity cut", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg->Sumw2();
	h_mgg->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	
	h_mgg_const = new TH2F("mgg_const", "No FCAL multiplicity cut", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const->Sumw2();
	h_mgg_const->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_const->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	
	// with minimum energy cuts:
	
	h_mgg_ecut = new TH2F("mgg_ecut", Form("E_{1,2} > %.1f",m_MIN_FCAL_ENERGY), 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_ecut->Sumw2();
	h_mgg_ecut->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_ecut->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	
	h_mgg_const_ecut = new TH2F("mgg_const_ecut", 
		Form("E_{1,2} > %.1f",m_MIN_FCAL_ENERGY), 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_ecut->Sumw2();
	h_mgg_const_ecut->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_const_ecut->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	
	// with 'good' multiplicity = 2:
	
	h_mgg_good_mult = new TH2F("mgg_good_mult", "2 Good FCAL showers", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_good_mult->Sumw2();
	h_mgg_good_mult->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_good_mult->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	
	h_mgg_const_good_mult = new TH2F("mgg_const_good_mult", "2 Good FCAL showers", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_good_mult->Sumw2();
	h_mgg_const_good_mult->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_const_good_mult->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	
	// with total multiplicity = 2:
	
	h_mgg_mult = new TH2F("mgg_mult", "2 FCAL showers", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_mult->Sumw2();
	h_mgg_mult->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_mult->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	
	h_mgg_const_mult = new TH2F("mgg_const_mult", "2 Good FCAL showers", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_mult->Sumw2();
	h_mgg_const_mult->GetXaxis()->SetTitle("#theta_{#eta} [#circ]");
	h_mgg_const_mult->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	
	// vary the fiducial cut:
	
	int loc_n_fiducial_cuts = m_fiducial_cuts.size();
	for(int icut=0; icut<loc_n_fiducial_cuts; icut++) {
		TH2F *loc_h_mgg = new TH2F(Form("mgg_fid_%02d", icut),
			Form("Inner layers removed by fiducial cut: %.1f", m_fiducial_cuts[icut]), 
			650, 0.0, 6.5, 600, 0.0, 1.2);
		loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{constr.} [GeV/c^{2}]");
		loc_h_mgg->Sumw2();
		h_mgg_fid.push_back(loc_h_mgg);
	}
	
	dir_primex_eta->cd("../");
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_primex_eta_analysis_FCAL::brun(JEventLoop *eventLoop, int32_t runnumber)
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
	// Set the target mass depending on run number:
	
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
jerror_t JEventProcessor_primex_eta_analysis_FCAL::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
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
		
		if(locThrownEnergy<m_MIN_BEAM_ENERGY) {
			japp->RootFillUnLock(this);
			return NOERROR;
		}
	}
	
	//-----------------------------------------------------//
	// RF Timing Histograms:
	
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
		if(fabs(loc_t) < m_FCAL_RF_CUT) {
			locFCALEnergySum += (*show)->getEnergy();
			locNFCALShowers++;
			if(!fcal_fiducial_cut(loc_pos, 2.0) && (*show)->getEnergy() > m_MIN_FCAL_ENERGY) {
				locNGoodFCALShowers++;
			}
		}
	}
	
	int    locNBCALShowers  = 0, locNBCALShowers_1ns = 0;
	double locBCALEnergySum = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	for(vector<const DBCALShower*>::const_iterator show = locDBCALShowers.begin(); 
		show != locDBCALShowers.end(); show++) {
		DVector3 loc_pos((*show)->x, (*show)->y, (*show)->z);
		loc_pos -= m_beamSpot;
		double loc_t = (*show)->t - (loc_pos.Mag()/m_c) - locRFTime;
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
		if(fabs(loc_t) < m_CCAL_RF_CUT) {
			locCCALEnergySum += (*show)->E;
			locNCCALShowers++;
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

void JEventProcessor_primex_eta_analysis_FCAL::eta_gg_analysis(
	vector<const DFCALShower*> fcal_showers, 
	vector<const DBeamPhoton*> beam_photons, 
	vector<const DBCALShower*> bcal_showers, 
	vector<const DTOFPoint*> tof_points, 
	vector<const DSCHit*> sc_hits, 
	int n_fcal_showers, int n_good_fcal_showers, 
	double bcal_energy_sum, int n_bcal_showers, int n_bcal_showers_1ns, double bcal_phi, double bcal_rfdt,
	double rfTime, bool is_mc, double thrown_beam_energy, double thrown_eta_angle)
{
	// Reject events with more than 1 shower in the BCAL:
	//if(bcal_energy_sum>m_MIN_BCAL_ENERGY || n_bcal_showers>1) return;
	
	int n_fcal_showers_total = static_cast<int>(fcal_showers.size());
	
	//----------------------------------------------------------------------------------//
	// Get FCAL shower multiplicities for various RF timing cuts:
	
	int loc_n_fiducial_cuts = m_fiducial_cuts.size();
	
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
		
		double px1 = e1*pos1.X()/pos1.Mag();
		double py1 = e1*pos1.Y()/pos1.Mag();
		double pz1 = e1*pos1.Z()/pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		check_TOF_match(pos1, rfTime, tof_points, tof_dx1, tof_dy1, tof_dt1, m_TOF_RF_CUT);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
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
			if((tof_dr1 < m_FCAL_TOF_CUT) && (tof_dr2 < m_FCAL_TOF_CUT)) continue;
			
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
			
			// opening angle:
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
			
			// invariant mass:
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			//-----------------------------------------------------//
			
			vector<int> loc_fiducial_cuts;
			for(int i=0; i<loc_n_fiducial_cuts; i++) loc_fiducial_cuts.push_back(0);
			
			if(n_fcal_showers==2 && 
				(e1 > m_MIN_FCAL_ENERGY) && (fabs(t1) < m_FCAL_RF_CUT) && 
				(e2 > m_MIN_FCAL_ENERGY) && (fabs(t2) < m_FCAL_RF_CUT))
			{
				for(int icut=0; icut<loc_n_fiducial_cuts; icut++) {
					int loc_fid_cut1 = fcal_fiducial_cut(pos1, m_fiducial_cuts[icut]);
					int loc_fid_cut2 = fcal_fiducial_cut(pos2, m_fiducial_cuts[icut]);
					if((loc_fid_cut1==0) && (loc_fid_cut2==0)) {
						loc_fiducial_cuts[icut] = 1;
					}
				}
			}
			
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
				else { continue; }
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				//double eeta = energy_after_recoil(eb, prod_th, m_eta, ParticleMass(m_Target));
				double eeta = energy_after_recoil(eb, prod_th, m_eta, m_Proton);
				
				// Apply a cut on the elasticity
				//  (ratio of measured energy of 2-photons, to the calculated energy above):
				
				bool elas_cut = false;
				double loc_elas_mean  = m_ELAS_CUT_MU_P0 + m_ELAS_CUT_MU_P1*prod_th;
				double loc_elas_width = m_ELAS_CUT_WIDTH * m_ELAS_CUT_SIGMA;
				if(fabs((Egg/eeta)-loc_elas_mean)<loc_elas_width) elas_cut = true;
				
				// apply the elasticity cut now:
				if(!elas_cut) continue;
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a coherently-produced eta meson:
				
				double sig1 = fcal_energy_res(e1);
				double sig2 = fcal_energy_res(e2);
				double sigr = pow(sig1/sig2,2.0);
				
				double e1c = e1/(1.+sigr) + (eeta-e2)/(1.+(1./sigr));
				double e2c = eeta - e1c;
				double invmass_const = sqrt(2.*e1c*e2c*(1.-cos12)); // energy-constrained invariant mass
				
				// re-compute the polar angle of the two-photon pair using these adjusted energies:
				/*
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
				*/
				//-----------------------------------------------------//
				// Vary FCAL cuts:
				
				// invariant mass (with and without energy-constrain) vs polar angle with no cuts:
				if((fabs(t1) < m_FCAL_RF_CUT) && (fabs(t2) < m_FCAL_RF_CUT)) {
					h_mgg->Fill(prod_th, invmass, fill_weight);
					h_mgg_const->Fill(prod_th, invmass, fill_weight);
					
					// next, apply minimum energy cut:
					if((e1 > m_MIN_FCAL_ENERGY) && (e2 > m_MIN_FCAL_ENERGY)) {
						h_mgg_ecut->Fill(prod_th, invmass, fill_weight);
						h_mgg_const_ecut->Fill(prod_th, invmass, fill_weight);
						
						// multiplicity cut (only 2 'good' showers):
						if(n_good_fcal_showers==2) {
							h_mgg_good_mult->Fill(prod_th, invmass, fill_weight);
							h_mgg_const_good_mult->Fill(prod_th, invmass, fill_weight);
							
							// multiplicity cut (only 2 showers total within timing cut):
							if(n_fcal_showers==2) {
								h_mgg_mult->Fill(prod_th, invmass, fill_weight);
								h_mgg_const_mult->Fill(prod_th, invmass, fill_weight);
							}
						}
					}
				}
				
				for(int icut=0; icut<loc_n_fiducial_cuts; icut++) {
					if(loc_fiducial_cuts[icut]==1) h_mgg_fid[icut]->Fill(prod_th, invmass_const, fill_weight);
				}
				
			} // loop over DBeamPhotons
		} // inner loop over DFCALShowers
	} // outer loop over DFCALShowers
	
	
	return;
}

int JEventProcessor_primex_eta_analysis_FCAL::fcal_fiducial_cut(DVector3 pos, double layer_cut) 
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

void JEventProcessor_primex_eta_analysis_FCAL::check_TOF_match(DVector3 pos1, double rfTime, 
	vector<const DTOFPoint*> tof_points, 
	double &dx_min, double &dy_min, double &dt_min, double rf_time_cut) {
	
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

double JEventProcessor_primex_eta_analysis_FCAL::energy_after_recoil(double eb, double theta, 
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

double JEventProcessor_primex_eta_analysis_FCAL::fcal_energy_res(double e)
{
	// hard-coded values for the FCAL energy resolution (taken from GlueX NIM paper)
	
	double a = 0.062, b = 0.047;
	double sig = (a*a)/e + (b*b);
	sig = sqrt(sig) * e;
	return sig;
}

double JEventProcessor_primex_eta_analysis_FCAL::get_acc_scaling_factor(double eb)
{
	if(eb > m_TAGMEnergyBoundHi)
		return m_HodoscopeHiFactor;
	else if(eb > m_TAGMEnergyBoundLo)
		return m_MicroscopeFactor;
	else
		return m_HodoscopeLoFactor;
}

