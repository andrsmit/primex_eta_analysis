// $Id$
//
//    File: JEventProcessor_eta_gg_tree.cc
// Created: Fri Aug 11 14:26:44 EDT 2023
// Creator: andrsmit (on Linux ifarm1802.jlab.org 3.10.0-1160.92.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_eta_gg_tree.h"

extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_eta_gg_tree());
	}
} // "C"

thread_local DTreeFillData JEventProcessor_eta_gg_tree::dTreeFillData;

//------------------
// JEventProcessor_eta_gg_tree (Constructor)
//------------------
JEventProcessor_eta_gg_tree::JEventProcessor_eta_gg_tree()
{
	/*
	Basic criteria for selecting events to write out:
	  There is at least 1 pair of FCAL showers that have an invariant mass greater than 0.3 within 3GeV of tagged photon
	
	For all events that satisfy that cut, write out:
		- All FCAL showers within RF timing cut
		- All BCAL showers within +/-12ns of RF
		- All Beam photons in main RF bunch and selected sidebands
		- All TOF points within +/-2ns of RF
		- All DMCThrown information
	*/
	
	// default values for the RF timing cuts for each sub-detector:
	m_FCAL_RF_CUT =  2.004; // [ns]
	m_BEAM_RF_CUT =  2.004; // [ns]
	m_TOF_RF_CUT  =  2.004; // [ns]
	m_BCAL_RF_CUT = 12.0;   // [ns]
	
	m_MIN_BEAM_ENERGY = 8.0; // energy of the tagged photon energy [GeV]
	m_DELTA_E_CUT     = 3.0; // energy difference between FCAL showers and coherent eta meson [GeV]
		
	// miscellaneous:
	m_USE_LOG_WEIGHT = 0; // use log-weighted FCAL position
	m_SAVE_MC_NOHITS = 1; // save MCThrown information even when an event doesn't pass the selection criteria
	
	//-------------------------------------------------------------------------------------//
	// allow for command-line overriding of the default values:
	
	gPARMS->SetDefaultParameter("eta_gg_tree:FCAL_RF_CUT",     m_FCAL_RF_CUT);
	gPARMS->SetDefaultParameter("eta_gg_tree:BCAL_RF_CUT",     m_BCAL_RF_CUT);
	gPARMS->SetDefaultParameter("eta_gg_tree:TOF_RF_CUT",      m_TOF_RF_CUT);
	gPARMS->SetDefaultParameter("eta_gg_tree:BEAM_RF_CUT",     m_BEAM_RF_CUT);
	gPARMS->SetDefaultParameter("eta_gg_tree:MIN_BEAM_ENERGY", m_MIN_BEAM_ENERGY);
	gPARMS->SetDefaultParameter("eta_gg_tree:DELTA_E_CUT",     m_DELTA_E_CUT);
	gPARMS->SetDefaultParameter("eta_gg_tree:USE_LOG_WEIGHT",  m_USE_LOG_WEIGHT);
	gPARMS->SetDefaultParameter("eta_gg_tree:SAVE_MC_NOHITS",  m_SAVE_MC_NOHITS);
}

//------------------
// init
//------------------
jerror_t JEventProcessor_eta_gg_tree::init(void)
{
	dTreeInterface = DTreeInterface::Create_DTreeInterface("eta_gg","eta_gg.root");
	
	DTreeBranchRegister locTreeBranchRegister;
	
	locTreeBranchRegister.Register_Single<Int_t>("eventNum");
	locTreeBranchRegister.Register_Single<Double_t>("rfTime");
	
	// Beam Photon:
	locTreeBranchRegister.Register_Single<Int_t>("nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_counter","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("tag_sys","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("beam_e","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("beam_t","nbeam");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("acc_scale_factor","nbeam");
	
	// FCAL Showers:
	locTreeBranchRegister.Register_Single<Int_t>("nfcal");
	// reconstructed quantities:
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_e","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_x","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_y","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_z","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_t","nfcal");
	// info about shower:
	locTreeBranchRegister.Register_FundamentalArray<Int_t>(   "fcal_nblocks","nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_e1e9",   "nfcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("fcal_e9e25",  "nfcal");
	
	// BCAL Showers:
	locTreeBranchRegister.Register_Single<Int_t>("nbcal");
	// reconstructed quantities:
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("bcal_e", "nbcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("bcal_x", "nbcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("bcal_y", "nbcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("bcal_z", "nbcal");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("bcal_t", "nbcal");
	
	// TOF Points:
	locTreeBranchRegister.Register_Single<Int_t>("ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_x","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_y","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_z","ntof");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("tof_t","ntof");
	
	// SC Hits:
	locTreeBranchRegister.Register_Single<Int_t>("nsc");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("sc_sector","nsc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("sc_phi","nsc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("sc_dE","nsc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("sc_t","nsc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("sc_pulse_height","nsc");
	
	// MC Thrown:
	locTreeBranchRegister.Register_Single<Int_t>("nmc");
	locTreeBranchRegister.Register_FundamentalArray<Int_t>("mc_pdgtype","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_x","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_y","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_z","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_t","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_e","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_p","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_theta","nmc");
	locTreeBranchRegister.Register_FundamentalArray<Double_t>("mc_phi","nmc");
	
	// MC Reaction:
	locTreeBranchRegister.Register_Single<Double_t>("thrownBeamEnergy");
	
	dTreeInterface->Create_Branches(locTreeBranchRegister);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_eta_gg_tree::brun(JEventLoop *eventLoop, int32_t runnumber)
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
jerror_t JEventProcessor_eta_gg_tree::evnt(JEventLoop *eventLoop, uint64_t eventnumber)
{
	vector<const DMCThrown*> locMCThrown;	
	eventLoop->Get(locMCThrown, "Primary");
	
	int locIsMC = 0;
	if(locMCThrown.size()) locIsMC = 1;
	
	const DMCReaction *locMCReaction = NULL;
	try {
		eventLoop->GetSingle(locMCReaction);
	} catch(...) {
		if(locIsMC) return NOERROR;
	}
	
	//-----------------------------------------------------//
	// Get RF Time
	
	const DEventRFBunch *locRFBunch = NULL;
	try {
		eventLoop->GetSingle(locRFBunch, "CalorimeterOnly");
	} catch(...) {
		if(locIsMC && m_SAVE_MC_NOHITS) {
			write_events(eventnumber, 0.0, locMCThrown, locMCReaction);
			dTreeInterface->Fill(dTreeFillData);
		}
		return NOERROR;
	}
	double locRFTime = locRFBunch->dTime;
	if(locRFBunch->dNumParticleVotes < 2) {
		if(locIsMC && m_SAVE_MC_NOHITS) {
			write_events(eventnumber, 0.0, locMCThrown, locMCReaction);
			dTreeInterface->Fill(dTreeFillData);
		}
		return NOERROR;
	}
	
	//-----------------------------------------------------//
	// Data objects
	
	vector<const DBeamPhoton*> locBeamPhotons;
	eventLoop->Get(locBeamPhotons);
	
	vector<const DFCALShower*> locFCALShowers;
	eventLoop->Get(locFCALShowers);
	
	vector<const DCCALShower*> locCCALShowers;
	eventLoop->Get(locCCALShowers);
	
	vector<const DBCALShower*> locBCALShowers;
	eventLoop->Get(locBCALShowers);
	
	vector<const DTOFPoint*> locTOFPoints;
	eventLoop->Get(locTOFPoints);
	
	vector<const DSCHit*> locSCHits;
	eventLoop->Get(locSCHits);
	
	//-----------------------------------------------------//
	// Look for pair of FCAL showers with invariant mass > 0.3 GeV/c^2:
	
	int locNFCALShowers = locFCALShowers.size();
	
	bool locEventSelector = false;
	
	for(int ishow=0; ishow<(locNFCALShowers-1); ishow++) {
		
		const DFCALShower *show1 = locFCALShowers[ishow];
		DVector3 pos1;
		if(m_USE_LOG_WEIGHT) {
			pos1 = show1->getPosition_log();
		} else {
			pos1 = show1->getPosition();
		}
		pos1 = pos1 - m_beamSpot + m_fcal_correction;
		
		double t1 = show1->getTime() - (pos1.Mag()/m_c) - locRFTime;
		if(fabs(t1) > m_FCAL_RF_CUT) continue;
		
		double  e1 = show1->getEnergy();
		double px1 = e1*pos1.X()/pos1.Mag();
		double py1 = e1*pos1.Y()/pos1.Mag();
		double pz1 = e1*pos1.Z()/pos1.Mag();
		
		for(int jshow=ishow+1; jshow<locNFCALShowers; jshow++) {
			
			const DFCALShower *show2 = locFCALShowers[jshow];
			DVector3 pos2;
			if(m_USE_LOG_WEIGHT) {
				pos2 = show2->getPosition_log();
			} else {
				pos2 = show2->getPosition();
			}
			pos2 = pos2 - m_beamSpot + m_fcal_correction;
			
			double t2 = show2->getTime() - (pos2.Mag()/m_c) - locRFTime;
			if(fabs(t2) > m_FCAL_RF_CUT) continue;
			
			double  e2 = show2->getEnergy();
			double px2 = e2*pos2.X()/pos2.Mag();
			double py2 = e2*pos2.Y()/pos2.Mag();
			double pz2 = e2*pos2.Z()/pos2.Mag();
			
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
			//double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			//if(invmass < 0.3) continue;
			//
			// 11.27.2024: 
			// I commented out the above two lines because cutting on the invariant mass 
			// and then plotting the energy-constrained invariant mass could lead to 
			// artificial structures to the left of the eta peak that would affect the background fit.
			
			for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin(); gam != locBeamPhotons.end(); gam++) {
				
				double eb = (*gam)->lorentzMomentum().E();
				double brfdt = (*gam)->time() - locRFTime;
				
				// remove beam photons below the minimum energy cut:
				if(eb < m_MIN_BEAM_ENERGY) continue;
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// Calculate the energy of the eta meson, assuming production on a free nucleon:
				double eeta = energy_after_recoil(eb, prod_th, m_eta, m_Proton);
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a photo-produced eta meson:
				
				double sig1 = fcal_energy_res(e1);
				double sig2 = fcal_energy_res(e2);
				double sigr = pow(sig1/sig2,2.0);
				
				// Energy-constrained invariant mass assuming production on free nucleon:
				double e1c = e1/(1.+sigr) + (eeta-e2)/(1.+(1./sigr));
				double e2c = eeta - e1c;
				double invmass_const = sqrt(2.*e1c*e2c*(1.-cos12));
				
				//-----------------------------------------------------//
				
				// only save events where there is a beam photon passing the cuts in the main RF bunch 
				// OR in one of the accidental sidebands that will be used for subtraction. 
				// NOTE: we need to use the exact same sidebands when analyzing the trees and doing the 
				// accidental subtraction that are used in these cuts.
				
				int bunch_val = 0;
				if(brfdt < (m_beam_bunches_main*2.004)) bunch_val = 1;
				else if(
					((m_beam_bunches_main+5.5)*4.008 < brfdt) &&
					(brfdt < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
				) bunch_val = 2;
				
				if(bunch_val==0) continue;
				
				double locDeltaE = Egg - eeta;
				if(fabs(locDeltaE) < 3.0 && invmass_const > 0.30) {
					locEventSelector = true;
				}
			}
		}
	}
	
	if(locEventSelector) {
		write_events(eventnumber, locRFTime, locBeamPhotons, locFCALShowers, locBCALShowers, locTOFPoints, locSCHits, 
			locMCThrown, locMCReaction);
		dTreeInterface->Fill(dTreeFillData);
	}
	else if(locIsMC && m_SAVE_MC_NOHITS) {
		write_events(eventnumber, locRFTime, locMCThrown, locMCReaction);
		dTreeInterface->Fill(dTreeFillData);
	}
	
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_eta_gg_tree::fini(void)
{
	delete dTreeInterface;
	return NOERROR;
}

double JEventProcessor_eta_gg_tree::energy_after_recoil(double eb, double theta, 
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

double JEventProcessor_eta_gg_tree::fcal_energy_res(double e)
{
	// hard-coded values for the FCAL energy resolution (taken from GlueX NIM paper)
	
	double a = 0.062, b = 0.047;
	double sig = (a*a)/e + (b*b);
	sig = sqrt(sig) * e;
	return sig;
}

double JEventProcessor_eta_gg_tree::get_acc_scaling_factor(double eb)
{
	if(eb > m_TAGMEnergyBoundHi)
		return m_HodoscopeHiFactor;
	else if(eb > m_TAGMEnergyBoundLo)
		return m_MicroscopeFactor;
	else
		return m_HodoscopeLoFactor;
}

void JEventProcessor_eta_gg_tree::write_events(uint64_t eventnumber, double rfTime, 
	vector<const DMCThrown*> mc_thrown, const DMCReaction *mc_reaction) {
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	dTreeFillData.Fill_Single<Int_t>("nbeam", 0);
	dTreeFillData.Fill_Single<Int_t>("nfcal", 0);
	dTreeFillData.Fill_Single<Int_t>("nbcal", 0);
	dTreeFillData.Fill_Single<Int_t>("ntof",  0);
	
	// MC Thrown:
	size_t n_mc_thrown = 0;
	for(vector<const DMCThrown*>::const_iterator mc = mc_thrown.begin(); mc != mc_thrown.end(); mc++) {
		
		dTreeFillData.Fill_Array<Int_t>("mc_pdgtype",    (*mc)->pdgtype, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_x",       (*mc)->position().X(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_y",       (*mc)->position().Y(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_z",       (*mc)->position().Z(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_t",       (*mc)->time(),                        n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_e",       (*mc)->energy(),                      n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_p",       (*mc)->momentum().Mag(),              n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_theta",   (*mc)->momentum().Theta()*180.0/M_PI, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_phi",     (*mc)->momentum().Phi()*180.0/M_PI,   n_mc_thrown);
		
		n_mc_thrown++;
	}
	dTreeFillData.Fill_Single<Int_t>("nmc", n_mc_thrown);
	
	dTreeFillData.Fill_Single<Double_t>("thrownBeamEnergy", mc_reaction->beam.energy());
	
	return;
}

void JEventProcessor_eta_gg_tree::write_events(uint64_t eventnumber, double rfTime, 
	vector<const DBeamPhoton*> beam_photons, 
	vector<const DFCALShower*> fcal_showers,
	vector<const DBCALShower*> bcal_showers,
	vector<const DTOFPoint*> tof_points,
	vector<const DSCHit*> sc_hits,
	vector<const DMCThrown*> mc_thrown, const DMCReaction* mc_reaction) {
	
	dTreeFillData.Fill_Single<Int_t>("eventNum", eventnumber);
	dTreeFillData.Fill_Single<Double_t>("rfTime", rfTime);
	
	// Beam Photons:
	size_t n_beam_photon = 0;
	for(vector<const DBeamPhoton*>::const_iterator gam = beam_photons.begin();
		gam != beam_photons.end(); gam++) {
		
		// Only write out beam photons in the prompt RF peak, and the sidebands used for accidental subtraction:
		//
		// NOTE: we need to use the exact same sidebands when analyzing the trees and doing the 
		// accidental subtraction that are used in these cuts.
		
		double brfdt = fabs((*gam)->time() - rfTime);
		
		int bunch_val = 0;
		if(brfdt < (m_beam_bunches_main*2.004)) bunch_val = 1;
		else if(
			((m_beam_bunches_main+5.5)*4.008 < brfdt) &&
			(brfdt < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
		) bunch_val = 2;
		
		if(bunch_val==0) continue;
		
		int loc_counter = (*gam)->dCounter;
		int loc_sys = -1;
		if((*gam)->dSystem == SYS_TAGH) loc_sys = 0;
		else if((*gam)->dSystem == SYS_TAGM) loc_sys = 1;
		
		double loc_eb = (*gam)->lorentzMomentum().E();
		dTreeFillData.Fill_Array<Int_t>("tag_counter", loc_counter,    n_beam_photon);
		dTreeFillData.Fill_Array<Int_t>("tag_sys",     loc_sys,        n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("beam_e",   loc_eb,         n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("beam_t",  (*gam)->time(),  n_beam_photon);
		dTreeFillData.Fill_Array<Double_t>("acc_scale_factor", get_acc_scaling_factor(loc_eb), n_beam_photon);
		
		n_beam_photon++;
	}
	dTreeFillData.Fill_Single<Int_t>("nbeam", n_beam_photon);
	
	// FCAL Showers:
	size_t n_fcal_shower = 0;
	for(vector<const DFCALShower*>::const_iterator show = fcal_showers.begin(); 
		show != fcal_showers.end(); show++) {
		
		dTreeFillData.Fill_Array<Double_t>("fcal_e", (*show)->getEnergy(), n_fcal_shower);
		DVector3 locPos = m_USE_LOG_WEIGHT ? (*show)->getPosition_log() : (*show)->getPosition();
		dTreeFillData.Fill_Array<Double_t>("fcal_x", locPos.X(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_y", locPos.Y(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_z", locPos.Z(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_t", (*show)->getTime(),  n_fcal_shower);
		
		dTreeFillData.Fill_Array<Int_t>(   "fcal_nblocks", (*show)->getNumBlocks(), n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_e1e9",    (*show)->getE1E9(),      n_fcal_shower);
		dTreeFillData.Fill_Array<Double_t>("fcal_e9e25",   (*show)->getE9E25(),     n_fcal_shower);
		
		n_fcal_shower++;
	}
	dTreeFillData.Fill_Single<Int_t>("nfcal", n_fcal_shower);
	
	// BCAL Showers:
	size_t n_bcal_shower = 0;
	for(vector<const DBCALShower*>::const_iterator show = bcal_showers.begin(); 
		show != bcal_showers.end(); show++) {
		
		dTreeFillData.Fill_Array<Double_t>("bcal_e", (*show)->E, n_bcal_shower);
		dTreeFillData.Fill_Array<Double_t>("bcal_x", (*show)->x, n_bcal_shower);
		dTreeFillData.Fill_Array<Double_t>("bcal_y", (*show)->y, n_bcal_shower);
		dTreeFillData.Fill_Array<Double_t>("bcal_z", (*show)->z, n_bcal_shower);
		dTreeFillData.Fill_Array<Double_t>("bcal_t", (*show)->t, n_bcal_shower);
		
		n_bcal_shower++;
	}
	dTreeFillData.Fill_Single<Int_t>("nbcal", n_bcal_shower);
	
	// TOF Points:
	size_t n_tof_points = 0;
	for(vector<const DTOFPoint*>::const_iterator tof = tof_points.begin(); 
		tof != tof_points.end(); tof++) {
		
		dTreeFillData.Fill_Array<Double_t>("tof_x", (*tof)->pos.X(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("tof_y", (*tof)->pos.Y(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("tof_z", (*tof)->pos.Z(), n_tof_points);
		dTreeFillData.Fill_Array<Double_t>("tof_t", (*tof)->t,       n_tof_points);
		
		n_tof_points++;
	}
	dTreeFillData.Fill_Single<Int_t>("ntof", n_tof_points);
	
	// SC Hits:
	size_t n_sc_hits = 0;
	for(vector<const DSCHit*>::const_iterator sc = sc_hits.begin(); sc != sc_hits.end(); sc++) {
		
		int sector = (*sc)->sector;
		double phi = m_sc_pos[sector-1][0].Phi() * (180./TMath::Pi());
		
		dTreeFillData.Fill_Array<Int_t>("sc_sector",          sector,              n_sc_hits);
		dTreeFillData.Fill_Array<Double_t>("sc_phi",          phi,                 n_sc_hits);
		dTreeFillData.Fill_Array<Double_t>("sc_dE",           (*sc)->dE,           n_sc_hits);
		dTreeFillData.Fill_Array<Double_t>("sc_t",            (*sc)->t,            n_sc_hits);
		dTreeFillData.Fill_Array<Double_t>("sc_pulse_height", (*sc)->pulse_height, n_sc_hits);
		
		n_sc_hits++;
	}
	dTreeFillData.Fill_Single<Int_t>("nsc", n_sc_hits);
	
	// MC Thrown:
	size_t n_mc_thrown = 0;
	for(vector<const DMCThrown*>::const_iterator mc = mc_thrown.begin(); mc != mc_thrown.end(); mc++) {
		
		dTreeFillData.Fill_Array<Int_t>("mc_pdgtype",    (*mc)->pdgtype, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_x",       (*mc)->position().X(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_y",       (*mc)->position().Y(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_z",       (*mc)->position().Z(),                n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_t",       (*mc)->time(),                        n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_e",       (*mc)->energy(),                      n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_p",       (*mc)->momentum().Mag(),              n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_theta",   (*mc)->momentum().Theta()*180.0/M_PI, n_mc_thrown);
		dTreeFillData.Fill_Array<Double_t>("mc_phi",     (*mc)->momentum().Phi()*180.0/M_PI,   n_mc_thrown);
		
		n_mc_thrown++;
	}
	dTreeFillData.Fill_Single<Int_t>("nmc", n_mc_thrown);
	
	if(n_mc_thrown>0) {
		dTreeFillData.Fill_Single<Double_t>("thrownBeamEnergy", mc_reaction->beam.energy());
	}
	
	return;
}
