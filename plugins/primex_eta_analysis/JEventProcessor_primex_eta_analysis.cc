// $Id$
//
//    File: JEventProcessor_primex_eta_analysis.cc
// Created: Tue Mar  4 10:14:54 AM EST 2025
// Creator: andrsmit (on Linux ifarm2401.jlab.org 5.14.0-503.19.1.el9_5.x86_64 x86_64)
//

/// For more information on the syntax changes between JANA1 and JANA2, visit: https://jeffersonlab.github.io/JANA2/#/jana1to2/jana1-to-jana2

#include "JEventProcessor_primex_eta_analysis.h"


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->Add(new JEventProcessor_primex_eta_analysis());
}
} // "C"

//------------------
// Init
//------------------
void JEventProcessor_primex_eta_analysis::Init()
{
	// default values for the RF timing cuts for each sub-detector:
	m_BeamRFCut =  2.004;
	m_FCALRFCut =  2.004;
	m_BCALRFCut = 12.028;
	m_CCALRFCut =  2.004;
	m_TOFRFCut  =  1.000;
	
	// default values for the minimum energy cuts:
	m_FCALEnergyCut      =  0.25; // minimum energy of each FCAL shower used in analysis
	m_FCALExtraEnergyCut =  0.05; // maximum unused energy in FCAL
	m_BCALEnergyCut      =  0.10; // currently un-implemented 
	m_minBeamEnergyCut   =  9.00;
	m_maxBeamEnergyCut   = 10.90;
	
	// default values for spacial cuts:
	m_FCALTOFCut      =  8.0; // distance between fcal shower and closest DTOFPoint [cm]
	m_BCALDeltaPhiCut = 30.0; // [deg]
	m_SCDeltaPhiCut   = 36.0; // [deg]
	
	// default value for elasticity cut:
	m_ElasWidth    = 0.031;
	m_ElasSigmaCut = 4.0;
	m_ElasMean_p0  = 1.0; // mu = p0 + p1*E_gamma
	m_ElasMean_p1  = 0.0;
	
	// miscellaneous:
	m_UseLogWeight  = 0; // force log-weighted FCAL position (should be automatically set in DFCALShower_factory)
	m_BypassTrigger = 0; // determines whether or not to check the trigger bits set for each event
	
	//-------------------------------------------------------------------------------------//
	// get JApplication for setting default parameters:
	
	auto app = GetApplication();
	
	//-------------------------------------------------------------------------------------//
	// allow for command-line overriding of the default values:
	
	app->SetDefaultParameter("primex_eta_analysis:BEAM_RF_CUT", m_BeamRFCut);
	app->SetDefaultParameter("primex_eta_analysis:FCAL_RF_CUT", m_FCALRFCut);
	app->SetDefaultParameter("primex_eta_analysis:BCAL_RF_CUT", m_BCALRFCut);
	app->SetDefaultParameter("primex_eta_analysis:CCAL_RF_CUT", m_CCALRFCut);
	app->SetDefaultParameter("primex_eta_analysis:TOF_RF_CUT",  m_TOFRFCut);
	
	app->SetDefaultParameter("primex_eta_analysis:MIN_FCAL_ENERGY",       m_FCALEnergyCut);
	app->SetDefaultParameter("primex_eta_analysis:MAX_FCAL_EXTRA_ENERGY", m_FCALExtraEnergyCut);
	app->SetDefaultParameter("primex_eta_analysis:MIN_BCAL_ENERGY",       m_BCALEnergyCut);
	app->SetDefaultParameter("primex_eta_analysis:MIN_BEAM_ENERGY",       m_minBeamEnergyCut);
	app->SetDefaultParameter("primex_eta_analysis:MAX_BEAM_ENERGY",       m_maxBeamEnergyCut);
	
	app->SetDefaultParameter("primex_eta_analysis:FCAL_TOF_CUT", m_FCALTOFCut);
	
	app->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_WIDTH",  m_ElasWidth);
	app->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_SIGMA",  m_ElasSigmaCut);
	app->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_MU_P0",  m_ElasMean_p0);
	app->SetDefaultParameter("primex_eta_analysis:ELAS_CUT_MU_P1",  m_ElasMean_p1);
	
	app->SetDefaultParameter("primex_eta_analysis:USE_LOG_WEIGHT",  m_UseLogWeight);
	app->SetDefaultParameter("primex_eta_analysis:BYPASS_TRIGGER",  m_BypassTrigger);
	
	//-------------------------------------------------------------------------------------//
	// initialize lock service:
	
	lockService = app->GetService<JLockService>();
	
	//-------------------------------------------------------------------------------------//
	// Initialize all histograms:
	
	InitializeHistograms();
}

//------------------
// BeginRun
//------------------
void JEventProcessor_primex_eta_analysis::BeginRun(const std::shared_ptr<const JEvent> &event)
{
	//--------------------------------------------------------------//
	// Get geometry information for each run from CCDB:
	
	DGeometry *locGeometry = DEvent::GetDGeometry(event);
	
	double locTargetZ = 65.0;
	double locX, locY, locZ;
	
	if(locGeometry==nullptr) {
		cerr << "No geometry accessbile to plugin." << endl;
		exit(1);
	}
	
	// Get target position:
	locGeometry->GetTargetZ(locTargetZ);
	
	// Get position of FCAL face:
	locGeometry->GetFCALPosition(locX, locY, locZ);
	m_fcalFace.SetXYZ(locX, locY, locZ);
	
	// Get position of CCAL face:
	locGeometry->GetCCALPosition(locX, locY, locZ);
	m_ccalFace.SetXYZ(locX, locY, locZ);
	
	// Get beam spot on center of target:
	std::map<string, float> locBeamSpot;
	DEvent::GetCalib(event, "PHOTON_BEAM/beam_spot", locBeamSpot);
	m_beamSpot.SetXYZ(locBeamSpot.at("x"), locBeamSpot.at("y"), locTargetZ);
	
	// Get start counter geomety:
	locGeometry->GetStartCounterGeom(m_scPos, m_scNorm);
	
	// Set the target according to the run number:
	int32_t locRunNumber = event->GetRunNumber();
	m_Target = GetTargetType(locRunNumber);
	
	//--------------------------------------------------------------//
	// Geometry corrections for PrimEx run periods:
	
	m_phaseVal = GetPrimExPhase(locRunNumber);
	
	if(m_phaseVal==1) {
		
		m_fcalCorrection.SetXYZ(0.455 - m_fcalFace.X(), -0.032 - m_fcalFace.Y(), 0.0);
		m_fcalFace += m_fcalCorrection;
		
		m_ccalCorrection.SetXYZ(-0.082 - m_ccalFace.X(), 0.061 - m_ccalFace.Y(), 0.0);
		if(locRunNumber>=61483) m_ccalCorrection.SetY(0.051 - m_ccalFace.Y());
		m_ccalFace += m_ccalCorrection;
		
		if(locRunNumber<61483) 
		{
			m_beamSpot.SetX( 0.027);
			m_beamSpot.SetY(-0.128);
		}
		else if(locRunNumber<61774) 
		{
			m_beamSpot.SetX( 0.001);
			m_beamSpot.SetY(-0.077);
		}
		else 
		{
			m_beamSpot.SetX( 0.038);
			m_beamSpot.SetY(-0.095);
		}
	}
	else {
		m_fcalCorrection.SetXYZ(0.0, 0.0, 0.0);
		m_ccalCorrection.SetXYZ(0.0, 0.0, 0.0);
	}
	
	//--------------------------------------------------------------//
	// Obtain the scaling factors for accidental beam bunches 
	
	std::map<string, float> locScalingFactors;
	if(DEvent::GetCalib(event, "ANALYSIS/accidental_scaling_factor", locScalingFactors)) {
		m_HodoscopeHiFactor    = locScalingFactors.at("HODOSCOPE_HI_FACTOR");
		m_HodoscopeHiFactorErr = locScalingFactors.at("HODOSCOPE_HI_FACTOR_ERR");
		m_MicroscopeFactor     = locScalingFactors.at("MICROSCOPE_FACTOR");
		m_MicroscopeFactorErr  = locScalingFactors.at("MICROSCOPE_FACTOR_ERR");
		m_HodoscopeLoFactor    = locScalingFactors.at("HODOSCOPE_LO_FACTOR");
		m_HodoscopeLoFactorErr = locScalingFactors.at("HODOSCOPE_LO_FACTOR_ERR");
		m_TAGMEnergyBoundHi    = locScalingFactors.at("MICROSCOPE_ENERGY_HI");
		m_TAGMEnergyBoundLo    = locScalingFactors.at("MICROSCOPE_ENERGY_LO");
	}
	else {
		// Set up defaults if constants weren't available:
		m_HodoscopeHiFactor    = 1.00;
		m_HodoscopeHiFactorErr = 0.01;
		m_MicroscopeFactor     = 1.00;
		m_MicroscopeFactorErr  = 0.01;
		m_HodoscopeLoFactor    = 1.00;
		m_HodoscopeLoFactorErr = 0.01;
		m_TAGMEnergyBoundHi    = 8.00;
		m_TAGMEnergyBoundLo    = 9.00;
	}
	//--------------------------------------------------------------//
}

//------------------
// Process
//------------------
void JEventProcessor_primex_eta_analysis::Process(const std::shared_ptr<const JEvent> &event)
{
	//----------------------------------------------------//
	// Get the event RF Time:
	
	auto locEventRFBunches = event->Get<DEventRFBunch>("CalorimeterOnly");
	
	// return if no RF bunch is available:
	if(locEventRFBunches.empty()) return;
	
	// skip events with less than 2 votes on RF bunch:
	if(locEventRFBunches[0]->dNumParticleVotes < 2) return;
	
	double locRFTime = locEventRFBunches[0]->dTime;
	
	//----------------------------------------------------//
	// Get Data Objects:
	
	auto locBeamPhotons = event->Get<DBeamPhoton>();
	auto locFCALShowers = event->Get<DFCALShower>();
	auto locBCALShowers = event->Get<DBCALShower>();
	auto locCCALShowers = event->Get<DCCALShower>();
	auto locTOFPoints   = event->Get<DTOFPoint>();
	auto locSCHits      = event->Get<DSCHit>();
	auto locMCThrown    = event->Get<DMCThrown>();
	
	//----------------------------------------------------//
	// Get Trigger Information:
	
	bool locTrigConditions[m_nTrigs];
	for(int itrig=0; itrig<m_nTrigs; itrig++) { locTrigConditions[itrig] = false; }
	
	if(locMCThrown.size() > 0) 
	{
		locTrigConditions[0] = true;
		locTrigConditions[1] = true;
	} else if(m_BypassTrigger) {
		locTrigConditions[0] = true;
		locTrigConditions[1] = true;
	} else {
		const DL1Trigger* trig = event->GetSingle<DL1Trigger>();
		if(trig == nullptr) { return; }
		
		uint32_t trigmask    = trig->trig_mask;	
		uint32_t fp_trigmask = trig->fp_trig_mask;
		
		if(!trigmask)   return;
		if(fp_trigmask) return;
		
		if(trigmask & (1 <<  1)) locTrigConditions[0] = true; // FCAL Energy Sum
		if(trigmask & (1 <<  2)) locTrigConditions[1] = true; // FCAL Energy Sum (low-threshold, phase 3 only)
		if(trigmask & (1 <<  3)) locTrigConditions[2] = true; // PS
		if(trigmask & (1 << 10)) locTrigConditions[3] = true; // CCAL Energy Sum
	}
	
	//-----------------------------------------------------//
	// Apply fill lock for multi-threaded running:
	
	lockService->RootFillLock(this);
	
	//-----------------------------------------------------//
	// Plot thrown distributions (if MC)
	
	bool   locIsMC         = false;
	double locThrownEnergy = 0.;
	double locThrownAngle  = 0.;
	
	if(locMCThrown.size()) {
		const DMCReaction* locDMCReaction = event->GetSingle<DMCReaction>();
		
		TLorentzVector locEtaMCP4(0, 0, 0, 0);
		vector<TLorentzVector> locPhotonsMCList; locPhotonsMCList.clear();
		vector<TLorentzVector> locPipsMCList; locPipsMCList.clear();
		vector<TLorentzVector> locPimsMCList; locPimsMCList.clear();
		//vector<TLorentzVector> locPsMCList; locPsMCList.clear();
		//vector<TLorentzVector> locNsMCList; locNsMCList.clear();
		
		for(unsigned int i = 0; i < locMCThrown.size(); i++) {
			const DMCThrown *mcthrown = locMCThrown[i];
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
		locThrownEnergy = locDMCReaction->beam.energy();
		locThrownAngle  = locEtaMCP4.Theta() * TMath::RadToDeg();
		
		for(int icut=0; icut<13; icut++) {
			double beamEnergyCut = 7.6 + 0.2*(double)(icut);
			if(locThrownEnergy>=beamEnergyCut) {
				h_ThrownAngle[icut]->Fill(locThrownAngle);
			}
		}
		if((locThrownEnergy<m_minBeamEnergyCut) || (locThrownEnergy>m_maxBeamEnergyCut)) {
			lockService->RootFillUnLock(this);
			return;
		}
	}
	
	//-----------------------------------------------------//
	// RF Timing Histograms:
	
	// Beam:
	for(vector<const DBeamPhoton*>::const_iterator gam = locBeamPhotons.begin(); 
		gam != locBeamPhotons.end(); gam++) {
		double locT = (*gam)->time() - locRFTime;
		if((*gam)->dSystem==SYS_TAGH) {
			for(int itrig=0; itrig<m_nTrigs; itrig++) {
				if(locTrigConditions[itrig]) h_taghRFdt[itrig]->Fill(locT);
			}
		} else {
			for(int itrig=0; itrig<m_nTrigs; itrig++) {
				if(locTrigConditions[itrig]) h_tagmRFdt[itrig]->Fill(locT);
			}
		}
	}
	
	// FCAL:
	int    locNFCALShowers  = 0, locNGoodFCALShowers = 0;
	double locFCALEnergySum = 0.;
	for(vector<const DFCALShower*>::const_iterator show = locFCALShowers.begin(); 
		show != locFCALShowers.end(); show++) {
		
		DVector3 locPos;
		if(m_UseLogWeight) {
			locPos = (*show)->getPosition_log();
		} else {
			locPos = (*show)->getPosition();
		}
		locPos = locPos - m_beamSpot + m_fcalCorrection;
		double locT = (*show)->getTime() - (locPos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<m_nTrigs; itrig++) {
			if(locTrigConditions[itrig]) h_fcalRFdt[itrig]->Fill(locT);
		}
		if(fabs(locT) < m_FCALRFCut) {
			locFCALEnergySum += (*show)->getEnergy();
			locNFCALShowers++;
			if(!FCALFiducialCut(locPos, 2.0) && (*show)->getEnergy() > m_FCALEnergyCut) {
				locNGoodFCALShowers++;
			}
		}
	}
	for(int itrig=0; itrig<m_nTrigs; itrig++) {
		if(locTrigConditions[itrig]) h_fcalEnergySum[itrig]->Fill(locFCALEnergySum);
	}
	
	// BCAL:
	int    locNBCALShowers  = 0, locNBCALShowers_1ns = 0;
	double locBCALEnergySum = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	for(vector<const DBCALShower*>::const_iterator show = locBCALShowers.begin(); 
		show != locBCALShowers.end(); show++) {
		DVector3 locPos((*show)->x, (*show)->y, (*show)->z);
		locPos -= m_beamSpot;
		double locT = (*show)->t - (locPos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<m_nTrigs; itrig++) {
			if(locTrigConditions[itrig]) h_bcalRFdt[itrig]->Fill(locT);
		}
		if(fabs(locT) < m_BCALRFCut) {
			locBCALEnergySum += (*show)->E;
			locNBCALShowers++;
			locBCALRFDT = locT;
			locBCALPhi  = locPos.Phi() * (180./TMath::Pi());
			if(fabs(locT) < 1.0) {
				locNBCALShowers_1ns++;
			}
		}
	}
	
	// CCAL:
	int    locNCCALShowers  = 0;
	double locCCALEnergySum = 0.;
	for(vector<const DCCALShower*>::const_iterator show = locCCALShowers.begin(); 
		show != locCCALShowers.end(); show++) {
		DVector3 locPos((*show)->x1, (*show)->y1, (*show)->z);
		locPos = locPos - m_beamSpot + m_ccalCorrection;
		double locT = (*show)->time - (locPos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<m_nTrigs; itrig++) {
			if(locTrigConditions[itrig]) h_ccalRFdt[itrig]->Fill(locT);
		}
		if(fabs(locT) < m_CCALRFCut) {
			locCCALEnergySum += (*show)->E;
			locNCCALShowers++;
		}
	}
	
	// TOF:
	for(vector<const DTOFPoint*>::const_iterator tof = locTOFPoints.begin(); 
		tof != locTOFPoints.end(); tof++) {
		DVector3 locPos = (*tof)->pos - m_beamSpot;
		double locT = (*tof)->t - (locPos.Mag()/m_c) - locRFTime;
		for(int itrig=0; itrig<m_nTrigs; itrig++) {
			if(locTrigConditions[itrig]) h_tofRFdt[itrig]->Fill(locT);
		}
	}
	
	// SC:
	for(vector<const DSCHit*>::const_iterator sc = locSCHits.begin(); 
		sc != locSCHits.end(); sc++) {
		double locT = (*sc)->t - locRFTime;
		for(int itrig=0; itrig<m_nTrigs; itrig++) {
			if(locTrigConditions[itrig]) h_scRFdt[itrig]->Fill(locT);
		}
	}
	
	//-----------------------------------------------------//
	// eta->2gamma Analysis with cuts:
	
	if(locNFCALShowers>1) {
		EtaGGAnalysis(
			locFCALShowers, locBeamPhotons, locBCALShowers, locTOFPoints, locSCHits,
			locNFCALShowers, locNGoodFCALShowers, locBCALEnergySum, locNBCALShowers, locNBCALShowers_1ns, locBCALPhi, locBCALRFDT,
			locRFTime, locIsMC, locThrownEnergy, locThrownAngle
		);
	}
	
	// Release fill lock and we're done:
	lockService->RootFillUnLock(this);
	
}


void JEventProcessor_primex_eta_analysis::EtaGGAnalysis(
	vector<const DFCALShower*> fcalShowers, 
	vector<const DBeamPhoton*> beamPhotons, 
	vector<const DBCALShower*> bcalShowers, 
	vector<const DTOFPoint*> tofPoints, 
	vector<const DSCHit*> scHits, 
	int nFCALShowers, int nGoodFCALShowers, 
	double bcalEnergySum, int nBCALShowers, int nBCALShowers_1ns, double bcalPhi, double bcalRFdt,
	double rfTime, bool isMC, double thrownBeamEnergy, double thrownEtaAngle
) {
	
	// Apply multiplicity cut on the number of FCAL showers:
	if((nFCALShowers!=2) || (nGoodFCALShowers!=2)) return;
	
	int nFCALShowersTotal = static_cast<int>(fcalShowers.size());
	
	//----------------------------------------------------------------------------------//
	// Loop over all possible combinations of pairs of FCAL showers that pass the cuts:
	//  (there should only be 1 pair at most, per event)
	
	for(int ishow=0; ishow<nFCALShowersTotal; ishow++) {
		
		const DFCALShower *show1 = fcalShowers[ishow];
		DVector3 pos1;
		if(m_UseLogWeight) {
			pos1 = show1->getPosition_log();
		} else {
			pos1 = show1->getPosition();
		}
		pos1 = pos1 - m_beamSpot + m_fcalCorrection;
		
		double t1  = show1->getTime() - (pos1.Mag()/m_c) - rfTime;
		double e1  = show1->getEnergy();
		
		// apply minimum energy and RF timing cuts:
		if(e1 < m_FCALEnergyCut || fabs(t1) >= m_FCALRFCut) continue;
		
		// apply fiducial cut to remove the innermost two FCAL layers:
		if(FCALFiducialCut(pos1, 2.0)) continue;
		
		double px1 = e1*pos1.X()/pos1.Mag();
		double py1 = e1*pos1.Y()/pos1.Mag();
		double pz1 = e1*pos1.Z()/pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		CheckTOFMatch(pos1, rfTime, tofPoints, tof_dx1, tof_dy1, tof_dt1, m_TOFRFCut);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		// plot FCAL-TOF matching distributions for monitoring:
		if(nFCALShowers==2) {
			h_fcalTOFdx->Fill(tof_dx1);
			h_fcalTOFdy->Fill(tof_dy1);
			h_fcalTOFdr->Fill(tof_dr1);
			h_fcalTOFdt->Fill(t1-tof_dt1);
			if(tof_dr1 < m_FCALTOFCut) {
				h_fcalTOFdt_cut->Fill(t1-tof_dt1);
			}
		}
		
		for(int jshow=ishow+1; jshow<nFCALShowersTotal; jshow++) {
			
			const DFCALShower *show2 = fcalShowers[jshow];
			DVector3 pos2;
			if(m_UseLogWeight) {
				pos2 = show2->getPosition_log();
			} else {
				pos2 = show2->getPosition();
			}
			pos2 = pos2 - m_beamSpot + m_fcalCorrection;
			
			double t2  = show2->getTime() - (pos2.Mag()/m_c) - rfTime;
			double e2  = show2->getEnergy();
			
			// apply minimum energy and RF timing cuts:
			if(e2 < m_FCALEnergyCut || fabs(t2) >= m_FCALRFCut) continue;
			
			// apply fiducial cut to remove the innermost two FCAL layers:
			if(FCALFiducialCut(pos2, 2.0)) continue;
			
			double px2 = e2*pos2.X()/pos2.Mag();
			double py2 = e2*pos2.Y()/pos2.Mag();
			double pz2 = e2*pos2.Z()/pos2.Mag();
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			CheckTOFMatch(pos2, rfTime, tofPoints, tof_dx2, tof_dy2, tof_dt2, m_TOFRFCut);
			double tof_dr2 = sqrt(pow(tof_dx2,2.0)+pow(tof_dy2,2.0));
			
			//-----------------------------------------------------//
			// TOF Veto
			
			// reject combinations of FCAL showers where both showers are near a TOF hit:
			bool isTOFVeto = false;
			if((tof_dr1 < m_FCALTOFCut) && (tof_dr2 < m_FCALTOFCut)) isTOFVeto = true;
			
			// count the number of FCAL-TOF matches for monitoring:
			int nTOFMatches = 0;
			if(tof_dr1 < m_FCALTOFCut) nTOFMatches++;
			if(tof_dr2 < m_FCALTOFCut) nTOFMatches++;
			h_fcalTOFMatches->Fill(nTOFMatches);
			
			if(isTOFVeto) continue;
			
			//-----------------------------------------------------//
			// Two-Photon kinematics:
			
			double Egg  =  e1 +  e2; // energy of 2-photon pair
			double pggx = px1 + px2; // momentum along x-axis
			double pggy = py1 + py2; // momentum along y-axis
			double pggz = pz1 + pz2; // momentum along z-axis
			
			// transverse momentum:
			double pggt = sqrt(pow(pggx,2.0) + pow(pggy,2.0));
			
			// polar angle:
			double prodTheta = (180./TMath::Pi()) * atan2(pggt,pggz);
			
			// azimuthal angle:
			double prodPhi = (180./TMath::Pi()) * atan2(pggy,pggx);
			
			// opening angle:
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
			
			// invariant mass:
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			//-----------------------------------------------------//
			// Different options for applying BCAL+SC Veto:
			
			// Check SC Matches:
			
			int nSCHits = 0;
			int nSCHits_coplanar = 0;
			for(vector<const DSCHit*>::const_iterator sc = scHits.begin(); sc != scHits.end(); sc++) {
				
				// only check hits between 1ns < (t_sc - t_RF) < 7ns 
				//    and with dE > 0.0002 (from DNeutralShower_factory)
				
				double locT  = (*sc)->t - rfTime;
				double locdE = (*sc)->dE;
				
				if((0.5<invmass) && (invmass<0.60)) h_scRFdt_cut->Fill(locT);
				
				if((1.0 < locT) && (locT < 9.0) && (locdE > 0.0002)) {
					nSCHits++;
					int sector = (*sc)->sector;
					double locDeltaPhi = prodPhi - (m_scPos[sector-1][0].Phi() * TMath::RadToDeg());
					if(IsCoplanarSC(locDeltaPhi)) nSCHits_coplanar++;
				}
			}
			
			vector<int> vetoOptions;
			for(int iveto=0; iveto<m_nVetos; iveto++) vetoOptions.push_back(0);
			
			// Option 0 (no veto): No Veto is applied:
			vetoOptions[0] = 1;
			
			// Option 1 (strict): Remove events with any BCAL shower within +/-12ns:
			if(nBCALShowers==0) vetoOptions[1] = 1;
			
			// Option 2 (loose): Remove events with any BCAL shower within +/-1ns:
			if(nBCALShowers_1ns==0) vetoOptions[2] = 1;
			
			// Options 3-6:
			//   If there are no showers within +/-12ns in the BCAL, accept event.
			//   If there is exactly 1 such shower, then:
			//     Require coplanarity cut (option 3)
			//     Require coplanarity cut + timing > 1ns (option 4)
			//     Require coplanarity cut + timing > 1ns + all SC hits should also be coplanar (option 5)
			//     Require coplanarity cut + timing > 1ns + all SC hits should also be coplanar + no more than 1 SC hit (option 6)
			
			if(nBCALShowers==0)
			{
				vetoOptions[3] = 1;
				vetoOptions[4] = 1;
				if(nSCHits==nSCHits_coplanar) {
					vetoOptions[5] = 1;
					if(nSCHits<=1) {
						vetoOptions[6] = 1;
					}
				}
			}
			else if(nBCALShowers==1)
			{
				if(IsCoplanarBCAL(prodPhi-bcalPhi)) {
					vetoOptions[3] = 1;
					if(bcalRFdt>1.0) {
						vetoOptions[4] = 1;
						if(nSCHits==nSCHits_coplanar) {
							vetoOptions[5] = 1;
							if(nSCHits<=1) {
								vetoOptions[6] = 1;
							}
						}
					}
				}
			}
			
			// Option 7: strict SC + strict BCAL vetos:
			if(nSCHits==0 && nBCALShowers==0) vetoOptions[7] = 1;
			
			//-----------------------------------------------------//
			// Loop over beam photons:
			
			for(vector<const DBeamPhoton*>::const_iterator gam = beamPhotons.begin(); 
				gam != beamPhotons.end(); gam++) {
				
				double eb    = (*gam)->lorentzMomentum().E(); // energy of beam photon
				double brfdt = (*gam)->time() - rfTime;
				
				// remove beam photons below the minimum energy cut:
				if((eb < m_minBeamEnergyCut) || (eb > m_maxBeamEnergyCut)) continue;
				
				// Accidental subtraction procedure:
				//   - Fill histograms with a weight of 1.0 for beam photons within main RF bunch
				//   - Select two side-bands to the left and two-sidebands to the right (4 in total)
				//   - Fill histograms with a weight of -1/4 for beam photons in these sidebands.
				//      - An extra scaling factor is applied to out-of-time beam photons to account
				//        for the non-uniformity of beam bunches.
				//      - This scaling factor comes from the /ANALYSIS/accidental_scaling_factor tabls in the CCDB.
				//      - reference: GlueX-doc-4122 (B. Zihlmann)
				
				double fillWeight = 0.0;
				if(fabs(brfdt) < m_BeamRFCut) 
				{
					fillWeight = 1.0;
				}
				else if((-30.060<=brfdt && brfdt<=-22.044) || (22.044<=brfdt && brfdt<= 30.060)) 
				{
					fillWeight = -0.25*GetAccScalingFactor(eb);
				}
				//else { continue; }
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				double etaEnergy_coh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
				
				// Calculate the energy of the eta meson, assuming production on a free nucleon:
				double etaEnergy = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(Proton));
				
				// Apply a cut on the elasticity
				//  (ratio of measured energy of 2-photons, to the calculated energy above):
				bool isElastic = IsElasticCut(Egg, etaEnergy, prodTheta);
				
				// set a variable to indicate if the two-photon mass is consistent with an eta meson:
				bool isEta = IsEtaCut(invmass);
				
				// Plot timing distribution of beam photons after elasticity cut to see the level of accidentals:
				if(isElastic && isEta && vetoOptions[6]) {
					h_beamRFdt_cut->Fill(brfdt);
				}
				
				// If the beam photon wasn't in the main RF bunch or selected sidebands, skip it:
				if(fillWeight==0.0) continue;
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a coherently-produced eta meson:
				
				double sig1 = GetFCALEnergyRes(e1);
				double sig2 = GetFCALEnergyRes(e2);
				double sigr = pow(sig1/sig2,2.0);
				
				// Energy-constrained invariant mass assuming production on free nucleon:
				double e1c = e1/(1.+sigr) + (etaEnergy-e2)/(1.+(1./sigr));
				double e2c = etaEnergy - e1c;
				double invmassConstr = sqrt(2.*e1c*e2c*(1.-cos12));
				
				// Energy-constrained invariant mass assuming coherent production on nucleus:
				double e1c_coh = e1/(1.+sigr) + (etaEnergy_coh-e2)/(1.+(1./sigr));
				double e2c_coh = etaEnergy_coh - e1c;
				double invmassConstr_coh = sqrt(2.*e1c_coh*e2c_coh*(1.-cos12)); // energy-constrained invariant mass
				
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
				// Missing Mass 
				//   - This is calculated assuming the reaction took place on a free nucleon
				
				double  mmSq = pow(ParticleMass(Proton),2.0) + pow(ParticleMass(Eta),2.0) 
					+ 2.0*ParticleMass(Proton)*eb 
					- 2.0*ParticleMass(Proton)*Egg 
					- 2.0*eb*(Egg - sqrt(pow(Egg,2.0)-pow(ParticleMass(Eta),2.0))*cos(prodTheta*TMath::DegToRad()));
				
				double mmSq_coh = pow(ParticleMass(m_Target),2.0) + pow(ParticleMass(Eta),2.0) 
					+ 2.0*ParticleMass(m_Target)*eb 
					- 2.0*ParticleMass(m_Target)*Egg 
					- 2.0*eb*(Egg - sqrt(pow(Egg,2.0)-pow(ParticleMass(Eta),2.0))*cos(prodTheta*TMath::DegToRad()));
				
				//-----------------------------------------------------//
				// Default Cuts:
				
				if(vetoOptions[6]) {
					
					// Elasticity vs. Invmass Ratio:
					h_elasVSmgg->Fill(invmass/ParticleMass(Eta), Egg/etaEnergy, fillWeight);
					
					// Elasticity vs. Polar Angle for events where the mass is close to the eta:
					if(isEta) {
						h_elas->Fill(prodTheta, Egg/eb, fillWeight); // measured energy over beam energy
						h_elasCorr->Fill(prodTheta, Egg/etaEnergy, fillWeight); // measured energy over eta energy
						h_elasCorr_coh->Fill(prodTheta, Egg/etaEnergy_coh, fillWeight); // measured energy over eta energy (coherent-nuclear)
						
						if(fillWeight==1.0) h_elasCorrMain->Fill(prodTheta, Egg/etaEnergy);
						else                h_elasCorrSide->Fill(prodTheta, Egg/etaEnergy, -1.0*fillWeight);
					}
					
					// Energy-constrained Invariant Mass vs. Polar Angle:
					if(isElastic) {
						// invariant mass vs. polar angle:
						h_mgg->Fill(prodTheta, invmass, fillWeight);
						
						// energy-constrained invariant mass vs. polar angle:
						h_mggConstr->Fill(prodTheta, invmassConstr, fillWeight);
						
						// energy-constrained invariant mass vs. polar angle (assuming coherent production):
						h_mggConstr_coh->Fill(prodTheta, invmassConstr_coh, fillWeight);
						
						// for monitoring purposes, plot the invariant mass separately for 
						//   beam photons in the main RF bunch and for beam photons in the accidental sidebands:
						
						if(fillWeight==1.0) {
							h_mggMain->Fill(prodTheta, invmass);
							h_mggConstrMain->Fill(prodTheta, invmassConstr);
						} else {
							h_mggSide->Fill(prodTheta, invmass, -1.0*fillWeight);
							h_mggConstrSide->Fill(prodTheta, invmassConstr, -1.0*fillWeight);
						}
						
						// for MC, plot invariant mass vs thrown information:
						if(isMC) {
							h_mggThrown->Fill(thrownEtaAngle, invmass, fillWeight);
							h_mggConstrThrown->Fill(thrownEtaAngle, invmassConstr, fillWeight);
							
							// plot reconstructed vs. thrown angle:
							if(0.5 < invmassConstr && invmassConstr < 0.6) {
								h_recVSthrown->Fill(thrownEtaAngle, prodTheta, fillWeight);
							}
						}
						
						if(isEta) {
							// plot x-y distribution of showers:
							h_xy_1->Fill(pos1.X(), pos1.Y(), fillWeight);
							h_xy_2->Fill(pos2.X(), pos2.Y(), fillWeight);
						}
					}
					
					// Missing Mass vs. Polar Angle:
					h_mm->Fill(prodTheta, mmSq-pow(ParticleMass(Proton),2.0), fillWeight);
					h_mm_coh->Fill(prodTheta, mmSq_coh-pow(ParticleMass(m_Target),2.0), fillWeight);
					if(isEta) {
						h_mm_etaCut->Fill(prodTheta, mmSq-pow(ParticleMass(Proton),2.0), fillWeight);
						h_mm_coh_etaCut->Fill(prodTheta, mmSq_coh-pow(ParticleMass(m_Target),2.0), fillWeight);
						if(isElastic) {
							h_mm_elasEtaCut->Fill(prodTheta, mmSq-pow(ParticleMass(Proton),2.0), fillWeight);
							h_mm_coh_elasEtaCut->Fill(prodTheta, mmSq_coh-pow(ParticleMass(m_Target),2.0), fillWeight);
						}
					}
				} // end cut requiring default BCAL+SC veto option
				
				//-----------------------------------------------------//
				// Plot results with different veto options:
				
				for(int iveto=0; iveto<m_nVetos; iveto++) {
					if(vetoOptions[iveto]) {
						if(isEta) h_elas_veto[iveto]->Fill(prodTheta, Egg/etaEnergy, fillWeight);
						if(isElastic) {
							h_mgg_veto[iveto]->Fill(prodTheta, invmass, fillWeight);
							h_mggConstr_veto[iveto]->Fill(prodTheta, invmassConstr, fillWeight);
							h_mggConstr_coh_veto[iveto]->Fill(prodTheta, invmassConstr_coh, fillWeight);
							if(isEta) {
								h_mm_veto[iveto]->Fill(prodTheta, mmSq-pow(ParticleMass(Proton),2.0), fillWeight);
								h_mm_coh_veto[iveto]->Fill(prodTheta, mmSq_coh-pow(ParticleMass(m_Target),2.0), fillWeight);
							}
						}
					}
				}
				
			}	// loop over DBeamPhotons
		} // inner loop over DFCALShowers
	} // outer loop over DFCALShowers
	
	return;
}

//------------------------------------------//

bool JEventProcessor_primex_eta_analysis::IsElasticCut(double Egg, double Eeta, double theta)
{
	double ElasPeakMean  = m_ElasMean_p0 + m_ElasMean_p1*theta;
	double ElasPeakWidth = m_ElasWidth * m_ElasSigmaCut;
	if(fabs((Egg/Eeta)-ElasPeakMean)<ElasPeakWidth) return true;
	else return false;
}

bool JEventProcessor_primex_eta_analysis::IsEtaCut(double invmass)
{
	if((0.5<=invmass) && (invmass<0.6)) return true;
	else return false;
}

bool JEventProcessor_primex_eta_analysis::IsCoplanarBCAL(double deltaPhi)
{
	// deltaPhi = #phi_gg - #phi_SC
	
	if(m_phaseVal==1) {
		// FIELD OFF:
		if(fabs(fabs(deltaPhi)-180.0) < m_BCALDeltaPhiCut)
			return true;
		else
			return false;
	} else {
		// FIELD ON:
		if(((140.0<deltaPhi) && (deltaPhi<240.0)) || ((-220.0<deltaPhi) && (deltaPhi<-120.0)))
			return true;
		else 
			return false;
	}
}

bool JEventProcessor_primex_eta_analysis::IsCoplanarSC(double deltaPhi)
{
	// deltaPhi = #phi_gg - #phi_SC
	
	if(fabs(fabs(deltaPhi)-180.0) < m_SCDeltaPhiCut)
		return true;
	else
		return false;
}

//------------------------------------------//

int JEventProcessor_primex_eta_analysis::FCALFiducialCut(TVector3 pos, double cutLayer) {
	
	int locFiducialCut = 0;
	
	double fcalInnerLayerCut = (1.5 + cutLayer) * 4.0157;
	
	double fcalFaceX = m_beamSpot.X() + (pos.X() * (m_fcalFace.Z() - m_beamSpot.Z())/pos.Z()) - m_fcalFace.X();
	double fcalFaceY = m_beamSpot.Y() + (pos.Y() * (m_fcalFace.Z() - m_beamSpot.Z())/pos.Z()) - m_fcalFace.Y();
	
	if((fabs(fcalFaceX) < fcalInnerLayerCut) && (fabs(fcalFaceY) < fcalInnerLayerCut)) locFiducialCut = 1;
	
	// exclude showers near outer layers:
	
	double fcalFaceR = sqrt(pow(fcalFaceX,2.0) + pow(fcalFaceY,2.0));
	if(fcalFaceR > 100.0) locFiducialCut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	/*
	if(m_phaseVal < 2) {
		if((-32.<fcalFaceY) && (fcalFaceY<-20.) && (-8.<fcalFaceX) && (fcalFaceX<4.)) {
			locFiducialCut = 1;
		}
	}
	*/
	return locFiducialCut;
}

double JEventProcessor_primex_eta_analysis::GetEnergyAfterRecoil(double eb, double theta, 
	double m0, double mp) 
{
	// convert angle to radians:
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

double JEventProcessor_primex_eta_analysis::GetFCALEnergyRes(double e)
{
	// Hard-coded values for the FCAL energy resolution (taken from GlueX NIM paper)
	
	double a = 0.062, b = 0.047;
	double sig = (a*a)/e + (b*b);
	sig = sqrt(sig) * e;
	return sig;
}

double JEventProcessor_primex_eta_analysis::GetAccScalingFactor(double eb)
{
	if(eb > m_TAGMEnergyBoundHi)
		return m_HodoscopeHiFactor;
	else if(eb > m_TAGMEnergyBoundLo)
		return m_MicroscopeFactor;
	else
		return m_HodoscopeLoFactor;
}

void JEventProcessor_primex_eta_analysis::CheckTOFMatch(DVector3 pos1, double rfTime, 
	vector<const DTOFPoint*> tofPoints, double &dxMin, double &dyMin, double &dtMin, double rfTimingCut) {
	
	dxMin = 1000.;
	dyMin = 1000.;
	dtMin = 1000.;
	
	for(vector<const DTOFPoint*>::const_iterator tof = tofPoints.begin(); 
		tof != tofPoints.end(); tof++) {
		
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
		
		if(fabs(dt) < rfTimingCut) {
			if((dx*dx + dy*dy) < (dxMin*dxMin + dyMin*dyMin)) {
				dxMin = dx;
				dyMin = dy;
				dtMin = dt;
			}
		}
	}
	return;
}

//------------------------------------------//

Particle_t JEventProcessor_primex_eta_analysis::GetTargetType(int32_t runNumber) {
	
	// Set Target type according to RunPeriod and RunNumber:
	
	if(runNumber<60000)
	{
		// GlueX
		return Proton;
	}
	else if(runNumber<70000)
	{
		// PrimEx-eta Phase 1
		if(runNumber<61355) return Be9;
		else return Helium;
	}
	else if(runNumber<80000)
	{
		// GlueX
		return Proton;
	}
	else if(runNumber<90000)
	{
		// PrimEx-eta Phase 2
		if(runNumber<81382) return Be9;
		else return Helium;
	}
	else if(runNumber<100000)
	{
		// SRC/CT
		if((90034<=runNumber) && (runNumber<=90200)) {
			return Helium;
		}
		else if((90207<=runNumber) && (runNumber<=90249)) {
			return Deuteron;
		}
		else if((90263<=runNumber) && (runNumber<=90536)) {
			return C12;
		}
		else if((90558<=runNumber) && (runNumber<=90601)) {
			return Deuteron;
		}
		else if((90607<=runNumber) && (runNumber<=90660)) {
			return Helium;
		}
		else {
			return Proton;
		}
	}
	else if(runNumber<110000)
	{
		// CPP/NPP
		return Proton; // placeholder
	}
	else if(runNumber<120000)
	{
		// PrimEx-eta Phase 3
		if(runNumber<110622) return Be9;
		else return Helium;
	}
	else
	{
		return Proton;
	}
}

int JEventProcessor_primex_eta_analysis::GetPrimExPhase(int32_t runNumber)
{
	if(( 60000<=runNumber) && (runNumber<= 69999)) return 1;
	if(( 80000<=runNumber) && (runNumber<= 89999)) return 2;
	if((110000<=runNumber) && (runNumber<=119999)) return 3;
	return 0;
}

//-------------------------------------------------------------------------------------//

void JEventProcessor_primex_eta_analysis::InitializeHistograms()
{
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	TDirectory *dirPrimExEta = new TDirectoryFile("primex_eta_analysis", "primex_eta_analysis");
	dirPrimExEta->cd();
	
	//------------------------------------------//
	// Thrown angular distributions for different energy ranges:
	
	TDirectory *dirThrown = new TDirectoryFile("thrown", "thrown");
	dirThrown->cd();
	for(int icut=0; icut<13; icut++) {
		double beamEnergyCut = 7.6 + 0.2*(double)(icut);
		h_ThrownAngle[icut] = new TH1F(Form("ThrownAngle_%02d", icut), 
			Form("Thrown Angle of #eta (E_{#gamma} > %.1f GeV)", beamEnergyCut), 
			nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin);
	}
	dirThrown->cd("../");
	
	//------------------------------------------//
	// Timing and energy sum for different trigger types:
	
	TDirectory *dirTiming = new TDirectoryFile("rfTiming", "rfTiming");
	dirTiming->cd();
	for(int itrig=0; itrig<m_nTrigs; itrig++) {
		h_fcalRFdt[itrig] = new TH1F(Form("fcalRFdt_%d",itrig), 
			Form("FCAL - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100., 100.);
		h_bcalRFdt[itrig] = new TH1F(Form("bcalRFdt_%d",itrig), 
			Form("BCAL - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100., 100.);
		h_ccalRFdt[itrig] = new TH1F(Form("ccalRFdt_%d",itrig), 
			Form("CCAL - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100., 100.);
		h_tofRFdt[itrig] = new TH1F(Form("tofRFdt_%d",itrig), 
			Form("TOF - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100., 100.);
		h_scRFdt[itrig] = new TH1F(Form("scRFdt_%d",itrig), 
			Form("SC - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100., 100.);
		h_taghRFdt[itrig] = new TH1F(Form("taghRFdt_%d",itrig), 
			Form("TAGH - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100.0, 100.0);
		h_tagmRFdt[itrig] = new TH1F(Form("tagmRFdt_%d",itrig), 
			Form("TAGM - RF Time (%s); [ns]", m_triggerNames[itrig].c_str()), 
			2000, -100.0, 100.0);
	}
	dirTiming->cd("../");
	
	TDirectory *dirFCALEnergy = new TDirectoryFile("fcalEnergySum", "fcalEnergySum");
	dirFCALEnergy->cd();
	for(int itrig=0; itrig<m_nTrigs; itrig++) {
		h_fcalEnergySum[itrig] = new TH1F(Form("fcalEnergySum_%d",itrig), 
			Form("FCAL Shower Energy Sum (%s); E_{FCAL} [GeV]", m_triggerNames[itrig].c_str()), 
			1200, 0., 12.);
	}
	dirFCALEnergy->cd("../");
	
	//------------------------------------------//
	// FCAL-TOF Matching:
	
	TDirectory *dirFCALTOF = new TDirectoryFile("fcalTOFMatch", "fcalTOFMatch");
	dirFCALTOF->cd();
	h_fcalTOFdx = new TH1F("fcalTOFdx", "x_{FCAL} - x_{TOF} (closest DTOFPoint); [cm]", 2000, -100., 100.);
	h_fcalTOFdy = new TH1F("fcalTOFdy", "y_{FCAL} - y_{TOF} (closest DTOFPoint); [cm]", 2000, -100., 100.);
	h_fcalTOFdr = new TH1F("fcalTOFdr", "Distance between FCAL Shower and closest DTOFPoint; [cm]", 1000, 0., 100.);
	
	h_fcalTOFdt      = new TH1F("fcalTOFdt",      "t_{FCAL} - t_{TOF}; [ns]", 2000, -100., 100.);
	h_fcalTOFdt_cut  = new TH1F("fcalTOFdt_cut",  "t_{FCAL} - t_{TOF}; [ns]", 2000, -100., 100.);
	h_fcalTOFMatches = new TH1F("fcalTOFMatches", "Number of TOF Matches per 2-#gamma Event", 3, -0.5, 2.5);
	dirFCALTOF->cd("../");
	
	TDirectory *dirGG = new TDirectoryFile("eta_gg", "eta_gg");
	dirGG->cd();
	
	//------------------------------------------//
	// RF timing after cuts are applies:
	
	h_beamRFdt_cut = new TH1F("beamRFdt_cut", "; t_{#gamma} - t_{RF} (ns); counts / 0.1 ns", 2000, -100., 100.);
	h_scRFdt_cut   = new TH1F("scRFdt_cut", "; t_{SC} - t_{RF} (ns); counts / 0.1 ns", 2000, -100., 100.);
	
	//------------------------------------------//
	// Elasticity:
	
	h_elas = new TH2F("elas", "Elasticity", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0., 2.);
	h_elas->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elas->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#gamma}");
	
	h_elasCorr = new TH2F("elasCorr", "Elasticity (corrected for nucleon recoil)", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0., 2.);
	h_elasCorr->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elasCorr->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	h_elasCorrMain = new TH2F("elasCorrMain", "Elasticity (corrected for nucleon recoil)", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0., 2.);
	h_elasCorrMain->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elasCorrMain->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	h_elasCorrSide = new TH2F("elasCorrSide", "Elasticity (corrected for nucleon recoil)", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0., 2.);
	h_elasCorrSide->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elasCorrSide->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	h_elasCorr_coh = new TH2F("elasCorr_coh", "Elasticity (corrected for helium recoil)", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0., 2.);
	h_elasCorr_coh->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_elasCorr_coh->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	
	//------------------------------------------//
	// Invariant Mass:
	
	h_mgg = new TH2F("mgg", "Two-Photon Invariant Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg->Sumw2();
	
	h_mggMain = new TH2F("mggMain", "Two-Photon Invariant Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mggMain->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mggMain->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mggMain->Sumw2();
	
	h_mggSide = new TH2F("mggSide", "Two-Photon Invariant Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mggSide->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mggSide->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mggSide->Sumw2();
	
	h_mggConstr = new TH2F("mggConstr", "Energy-Constrained Inv Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mggConstr->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mggConstr->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mggConstr->Sumw2();
	
	h_mggConstrMain = new TH2F("mggConstrMain", "Energy-Constrained Inv Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mggConstrMain->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mggConstrMain->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mggConstrMain->Sumw2();
	
	h_mggConstrSide = new TH2F("mggConstrSide", "Energy-Constrained Inv Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mggConstrSide->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mggConstrSide->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mggConstrSide->Sumw2();
	
	h_mggConstr_coh = new TH2F("mggConstr_coh", "Energy-Constrained Inv Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mggConstr_coh->GetXaxis()->SetTitle("#theta_{#gamma#gamma}^{constr} [#circ]");
	h_mggConstr_coh->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mggConstr_coh->Sumw2();
	
	//------------------------------------------//
	// Missing Mass:
	
	h_mm = new TH2F("mm", "Squared Missing Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
	h_mm->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm->Sumw2();
	
	h_mm_etaCut = new TH2F("mm_etaCut", "Squared Missing Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
	h_mm_etaCut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_etaCut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_etaCut->Sumw2();
	
	h_mm_elasEtaCut = new TH2F("mm_elasEtaCut", "Squared Missing Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
	h_mm_elasEtaCut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_elasEtaCut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_elasEtaCut->Sumw2();
	
	h_mm_coh = new TH2F("mm_coh", "Squared Missing Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
	h_mm_coh->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_coh->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_coh->Sumw2();
	
	h_mm_coh_etaCut = new TH2F("mm_coh_etaCut", "Squared Missing Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
	h_mm_coh_etaCut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_coh_etaCut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_coh_etaCut->Sumw2();
	
	h_mm_coh_elasEtaCut = new TH2F("mm_coh_elasEtaCut", "Squared Missing Mass", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
	h_mm_coh_elasEtaCut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_coh_elasEtaCut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_coh_elasEtaCut->Sumw2();
	
	//------------------------------------------//
	// Miscellaneous:
	
	h_elasVSmgg = new TH2F("elasVSmgg", "Elasticity vs. Invmass Ratio", 500, 0.0, 2.0, 500, 0.0, 2.0);
	h_elasVSmgg->GetXaxis()->SetTitle("m_{#gamma#gamma} / m_{#eta,PDG}");
	h_elasVSmgg->GetXaxis()->SetTitleSize(0.05);
	h_elasVSmgg->GetXaxis()->CenterTitle(true);
	h_elasVSmgg->GetYaxis()->SetTitle("E_{#gamma#gamma} / E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
	h_elasVSmgg->GetYaxis()->SetTitleSize(0.05);
	h_elasVSmgg->GetYaxis()->CenterTitle(true);
	
	h_xy_1 = new TH2F("xy_1", "Postion of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -100., 100., 500, -100., 100.);
	h_xy_2 = new TH2F("xy_2", "Postion of Shower 2; x_{2} [cm]; y_{2} [cm]", 500, -100., 100., 500, -100., 100.);
	
	// For MC:
	
	h_recVSthrown = new TH2F("recVSthrown", "Reconstructed Angle vs. Thrown Angle; #theta_{thrown} [#circ]; #theta_{rec} [#circ]",
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin,
		nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin);
	h_recVSthrown->Sumw2();
	
	h_mggThrown = new TH2F("mggThrown", "Two-Photon Invariant Mass", 
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nInvmassBins,     m_minInvmassBin,     m_maxInvmassBin);
	h_mggThrown->GetXaxis()->SetTitle("#theta_{thrown} [#circ]");
	h_mggThrown->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV/c^{2}]");
	h_mggThrown->Sumw2();
	
	h_mggConstrThrown = new TH2F("mggConstrThrown", "Two-Photon Invariant Mass", 
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nInvmassBins,     m_minInvmassBin,     m_maxInvmassBin);
	h_mggConstrThrown->GetXaxis()->SetTitle("#theta_{thrown} [#circ]");
	h_mggConstrThrown->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV/c^{2}]");
	h_mggConstrThrown->Sumw2();
	
	//------------------------------------------//
	// Different Veto Options:
	
	for(int ihist=0; ihist<m_nVetos; ihist++) {
		
		h_elas_veto[ihist] = new TH2F(Form("elas_veto_%d",ihist), 
			Form("Elasticity (Veto Option %d)", ihist+1), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
		h_elas_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_elas_veto[ihist]->GetYaxis()->SetTitle(
			"E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
		
		h_mgg_veto[ihist] = new TH2F(Form("mgg_veto_%d",ihist), 
			Form("Two-Photon Invariant Mass (Veto Option %d)", ihist+1), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mgg_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg_veto[ihist]->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
		h_mgg_veto[ihist]->Sumw2();
		
		h_mggConstr_veto[ihist] = new TH2F(Form("mggConstr_veto_%d",ihist), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)", ihist+1), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_veto[ihist]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
		h_mggConstr_veto[ihist]->Sumw2();
		
		h_mggConstr_coh_veto[ihist] = new TH2F(Form("mggConstr_coh_veto_%d",ihist), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)", ihist+1), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr_coh_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_coh_veto[ihist]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
		h_mggConstr_coh_veto[ihist]->Sumw2();
		
		h_mm_veto[ihist] = new TH2F(Form("mm_veto_%d",ihist),
			Form("Missing Mass Squared (Veto Option %d)", ihist+1), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
		h_mm_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_veto[ihist]->GetYaxis()->SetTitle("M_{miss}^{2} [GeV^{2}/c^{4}]");
		h_mm_veto[ihist]->Sumw2();
		
		h_mm_coh_veto[ihist] = new TH2F(Form("mm_coh_veto_%d",ihist),
			Form("Missing Mass Squared (Energy-constrained) (Veto Option %d)", ihist+1), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
		h_mm_coh_veto[ihist]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_coh_veto[ihist]->GetYaxis()->SetTitle("M_{miss}^{2} [GeV^{2}/c^{4}]");
		h_mm_coh_veto[ihist]->Sumw2();
	}
	
	dirGG->cd("../");
	
	dirPrimExEta->cd("../");
	
	return;
}
