#include "EtaAna.h"
#include "MyReadConfig.h"

// Default Constructor:

EtaAna::EtaAna() {
	
	// set default run number to 61691:
	
	m_runNumber = 61691;
	
	m_random = new TRandom3(0);
	
	// set defaults for cut values:
	
	m_FCALRFCut =  2.004;
	m_BCALRFCut = 12.028;
	m_BeamRFCut =  2.004;
	m_TOFRFCut  =  1.000;
	
	m_FCALEnergyCut      = 0.25; // threshold for photons to be used in this analysis
	m_FCALExtraEnergyCut = 0.25; // threshold for 'extra' showers in the FCAL
		// There should only be two photons with energy > m_FCALExtraEnergyCut. However, we
		// can vary m_FCALEnergyCut separately from m_FCALExtraEnergyCut in order to
		// change the energy threshold of the showers which are used for reconstructing etas,
		// without having to include "extra" higher energy showers.
	
	m_BCALEnergyCut      = 0.1;
	
	m_minBeamEnergyCut =  8.0;
	m_maxBeamEnergyCut = 10.9;
	
	m_FCALTOFCut      =  8.0; // [cm]
	m_BCALDeltaPhiCut = 30.0; // [deg]
	m_SCDeltaPhiCut   = 36.0; // [deg]
	
	m_ElasMean_p0  = 1.0;
	m_ElasMean_p1  = 0.0;
	m_ElasWidth    = 0.031;
	m_ElasSigmaCut = 4.0;
	
	// Unless otherwise specified, assume phase=1, target=Helium:
	
	m_phaseVal = 1;
	
	m_Target        = Helium;
	m_targetLength  = 29.535;
	m_targetDensity = 0.1217;
	m_targetAtten   = 0.00821;
	
	// Set defaults for Geometry from CCDB:
	
	m_fcalFace.SetXYZ(0.529, -0.002, 624.906);
	m_fcalCorrection.SetXYZ(0., 0., 0);
	
	m_vertex.SetXYZ(0.1914, -0.0769, 65.);
	
	// Set event number to zero on initialization:
	
	m_event = 0;
	
	// initialize as NULL pointers:
	
	h_thrown        = NULL;
	h_invmassMatrix = NULL;
}

//----------------------------------------------------------//
//----------------------------------------------------------//
//---                                                    ---//
//---       Private Member Function Definitions:         ---//
//---                                                    ---//
//----------------------------------------------------------//
//----------------------------------------------------------//

int EtaAna::SetGeometry() {
	
	if(m_phaseVal==1) {
		
		if(m_runNumber<61355) {
			
			// Beryllium Target runs:
			
			m_Target        = Be9;
			m_targetLength  = 1.7755;
			m_targetDensity = 1.848;
			m_targetAtten   = 0.01172;
			m_vertex.SetZ(64.935);
			
		} else {
			
			// Helium Target runs:
			
			m_Target        = Helium;
			m_targetLength  = 29.535;
			m_targetDensity = 0.1217;
			m_targetAtten   = 0.00821;
			m_vertex.SetZ(65.0);
		}
		
		// Correction to alignment after Simon updated beam spot:
		
		m_fcalFace.SetXYZ(0.455, -0.032, 624.906);
		m_fcalCorrection.SetXYZ(0.455-0.529, -0.032+0.002, 0.0);
		
		if(m_runNumber<61483) {
			m_vertex.SetX( 0.027);
			m_vertex.SetY(-0.128);
		} else if(m_runNumber<61774) {
			m_vertex.SetX( 0.001);
			m_vertex.SetY(-0.077);
		} else {
			m_vertex.SetX( 0.038);
			m_vertex.SetY(-0.095);
		}
	}
	else if(m_phaseVal==3) {
		
		if(m_runNumber<110622) {
			
			// Beryllium Target runs:
			
			m_Target        = Be9;
			m_targetLength  = 1.7755;
			m_targetDensity = 1.848;
			m_targetAtten   = 0.01172;
			m_vertex.SetZ(64.935);
			
		} else {
			
			// Helium Target runs:
			
			m_Target        = Helium;
			m_targetLength  = 29.535;
			m_targetDensity = 0.1217;
			m_targetAtten   = 0.00821;
			m_vertex.SetZ(65.0);
		}
		
		// Correction to alignment after Simon updated beam spot:
		
		m_fcalFace.SetXYZ(0.189, 0.022, 624.32);
		m_fcalCorrection.SetXYZ(0.0, 0.0, 0.0);
	}
	else {
		std::cout << "Unsupported run period provided. Skipping Run." << std::endl;
		return 1;
	}
	
	return 0;
}

int EtaAna::LoadTree() {
	
	m_tree = (TTree*)m_inputFile->Get("eta_gg");
	if(m_tree==NULL) return 0;
	
	// Reset event count to zero when laoding a new Tree:
	m_event = 0;
	
	return m_tree->GetEntries();
}

int EtaAna::CheckEventMultiplicities() {
	
	int nSystemsOver = 0;
	if(m_nfcal>MAX_FCAL) {
		nSystemsOver++;
		printf("    Too many FCAL showers reconstructed in event %d\n",m_event);
	}
	if(m_nbcal>MAX_BCAL) {
		nSystemsOver++;
		printf("    Too many BCAL showers reconstructed in event %d\n",m_event);
	}
	if(m_nbeam>MAX_BEAM) {
		nSystemsOver++;
		printf("    Too many Beam photons reconstructed in event %d\n",m_event);
	}
	if(m_ntof >MAX_TOF ) {
		nSystemsOver++;
		printf("    Too many TOF points reconstructed in event %d\n",m_event);
	}
	if(m_nsc  >MAX_SC  ) {
		nSystemsOver++;
		printf("    Too many SC hits reconstructed in event %d\n",m_event);
	}
	return nSystemsOver;
}

int EtaAna::AcceptRejectEvent() {
	//
	// Accept-reject filter to create a realistic z-vertex distribution for simulated events
	//
	if(m_nmc==0) return 0;
	
	int reject = 0;
	
	double vertexZ = m_mcZ[0];
	
	// shift coordinate system so that upstream entrance of target is at z=0:
	
	double locZ = vertexZ - m_vertex.Z() + (m_targetLength/2.0);
	
	// use attenuation length from XCOM database to calculate probability of photon 
	// absorption.
	
	double vertexWeight;
	if(locZ < 0.) {
		vertexWeight = 1.0;
	} else if(locZ > m_targetLength) {
		vertexWeight = TMath::Exp(-m_targetAtten * m_targetDensity * m_targetLength);
	} else {
		vertexWeight = TMath::Exp(-m_targetAtten * m_targetDensity * locZ);
	}
	
	if(vertexWeight < m_random->Uniform()) {
		reject = 1;
	}
	
	return reject;
}

void EtaAna::ReadEvent() {
	
	if(m_event == 0) {
		
		// Set Branch Address on first Event:
		
		m_tree->SetBranchAddress("rfTime",             &m_rfTime);
		m_tree->SetBranchAddress("nbeam",              &m_nbeam);
		m_tree->SetBranchAddress("tag_counter",        &m_tagCounter);
		m_tree->SetBranchAddress("tag_sys",            &m_tagSystem);
		m_tree->SetBranchAddress("beam_e",             &m_beamE);
		m_tree->SetBranchAddress("beam_t",             &m_beamT);
		m_tree->SetBranchAddress("acc_scale_factor",   &m_accScaleFactor);
		
		m_tree->SetBranchAddress("nfcal",              &m_nfcal);
		m_tree->SetBranchAddress("fcal_e",             &m_fcalE);
		m_tree->SetBranchAddress("fcal_x",             &m_fcalX);
		m_tree->SetBranchAddress("fcal_y",             &m_fcalY);
		m_tree->SetBranchAddress("fcal_z",             &m_fcalZ);
		m_tree->SetBranchAddress("fcal_t",             &m_fcalT);
		m_tree->SetBranchAddress("fcal_nblocks",       &m_fcalNblocks);
		
		m_tree->SetBranchAddress("nbcal",              &m_nbcal);
		m_tree->SetBranchAddress("bcal_e",             &m_bcalE);
		m_tree->SetBranchAddress("bcal_x",             &m_bcalX);
		m_tree->SetBranchAddress("bcal_y",             &m_bcalY);
		m_tree->SetBranchAddress("bcal_z",             &m_bcalZ);
		m_tree->SetBranchAddress("bcal_t",             &m_bcalT);
		
		m_tree->SetBranchAddress("ntof",               &m_ntof);
		m_tree->SetBranchAddress("tof_x",              &m_tofX);
		m_tree->SetBranchAddress("tof_y",              &m_tofY);
		m_tree->SetBranchAddress("tof_z",              &m_tofZ);
		m_tree->SetBranchAddress("tof_t",              &m_tofT);
		
		m_tree->SetBranchAddress("nsc",                &m_nsc);
		m_tree->SetBranchAddress("sc_sector",          &m_scSector);
		m_tree->SetBranchAddress("sc_phi",             &m_scPhi);
		m_tree->SetBranchAddress("sc_dE",              &m_scdE);
		m_tree->SetBranchAddress("sc_t",               &m_scT);
		m_tree->SetBranchAddress("sc_pulse_height",    &m_scPulseHeight);
		
		m_tree->SetBranchAddress("nmc",                &m_nmc);
		m_tree->SetBranchAddress("mc_pdgtype",         &m_mcPDGType);
		m_tree->SetBranchAddress("mc_x",               &m_mcX);
		m_tree->SetBranchAddress("mc_y",               &m_mcY);
		m_tree->SetBranchAddress("mc_z",               &m_mcZ);
		m_tree->SetBranchAddress("mc_t",               &m_mcT);
		m_tree->SetBranchAddress("mc_e",               &m_mcE);
		m_tree->SetBranchAddress("mc_p",               &m_mcP);
		m_tree->SetBranchAddress("mc_theta",           &m_mcTheta);
		m_tree->SetBranchAddress("mc_phi",             &m_mcPhi);
	}
	
	m_tree->GetEvent(m_event);
	
	return;
}

void EtaAna::PlotThrown(double energy, double angle) {
	
	if(m_nmc==0) return;
	if(h_thrown==NULL) {
		h_thrown = new TH2F("thrown", 
			"Thrown Angle vs. Thrown Beam Energy; E_{#gamma}(thrown) [GeV]; #theta(thrown) [#circ]",
			70, 5.0, 12.0, 500, 0.0, 5.0);
		h_thrown->SetDirectory(0);
	}
	h_thrown->Fill(energy, angle);
	return;
}

void EtaAna::FillAngularMatrix(int vetoOption, double thrownEnergy, double thrownAngle, 
	double recAngle, double weight) {
	
	double minBeamEnergy      =  7.0;
	double maxBeamEnergy      = 12.0;
	double beamEnergyBinSize  =  0.1;
	
	double minRecAngle        = 0.0;
	double maxRecAngle        = 5.5;
	double recAngleBinSize    = 0.01;
	
	double minThrownAngle     = 0.0;
	double maxThrownAngle     = 5.0;
	double thrownAngleBinSize = 0.01;
	
	if((minBeamEnergy>thrownEnergy) || (thrownEnergy>maxBeamEnergy)) return;
	
	// check if matrix has already been initialized:
	if(h_AngularMatrix.size()==0) {
		
		int nBeamEnergyBins  = (int)((maxBeamEnergy-minBeamEnergy)/beamEnergyBinSize);
		int nRecAngleBins    = (int)((maxRecAngle-minRecAngle)/recAngleBinSize);
		int nThrownAngleBins = (int)((maxThrownAngle-minThrownAngle)/thrownAngleBinSize);
		
		for(int iveto=0; iveto<m_nVetos; iveto++) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_veto_%d",iveto), 
				Form("Veto Option %d; #theta(thrown) [#circ]; #theta(rec) [#circ]", iveto),
				nThrownAngleBins, minThrownAngle, maxThrownAngle, 
				nRecAngleBins,    minRecAngle,    maxRecAngle,
				nBeamEnergyBins,  minBeamEnergy,  maxBeamEnergy);
			hMatrix->SetDirectory(0);
			h_AngularMatrix.push_back(hMatrix);
		}
	}
	
	// find the index associated with this beam energy:
	h_AngularMatrix[vetoOption]->Fill(thrownAngle, recAngle, thrownEnergy, weight);
	
	return;
}

void EtaAna::FillInvmassMatrix(double theta, double mgg, double beamEnergy, double weight) {
	
	// check if matrix has already been initialized:
	if(h_invmassMatrix==NULL) {
		h_invmassMatrix = new TH3F("invmassMatrix", 
			"Invariant Mass Matrix; #theta(rec) [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]; E_{#gamma} [GeV]",
			650, 0.0, 6.5, 600, 0.0, 1.2, 50, 7.0, 12.0);
		h_invmassMatrix->SetDirectory(0);
		
		h_invmassMatrix_acc = new TH3F("invmassMatrix_acc", 
			"Invariant Mass Matrix; #theta(rec) [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]; E_{#gamma} [GeV]",
			650, 0.0, 6.5, 600, 0.0, 1.2, 50, 7.0, 12.0);
		h_invmassMatrix_acc->SetDirectory(0);
	}
	
	if(weight<0.0) h_invmassMatrix_acc->Fill(theta, mgg, beamEnergy);
	else h_invmassMatrix->Fill(theta, mgg, beamEnergy);
	
	return;
}

int EtaAna::GetFCALShowerList(vector<int> &goodShowers, int &nGoodFCALShowers, 
	double energyCut, double extraEnergyCut, double fiducialCut, double timingCut) {
	
	int nFCALShowers = 0;
	nGoodFCALShowers = 0;
	for(int ishow=0; ishow<m_nfcal; ishow++) {
		
		TVector3 locPos = GetFCALPosition(ishow);
		double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		
		// to be consistent between Data and MC, we shoud remove showers from the dead region as soon as possible: 
		if(m_phaseVal==1) {
			
			double locFCALFaceX = m_vertex.X() - m_fcalFace.X() 
				+ (locPos.X() * (m_fcalFace.Z() - m_vertex.Z())/locPos.Z());
			double locFCALFaceY = m_vertex.Y() - m_fcalFace.Y() 
				+ (locPos.Y() * (m_fcalFace.Z() - m_vertex.Z())/locPos.Z());
			
			if((-32. < locFCALFaceY && locFCALFaceY < -20.) 
				&& (-8. < locFCALFaceX && locFCALFaceX < 4.)) continue;
		}
		
		if(fabs(locT) < timingCut) {
			nFCALShowers++;
			if(m_fcalE[ishow] > extraEnergyCut) nGoodFCALShowers++;
			
			if((m_fcalE[ishow] > energyCut) && !FCALFiducialCut(locPos, fiducialCut)) {
				goodShowers.push_back(ishow);
			}
		}
	}
	
	return nFCALShowers;
}

int EtaAna::GetBCALShowerList(vector<int> &goodShowers, double energyCut, double timingCut) {
	
	int nBCALShowers = 0;
	for(int ishow=0; ishow<m_nbcal; ishow++) {
		
		TVector3 locPos = GetBCALPosition(ishow);
		double locT = m_bcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		
		if(fabs(locT) < timingCut) {
			nBCALShowers++;
			if(m_bcalE[ishow] > energyCut) {
				goodShowers.push_back(ishow);
			}
		}
	}
	
	return nBCALShowers;
}

int EtaAna::GetBeamPhotonList(vector<pair<int,double>> &goodPhotons, double minEnergyCut, double maxEnergyCut) {
	
	int nBeamPhotons = 0;
	for(int igam=0; igam<m_nbeam; igam++) {
		
		double locRFdt   = m_beamT[igam] - m_rfTime;
		double locWeight = 0.0;
		
		double locBeamCut = beamBunchesMain*m_BeamRFCut;
		
		if(fabs(locRFdt) < locBeamCut) locWeight = 1.0;
		else if(
			((beamBunchesMain+5.5)*4.008 < fabs(locRFdt)) && 
			(fabs(locRFdt) < (beamBunchesMain+5.5+beamBunchesAcc)*4.008)
		) locWeight = -1.0/(2.0*beamBunchesAcc);
		else continue;
		
		if(locWeight < 0.0) locWeight *= m_accScaleFactor[igam];
		if((m_beamE[igam] > minEnergyCut) && (m_beamE[igam] < maxEnergyCut)) {
			nBeamPhotons++;
			goodPhotons.push_back({igam,locWeight});
		}
	}
	return nBeamPhotons;
}

TVector3 EtaAna::GetFCALPosition(int index) {
	
	TVector3 pos(m_fcalX[index], m_fcalY[index], m_fcalZ[index]);
	pos = pos + m_fcalCorrection - m_vertex;
	
	return pos;
}
TVector3 EtaAna::GetBCALPosition(int index) {
	
	TVector3 pos(m_bcalX[index], m_bcalY[index], m_bcalZ[index]);
	pos = pos - m_vertex;
	
	return pos;
}

int EtaAna::FCALFiducialCut(TVector3 pos, double cutLayer) {
	
	int locFiducialCut = 0;
	
	double fcalInnerLayerCut = (1.5 + cutLayer) * m_fcalBlockSize;
	
	double fcalFaceX = m_vertex.X() + (pos.X() * (m_fcalFace.Z() - m_vertex.Z())/pos.Z()) - m_fcalFace.X();
	double fcalFaceY = m_vertex.Y() + (pos.Y() * (m_fcalFace.Z() - m_vertex.Z())/pos.Z()) - m_fcalFace.Y();
	
	if((fabs(fcalFaceX) < fcalInnerLayerCut) && (fabs(fcalFaceY) < fcalInnerLayerCut)) locFiducialCut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	
	if(m_phaseVal < 2) {
		if((-32.<fcalFaceY) && (fcalFaceY<-20.) && (-8.<fcalFaceX) && (fcalFaceX<4.)) {
			locFiducialCut = 1;
		}
	}
	
	return locFiducialCut;
}

void EtaAna::CheckTOFMatch(TVector3 pos, double &dxMin, double &dyMin, double &dtMin, double rfTimingCut) {
	
	dxMin = 1000.;
	dyMin = 1000.;
	dtMin = 1000.;
	
	for(int itof=0; itof<m_ntof; itof++) {
		
		double xt = m_tofX[itof] - m_vertex.X();
		double yt = m_tofY[itof] - m_vertex.Y();
		double zt = m_tofZ[itof] - m_vertex.Z();
		double rt = sqrt(xt*xt + yt*yt + zt*zt);
		double dt = m_tofT[itof] - (rt/m_c) - m_rfTime;
		xt *= pos.Z() / zt;
		yt *= pos.Z() / zt;
		double dx = pos.X() - xt;
		double dy = pos.Y() - yt;
		
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

double EtaAna::GetEnergyAfterRecoil(double eb, double theta, double m0, double mp) {
	
	theta *= TMath::DegToRad();
	
	double t1 = eb*cos(theta);
	double t2 = mp+eb;
	double t3 = mp*eb + m0*m0*0.5;
	
	double a = t1*t1-t2*t2;
	double b = 2.*t2*t3;
	double c = -m0*m0*t1*t1-t3*t3;
	double d = b*b - 4.*a*c;
	
	if((d < 0.) || (a == 0.)) {
		cout << "IMAGINARY ETA ENERGY!!!" << endl;
		return 0.;
	}
	
	double energy = (-b-sqrt(d))/2./a;
	return energy;
}

double EtaAna::GetFCALEnergyResolution(double e) {
	
	// hard-coded values for the FCAL energy resolution (taken from GlueX NIM paper)
	
	double a = 0.062, b = 0.047;
	double sig = (a*a)/e + (b*b);
	sig = sqrt(sig) * e;
	return sig;
}

void EtaAna::GetThrownEnergyAndAngle(double &thrownEnergy, double &thrownAngle) {
	if(m_nmc!=2) return;
	
	thrownEnergy = -1.0*ParticleMass(Helium);
	thrownAngle  =  0.0;
	for(int imc=0; imc<m_nmc; imc++) {
		thrownEnergy += m_mcE[imc];
		if(m_mcPDGType[imc]==PDGtype(Eta)) thrownAngle = m_mcTheta[imc];
	}
	return;
}

bool EtaAna::IsElasticCut(double Egg, double Eeta, double theta) {
	double ElasPeakMean  = m_ElasMean_p0 + m_ElasMean_p1*theta;
	double ElasPeakWidth = m_ElasWidth * m_ElasSigmaCut;
	if(fabs((Egg/Eeta)-ElasPeakMean)<ElasPeakWidth) return true;
	else return false;
}

bool EtaAna::IsEtaCut(double invmass) {
	if((invmass<=0.5) && (invmass<0.6)) return true;
	else return false;
}

//----------------------------------------------------------//
//----------------------------------------------------------//
//---                                                    ---//
//---       Public Member Function Definitions:          ---//
//---                                                    ---//
//----------------------------------------------------------//
//----------------------------------------------------------//

int EtaAna::GetPrimexPhase(int runNumber) {
	
	int locPhase = 0;
	if(runNumber < 60000) {
		return 0;
	} else if((runNumber >= 60000) && (runNumber <= 69999)) {
		return 1;
	} else if((runNumber >= 80000) && (runNumber <= 89999)) {
		return 2;
	} else if((runNumber >=110000) && (runNumber <=119999)) {
		return 3;
	} else {
		return 0;
	}
}

int EtaAna::SetRunNumber(int runNumber) { 
	
	m_runNumber = runNumber;
	m_phaseVal  = GetPrimexPhase(runNumber);
	return SetGeometry();
}

int EtaAna::SetCuts(TString configFileName) {
	
	if(gSystem->AccessPathName(configFileName.Data())) {
		printf("Problem accessing config file. Using default cut values.\n");
		return 1;
	} else {
		printf("Reading config file from: %s\n",configFileName.Data());
	}
	
	MyReadConfig *ReadFile = new MyReadConfig();
	ReadFile->ReadConfigFile(configFileName.Data());
	
	if(ReadFile->GetConfigName("FCALRFDT") != "") {
		m_FCALRFCut = ReadFile->GetConfig1Par("FCALRFDT")[0];
	}
	if(ReadFile->GetConfigName("BCALRFDT") != "") {
		m_BCALRFCut = ReadFile->GetConfig1Par("BCALRFDT")[0];
	}
	if(ReadFile->GetConfigName("BEAMRFDT") != "") {
		m_BeamRFCut = ReadFile->GetConfig1Par("BEAMRFDT")[0];
	}
	if(ReadFile->GetConfigName("TOFRFDT") != "") {
		m_TOFRFCut = ReadFile->GetConfig1Par("TOFRFDT")[0];
	}
	
	if(ReadFile->GetConfigName("FCALENERGY") != "") {
		m_FCALEnergyCut = ReadFile->GetConfig1Par("FCALENERGY")[0];
	}
	if(ReadFile->GetConfigName("EXTRAFCALENERGY") != "") {
		m_FCALExtraEnergyCut = ReadFile->GetConfig1Par("EXTRAFCALENERGY")[0];
	}
	if(ReadFile->GetConfigName("BCALENERGY") != "") {
		m_BCALEnergyCut = ReadFile->GetConfig1Par("BCALENERGY")[0];
	}
	if(ReadFile->GetConfigName("MINBEAMENERGY") != "") {
		m_minBeamEnergyCut = ReadFile->GetConfig1Par("MINBEAMENERGY")[0];
	}
	if(ReadFile->GetConfigName("MAXBEAMENERGY") != "") {
		m_maxBeamEnergyCut = ReadFile->GetConfig1Par("MAXBEAMENERGY")[0];
	}
	
	if(ReadFile->GetConfigName("FCALTOFCUT") != "") {
		m_FCALTOFCut = ReadFile->GetConfig1Par("FCALTOFCUT")[0];
	}
	
	if(ReadFile->GetConfigName("BCALDELTAPHI") != "") {
		m_BCALDeltaPhiCut = ReadFile->GetConfig1Par("BCALDELTAPHI")[0];
	}
	if(ReadFile->GetConfigName("SCDELTAPHI") != "") {
		m_SCDeltaPhiCut = ReadFile->GetConfig1Par("SCDELTAPHI")[0];
	}
	
	if(ReadFile->GetConfigName("ELASMEAN_P0") != "") {
		m_ElasMean_p0 = ReadFile->GetConfig1Par("ELASMEAN_P0")[0];
	}
	if(ReadFile->GetConfigName("ELASMEAN_P1") != "") {
		m_ElasMean_p1 = ReadFile->GetConfig1Par("ELASMEAN_P1")[0];
	}
	if(ReadFile->GetConfigName("ELASWIDTH") != "") {
		m_ElasWidth = ReadFile->GetConfig1Par("ELASWIDTH")[0];
	}
	if(ReadFile->GetConfigName("ELASSIGMA") != "") {
		m_ElasSigmaCut = ReadFile->GetConfig1Par("ELASSIGMA")[0];
	}
	
	delete ReadFile;
	return 0;
}

void EtaAna::DumpCuts() {
	
	printf("\n===================================\n");
	printf("CUTS:\n\n");
	printf("FCAL-RF Timing Cut: %.3f ns\n", m_FCALRFCut);
	printf("BCAL-RF Timing Cut: %.3f ns\n", m_BCALRFCut);
	printf("Beam-RF Timing Cut: %.3f ns\n", m_BeamRFCut);
	printf(" TOF-RF Timing Cut: %.3f ns\n", m_TOFRFCut );
	printf("\n");
	printf("FCAL Energy Cut: %.2f GeV\n", m_FCALEnergyCut);
	printf("FCAL Energy Cut for extra showers: %.2f GeV\n", m_FCALExtraEnergyCut);
	printf("BCAL Energy Cut: %.2f GeV\n", m_BCALEnergyCut);
	printf("Beam Energy Cut: %.2f GeV - %.2f GeV\n", m_minBeamEnergyCut, m_maxBeamEnergyCut);
	printf("\n");
	printf("FCAL-TOF distance cut: %.2f cm\n", m_FCALTOFCut);
	printf("BCAL DeltaPhi Cut: %.1f deg.\n", m_BCALDeltaPhiCut);
	printf("  SC DeltaPhi Cut: %.1f deg.\n", m_SCDeltaPhiCut  );
	printf("\n");
	printf("===================================\n\n\n");
	
	return;
}

void EtaAna::RunAnalysis(TString inputFileName, int analysisOption) {
	
	m_inputFile = new TFile(inputFileName.Data(), "READ");
	int nTotalEvents = LoadTree();
	
	while(m_event < nTotalEvents) {
		ReadEvent();
		switch(analysisOption) {
			case 0:
				EtaggAnalysis();
				break;
			case 1:
				EtaggAnalysis_FCAL();
				break;
			default:
				EtaggAnalysis();
				break;
		}
		m_event++;
	}
	m_inputFile->Close();
	
	return;
}

void EtaAna::InitHistograms(int analysisOption) {
	
	switch(analysisOption) {
		case 0:
		{
			// DEFAULT ANALYSIS HISTOGRAMS:
			
			h_mcVertex         = new TH1F("vertex",          "Vertex Z Position (unweighted)", 1000, 0., 100.);
			h_mcVertexAccepted = new TH1F("vertex_accepted", "Vertex Z Position (weighted)",   1000, 0., 100.);
			
			h_fcalRFdt     = new TH1F("fcal_rf_dt",     "t_{FCAL} - t_{RF}; [ns]", 10000, -100., 100.);
			h_bcalRFdt     = new TH1F("bcal_rf_dt",     "t_{BCAL} - t_{RF}; [ns]", 10000, -100., 100.);
			h_tofRFdt      = new TH1F( "tof_rf_dt",      "t_{TOF} - t_{RF}; [ns]", 10000, -100., 100.);
			h_scRFdt       = new TH1F(  "sc_rf_dt",       "t_{SC} - t_{RF}; [ns]", 10000, -100., 100.);
			
			h_beamRFdt     = new TH1F("beam_rf_dt",     "t_{CCAL} - t_{RF}; [ns]", 10000, -100., 100.);
			h_beamRFdt_cut = new TH1F("beam_rf_dt_cut", "t_{Beam} - t_{RF}; [ns]", 10000, -100., 100.);
			
			for(int iveto=0; iveto<m_nVetos; iveto++) {
				h_elasticity[iveto] = new TH2F(Form("elasticity_veto_%d",iveto),
					Form("Elasticity (Veto Option %d)",iveto), 650, 0.0, 6.5, 1000, 0.0, 2.0);
				h_elasticity[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_elasticity[iveto]->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
				
				h_elasticityConstr[iveto] = new TH2F(Form("elasticity_const_veto_%d",iveto),
					Form("Mass-Constrained Elasticity (Veto Option %d)",iveto), 650, 0.0, 6.5, 1000, 0.0, 2.0);
				h_elasticityConstr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_elasticityConstr[iveto]->GetYaxis()->SetTitle(
					"E_{#gamma#gamma}^{Constr}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
				
				h_mgg[iveto] = new TH2F(Form("mgg_veto_%d",iveto), 
					Form("Two-Photon Invariant Mass (Veto Option %d)",iveto), 550, 0., 5.5, 400, 0.3, 1.1);
				h_mgg[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_mgg[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV/c^{2}]");
				h_mgg[iveto]->Sumw2();
				
				h_mggConstr[iveto] = new TH2F(Form("mgg_const_veto_%d",iveto), 
					Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 550, 0., 5.5, 400, 0.3, 1.1);
				h_mggConstr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_mggConstr[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
				h_mggConstr[iveto]->Sumw2();
				
				h_mggConstr_coh[iveto] = new TH2F(Form("mgg_const_coh_veto_%d",iveto), 
					Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 550, 0., 5.5, 400, 0.3, 1.1);
				h_mggConstr_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_mggConstr_coh[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
				h_mggConstr_coh[iveto]->Sumw2();
				
				h_pt[iveto] = new TH2F(Form("pt_veto_%d",iveto), 
					Form("Transverse Momentum Difference (Veto Option %d)",iveto), 550, 0., 5.5, 1000, -0.5, 0.5);
				h_pt[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_pt[iveto]->GetYaxis()->SetTitle("p_{T} - p_{T}^{Calc} [GeV/c]");
				h_pt[iveto]->Sumw2();
				
				h_ptCoh[iveto] = new TH2F(Form("ptConstr_veto_%d",iveto), 
					Form("Transverse Momentum Differnce (Veto Option %d)",iveto), 550, 0., 5.5, 1000, -0.5, 0.5);
				h_ptCoh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_ptCoh[iveto]->GetYaxis()->SetTitle("p_{T} - p_{T}^{Calc} [GeV/c]");
				h_ptCoh[iveto]->Sumw2();
				
				h_mm[iveto] = new TH2F(Form("mm_veto_%d",iveto), 
					Form("Missing Mass (Veto Option %d)",iveto), 550, 0., 5.5, 1000, -20.0, 20.0);
				h_mm[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_mm[iveto]->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
				h_mm[iveto]->Sumw2();
				
				h_mm_coh[iveto] = new TH2F(Form("mm_coh_veto_%d",iveto), 
					Form("Missing Mass (Veto Option %d)",iveto), 550, 0., 5.5, 1000, -20.0, 20.0);
				h_mm_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				h_mm_coh[iveto]->GetYaxis()->SetTitle("#Deltam^{2} - m_{He}^{2} [GeV^{2}/c^{4}]");
				h_mm_coh[iveto]->Sumw2();
			}
			
			/*
			h_pT_vs_elas     = new TH2F("pT_vs_elas", 
				"; E_{#gamma#gamma} / E_{Calc}#left(E_{#gamma},#theta_{#gamma#gamma}#right); p_{T} / p_{T}^{Calc}",
				1000, 0.0, 2.0, 1000, 0.0, 2.0);
			h_pT_vs_elas_cut = new TH2F("pT_vs_elas_cut", 
				"; E_{#gamma#gamma} / E_{Calc}#left(E_{#gamma},#theta_{#gamma#gamma}#right); p_{T} / p_{T}^{Calc}",
				1000, 0.0, 2.0, 1000, 0.0, 2.0);
			*/
			
			h_scDeltaPhi = new TH2F("scDeltaPhi", 
				"#phi_{SC} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
				650, 0.0, 6.5, 2000, -360.0, 360.0);
			
			h_bcalDeltaPhi = new TH2F("bcalDeltaPhi", 
				"#phi_{BCAL} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
				650, 0.0, 6.5, 2000, -360.0, 360.0);
			
			h_xy1 = new TH2F("xy1", "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -100., 100., 500, -100., 100.);
			h_xy2 = new TH2F("xy2", "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 500, -100., 100., 500, -100., 100.);
			
			break;
		}
		case 1:
		{
			// VARY FCAL CUTS:
			
			h_FCAL_mgg = new TH2F("FCAL_mgg", "No Multiplicity Cut", 550, 0., 5.5, 400, 0.3, 1.1);
			h_FCAL_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
			h_FCAL_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
			h_FCAL_mgg->Sumw2();
			
			// with minimum energy cuts:
			
			h_FCAL_mggECut = new TH2F("FCAL_mgg_ecut", Form("E_{1,2} > %.2f GeV", m_FCALEnergyCut), 550, 0., 5.5, 400, 0.3, 1.1);
			h_FCAL_mggECut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
			h_FCAL_mggECut->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
			h_FCAL_mggECut->Sumw2();
			
			// with fiducial cuts:
			
			h_FCAL_mggFidCut = new TH2F("FCAL_mgg_fidcut", "Fiducial Cut Applied", 550, 0., 5.5, 400, 0.3, 1.1);
			h_FCAL_mggFidCut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
			h_FCAL_mggFidCut->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
			h_FCAL_mggFidCut->Sumw2();
			
			// with both energy and fiducial cuts:
			
			h_FCAL_mggCuts = new TH2F("FCAL_mgg_cuts", "Fiducial+Energy Cuts Applied", 550, 0., 5.5, 400, 0.3, 1.1);
			h_FCAL_mggCuts->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
			h_FCAL_mggCuts->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
			h_FCAL_mggCuts->Sumw2();
			
			// with 'good' multiplicity = 2:
			
			h_FCAL_mggGoodMult = new TH2F("FCAL_mgg_good_mult", "2 Good FCAL Showers", 550, 0., 5.5, 400, 0.3, 1.1);
			h_FCAL_mggGoodMult->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
			h_FCAL_mggGoodMult->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
			h_FCAL_mggGoodMult->Sumw2();
			
			// with total multiplicity = 2:
			
			h_FCAL_mggMult = new TH2F("FCAL_mgg_mult", "2 FCAL Showers", 550, 0., 5.5, 400, 0.3, 1.1);
			h_FCAL_mggMult->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
			h_FCAL_mggMult->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
			h_FCAL_mggMult->Sumw2();
			
			// vary the size of the fiducial cut:
			
			m_fcalFiducialCuts.clear();
			for(int icut=0; icut<13; icut++) {
				double locCut = 0.0 + 0.25*(double)(icut);
				m_fcalFiducialCuts.push_back(locCut);
			}
			
			for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
				TH2F *loc_h_mgg = new TH2F(Form("FCAL_mgg_fid_%02d", icut),
					Form("Inner layers removed by fiducial cut: %.1f", m_fcalFiducialCuts[icut]), 
					550, 0.0, 5.5, 400, 0.3, 1.1);
				loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
				loc_h_mgg->Sumw2();
				h_FCAL_mggFidCutVec.push_back(loc_h_mgg);
			}
			
			// vary the minimum energy cut:
			
			m_fcalEnergyCuts.clear();
			for(int icut=0; icut<20; icut++) {
				double locCut = 0.0 + 0.05*(double)(icut);
				m_fcalEnergyCuts.push_back(locCut);
			}
			
			for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
				TH2F *loc_h_mgg = new TH2F(Form("FCAL_mgg_ecut_%02d", icut),
					Form("Minimum Energy Cut: %.2f GeV", m_fcalEnergyCuts[icut]), 
					550, 0.0, 5.5, 400, 0.3, 1.1);
				loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
				loc_h_mgg->Sumw2();
				h_FCAL_mggECutVec.push_back(loc_h_mgg);
			}
			
			for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
				TH2F *loc_h_mgg = new TH2F(Form("FCAL_mgg_extra_ecut_%02d", icut),
					Form("Minimum Energy Cut: %.2f GeV", m_fcalEnergyCuts[icut]), 
					550, 0.0, 5.5, 400, 0.3, 1.1);
				loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
				loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
				loc_h_mgg->Sumw2();
				h_FCAL_mggExtraECutVec.push_back(loc_h_mgg);
			}
			
			break;
		}
	}
	
	return;
}

void EtaAna::ResetHistograms(int analysisOption) {
	
	switch(analysisOption) {
		case 0:
		{
			h_mcVertex->Reset();
			h_mcVertexAccepted->Reset();
			
			h_fcalRFdt->Reset();
			h_bcalRFdt->Reset();
			h_tofRFdt->Reset();
			h_scRFdt->Reset();
			h_beamRFdt->Reset();
			h_beamRFdt_cut->Reset();
			
			for(int iveto=0; iveto<m_nVetos; iveto++) {
				h_elasticity[iveto]->Reset();
				h_elasticityConstr[iveto]->Reset();
				h_mgg[iveto]->Reset();
				h_mggConstr[iveto]->Reset();
				h_mggConstr_coh[iveto]->Reset();
				h_pt[iveto]->Reset();
				h_ptCoh[iveto]->Reset();
				h_mm[iveto]->Reset();
				h_mm_coh[iveto]->Reset();
			}
			/*
			h_pT_vs_elas->Reset();
			h_pT_vs_elas_cut->Reset();
			*/
			h_scDeltaPhi->Reset();
			h_bcalDeltaPhi->Reset();
			
			h_xy1->Reset();
			h_xy2->Reset();
			
			break;
		}
		case 1:
		{
			h_FCAL_mgg->Reset();
			h_FCAL_mggECut->Reset();
			h_FCAL_mggFidCut->Reset();
			h_FCAL_mggCuts->Reset();
			h_FCAL_mggGoodMult->Reset();
			h_FCAL_mggMult->Reset();
			for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) h_FCAL_mggFidCutVec[icut]->Reset();
			for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
				h_FCAL_mggECutVec[icut]->Reset();
				h_FCAL_mggExtraECutVec[icut]->Reset();
			}
			
			if(h_AngularMatrix_FCALECut.size()) {
				for(int i=0; i<h_AngularMatrix_FCALECut.size(); i++) {
					h_AngularMatrix_FCALECut[i]->Reset();
				}
			}
			if(h_AngularMatrix_FCALExtraECut.size()) {
				for(int i=0; i<h_AngularMatrix_FCALExtraECut.size(); i++) {
					h_AngularMatrix_FCALExtraECut[i]->Reset();
				}
			}
			if(h_AngularMatrix_FCALFidCut.size()) {
				for(int i=0; i<h_AngularMatrix_FCALFidCut.size(); i++) {
					h_AngularMatrix_FCALFidCut[i]->Reset();
				}
			}
			
			break;
		}
	}
	
	if(h_thrown!=NULL) h_thrown->Reset();
	
	if(h_AngularMatrix.size()) {
		for(int i=0; i<h_AngularMatrix.size(); i++) {
			h_AngularMatrix[i]->Reset();
		}
	}
	
	if(h_invmassMatrix!=NULL) {
		h_invmassMatrix->Reset();
		h_invmassMatrix_acc->Reset();
	}
	
	return;
}

void EtaAna::WriteHistograms(int analysisOption) {
	
	cout << "writing histograms to: " << m_outputFileName << "..." << endl;
	TFile *fOut = new TFile(m_outputFileName.c_str(), "RECREATE");
	fOut->cd();
	
	switch(analysisOption) {
		case 0:
		{
			h_mcVertex->Write();
			h_mcVertexAccepted->Write();
			
			h_fcalRFdt->Write();
			h_bcalRFdt->Write();
			h_tofRFdt->Write();
			h_scRFdt->Write();
			h_beamRFdt->Write();
			h_beamRFdt_cut->Write();
			
			for(int iveto=0; iveto<m_nVetos; iveto++) {
				h_elasticity[iveto]->Write();
				h_elasticityConstr[iveto]->Write();
				h_mgg[iveto]->Write();
				h_mggConstr[iveto]->Write();
				h_mggConstr_coh[iveto]->Write();
				h_pt[iveto]->Write();
				h_ptCoh[iveto]->Write();
				h_mm[iveto]->Write();
				h_mm_coh[iveto]->Write();
			}
			/*
			h_pT_vs_elas->Write();
			h_pT_vs_elas_cut->Write();
			*/
			h_scDeltaPhi->Write();
			h_bcalDeltaPhi->Write();
			
			h_xy1->Write();
			h_xy2->Write();
			break;
		}
		case 1:
		{
			h_FCAL_mgg->Write();
			h_FCAL_mggECut->Write();
			h_FCAL_mggFidCut->Write();
			h_FCAL_mggCuts->Write();
			h_FCAL_mggGoodMult->Write();
			h_FCAL_mggMult->Write();
			
			TDirectory *dirECut = new TDirectoryFile("ECut", "ECut");
			dirECut->cd();
			for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
				h_FCAL_mggECutVec[icut]->Write();
			}
			if(h_AngularMatrix_FCALECut.size()) {
				for(int i=0; i<h_AngularMatrix_FCALECut.size(); i++) {
					h_AngularMatrix_FCALECut[i]->Write();
				}
			}
			dirECut->cd("../");
			
			TDirectory *dirExtraECut = new TDirectoryFile("ExtraECut", "ExtraECut");
			dirExtraECut->cd();
			for(int icut=0; icut<h_FCAL_mggExtraECutVec.size(); icut++) {
				h_FCAL_mggExtraECutVec[icut]->Write();
			}
			if(h_AngularMatrix_FCALExtraECut.size()) {
				for(int i=0; i<h_AngularMatrix_FCALExtraECut.size(); i++) {
					h_AngularMatrix_FCALExtraECut[i]->Write();
				}
			}
			dirExtraECut->cd("../");
			
			TDirectory *dirFidCut = new TDirectoryFile("FidCut", "FidCut");
			dirFidCut->cd();
			for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
				h_FCAL_mggFidCutVec[icut]->Write();
			}
			if(h_AngularMatrix_FCALFidCut.size()) {
				for(int i=0; i<h_AngularMatrix_FCALFidCut.size(); i++) {
					h_AngularMatrix_FCALFidCut[i]->Write();
				}
			}
			dirFidCut->cd("../");
			
			break;
		}
	}
	
	if(h_thrown!=NULL) h_thrown->Write();
	
	if(h_AngularMatrix.size()) {
		TDirectory *dirMatrix = new TDirectoryFile("AngularMatrix","AngularMatrix");
		dirMatrix->cd();
		for(int i=0; i<h_AngularMatrix.size(); i++) {
			h_AngularMatrix[i]->Write();
		}
		dirMatrix->cd("../");
	}
	
	if(h_invmassMatrix!=NULL) {
		h_invmassMatrix->Write();
		h_invmassMatrix_acc->Write();
	}
	
	fOut->Write();
	
	return;
}
