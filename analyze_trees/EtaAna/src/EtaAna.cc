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
	m_FCALExtraEnergyCut = 0.05; // threshold for 'extra' showers in the FCAL
		// There should only be two photons with energy > m_FCALExtraEnergyCut. However, we
		// can vary m_FCALEnergyCut separately from m_FCALExtraEnergyCut in order to
		// change the energy threshold of the showers which are used for reconstructing etas,
		// without having to include "extra" higher energy showers.
	
	m_BCALEnergyCut = 0.10;
	
	m_minBeamEnergyCut =  9.0;
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
	
	// for event-by-event acceptance correction:
	
	h_acceptance = nullptr;
	
	// initialize as NULL pointers:
	
	h_thrown                  = nullptr;
	
	h_invmassMatrix           = nullptr;
	h_invmassMatrix_prompt    = nullptr;
	h_invmassMatrix_acc       = nullptr;
	h_AngularMatrix           = nullptr;
	
	h_AngularMatrix_noTOF     = nullptr;
	h_AngularMatrix_TOF       = nullptr;
	h_AngularMatrix_singleTOF = nullptr;
	
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
	else if(m_phaseVal==2) {
		
		if(m_runNumber<81396) {
			
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
		
		if(m_tree->GetListOfBranches()->FindObject("thrownBeamEnergy")) {
			m_tree->SetBranchAddress("thrownBeamEnergy",   &m_thrownBeamEnergy);
		}
		else {
			m_thrownBeamEnergy = 0.0;
		}
	}
	
	m_tree->GetEvent(m_event);
	
	return;
}

void EtaAna::PlotThrown(double energy, double angle) {
	
	if(m_nmc==0) return;
	if(h_thrown==nullptr) {
		
		int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
		int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
		
		h_thrown = new TH2F("thrown", 
			"Thrown Angle vs. Thrown Beam Energy; E_{#gamma}(thrown) [GeV]; #theta(thrown) [#circ]",
			nBeamEnergyBins, m_minBeamEnergyBin, m_maxBeamEnergyBin, 
			nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin);
		h_thrown->SetDirectory(0);
	}
	h_thrown->Fill(energy, angle);
	return;
}

void EtaAna::FillAngularMatrix(double thrownEnergy, double thrownAngle, 
	double recAngle, double weight) {
	
	if((m_minBeamEnergyBin>thrownEnergy) || (thrownEnergy>m_maxBeamEnergyBin)) return;
	
	// check if matrix has been initialized:
	if(h_AngularMatrix==nullptr) {
		return;
	}
	
	// find the index associated with this beam energy:
	h_AngularMatrix->Fill(thrownAngle, recAngle, thrownEnergy, weight);
	
	return;
}

void EtaAna::FillAngularMatrix_vetos(int vetoOption, double thrownEnergy, double thrownAngle, 
	double recAngle, double weight) {
	
	if((m_minBeamEnergyBin>thrownEnergy) || (thrownEnergy>m_maxBeamEnergyBin)) return;
	
	// check if matrix has been initialized:
	if(h_AngularMatrix_vetos.size()==0) {
		return;
	}
	
	// find the index associated with this beam energy:
	h_AngularMatrix_vetos[vetoOption]->Fill(thrownAngle, recAngle, thrownEnergy, weight);
	
	return;
}

void EtaAna::FillInvmassMatrix(double theta, double mgg, double beamEnergy, double weight) {
	
	h_invmassMatrix->Fill(theta, mgg, beamEnergy, weight);
	
	if(weight<0.0) h_invmassMatrix_acc->Fill(theta, mgg, beamEnergy);
	else h_invmassMatrix_prompt->Fill(theta, mgg, beamEnergy);
	
	return;
}

int EtaAna::GetFCALShowerList(vector<int> &goodShowers, int &nGoodFCALShowers, 
	double energyCut, double extraEnergyCut, double fiducialCut, double timingCut) {
	
	int nFCALShowers = 0;
	nGoodFCALShowers = 0;
	for(int ishow=0; ishow<m_nfcal; ishow++) {
		
		TVector3 locPos = GetFCALPosition(ishow);
		double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		
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

int EtaAna::GetSCHitList(vector<int> &goodHits) {
	
	// Timing and dE cuts are hard coded for now
	
	int nSCHits = 0;
	for(int ihit=0; ihit<m_nsc; ihit++) {
		
		// only check hits between 1ns < (t_sc - t_RF) < 7ns 
		//    and with dE > 0.0002 (from DNeutralShower_factory)
		
		double locT  = m_scT[ihit] - m_rfTime;
		double locdE = m_scdE[ihit];
		
		if((1.0 < locT) && (locT < 9.0) && (locdE > 0.0002)) {
			nSCHits++;
			goodHits.push_back(ihit);
		}
	}
	
	return nSCHits;
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
		
		if((locWeight < 0.0) && (m_nmc==0)) {
			locWeight *= m_accScaleFactor[igam];
		}
		if((minEnergyCut < m_beamE[igam]) && (m_beamE[igam] < maxEnergyCut)) {
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
	//if(m_nmc!=2) return;
	
	/*
	The following code calculates the thrown beam energy assuming an eta 
	was produced coherently on a He-4 nucleus at rest. However, this is only accurate
	for the coherent grid simulation. Not the bggen or other background simulations.
	That's why the last line of this function reads the thrown beam energy from a dedicated 
	branch. This is how it should be done for all simulations, but when writing the trees initially,
	this branch was not included. 
	
	In the future, when the trees are updated to include this branch, we should only use the last
	line for getting the thrown beam energy.
	*/
	
	thrownEnergy = -1.0*ParticleMass(Helium);
	thrownAngle  =  0.0;
	for(int imc=0; imc<m_nmc; imc++) {
		thrownEnergy += m_mcE[imc];
		if(m_mcPDGType[imc]==PDGtype(Eta)) thrownAngle = m_mcTheta[imc];
	}
	
	if(m_thrownBeamEnergy>1.0) thrownEnergy = m_thrownBeamEnergy;
	
	return;
}

void EtaAna::GetThrownEnergyAndAngleBGGEN(double &thrownEnergy, double &thrownAngle) {
	
	// in principle all simulations should be analyzed in this way, but at the 
	// time of writing this (2/18/25) only the BGGEN mc trees have been updated to 
	// include the "thrownBeamEnergy" branch.
	
	thrownAngle  =  0.0;
	for(int imc=0; imc<m_nmc; imc++) {
		if(m_mcPDGType[imc]==PDGtype(Eta)) thrownAngle = m_mcTheta[imc];
	}
	thrownEnergy = m_thrownBeamEnergy;
	return;
}

bool EtaAna::IsElasticCut(double Egg, double Eeta, double theta) {
	double ElasPeakMean  = m_ElasMean_p0 + m_ElasMean_p1*theta;
	double ElasPeakWidth = m_ElasWidth * m_ElasSigmaCut;
	if(fabs((Egg/Eeta)-ElasPeakMean)<ElasPeakWidth) return true;
	else return false;
}

bool EtaAna::IsEtaCut(double invmass) {
	if((0.5<=invmass) && (invmass<0.6)) return true;
	else return false;
}

bool EtaAna::IsCoplanarBCAL(double deltaPhi) {
	
	// deltaPhi = #phi_gg - #phi_BCAL
	
	if(m_phaseVal>1) {
		if(((140.0<deltaPhi) && (deltaPhi<240.0)) || ((-220.0<deltaPhi) && (deltaPhi<-120.0)))
			return true;
		else 
			return false;
	}
	else {
		if(fabs(fabs(deltaPhi)-180.0) < m_BCALDeltaPhiCut)
			return true;
		else
			return false;
	}
}

bool EtaAna::IsCoplanarSC(double deltaPhi) {
	
	// deltaPhi = #phi_gg - #phi_SC
	
	if(fabs(fabs(deltaPhi)-180.0) < m_SCDeltaPhiCut)
		return true;
	else
		return false;
}

void EtaAna::SmearShowerEnergy(double &e) {
	// smear the cluster energy by 10%:
	double locSig = 0.1*e;
	e += m_random->Gaus(0.0, locSig);
	return;
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

int EtaAna::GetAcceptanceHistogram() {
	
	if(h_acceptance!=nullptr) return 0;
	
	TString accFileName = "/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles/grid_acceptance.root";
	if(gSystem->AccessPathName(accFileName.Data())) return 1;
	
	TFile *fAcc = new TFile(accFileName.Data(), "READ");
	h_acceptance = (TH2F*)fAcc->Get("grid_acceptance");
	h_acceptance->SetDirectory(0);
	fAcc->Close();
	
	return 0;
}

void EtaAna::RunAnalysis(TString inputFileName, int analysisOption) {
	
	m_inputFile = new TFile(inputFileName.Data(), "READ");
	int nTotalEvents = LoadTree();
	
	/*
	// try to read in grid-acceptance:
	if(GetAcceptanceHistogram()) {
		std::cout << "Unable to get grid acceptance from ROOT file" << std::endl;
	}
	*/
	
	while(m_event < nTotalEvents) {
		ReadEvent();
		if(CheckEventMultiplicities()) {
			printf("    Skipping event %d\n", m_event);
			m_event++;
			continue;
		}
		switch(analysisOption) {
			case 0:
				EtaggAnalysis();
				break;
			case 1:
				EtaggAnalysis_matrix();
				break;
			case 2:
				EtaggAnalysis_FCAL();
				break;
			case 5:
				EtaggAnalysis_TOF();
				break;
			case 6:
				EtaggAnalysis_bggen();
				break;
			case 7:
				Omega3gAnalysis();
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
	
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	switch(analysisOption) {
		case 0:
		{
			InitializeDefaultHists();
			break;
		}
		case 1:
		{
			InitializeMatrixHists();
			break;
		}
		case 2:
		{
			InitializeFCALHists();
			break;
		}
		case 5:
		{
			InitializeTOFHists();
			break;
		}
		case 6:
		{
			InitializeReactionTypes();
			InitializeBGGENHists();
			break;
		}
		case 7:
		{
			InitializeOmegaHists();
			break;
		}
	}
	
	return;
}

void EtaAna::ResetHistograms(int analysisOption) {
	
	switch(analysisOption) {
		case 0:
		{
			ResetDefaultHists();
			break;
		}
		case 1:
		{
			ResetMatrixHists();
			break;
		}
		case 2:
		{
			ResetFCALHists();
			break;
		}
		case 5:
		{
			ResetTOFHists();
			break;
		}
		case 6:
		{
			ResetBGGENHists();
			break;
		}
		case 7:
		{
			ResetOmegaHists();
			break;
		}
	}
	
	if(h_thrown!=nullptr) h_thrown->Reset();
	
	return;
}

void EtaAna::WriteHistograms(int analysisOption) {
	
	cout << "writing histograms to: " << m_outputFileName << "..." << endl;
	TFile *fOut = new TFile(m_outputFileName.c_str(), "RECREATE");
	fOut->cd();
	
	switch(analysisOption) {
		case 0:
		{
			WriteDefaultHists();
			break;
		}
		case 1:
		{
			WriteMatrixHists();
			break;
		}
		case 2:
		{
			WriteFCALHists();
			break;
		}
		case 5:
		{
			WriteTOFHists();
			break;
		}
		case 6:
		{
			WriteBGGENHists();
			break;
		}
		case 7:
		{
			WriteOmegaHists();
			break;
		}
	}
	
	if(h_thrown!=nullptr) {
		printf("\n  Writing thrown histogram...\n");
		h_thrown->Write();
		printf("  Done.\n");
	}
	
	printf("\n  Writing ROOT file...\n");
	fOut->Write();
	printf("  Done.\n");
	
	return;
}
