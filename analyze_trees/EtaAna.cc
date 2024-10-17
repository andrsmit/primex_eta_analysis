#include "EtaAna.h"

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
	
	m_FCALEnergyCut = 0.5;
	m_BCALEnergyCut = 0.1;
	
	m_minBeamEnergyCut =  8.0;
	m_maxBeamEnergyCut = 10.9;
	
	m_FCALTOFCut    = 8.0;
	m_SCDeltaPhiCut = 36.0;
	
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
	
	// set up FCAL channel number array:
	int num_active_blocks = 0;
	for(int row = 0; row < 59; row++){
		for(int col = 0; col < 59; col++){
			
			// transform to beam axis
			m_fcalPositionOnFace[row][col] = TVector2((col - 29) * 4.0157, (row - 29) * 4.0157);
			
			double thisRadius = m_fcalPositionOnFace[row][col].Mod();
			
			if((thisRadius < 120.471) && (thisRadius > 5.73585)){
				
				// build the "channel map"
				m_fcalChannelNumber[row][col]   = num_active_blocks;
				m_fcalRow[num_active_blocks]    = row;
				m_fcalColumn[num_active_blocks] = col;
				
				num_active_blocks++;
			}
		}
	}
}

void EtaAna::RunAnalysis(TString inputFileName) {
	
	m_inputFile = new TFile(inputFileName.Data(), "READ");
	int nTotalEvents = LoadTree();
	
	while(m_event < nTotalEvents) {
		ReadEvent();
		EtaggAnalysis();
		m_event++;
	}
	m_inputFile->Close();
	
	return;
}

int EtaAna::CheckEventMultiplicities() {
	
	int nSystemsOver = 0;
	if(m_nfcal>MAX_FCAL) nSystemsOver++;
	if(m_nbcal>MAX_BCAL) nSystemsOver++;
	if(m_nbeam>MAX_BEAM) nSystemsOver++;
	if(m_ntof >MAX_TOF ) nSystemsOver++;
	if(m_nsc  >MAX_SC  ) nSystemsOver++;
	return nSystemsOver;
}

void EtaAna::EtaggAnalysis() {
	
	if(m_nmc>0) {
		h_mcVertex->Fill(m_mcZ[0]);
		if(AcceptRejectEvent()) return;
		h_mcVertexAccepted->Fill(m_mcZ[0]);
	}
	
	if(CheckEventMultiplicities()) {
		printf("Skipping event %d (error code: %d)\n", m_event, CheckEventMultiplicities());
		return;
	}
	
	//
	// Make a list of good FCAL showers to use in analysis:
	//
	vector<int> locGoodFCALShowers; locGoodFCALShowers.clear();
	int locNFCALShowers = 0, locNGoodFCALShowers = 0;
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
		
		h_fcalRFdt->Fill(locT);
		
		if(fabs(locT) < m_FCALRFCut) {
			locNFCALShowers++;
			if((m_fcalE[ishow] > m_FCALEnergyCut) && !FCALFiducialCut(locPos, 2.0)) {
				locNGoodFCALShowers++;
				locGoodFCALShowers.push_back(ishow);
			}
		}
	}
	
	//
	// Make a list of good BCAL showers to use in analysis:
	//
	vector<int> locGoodBCALShowers; locGoodBCALShowers.clear();
	int locNBCALShowers = 0, locNBCALShowers_1ns = 0;
	double locBCALEnergySum = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	for(int ishow=0; ishow<m_nbcal; ishow++) {
		
		TVector3 locPos = GetBCALPosition(ishow);
		double locT = m_bcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		
		h_bcalRFdt->Fill(locT);
		
		if(fabs(locT) < m_BCALRFCut) {
			locBCALEnergySum += m_bcalE[ishow];
			locNBCALShowers++;
			locBCALRFDT = locT;
			locBCALPhi  = locPos.Phi() * TMath::RadToDeg();
			locGoodBCALShowers.push_back(ishow);
			if(fabs(locT) < 1.0) {
				locNBCALShowers_1ns++;
			}
		}
	}
	
	// Make a list of prompt and selected-sideband beam photons to use in analysis:
	
	vector<pair<int,double>> locGoodBeamPhotons;
	locGoodBeamPhotons.clear();
	for(int igam=0; igam<m_nbeam; igam++) {
		
		double locRFdt   = m_beamT[igam] - m_rfTime;
		double locWeight = 0.0;
		
		h_beamRFdt->Fill(locRFdt);
		
		double locBeamCut = beamBunchesMain*m_BeamRFCut;
		
		if(fabs(locRFdt) < locBeamCut) locWeight = 1.0;
		else if(
			((beamBunchesMain+5.5)*4.008 < fabs(locRFdt)) && 
			(fabs(locRFdt) < (beamBunchesMain+5.5+beamBunchesAcc)*4.008)
		) locWeight = -1.0/(2.0*beamBunchesAcc);
		else continue;
		
		if(locWeight < 0.0) locWeight *= m_accScaleFactor[igam];
		if((m_beamE[igam] > m_minBeamEnergyCut) && (m_beamE[igam] < m_maxBeamEnergyCut)) {
			locGoodBeamPhotons.push_back({igam,locWeight});
		}
	}
	
	if(locNFCALShowers < 2) return;
	
	//=====================================================================================//
	
	// Apply multiplicity cut on the number of FCAL showers: 
	if((locNFCALShowers!=2) || (locNGoodFCALShowers!=2)) return;
	
	for(int ishow=0; ishow<(locNGoodFCALShowers-1); ishow++) {
		
		int show1 = locGoodFCALShowers[ishow];
		TVector3 pos1 = GetFCALPosition(show1);
		
		double t1 = m_fcalT[show1] - (pos1.Mag()/m_c) - m_rfTime;
		double e1 = m_fcalE[show1];
		
		// apply minimum energy and timing cuts (this is redundant):
		if((e1 < m_FCALEnergyCut) || (fabs(t1) >= m_FCALRFCut)) continue;
		
		// apply fiducial cut (redundant):
		if(FCALFiducialCut(pos1, 2.0)) continue;
		
		double px1 = e1*pos1.X() / pos1.Mag();
		double py1 = e1*pos1.Y() / pos1.Mag();
		double pz1 = e1*pos1.Z() / pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		CheckTOFMatch(pos1, tof_dx1, tof_dy1, tof_dt1, m_TOFRFCut);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		for(int jshow=(ishow+1); jshow<locNGoodFCALShowers; jshow++) {
			
			int show2 = locGoodFCALShowers[jshow];
			TVector3 pos2 = GetFCALPosition(show2);
			
			double t2 = m_fcalT[show2] - (pos2.Mag()/m_c) - m_rfTime;
			double e2 = m_fcalE[show2];
			
			// apply minimum energy and timing cuts (this is redundant):
			if((e2 < m_FCALEnergyCut) || (fabs(t2) >= m_FCALRFCut)) continue;
			
			// apply fiducial cut (redundant):
			if(FCALFiducialCut(pos2, 2.0)) continue;
			
			double px2 = e2*pos2.X() / pos2.Mag();
			double py2 = e2*pos2.Y() / pos2.Mag();
			double pz2 = e2*pos2.Z() / pos2.Mag();
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			CheckTOFMatch(pos2, tof_dx2, tof_dy2, tof_dt2, m_TOFRFCut);
			double tof_dr2 = sqrt(pow(tof_dx2,2.0)+pow(tof_dy2,2.0));
			
			//-----------------------------------------------------//
			// TOF Veto
			
			// reject combinations of FCAL showers where both showers are near a TOF hit:
			bool isTOFVeto = false;
			if((tof_dr1 < m_FCALTOFCut) && (tof_dr2 < m_FCALTOFCut)) isTOFVeto = true;
			
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
			double prodTheta = atan2(pggt,pggz) * TMath::RadToDeg();
			
			// azimuthal angle:
			double prodPhi = atan2(pggy,pggx) * TMath::RadToDeg();
			
			// opening angle:
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
			
			// invariant mass:
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			//-----------------------------------------------------//
			// check different veto options:
			
			// Check SC Matches:
			
			int locNSCHits = 0;
			int locNSCHits_coplanar = 0;
			for(int isc = 0; isc < m_nsc; isc++) {
				
				// only check hits between 1ns < (t_sc - t_RF) < 7ns 
				//    and with dE > 0.0002 (from DNeutralShower_factory)
				
				double locT  = m_scT[isc] - m_rfTime;
				double locdE = m_scdE[isc];
				
				if((1.0 < locT) && (locT < 7.0) && (locdE > 0.0002)) {
					locNSCHits++;
					if(fabs(fabs(m_scPhi[isc]-prodPhi)-180.0) < m_SCDeltaPhiCut) locNSCHits_coplanar++;
				}
			}
			
			vector<int> locVetoOptions; locVetoOptions.clear();
			for(int iveto=0; iveto<m_nVetos; iveto++) locVetoOptions.push_back(0);
			
			// Option 0 (no veto): No Veto is applied:
			locVetoOptions[0] = 1;
			
			// Option 1 (strict): Remove events with any BCAL shower within +/-12ns:
			if(locNBCALShowers==0) locVetoOptions[1] = 1;
			
			// Option 2 (loose): Remove events with any BCAL shower within +/-1ns:
			if(locNBCALShowers_1ns==0) locVetoOptions[2] = 1;
			
			// Option 3 (looser): Keep events where there is EITHER (i) no BCAL shower within +/-12ns, 
			//        OR (ii) 1 BCAL shower that has opposite phi angle to two-photon pair in FCAL:
			if((locNBCALShowers==0) ||
				((locNBCALShowers==1) && (fabs(fabs(locBCALPhi-prodPhi)-180.0) < m_BCALDeltaPhiCut))) locVetoOptions[3] = 1;
			
			// Option 4 (improvement on option 3): Keep events where there is EITHER (i) no BCAL shower within +/-12ns, 
			//        OR (ii) 1 BCAL shower that has opposite phi angle to two-photon pair in FCAL 
			//        and is more than 1ns removed from RF time:
			if((locNBCALShowers==0) || 
				((locNBCALShowers==1) && (fabs(fabs(locBCALPhi-prodPhi)-180.0) < m_BCALDeltaPhiCut) && (locBCALRFDT>1.0))) {
				locVetoOptions[4] = 1;
				
				// Option 5 (add in SC Veto): Remove events where there is a hit in the SC outside of the range:
				//          150 < |phi_SC - phi_FCAL| < 210:
				if(locNSCHits_coplanar==locNSCHits) locVetoOptions[5] = 1;
			}
			
			// Option 6: Use tight veto on SC only:
			if(locNSCHits==0) locVetoOptions[6] = 1;
			
			// Option 7: tight SC + tight BCAL vetos:
			if(locNSCHits==0 && locNBCALShowers==0) locVetoOptions[7] = 1;
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(int igam=0; igam<(int)locGoodBeamPhotons.size(); igam++) {
				
				int ibeam = locGoodBeamPhotons[igam].first;
				double fillWeight = locGoodBeamPhotons[igam].second;
				
				double eb    = m_beamE[ibeam];
				double brfdt = m_beamT[ibeam] - m_rfTime;
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				double etaEnergy_coh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
				
				// Calculate the energy of the eta meson, assuming production on a free nucleon:
				double etaEnergy = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(Proton));
				
				// Apply a cut on the elasticity
				//  (ratio of measured energy of 2-photons, to the calculated energy above):
				
				bool isElastic = false;
				double locElasMean  = m_ElasMean_p0 + m_ElasMean_p1*prodTheta;
				double locElasWidth = m_ElasWidth * m_ElasSigmaCut;
				if(fabs((Egg/etaEnergy)-locElasMean)<locElasWidth) isElastic = true;
				
				// set a variable to indicate if the two-photon mass is consistent with an eta meson:
				bool isEta = false;
				if(0.497862<invmass && invmass<0.597862) isEta = true;
				
				// Plot timing distribution of beam photons after elasticity cut to see the level of accidentals:
				if(isElastic && locVetoOptions[5]) {
					h_beamRFdt_cut->Fill(brfdt);
				}
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a coherently-produced eta meson:
				
				double sig1 = GetFCALEnergyResolution(e1);
				double sig2 = GetFCALEnergyResolution(e2);
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
				double prod_th_const = atan2(pggtc,pggzc) * TMath::RadToDeg();
				
				if(isEta && isElastic && locVetoOptions[1]) {
					h_xy1->Fill(pos1.X(),pos1.Y());
					h_xy2->Fill(pos2.X(),pos2.Y());
				}
				
				for(int iveto=0; iveto<m_nVetos; iveto++) {
					if(locVetoOptions[iveto]==0) continue;
					
					if(isEta) h_elasticity[iveto]->Fill(prodTheta, Egg/etaEnergy, fillWeight);
					if(isElastic) {
						h_mgg[iveto]->Fill(prodTheta, invmass, fillWeight);
						h_mggConstr[iveto]->Fill(prodTheta, invmassConstr, fillWeight);
						h_mggConstr_coh[iveto]->Fill(prodTheta, invmassConstr_coh, fillWeight);
						if(isEta) {
							h_mggConstr_etaCut[iveto]->Fill(prodTheta, invmassConstr, fillWeight);
						}
					}
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
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
	else {
		std::cout << "Unsupported run period provided. Skipping Run." << std::endl;
		return 1;
	}
	
	return 0;
}

void EtaAna::SetOutputFileName(string fileName) { m_outputFileName = fileName; return; }

int EtaAna::LoadTree() {
	
	m_tree = (TTree*)m_inputFile->Get("eta_gg");
	if(m_tree==NULL) return 0;
	
	// Reset event count to zero when laoding a new Tree:
	m_event = 0;
	
	return m_tree->GetEntries();
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

void EtaAna::InitHistograms() {
	
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
		
		h_mgg[iveto] = new TH2F(Form("mgg_veto_%d",iveto), 
			Form("Two-Photon Invariant Mass (Veto Option %d)",iveto), 650, 0., 6.5, 600, 0., 1.2);
		h_mgg[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg[iveto]->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
		h_mgg[iveto]->Sumw2();
		
		h_mggConstr[iveto] = new TH2F(Form("mgg_const_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 650, 0., 6.5, 600, 0., 1.2);
		h_mggConstr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr[iveto]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr[iveto]->Sumw2();
		
		h_mggConstr_coh[iveto] = new TH2F(Form("mgg_const_coh_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 650, 0., 6.5, 600, 0., 1.2);
		h_mggConstr_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_coh[iveto]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr_coh[iveto]->Sumw2();
		
		h_mggConstr_etaCut[iveto] = new TH2F(Form("mgg_const_etacut_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 650, 0., 6.5, 600, 0., 1.2);
		h_mggConstr_etaCut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_etaCut[iveto]->GetYaxis()->SetTitle("M_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr_etaCut[iveto]->Sumw2();
	}
	
	h_xy1 = new TH2F("xy1", "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -100., 100., 500, -100., 100.);
	h_xy2 = new TH2F("xy2", "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 500, -100., 100., 500, -100., 100.);
	
	return;
}

void EtaAna::ResetHistograms() {
	
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
		h_mgg[iveto]->Reset();
		h_mggConstr[iveto]->Reset();
		h_mggConstr_coh[iveto]->Reset();
		h_mggConstr_etaCut[iveto]->Reset();
	}
	
	h_xy1->Reset();
	h_xy2->Reset();
	
	return;
}

void EtaAna::WriteHistograms() {
	
	cout << "writing histograms to: " << m_outputFileName << "..." << endl;
	TFile *fOut = new TFile(m_outputFileName.c_str(), "RECREATE");
	fOut->cd();
	
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
		h_mgg[iveto]->Write();
		h_mggConstr[iveto]->Write();
		h_mggConstr_coh[iveto]->Write();
		h_mggConstr_etaCut[iveto]->Write();
	}
	
	h_xy1->Write();
	h_xy2->Write();
	
	fOut->Write();
	
	return;
}
