#include "EtaAna.h"

void EtaAna::EtaggAnalysis_FCAL() {
	
	double locThrownBeamEnergy, locThrownAngle;
	if(m_nmc>0) {
		if(AcceptRejectEvent()) return;
		GetThrownEnergyAndAngle(locThrownBeamEnergy, locThrownAngle);
		PlotThrown(locThrownBeamEnergy, locThrownAngle);
		if(h_AngularMatrix_FCALECut.size()==0) {
			InitializeAngularMatrices_FCAL();
		}
	}
	if(CheckEventMultiplicities()) {
		printf("    Skipping event %d\n", m_event);
		return;
	}
	
	vector<pair<int,double>> locGoodBeamPhotons; locGoodBeamPhotons.clear();
	int locNBeamPhotons = GetBeamPhotonList(locGoodBeamPhotons, m_minBeamEnergyCut, m_maxBeamEnergyCut);
	
	vector<int> locGoodFCALShowers; locGoodFCALShowers.clear();
	int locNFCALShowers_EnergyCut;
	int locNFCALShowers = GetFCALShowerList(locGoodFCALShowers, locNFCALShowers_EnergyCut, 0.0, m_FCALExtraEnergyCut, 0.0, m_FCALRFCut);
	
	vector<int> locGoodBCALShowers; locGoodBCALShowers.clear();
	int locNBCALShowers = GetBCALShowerList(locGoodBCALShowers, 0.0, m_BCALRFCut);
	
	int locNBCALShowers_1ns = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	
	for(int ishow=0; ishow<locGoodBCALShowers.size(); ishow++) {
		int showIndex = locGoodBCALShowers[ishow];
		TVector3 locPos = GetBCALPosition(showIndex);
		double locT = m_bcalT[showIndex] - (locPos.Mag()/m_c) - m_rfTime;
		
		locBCALRFDT = locT;
		locBCALPhi  = locPos.Phi() * TMath::RadToDeg();
		if(fabs(locT) < 1.0) {
			locNBCALShowers_1ns++;
		}
	}
	
	/*
	The vector 'locGoodFCALShowers' stores all FCAL showers within +/-2ns of the RF time.
	In the following code, we loop over different values of minimum energy cuts and determine whether or not 
	there are exactly 2 showers (no more) with that cut value:
	*/
	vector<int> locEnergyCuts; locEnergyCuts.clear();
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		
		int locGoodFCALShowers_cut = 0;
		for(int ishow=0; ishow<locGoodFCALShowers.size(); ishow++) {
			if((m_fcalE[locGoodFCALShowers[ishow]] > m_fcalEnergyCuts[icut]))
				locGoodFCALShowers_cut++;
		}
		if(locGoodFCALShowers_cut==2) locEnergyCuts.push_back(1);
		else locEnergyCuts.push_back(0);
	}
	
	//=====================================================================================//
	
	int locNGoodFCALShowers = (int)locGoodFCALShowers.size();
	if(locNGoodFCALShowers<2) return;
	
	for(int ishow=0; ishow<(locNGoodFCALShowers-1); ishow++) {
		
		int show1 = locGoodFCALShowers[ishow];
		TVector3 pos1 = GetFCALPosition(show1);
		
		double t1 = m_fcalT[show1] - (pos1.Mag()/m_c) - m_rfTime;
		double e1 = m_fcalE[show1];
		
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
			
			bool isHadronicVeto = true;
			
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
			
			if((locNBCALShowers==0) || 
				((locNBCALShowers==1) && (fabs(fabs(locBCALPhi-prodPhi)-180.0) < m_BCALDeltaPhiCut) && (locBCALRFDT>1.0))) {
				
				if(locNSCHits_coplanar==locNSCHits) isHadronicVeto = false;
			}
			if(isHadronicVeto) continue;
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(int igam=0; igam<locNBeamPhotons; igam++) {
				
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
				bool isElastic = IsElasticCut(Egg, etaEnergy, prodTheta);
				
				// set a variable to indicate if the two-photon mass is consistent with an eta meson:
				bool isEta = IsEtaCut(invmass);
				
				if(!isElastic) continue;
				
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
				
				//-----------------------------------------------------//
				
				h_FCAL_mgg->Fill(prodTheta, invmassConstr, fillWeight);
				if((e1 > m_FCALEnergyCut) && (e2 > m_FCALEnergyCut)) {
					h_FCAL_mggECut->Fill(prodTheta, invmassConstr, fillWeight);
				}
				if(!FCALFiducialCut(pos1, 2.0) && !FCALFiducialCut(pos2, 2.0)) {
					h_FCAL_mggFidCut->Fill(prodTheta, invmassConstr, fillWeight);
					if((e1 > m_FCALEnergyCut) && (e2 > m_FCALEnergyCut)) {
						h_FCAL_mggCuts->Fill(prodTheta, invmassConstr, fillWeight);
						if(locNFCALShowers_EnergyCut==2) {
							h_FCAL_mggGoodMult->Fill(prodTheta, invmassConstr, fillWeight);
							if(locNFCALShowers==2) {
								h_FCAL_mggMult->Fill(prodTheta, invmassConstr, fillWeight);
							}
						}
					}
				}
				
				if(locNFCALShowers_EnergyCut==2) {
					// Vary the size of the fiducial cut used:
					if((e1 > m_FCALEnergyCut) && (e2 > m_FCALEnergyCut)) {
						for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
							if(!FCALFiducialCut(pos1, m_fcalFiducialCuts[icut]) && !FCALFiducialCut(pos2, m_fcalFiducialCuts[icut])) {
								h_FCAL_mggFidCutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
								if(isEta && (m_nmc>0)) {
									h_AngularMatrix_FCALFidCut[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
								}
							}
						}
					}
					
					// Vary the minimum energy cut used:
					if(!FCALFiducialCut(pos1, 2.0) && !FCALFiducialCut(pos2, 2.0)) {
						for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
							if(e1>m_fcalEnergyCuts[icut] && e2>m_fcalEnergyCuts[icut]) {
								h_FCAL_mggECutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
								if(isEta && (m_nmc>0)) {
									h_AngularMatrix_FCALECut[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
								}
							}
						}
					}
				}
				
				// Vary the minimum energy cut used for looking for extra showers:
				if(!FCALFiducialCut(pos1,2.0) && !FCALFiducialCut(pos2,2.0) && (e1>m_FCALEnergyCut) && (e2>m_FCALEnergyCut)) {
					for(int icut=0; icut<locEnergyCuts.size(); icut++) {
						if(locEnergyCuts[icut]==1) {
							h_FCAL_mggExtraECutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
							if(isEta && (m_nmc>0)) {
								h_AngularMatrix_FCALExtraECut[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
							}
						}
					}
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

void EtaAna::InitializeAngularMatrices_FCAL() {
	
	double minBeamEnergy      =  7.0;
	double maxBeamEnergy      = 12.0;
	double beamEnergyBinSize  =  0.1;
	
	double minRecAngle        = 0.0;
	double maxRecAngle        = 5.5;
	double recAngleBinSize    = 0.01;
	
	double minThrownAngle     = 0.0;
	double maxThrownAngle     = 5.0;
	double thrownAngleBinSize = 0.01;
	
	int nBeamEnergyBins  = (int)((maxBeamEnergy-minBeamEnergy)/beamEnergyBinSize);
	int nRecAngleBins    = (int)((maxRecAngle-minRecAngle)/recAngleBinSize);
	int nThrownAngleBins = (int)((maxThrownAngle-minThrownAngle)/thrownAngleBinSize);
	
	for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
		TH3F *hMatrix = new TH3F(Form("AngularMatrix_FCALFidCut_%02d",icut), 
			Form("Inner layers removed by fiducial cut: %.1f; #theta(thrown) [#circ]; #theta(rec) [#circ]", m_fcalFiducialCuts[icut]),
			nThrownAngleBins, minThrownAngle, maxThrownAngle, 
			nRecAngleBins,    minRecAngle,    maxRecAngle,
			nBeamEnergyBins,  minBeamEnergy,  maxBeamEnergy);
		hMatrix->SetDirectory(0);
		h_AngularMatrix_FCALFidCut.push_back(hMatrix);
	}
	
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		TH3F *hMatrix = new TH3F(Form("AngularMatrix_FCALECut_%02d",icut), 
			Form("Minimum Energy Cut: %.2f GeV; #theta(thrown) [#circ]; #theta(rec) [#circ]", m_fcalEnergyCuts[icut]),
			nThrownAngleBins, minThrownAngle, maxThrownAngle, 
			nRecAngleBins,    minRecAngle,    maxRecAngle,
			nBeamEnergyBins,  minBeamEnergy,  maxBeamEnergy);
		hMatrix->SetDirectory(0);
		h_AngularMatrix_FCALECut.push_back(hMatrix);
	}
	
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		TH3F *hMatrix = new TH3F(Form("AngularMatrix_FCALExtraECut_%02d",icut), 
			Form("Extra Energy Cut: %.2f GeV; #theta(thrown) [#circ]; #theta(rec) [#circ]", m_fcalEnergyCuts[icut]),
			nThrownAngleBins, minThrownAngle, maxThrownAngle, 
			nRecAngleBins,    minRecAngle,    maxRecAngle,
			nBeamEnergyBins,  minBeamEnergy,  maxBeamEnergy);
		hMatrix->SetDirectory(0);
		h_AngularMatrix_FCALExtraECut.push_back(hMatrix);
	}
	
	return;
}
