#include "EtaAna.h"

void EtaAna::EtaggAnalysis_TOF() {
	
	if(m_nmc>0) {
		if(AcceptRejectEvent()) return;
	}
	
	double locThrownBeamEnergy = 0.0, locThrownAngle = 0.0;
	if(m_FillThrown) {
		GetThrownEnergyAndAngle(locThrownBeamEnergy, locThrownAngle);
		//if((locThrownBeamEnergy < m_minBeamEnergyCut) || (locThrownBeamEnergy >= m_maxBeamEnergyCut)) return;
		PlotThrown(locThrownBeamEnergy, locThrownAngle);
	}
	
	vector<pair<int,double>> locGoodBeamPhotons; locGoodBeamPhotons.clear();
	int locNBeamPhotons = GetBeamPhotonList(locGoodBeamPhotons, m_minBeamEnergyCut, m_maxBeamEnergyCut);
	
	vector<int> locGoodFCALShowers; locGoodFCALShowers.clear();
	int locNFCALShowers_EnergyCut;
	int locNFCALShowers = GetFCALShowerList(locGoodFCALShowers, locNFCALShowers_EnergyCut, m_FCALEnergyCut, m_FCALExtraEnergyCut, 2.0, m_FCALRFCut);
	
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
	
	//=====================================================================================//
	
	// Apply multiplicity cut on the number of FCAL showers: 
	int locNGoodFCALShowers = (int)locGoodFCALShowers.size();
	if((locNFCALShowers_EnergyCut!=2) || (locNGoodFCALShowers!=2)) return;
	
	for(int ishow=0; ishow<(locNGoodFCALShowers-1); ishow++) {
		
		int show1 = locGoodFCALShowers[ishow];
		TVector3 pos1 = GetFCALPosition(show1);
		
		double  e1 = m_fcalE[show1];
		
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
			
			// standard TOF veto:
			
			bool locTOFVeto = false;
			if((tof_dr1 < m_FCALTOFCut) && (tof_dr2 < m_FCALTOFCut)) locTOFVeto = true;
			
			// only require one hit to be matched with TOF to veto:
			
			bool locTOFVeto_single = false;
			if((tof_dr1 < m_FCALTOFCut) || (tof_dr2 < m_FCALTOFCut)) locTOFVeto_single = true;
			
			// vary the timing cut used for TOF veto:
			
			vector<bool> locTOFVeto_dT;
			for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
				locTOFVeto_dT.push_back(false);
			}
			double locTOFdx1, locTOFdy1, locTOFdt1;
			double locTOFdx2, locTOFdy2, locTOFdt2;
			for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
				CheckTOFMatch(pos1, locTOFdx1, locTOFdy1, locTOFdt1, m_TOFTimingCuts[icut]);
				double locTOFdr1 = sqrt(pow(locTOFdx1,2.0)+pow(locTOFdy1,2.0));
				CheckTOFMatch(pos2, locTOFdx2, locTOFdy2, locTOFdt2, m_TOFTimingCuts[icut]);
				double locTOFdr2 = sqrt(pow(locTOFdx2,2.0)+pow(locTOFdy2,2.0));
				if((locTOFdr1 < m_FCALTOFCut) && (locTOFdr2 < m_FCALTOFCut)) {
					locTOFVeto_dT[icut] = true;
				}
			}
			
			// vary the distance cut used for TOF veto:
			
			vector<bool> locTOFVeto_dR;
			for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
				locTOFVeto_dR.push_back(false);
			}
			for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
				if((tof_dr1 < m_TOFDistanceCuts[icut]) && (tof_dr2 < m_TOFDistanceCuts[icut])) {
					locTOFVeto_dR[icut] = true;
				}
			}
			
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
				
				h_mgg_noTOF->Fill(prodTheta, invmassConstr, fillWeight);
				if(!locTOFVeto) h_mgg_TOF->Fill(prodTheta, invmassConstr, fillWeight);
				if(!locTOFVeto_single) h_mgg_singleTOF->Fill(prodTheta, invmassConstr, fillWeight);
				if(m_FillThrown && isEta) {
					h_AngularMatrix_noTOF->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
					if(!locTOFVeto) h_AngularMatrix_TOF->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
					if(!locTOFVeto_single) h_AngularMatrix_singleTOF->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
				}
				for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
					if(!locTOFVeto_dT[icut]) {
						h_mgg_TOFTimingCutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
						if(m_FillThrown && isEta) {
							h_AngularMatrix_TOFTimingCutVec[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
						}
					}
				}
				for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
					if(!locTOFVeto_dR[icut]) {
						h_mgg_TOFDistanceCutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
						if(m_FillThrown && isEta) {
							h_AngularMatrix_TOFDistanceCutVec[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
						}
					}
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

void EtaAna::InitializeAngularMatrices_TOF() {
	
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	h_AngularMatrix_noTOF = new TH3F("AngularMatrix_noTOF", 
		"No TOF Veto; #theta(thrown) [#circ]; #theta(rec) [#circ]; E_{#gamma}(thrown) [GeV]",
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
		nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
	h_AngularMatrix_noTOF->SetDirectory(0);
	
	h_AngularMatrix_TOF = new TH3F("AngularMatrix_TOF", 
		"Standard TOF Veto; #theta(thrown) [#circ]; #theta(rec) [#circ]; E_{#gamma}(thrown) [GeV]",
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
		nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
	h_AngularMatrix_TOF->SetDirectory(0);
	
	h_AngularMatrix_singleTOF = new TH3F("AngularMatrix_singleTOF", 
		"Veto on single FCAL-TOF Match; #theta(thrown) [#circ]; #theta(rec) [#circ]; E_{#gamma}(thrown) [GeV]",
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
		nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
	h_AngularMatrix_singleTOF->SetDirectory(0);
	
	for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
		TH3F *hMatrix = new TH3F(Form("AngularMatrix_TOFTimingCut_%02d",icut), 
			Form("Timing Cut for TOF Veto: #left|t_{TOF} - t_{RF}#right| < %.2f ns", m_TOFTimingCuts[icut]),
			nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
			nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
			nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
		hMatrix->SetDirectory(0);
		hMatrix->GetXaxis()->SetTitle("#theta(thrown) [#circ]");
		hMatrix->GetYaxis()->SetTitle("#theta(rec) [#circ]");
		hMatrix->GetZaxis()->SetTitle("E_{#gamma}(thrown) [#circ]");
		h_AngularMatrix_TOFTimingCutVec.push_back(hMatrix);
	}
	
	for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
		TH3F *hMatrix = new TH3F(Form("AngularMatrix_TOFDistanceCut_%02d",icut), 
			Form("Distance Cut for TOF Veto: #DeltaR_{FCAL-TOF} < %.1f cm", m_TOFDistanceCuts[icut]),
			nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
			nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
			nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
		hMatrix->SetDirectory(0);
		hMatrix->GetXaxis()->SetTitle("#theta(thrown) [#circ]");
		hMatrix->GetYaxis()->SetTitle("#theta(rec) [#circ]");
		hMatrix->GetZaxis()->SetTitle("E_{#gamma}(thrown) [#circ]");
		h_AngularMatrix_TOFDistanceCutVec.push_back(hMatrix);
	}
	return;
}
