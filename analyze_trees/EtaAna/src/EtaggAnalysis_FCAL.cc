#include "EtaAna.h"

void EtaAna::EtaggAnalysis_FCAL() {
	
	if(m_nmc>0) {
		if(AcceptRejectEvent()) return;
	}
	
	double locThrownBeamEnergy = 0.0, locThrownAngle = 0.0;
	if(m_FillThrown) {
		GetThrownEnergyAndAngle(locThrownBeamEnergy, locThrownAngle);
		if((locThrownBeamEnergy < m_minBeamEnergyCut) || (locThrownBeamEnergy >= m_maxBeamEnergyCut)) return;
		PlotThrown(locThrownBeamEnergy, locThrownAngle);
	}
	
	//-------------------------------------------//
	// Get list of selected beam photons:
	
	vector<pair<int,double>> locGoodBeamPhotons; locGoodBeamPhotons.clear();
	int locNBeamPhotons = GetBeamPhotonList(locGoodBeamPhotons, m_minBeamEnergyCut, m_maxBeamEnergyCut);
	
	//-------------------------------------------//
	// Get list of 'good' FCAL showers:
	
	vector<int> locGoodFCALShowers; locGoodFCALShowers.clear();
	int locNFCALShowersEnergyCut;
	int locNFCALShowersTotal = GetFCALShowerList(locGoodFCALShowers, locNFCALShowersEnergyCut, 
		0.0, m_FCALExtraEnergyCut, 0.0, m_FCALRFCut);
	
	//-------------------------------------------//
	// Get list of 'good' BCAL showers:
	
	vector<int> locGoodBCALShowers; locGoodBCALShowers.clear();
	int locNBCALShowers = GetBCALShowerList(locGoodBCALShowers, 0.0, m_BCALRFCut);
	
	//-------------------------------------------//
	// Get list of 'good' SC hits:
	
	vector<int> locGoodSCHits; locGoodSCHits.clear();
	int locNSCHits = GetSCHitList(locGoodSCHits);
	
	//-------------------------------------------//
	
	int locNBCALShowers_1ns = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0., locBCALTheta = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	
	for(int ishow=0; ishow<locGoodBCALShowers.size(); ishow++) {
		int showIndex = locGoodBCALShowers[ishow];
		TVector3 locPos = GetBCALPosition(showIndex);
		double locT = m_bcalT[showIndex] - (locPos.Mag()/m_c) - m_rfTime;
		
		locBCALRFDT  = locT;
		locBCALPhi   = locPos.Phi()   * TMath::RadToDeg();
		locBCALTheta = locPos.Theta() * TMath::RadToDeg();
		if(fabs(locT) < 1.0) {
			locNBCALShowers_1ns++;
		}
	}
	
	//-------------------------------------------//
	/*
	The vector 'locGoodFCALShowers' stores the indices for all FCAL showers within +/-2ns of the RF time.
	In the following code, we loop over different values of minimum energy cuts and determine whether or not 
	there are exactly 2 showers (no more) with that cut value:
	*/
	vector<int> locEnergyCuts; locEnergyCuts.clear();
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		
		int locGoodFCALShowers_cut = 0;
		for(int ishow=0; ishow<m_nfcal; ishow++) {
			
			TVector3 locPos = GetFCALPosition(ishow);
			double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
			if((fabs(locT) < m_FCALRFCut) && (m_fcalE[ishow] > m_fcalEnergyCuts[icut])) {
				locGoodFCALShowers_cut++;
			}
		}
		if(locGoodFCALShowers_cut==2) locEnergyCuts.push_back(1);
		else locEnergyCuts.push_back(0);
	}
	
	//=====================================================================================//
	
	int locNFCALShowersGood = (int)locGoodFCALShowers.size();
	if(locNFCALShowersGood<2) return;
	
	for(int ishow=0; ishow<(locNFCALShowersGood-1); ishow++) {
		
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
		
		for(int jshow=(ishow+1); jshow<locNFCALShowersGood; jshow++) {
			
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
			//if((tof_dr1 < m_FCALTOFCut) || (tof_dr2 < m_FCALTOFCut)) isTOFVeto = true;
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
			// BCAL+SC Veto:
			
			int locNSCHits_coplanar = 0;
			for(int isc=0; isc<locGoodSCHits.size(); isc++) {
				int hitIndex = locGoodSCHits[isc];
				double locDeltaPhi = prodPhi - m_scPhi[hitIndex];
				if(IsCoplanarSC(locDeltaPhi)) locNSCHits_coplanar++;
			}
			
			if(
				IsHadronicVeto(m_vetoOption, locNBCALShowers, locNBCALShowers_1ns, locNSCHits, locNSCHits_coplanar, 
					prodPhi, locBCALPhi, locBCALTheta, locBCALRFDT)
			) continue;
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(int igam=0; igam<locNBeamPhotons; igam++) {
				
				int ibeam = locGoodBeamPhotons[igam].first;
				double fillWeight = locGoodBeamPhotons[igam].second;
				
				double eb    = m_beamE[ibeam];
				double brfdt = m_beamT[ibeam] - m_rfTime;
				
				// Calculate the energy of the eta meson, assuming production on a free nucleon:
				double etaEnergy = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(Proton));
				
				// Calculate the energy of the eta meson, assuming production on a helium-4 nucleus::
				double etaEnergyCoh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
				
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
				double constraintEnergy = m_IsCohMC ? etaEnergyCoh : etaEnergy;
				
				double e1c = e1/(1.+sigr) + (constraintEnergy-e2)/(1.+(1./sigr));
				double e2c = constraintEnergy - e1c;
				double invmassConstr = sqrt(2.*e1c*e2c*(1.-cos12));
				
				//-----------------------------------------------------//
				
				h_mgg_FCAL->Fill(prodTheta, invmassConstr, fillWeight);
				if((e1 > m_FCALEnergyCut) && (e2 > m_FCALEnergyCut)) {
					h_mgg_FCALECut->Fill(prodTheta, invmassConstr, fillWeight);
				}
				if(!FCALFiducialCut(pos1, 2.0) && !FCALFiducialCut(pos2, 2.0)) {
					h_mgg_FCALFidCut->Fill(prodTheta, invmassConstr, fillWeight);
					if((e1 > m_FCALEnergyCut) && (e2 > m_FCALEnergyCut)) {
						h_mgg_FCALCuts->Fill(prodTheta, invmassConstr, fillWeight);
						if(locNFCALShowersEnergyCut==2) {
							h_mgg_FCALGoodMult->Fill(prodTheta, invmassConstr, fillWeight);
							if(locNFCALShowersTotal==2) {
								h_mgg_FCALMult->Fill(prodTheta, invmassConstr, fillWeight);
							}
						}
					}
				}
				
				if(locNFCALShowersEnergyCut==2) {
					// Vary the size of the fiducial cut used:
					if((e1 > m_FCALEnergyCut) && (e2 > m_FCALEnergyCut)) {
						for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
							if(!FCALFiducialCut(pos1, m_fcalFiducialCuts[icut]) && !FCALFiducialCut(pos2, m_fcalFiducialCuts[icut])) {
								h_mgg_FCALFidCutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
								if(m_FillThrown && IsEtaCut(invmassConstr)) {
									h_AngularMatrix_FCALFidCutVec[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
								}
							}
						}
					}
					
					// Vary the minimum energy cut used:
					if(!FCALFiducialCut(pos1, 2.0) && !FCALFiducialCut(pos2, 2.0)) {
						for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
							if(e1>m_fcalEnergyCuts[icut] && e2>m_fcalEnergyCuts[icut]) {
								h_mgg_FCALECutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
								if(m_FillThrown && IsEtaCut(invmassConstr)) {
									h_AngularMatrix_FCALECutVec[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
								}
							}
						}
					}
				}
				
				// Vary the minimum energy cut used for looking for extra showers:
				if(!FCALFiducialCut(pos1,2.0) && !FCALFiducialCut(pos2,2.0) && (e1>m_FCALEnergyCut) && (e2>m_FCALEnergyCut)) {
					for(int icut=0; icut<locEnergyCuts.size(); icut++) {
						if(locEnergyCuts[icut]==1) {
							h_mgg_FCALExtraECutVec[icut]->Fill(prodTheta, invmassConstr, fillWeight);
							if(m_FillThrown && IsEtaCut(invmassConstr)) {
								h_AngularMatrix_FCALExtraECutVec[icut]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
							}
						}
					}
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

void EtaAna::InitializeFCALHists()
{
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	// VARY FCAL CUTS:
	
	h_mgg_FCAL = new TH2F("mgg_FCAL", "No Multiplicity Cut", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_FCAL->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_FCAL->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_FCAL->Sumw2();
	
	// with minimum energy cuts:
	
	h_mgg_FCALECut = new TH2F("mgg_FCAL_ecut", Form("E_{1,2} > %.2f GeV", m_FCALEnergyCut), 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_FCALECut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_FCALECut->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_FCALECut->Sumw2();
	
	// with fiducial cuts:
	
	h_mgg_FCALFidCut = new TH2F("mgg_FCAL_fidcut", "Fiducial Cut Applied", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_FCALFidCut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_FCALFidCut->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_FCALFidCut->Sumw2();
	
	// with both energy and fiducial cuts:
	
	h_mgg_FCALCuts = new TH2F("mgg_FCAL_cuts", "Fiducial+Energy Cuts Applied", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_FCALCuts->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_FCALCuts->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_FCALCuts->Sumw2();
	
	// with 'good' multiplicity = 2:
	
	h_mgg_FCALGoodMult = new TH2F("mgg_FCAL_good_mult", "2 Good FCAL Showers", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_FCALGoodMult->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_FCALGoodMult->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_FCALGoodMult->Sumw2();
	
	// with total multiplicity = 2:
	
	h_mgg_FCALMult = new TH2F("mgg_FCAL_mult", "2 FCAL Showers", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_FCALMult->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_FCALMult->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_FCALMult->Sumw2();
	
	// vary the size of the fiducial cut:
	
	m_fcalFiducialCuts.clear();
	for(int icut=0; icut<17; icut++) {
		double locCut = 0.0 + 0.25*(double)(icut);
		m_fcalFiducialCuts.push_back(locCut);
	}
	
	for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
		TH2F *loc_h_mgg = new TH2F(Form("mgg_FCAL_fid_%02d", icut),
			Form("Inner layers removed by fiducial cut: %.1f", m_fcalFiducialCuts[icut]), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		loc_h_mgg->Sumw2();
		h_mgg_FCALFidCutVec.push_back(loc_h_mgg);
	}
	
	// vary the minimum energy cut:
	
	m_fcalEnergyCuts.clear();
	for(int icut=0; icut<20; icut++) {
		double locCut = 0.0 + 0.05*(double)(icut);
		m_fcalEnergyCuts.push_back(locCut);
	}
	
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		TH2F *loc_h_mgg = new TH2F(Form("mgg_FCAL_ecut_%02d", icut),
			Form("Minimum Energy Cut: %.2f GeV", m_fcalEnergyCuts[icut]), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		loc_h_mgg->Sumw2();
		h_mgg_FCALECutVec.push_back(loc_h_mgg);
	}
	
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		TH2F *loc_h_mgg = new TH2F(Form("mgg_FCAL_extra_ecut_%02d", icut),
			Form("Minimum Energy Cut: %.2f GeV", m_fcalEnergyCuts[icut]), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		loc_h_mgg->Sumw2();
		h_mgg_FCALExtraECutVec.push_back(loc_h_mgg);
	}
	
	if(m_FillThrown) 
	{
		for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_FCALFidCut_%02d",icut), 
				Form("Inner layers removed by fiducial cut: %.1f; #theta(thrown) [#circ]; #theta(rec) [#circ]", m_fcalFiducialCuts[icut]),
				nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
				nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
				nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
			hMatrix->SetDirectory(0);
			h_AngularMatrix_FCALFidCutVec.push_back(hMatrix);
		}
		
		for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_FCALECut_%02d",icut), 
				Form("Minimum Energy Cut: %.2f GeV; #theta(thrown) [#circ]; #theta(rec) [#circ]", m_fcalEnergyCuts[icut]),
				nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
				nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
				nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
			hMatrix->SetDirectory(0);
			h_AngularMatrix_FCALECutVec.push_back(hMatrix);
		}
		
		for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_FCALExtraECut_%02d",icut), 
				Form("Extra Energy Cut: %.2f GeV; #theta(thrown) [#circ]; #theta(rec) [#circ]", m_fcalEnergyCuts[icut]),
				nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
				nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
				nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
			hMatrix->SetDirectory(0);
			h_AngularMatrix_FCALExtraECutVec.push_back(hMatrix);
		}
	}
	
	return;
}

void EtaAna::ResetFCALHists()
{
	h_mgg_FCAL->Reset();
	h_mgg_FCALECut->Reset();
	h_mgg_FCALFidCut->Reset();
	h_mgg_FCALCuts->Reset();
	h_mgg_FCALGoodMult->Reset();
	h_mgg_FCALMult->Reset();
	for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) h_mgg_FCALFidCutVec[icut]->Reset();
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		h_mgg_FCALECutVec[icut]->Reset();
		h_mgg_FCALExtraECutVec[icut]->Reset();
	}
	
	if(h_AngularMatrix_FCALECutVec.size()) {
		for(int i=0; i<h_AngularMatrix_FCALECutVec.size(); i++) {
			h_AngularMatrix_FCALECutVec[i]->Reset();
		}
	}
	if(h_AngularMatrix_FCALExtraECutVec.size()) {
		for(int i=0; i<h_AngularMatrix_FCALExtraECutVec.size(); i++) {
			h_AngularMatrix_FCALExtraECutVec[i]->Reset();
		}
	}
	if(h_AngularMatrix_FCALFidCutVec.size()) {
		for(int i=0; i<h_AngularMatrix_FCALFidCutVec.size(); i++) {
			h_AngularMatrix_FCALFidCutVec[i]->Reset();
		}
	}
	return;
}

void EtaAna::WriteFCALHists()
{
	printf("\n  Writing FCAL histograms...\n");
	
	h_mgg_FCAL->Write();
	h_mgg_FCALECut->Write();
	h_mgg_FCALFidCut->Write();
	h_mgg_FCALCuts->Write();
	h_mgg_FCALGoodMult->Write();
	h_mgg_FCALMult->Write();
	
	TDirectory *dirECut = new TDirectoryFile("ECut", "ECut");
	dirECut->cd();
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		h_mgg_FCALECutVec[icut]->Write();
	}
	if(h_AngularMatrix_FCALECutVec.size()) {
		for(int i=0; i<h_AngularMatrix_FCALECutVec.size(); i++) {
			h_AngularMatrix_FCALECutVec[i]->Write();
		}
	}
	dirECut->cd("../");
	
	TDirectory *dirExtraECut = new TDirectoryFile("ExtraECut", "ExtraECut");
	dirExtraECut->cd();
	for(int icut=0; icut<m_fcalEnergyCuts.size(); icut++) {
		h_mgg_FCALExtraECutVec[icut]->Write();
	}
	if(h_AngularMatrix_FCALExtraECutVec.size()) {
		for(int i=0; i<h_AngularMatrix_FCALExtraECutVec.size(); i++) {
			h_AngularMatrix_FCALExtraECutVec[i]->Write();
		}
	}
	dirExtraECut->cd("../");
	
	TDirectory *dirFidCut = new TDirectoryFile("FidCut", "FidCut");
	dirFidCut->cd();
	for(int icut=0; icut<m_fcalFiducialCuts.size(); icut++) {
		h_mgg_FCALFidCutVec[icut]->Write();
	}
	if(h_AngularMatrix_FCALFidCutVec.size()) {
		for(int i=0; i<h_AngularMatrix_FCALFidCutVec.size(); i++) {
			h_AngularMatrix_FCALFidCutVec[i]->Write();
		}
	}
	dirFidCut->cd("../");
	
	printf("  Done.\n");
	return;
}

