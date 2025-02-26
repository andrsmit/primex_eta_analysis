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
	
	//-------------------------------------------//
	// Get list of selected beam photons:
	
	vector<pair<int,double>> locGoodBeamPhotons; locGoodBeamPhotons.clear();
	int locNBeamPhotons = GetBeamPhotonList(locGoodBeamPhotons, m_minBeamEnergyCut, m_maxBeamEnergyCut);
	
	//-------------------------------------------//
	// Get list of 'good' FCAL showers:
	
	vector<int> locGoodFCALShowers; locGoodFCALShowers.clear();
	int locNFCALShowers_EnergyCut;
	int locNFCALShowers = GetFCALShowerList(locGoodFCALShowers, locNFCALShowers_EnergyCut, m_FCALEnergyCut, m_FCALExtraEnergyCut, 2.0, m_FCALRFCut);
	
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
			// BCAL+SC Veto:
			
			bool isHadronicVeto = true;
			
			int locNSCHits = 0;
			int locNSCHits_coplanar = 0;
			for(int isc=0; isc<locGoodSCHits.size(); isc++) {
				int hitIndex = locGoodSCHits[isc];
				double locDeltaPhi = prodPhi - m_scPhi[hitIndex];
				if(IsCoplanarSC(locDeltaPhi)) locNSCHits_coplanar++;
			}
			
			if(locNBCALShowers==0) 
			{
				isHadronicVeto = false;
			}
			else if(locNBCALShowers==1) {
				if(IsCoplanarBCAL(prodPhi-locBCALPhi) && (locBCALRFDT>1.0)) {
					if(locNSCHits_coplanar==locNSCHits) {
						if(locNSCHits<=1) {
							isHadronicVeto = false;
						}
					}
				}
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

void EtaAna::InitializeTOFHists()
{
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	// VARY TOF CUTS:
	
	h_mgg_noTOF = new TH2F("mgg_noTOF", "No TOF Veto", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_noTOF->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_noTOF->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_noTOF->Sumw2();
	
	// with standard veto:
	
	h_mgg_TOF = new TH2F("mgg_TOF", "Standard TOF Veto", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_TOF->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_TOF->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_TOF->Sumw2();
	
	// with single-shower TOF veto:
	
	h_mgg_singleTOF = new TH2F("mgg_singleTOF", "Veto on single-shower TOF match", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_singleTOF->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_singleTOF->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
	h_mgg_singleTOF->Sumw2();
	
	// vary the timing cut used for TOF veto:
	
	m_TOFTimingCuts.clear();
	for(int icut=0; icut<7; icut++) {
		double locCut = 0.5 + 0.25*(double)(icut);
		m_TOFTimingCuts.push_back(locCut);
	}
	m_TOFTimingCuts.push_back(6.0);
	
	for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
		TH2F *loc_h_mgg = new TH2F(Form("mgg_TOFTiming_%02d", icut),
			Form("TOF Timing Cut: #left|t_{TOF} - t_{RF}#right| < %.2f ns", m_TOFTimingCuts[icut]), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		loc_h_mgg->Sumw2();
		h_mgg_TOFTimingCutVec.push_back(loc_h_mgg);
	}
	
	// vary the distance cut used for TOF veto:
	
	m_TOFDistanceCuts.clear();
	for(int icut=0; icut<15; icut++) {
		double locCut = 5.0 + 0.5*(double)(icut);
		m_TOFDistanceCuts.push_back(locCut);
	}
	m_TOFDistanceCuts.push_back(500.0);
	
	for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
		TH2F *loc_h_mgg = new TH2F(Form("mgg_TOFDistance_%02d", icut),
			Form("TOF Distance Cut: #DeltaR_{FCAL-TOF} < %.1f cm", m_TOFDistanceCuts[icut]), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		loc_h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		loc_h_mgg->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		loc_h_mgg->Sumw2();
		h_mgg_TOFDistanceCutVec.push_back(loc_h_mgg);
	}
	
	if(m_FillThrown) 
	{
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
	}
	return;
}

void EtaAna::ResetTOFHists()
{
	h_mgg_noTOF->Reset();
	if(h_AngularMatrix_noTOF!=nullptr) h_AngularMatrix_noTOF->Reset();
	
	h_mgg_TOF->Reset();
	if(h_AngularMatrix_TOF!=nullptr) h_AngularMatrix_TOF->Reset();
	
	h_mgg_singleTOF->Reset();
	if(h_AngularMatrix_singleTOF!=nullptr) h_AngularMatrix_singleTOF->Reset();
	
	for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
		h_mgg_TOFTimingCutVec[icut]->Reset();
	}
	for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
		h_mgg_TOFDistanceCutVec[icut]->Reset();
	}
	if(h_AngularMatrix_TOFTimingCutVec.size()) {
		for(int i=0; i<h_AngularMatrix_TOFTimingCutVec.size(); i++) {
			h_AngularMatrix_TOFTimingCutVec[i]->Reset();
		}
	}
	if(h_AngularMatrix_TOFDistanceCutVec.size()) {
		for(int i=0; i<h_AngularMatrix_TOFDistanceCutVec.size(); i++) {
			h_AngularMatrix_TOFDistanceCutVec[i]->Reset();
		}
	}
	return;
}

void EtaAna::WriteTOFHists()
{
	printf("\n  Writing TOF histograms...\n");
	
	h_mgg_noTOF->Write();
	if(h_AngularMatrix_noTOF!=nullptr) h_AngularMatrix_noTOF->Write();
	
	h_mgg_TOF->Write();
	if(h_AngularMatrix_TOF!=nullptr) h_AngularMatrix_TOF->Write();
	
	h_mgg_singleTOF->Write();
	if(h_AngularMatrix_singleTOF!=nullptr) h_AngularMatrix_singleTOF->Write();
	
	TDirectory *dirTOFTimingCut = new TDirectoryFile("TOFTimingCut", "TOFTimingCut");
	dirTOFTimingCut->cd();
	for(int icut=0; icut<m_TOFTimingCuts.size(); icut++) {
		h_mgg_TOFTimingCutVec[icut]->Write();
	}
	if(h_AngularMatrix_TOFTimingCutVec.size()) {
		for(int i=0; i<h_AngularMatrix_TOFTimingCutVec.size(); i++) {
			h_AngularMatrix_TOFTimingCutVec[i]->Write();
		}
	}
	dirTOFTimingCut->cd("../");
	
	TDirectory *dirTOFDistanceCut = new TDirectoryFile("TOFDistanceCut", "TOFDistanceCut");
	dirTOFDistanceCut->cd();
	for(int icut=0; icut<m_TOFDistanceCuts.size(); icut++) {
		h_mgg_TOFDistanceCutVec[icut]->Write();
	}
	if(h_AngularMatrix_TOFDistanceCutVec.size()) {
		for(int i=0; i<h_AngularMatrix_TOFDistanceCutVec.size(); i++) {
			h_AngularMatrix_TOFDistanceCutVec[i]->Write();
		}
	}
	dirTOFDistanceCut->cd("../");
	
	printf("  Done.\n");
	return;
}
