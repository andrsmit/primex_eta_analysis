#include "EtaAna.h"

void EtaAna::EtaggAnalysis() {
	
	if(m_nmc>0) {
		h_mcVertex->Fill(m_mcZ[0]);
		if(AcceptRejectEvent()) return;
		h_mcVertexAccepted->Fill(m_mcZ[0]);
	}
	
	double locThrownBeamEnergy = 0.0, locThrownAngle = 0.0;
	if(m_FillThrown) {
		GetThrownEnergyAndAngle(locThrownBeamEnergy, locThrownAngle);
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
	
	//-------------------------------------------//
	// RF Timing distributions:
	
	// FCAL:
	/*
	for(int ishow = 0; ishow < m_nfcal; ishow++) {
		TVector3 locPos = GetFCALPosition(ishow);
		double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		h_fcalRFdt->Fill(locT);
	}
	*/
	// BCAL:
	/*
	for(int ishow = 0; ishow < m_nbcal; ishow++) {
		TVector3 locPos = GetBCALPosition(ishow);
		double locT = m_bcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		h_bcalRFdt->Fill(locT);
	}
	*/
	
	// SC:
	for(int isc = 0; isc < m_nsc; isc++) {
		double locT  = m_scT[isc] - m_rfTime;
		double locdE = m_scdE[isc];
		if(locdE > 0.0002) {
			h_scRFdt->Fill(locT);
		}
	}
	
	// TOF:
	for(int itof = 0; itof < m_ntof; itof++) {
		TVector3 locPos(m_tofX[itof], m_tofY[itof], m_tofZ[itof]);
		locPos -= m_vertex;
		double locT = m_tofT[itof] - (locPos.Mag()/m_c) - m_rfTime;
		h_tofRFdt->Fill(locT);
	}
	
	// Beam:
	for(int igam = 0; igam < locNBeamPhotons; igam++) {
		double locT = m_beamT[locGoodBeamPhotons[igam].first] - m_rfTime;
		h_beamRFdt->Fill(locT);
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
			
			int tof_match1 = tof_dr1 < m_FCALTOFCut ? 1.0 : 0.0;
			int tof_match2 = tof_dr2 < m_FCALTOFCut ? 1.0 : 0.0;
			
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
			
			int locNSCHits_coplanar = 0;
			for(int isc=0; isc<locGoodSCHits.size(); isc++) {
				int hitIndex = locGoodSCHits[isc];
				double locDeltaPhi = prodPhi - m_scPhi[hitIndex];
				if(IsCoplanarSC(locDeltaPhi)) locNSCHits_coplanar++;
			}
			
			vector<int> locVetoOptions; locVetoOptions.clear();
			for(int iveto=0; iveto<m_nVetos; iveto++) locVetoOptions.push_back(0);
			
			// Option 0 (no veto): No Veto is applied:
			locVetoOptions[0] = 1;
			
			// Option 1 (strict): Remove events with any BCAL shower within +/-12ns:
			if(locNBCALShowers==0) locVetoOptions[1] = 1;
			
			// Option 2 (loose): Remove events with any BCAL shower within +/-1ns:
			if(locNBCALShowers_1ns==0) locVetoOptions[2] = 1;
			
			
			// Options 3-6:
			//   If there are no showers within +/-12ns in the BCAL, accept event
			//   If there is exactly 1 shower, then:
			//      Require coplanarity cut (option 3)
			//      Require coplanarity cut + timing > 1ns (option 4)
			//      Require coplanarity cut + timing > 1ns + all SC hits are coplanar (option 5)
			//      Require coplanarity cut + timing > 1ns + all SC hits are coplnar + less than 2 SC hits (option 6)
			
			if(locNBCALShowers==0) 
			{
				locVetoOptions[3] = 1;
				locVetoOptions[4] = 1;
				locVetoOptions[5] = 1;
				locVetoOptions[6] = 1;
			}
			else if(locNBCALShowers==1) 
			{
				if(IsCoplanarBCAL(prodPhi-locBCALPhi)) {
					locVetoOptions[3] = 1;
					if(locBCALRFDT>1.0) {
						locVetoOptions[4] = 1;
						if(locNSCHits_coplanar==locNSCHits) {
							locVetoOptions[5] = 1;
							if(locNSCHits<=1) {
								locVetoOptions[6] = 1;
							}
						}
					}
				}
			}
			
			// Option 7: strict SC + strict BCAL vetos:
			if(locNSCHits==0 && locNBCALShowers==0) locVetoOptions[7] = 1;
			
			// Option 8: BCAL and SC are completely decoupled, but both need to be satisfied:
			
			bool isBCALVeto = true;
			if(locNBCALShowers==0)
			{
				isBCALVeto = false;
			}
			else if(locNBCALShowers==1)
			{
				if(IsCoplanarBCAL(prodPhi-locBCALPhi)) {
					if(locBCALRFDT>1.0) {
						isBCALVeto = false;
					}
				}
			}
			
			bool isSCVeto = true;
			if(locNSCHits<=1)
			{
				if(locNSCHits==locNSCHits_coplanar) {
					isSCVeto = false;
				}
			}
			
			if(!isSCVeto && !isBCALVeto) {
				locVetoOptions[8] = 1;
			}
			
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
				
				// Plot timing distribution of beam photons after elasticity cut to see the level of accidentals:
				if(isElastic && locVetoOptions[6]) {
					h_beamRFdt_cut->Fill(brfdt);
				}
				
				TLorentzVector k_beam(0.0, 0.0, eb, eb);
				TLorentzVector p_eta(pggx, pggy, pggz, Egg);
				double t    = (k_beam-p_eta)*(k_beam-p_eta);
				double logt = log10(-t);
				
				// acceptance correction:
				/*
				double locAcc = 1.0;
				if(h_acceptance != nullptr) {
					locAcc = h_acceptance->GetBinContent(h_acceptance->GetXaxis()->FindBin(eb), 
						h_acceptance->GetYaxis()->FindBin(prodTheta));
				}
				*/
				
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
				double e2c_coh = etaEnergy_coh - e1c_coh;
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
				
				//-----------------------------------------------------//
				// Mass constraint
				
				// adjust the measured energies of the two-photons so that the invariant mass exactly equals the 
				// PDG mass of an eta:
				
				double ecorr = -1.0*e1 - sigr*e2 
					+ sqrt(pow(e1 + sigr*e2,2.0) - 4.0*sigr*(e1*e2 - pow(ParticleMass(Eta),2.0)/(2.0*(1.0-cos12))));
				double e1mc = e1 + 0.5*ecorr;
				double e2mc = e2 + 0.5*ecorr/sigr;
				
				double elasConstr = (e1mc+e2mc)/etaEnergy_coh;
				
				//-----------------------------------------------------//
				// Compare measured transverse momentum of eta to calculated value:
				
				double pTCalc    = sqrt(pow(etaEnergy,2.0)     - pow(ParticleMass(Eta),2.0))*sin(prodTheta*TMath::DegToRad());
				double pTCalcCoh = sqrt(pow(etaEnergy_coh,2.0) - pow(ParticleMass(Eta),2.0))*sin(prodTheta*TMath::DegToRad());
				
				//-----------------------------------------------------//
				// Missing Mass
				
				double mmSq = pow(ParticleMass(Proton),2.0) + pow(ParticleMass(Eta),2.0) 
					+ 2.0*ParticleMass(Proton)*eb 
					- 2.0*ParticleMass(Proton)*Egg 
					- 2.0*eb*(Egg - sqrt(pow(Egg,2.0)-pow(ParticleMass(Eta),2.0))*cos(prodTheta*TMath::DegToRad()));
				
				double mmSqCoh = pow(ParticleMass(m_Target),2.0) + pow(ParticleMass(Eta),2.0) 
					+ 2.0*ParticleMass(m_Target)*eb 
					- 2.0*ParticleMass(m_Target)*Egg 
					- 2.0*eb*(Egg - sqrt(pow(Egg,2.0)-pow(ParticleMass(Eta),2.0))*cos(prodTheta*TMath::DegToRad()));
				
				//-----------------------------------------------------//
				
				/*
				TLorentzVector locMesonP4(pggx, pggy, pggz, Egg);
				
				TVector3 locBoostVector_MesonCM = -1.0*(locMesonP4.BoostVector());
				TLorentzVector locProduct1_MesonCM(px1, py1, pz1, e1);
				TLorentzVector locProduct2_MesonCM(px2, py2, pz2, e2);
				locProduct1_MesonCM.Boost(locBoostVector_MesonCM);
				locProduct2_MesonCM.Boost(locBoostVector_MesonCM);
				
				h_deltaPhi_CM->Fill((locProduct1_MesonCM.Phi()-locProduct2_MesonCM.Phi())*TMath::RadToDeg());
				
				// plot energy vs angle in CM frame:
				
				if(((e1/Egg)<0.3) || ((e2/Egg)<0.3)) continue;
				
				if(isEta && isElastic && locVetoOptions[6]) {
					h_e_vs_theta->Fill(acos((px1*pggx + py1*pggy + pz1*pggz)/(e1*sqrt(pow(pggt,2.0)+pow(pggz,2.0)))), e1/Egg);
					h_e_vs_theta->Fill(acos((px2*pggx + py2*pggy + pz2*pggz)/(e2*sqrt(pow(pggt,2.0)+pow(pggz,2.0)))), e2/Egg);
				}
				*/
				
				if(isEta && isElastic) {
					
					if(locVetoOptions[6]) {
						h_xy1->Fill(pos1.X(),pos1.Y());
						h_xy2->Fill(pos2.X(),pos2.Y());
						h_t_vs_theta->Fill(prodTheta, logt, fillWeight);
						
						// plot likelihood of TOF match as a function of radius:
						
						double rho1 = sqrt(pow(pos1.X(),2.0) + pow(pos1.Y(),2.0));
						double rho2 = sqrt(pow(pos2.X(),2.0) + pow(pos2.Y(),2.0));
						h_tofMatch->Fill(rho1, tof_match1, fillWeight);
						h_tofMatch->Fill(rho2, tof_match2, fillWeight);
					}
					
					// plot SC-DeltaPhi:
					
					for(int isc=0; isc<locGoodSCHits.size(); isc++) {
						int hitIndex = locGoodSCHits[isc];
						double locDeltaPhi = prodPhi - m_scPhi[hitIndex];
						h_scDeltaPhi->Fill(prodTheta, locDeltaPhi, fillWeight);
						if(locGoodSCHits.size()==1) {
							h_scDeltaPhi_singleHit->Fill(prodTheta, locDeltaPhi, fillWeight);
						}
					}
					
					// plot BCAL-DeltaPhi:
					
					for(int ibcal = 0; ibcal < locGoodBCALShowers.size(); ibcal++) {
						double locDeltaPhi = prodPhi - GetBCALPosition(locGoodBCALShowers[ibcal]).Phi() * TMath::RadToDeg();
						h_bcalDeltaPhi->Fill(prodTheta, locDeltaPhi, fillWeight);
						if(locGoodBCALShowers.size()==1) {
							h_bcalDeltaPhi_singleHit->Fill(prodTheta, locDeltaPhi, fillWeight);
						}
					}
					
					// plot FCAL-RF dt for accepted events:
					
					for(int ishow = 0; ishow < m_nfcal; ishow++) {
						TVector3 locPos = GetFCALPosition(ishow);
						double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
						h_fcalRFdt->Fill(locT, fillWeight);
					}
					
					// plot BCAL-RF dt for accepted events:
					
					for(int ishow = 0; ishow < m_nbcal; ishow++) {
						TVector3 locPos = GetBCALPosition(ishow);
						double locT = m_bcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
						h_bcalRFdt->Fill(locT, fillWeight);
					}
					
					// plot number of SC hits:
					
					h_nSC->Fill(locNSCHits, fillWeight);
					if(locNBCALShowers==0) h_nSC_nobcal->Fill(locNSCHits, fillWeight);
					if(locNBCALShowers==1) h_nSC_onebcal->Fill(locNSCHits, fillWeight);
					if(locVetoOptions[4]) {
						h_nSC_bcalVeto->Fill(locNSCHits, fillWeight);
						h_nSC_extra->Fill(locNSCHits-locNSCHits_coplanar, fillWeight);
					}
				}
				
				for(int iveto=0; iveto<m_nVetos; iveto++) {
					
					if(locVetoOptions[iveto]==0) continue;
					if(isEta) {
						h_elasticity[iveto]->Fill(prodTheta, Egg/etaEnergy, fillWeight);
						h_elasticityConstr[iveto]->Fill(prodTheta, elasConstr, fillWeight);
						
						h_mm[iveto]->Fill(prodTheta, mmSq, fillWeight);
						h_mm_coh[iveto]->Fill(prodTheta, (mmSqCoh - pow(ParticleMass(Helium),2.0)), fillWeight);
						if(isElastic) {
							h_mm_elas[iveto]->Fill(prodTheta, mmSq, fillWeight);
							h_mm_elas_coh[iveto]->Fill(prodTheta, (mmSqCoh - pow(ParticleMass(Helium),2.0)), fillWeight);
						}
						
						h_pt[iveto]->Fill(prodTheta, pggt/pTCalc, fillWeight);
						h_ptCoh[iveto]->Fill(prodTheta, pggt/pTCalcCoh, fillWeight);
					}
					
					if(isElastic) {
						h_mgg[iveto]->Fill(prodTheta, invmass, fillWeight);
						h_mggConstr[iveto]->Fill(prodTheta, invmassConstr, fillWeight);
						h_mggConstr_coh[iveto]->Fill(prodTheta, invmassConstr_coh, fillWeight);
					}
					
					if(m_FillThrown && isEta && isElastic) {
						FillAngularMatrix_vetos(iveto, locThrownBeamEnergy, locThrownAngle, prodTheta, fillWeight);
					}
				}
				
				if((m_nmc>0) && locVetoOptions[6] && isElastic) {
					h_mgg_vs_vertexZ->Fill(m_mcZ[0], invmassConstr, fillWeight);
					h_mgg_vs_vertexR->Fill(sqrt(pow(m_mcX[0],2.0) + pow(m_mcY[0],2.0)), invmassConstr, fillWeight);
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}


void EtaAna::InitializeDefaultHists() 
{
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	// DEFAULT ANALYSIS HISTOGRAMS:
	
	h_mcVertex         = new TH1F("vertex",          "Vertex Z Position (thrown)",   1000, 0., 600.);
	h_mcVertexAccepted = new TH1F("vertex_accepted", "Vertex Z Position (filtered)", 1000, 0., 600.);
	
	// Timing histograms (for monitoring of calibrations):
	
	h_fcalRFdt     = new TH1F("fcal_rf_dt",     "t_{FCAL} - t_{RF}; [ns]", 10000, -100., 100.);
	h_bcalRFdt     = new TH1F("bcal_rf_dt",     "t_{BCAL} - t_{RF}; [ns]", 10000, -100., 100.);
	h_tofRFdt      = new TH1F( "tof_rf_dt",      "t_{TOF} - t_{RF}; [ns]", 10000, -100., 100.);
	h_scRFdt       = new TH1F(  "sc_rf_dt",       "t_{SC} - t_{RF}; [ns]", 10000, -100., 100.);
	h_beamRFdt     = new TH1F("beam_rf_dt",     "t_{CCAL} - t_{RF}; [ns]", 10000, -100., 100.);
	h_beamRFdt_cut = new TH1F("beam_rf_dt_cut", "t_{Beam} - t_{RF}; [ns]", 10000, -100., 100.);
	
	// To study veto likelihood versus polar angle:
	
	h_tofMatch = new TH2F("tofMatch", "TOF Match vs. Radial Distance; r_{FCAL} [cm]", 1500, 0.0, 150.0, 2, -0.5, 1.5);
	
	// Invariant mass, elasticity, and missing mass distributions for different BCAL/SC veto options:
	
	for(int iveto=0; iveto<m_nVetos; iveto++) {
		
		h_elasticity[iveto] = new TH2F(Form("elasticity_veto_%d",iveto),
			Form("Elasticity (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
		h_elasticity[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_elasticity[iveto]->GetYaxis()->SetTitle("E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
		
		h_elasticityConstr[iveto] = new TH2F(Form("elasticity_const_veto_%d",iveto),
			Form("Mass-Constrained Elasticity (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
		h_elasticityConstr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_elasticityConstr[iveto]->GetYaxis()->SetTitle(
			"E_{#gamma#gamma}^{Constr}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)");
		
		h_mgg[iveto] = new TH2F(Form("mgg_veto_%d",iveto), 
			Form("Two-Photon Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mgg[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV/c^{2}]");
		h_mgg[iveto]->Sumw2();
		
		h_mggConstr[iveto] = new TH2F(Form("mgg_const_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr[iveto]->Sumw2();
		
		h_mggConstr_coh[iveto] = new TH2F(Form("mgg_const_coh_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_coh[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr_coh[iveto]->Sumw2();
		
		h_mm[iveto] = new TH2F(Form("mm_veto_%d",iveto), 
			Form("Missing Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
		h_mm[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm[iveto]->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
		h_mm[iveto]->Sumw2();
		
		h_mm_coh[iveto] = new TH2F(Form("mm_coh_veto_%d",iveto), 
			Form("Missing Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
		h_mm_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_coh[iveto]->GetYaxis()->SetTitle("#Deltam^{2} - m_{He}^{2} [GeV^{2}/c^{4}]");
		h_mm_coh[iveto]->Sumw2();
		
		h_mm_elas[iveto] = new TH2F(Form("mm_elas_veto_%d",iveto), 
			Form("Missing Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
		h_mm_elas[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_elas[iveto]->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
		h_mm_elas[iveto]->Sumw2();
		
		h_mm_elas_coh[iveto] = new TH2F(Form("mm_elas_coh_veto_%d",iveto), 
			Form("Missing Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -20.0, 20.0);
		h_mm_elas_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mm_elas_coh[iveto]->GetYaxis()->SetTitle("#Deltam^{2} - m_{He}^{2} [GeV^{2}/c^{4}]");
		h_mm_elas_coh[iveto]->Sumw2();
		
		h_pt[iveto] = new TH2F(Form("pt_veto_%d",iveto), 
			Form("Transverse Momentum Ratio (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
		h_pt[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_pt[iveto]->GetYaxis()->SetTitle("p_{T} / p_{T}^{Calc} [GeV/c]");
		h_pt[iveto]->Sumw2();
		
		h_ptCoh[iveto] = new TH2F(Form("ptConstr_veto_%d",iveto), 
			Form("Transverse Momentum Ratio (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
		h_ptCoh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_ptCoh[iveto]->GetYaxis()->SetTitle("p_{T} / p_{T}^{Calc} [GeV/c]");
		h_ptCoh[iveto]->Sumw2();
	}
	
	// Kinematics:
	
	h_t_vs_theta = new TH2F("t_vs_theta", "; #theta_{#gamma#gamma} [#circ]; -t [GeV/c^{2}]", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 550, -4.0, 1.0);
	
	// in CM frame of eta:
	
	//h_e_vs_theta  = new TH2F("h_e_vs_theta", "; #theta_{cm} [#circ]; E_{CM} [GeV]", 500, 0.0, 1.0, 500, 0.0, 1.0);
	//h_deltaPhi_CM = new TH1F("deltaPhi_CM", "; #Delta#phi_{cm} [#circ]", 3600, -360.0, 360.0);
	
	// Variation of invariant mass vs. thrown z position:
	
	h_mgg_vs_vertexZ = new TH2F("mgg_vs_vertexZ", "; z_{thrown} [cm]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]",
		500, 0.0, 500.0, 300, 0.0, 1.2);
	
	// Variation of invariant mass vs. thrown radial distance:
	
	h_mgg_vs_vertexR = new TH2F("mgg_vs_vertexR", "; r_{thrown} [cm]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]",
		500, 0.0, 5.0, 300, 0.0, 1.2);
	
	// DeltaPhi distributions from BCAL and SC:
	
	h_scDeltaPhi             = new TH2F("scDeltaPhi", 
		"#phi_{SC} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		650, 0.0, 6.5, 2000, -360.0, 360.0);
	
	h_scDeltaPhi_singleHit   = new TH2F("scDeltaPhi_singleHit", 
		"#phi_{SC} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		650, 0.0, 6.5, 2000, -360.0, 360.0);
	
	h_bcalDeltaPhi           = new TH2F("bcalDeltaPhi", 
		"#phi_{BCAL} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		650, 0.0, 6.5, 2000, -360.0, 360.0);
	
	h_bcalDeltaPhi_singleHit = new TH2F("bcalDeltaPhi_singleHit", 
		"#phi_{BCAL} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		650, 0.0, 6.5, 2000, -360.0, 360.0);
	
	// SC Multiplicities:
	
	h_nSC          = new TH1F("nSC",          "Number of 'good' SC hits", 10, -0.5, 9.5);
	h_nSC_nobcal   = new TH1F("nSC_nobcal",   "Number of 'good' SC hits", 10, -0.5, 9.5);
	h_nSC_onebcal  = new TH1F("nSC_onebcal",  "Number of 'good' SC hits", 10, -0.5, 9.5);
	h_nSC_bcalVeto = new TH1F("nSC_bcalVeto", "Number of 'good' SC hits", 10, -0.5, 9.5);
	h_nSC_extra    = new TH1F("nSC_extra", "Number of non-coplanar SC hits after applying BCAL veto", 10, -0.5, 9.5);
	
	// Distribution of accepted events in FCAL:
	
	h_xy1 = new TH2F("xy1", "Position of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -100., 100., 500, -100., 100.);
	h_xy2 = new TH2F("xy2", "Position of Shower 2; x_{2} [cm]; y_{2} [cm]", 500, -100., 100., 500, -100., 100.);
	
	// Angular matrix needed for acceptance calculation:
	
	if(m_FillThrown) {
		for(int iveto=0; iveto<m_nVetos; iveto++) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_veto_%d",iveto), 
				Form("Veto Option %d; #theta(thrown) [#circ]; #theta(rec) [#circ]", iveto),
				nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
				nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
				nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
			hMatrix->SetDirectory(0);
			h_AngularMatrix_vetos.push_back(hMatrix);
		}
	}
	
	return;
}

void EtaAna::ResetDefaultHists()
{
	h_mcVertex->Reset();
	h_mcVertexAccepted->Reset();
	
	h_fcalRFdt->Reset();
	h_bcalRFdt->Reset();
	h_tofRFdt->Reset();
	h_scRFdt->Reset();
	h_beamRFdt->Reset();
	h_beamRFdt_cut->Reset();
	
	h_tofMatch->Reset();
	
	for(int iveto=0; iveto<m_nVetos; iveto++) {
		h_elasticity[iveto]->Reset();
		h_elasticityConstr[iveto]->Reset();
		h_mgg[iveto]->Reset();
		h_mggConstr[iveto]->Reset();
		h_mggConstr_coh[iveto]->Reset();
		h_mm[iveto]->Reset();
		h_mm_coh[iveto]->Reset();
		h_mm_elas[iveto]->Reset();
		h_mm_elas_coh[iveto]->Reset();
		h_pt[iveto]->Reset();
		h_ptCoh[iveto]->Reset();
	}
	h_t_vs_theta->Reset();
	
	h_mgg_vs_vertexZ->Reset();
	h_mgg_vs_vertexR->Reset();
	
	h_scDeltaPhi->Reset();
	h_scDeltaPhi_singleHit->Reset();
	h_bcalDeltaPhi->Reset();
	h_bcalDeltaPhi_singleHit->Reset();
	
	h_nSC->Reset();
	h_nSC_nobcal->Reset();
	h_nSC_onebcal->Reset();
	h_nSC_bcalVeto->Reset();
	h_nSC_extra->Reset();
	
	h_xy1->Reset();
	h_xy2->Reset();
	
	if(h_AngularMatrix_vetos.size()) {
		for(int i=0; i<h_AngularMatrix_vetos.size(); i++) {
			h_AngularMatrix_vetos[i]->Reset();
		}
	}
	
	return;
}

void EtaAna::WriteDefaultHists()
{
	printf("  Writing default analysis option histograms...\n");
	
	h_mcVertex->Write();
	h_mcVertexAccepted->Write();
	
	h_fcalRFdt->Write();
	h_bcalRFdt->Write();
	h_tofRFdt->Write();
	h_scRFdt->Write();
	h_beamRFdt->Write();
	h_beamRFdt_cut->Write();
	
	h_tofMatch->Write();
	
	for(int iveto=0; iveto<m_nVetos; iveto++) {
		
		TDirectoryFile *dirVeto = new TDirectoryFile(Form("VetoOption%d",iveto), Form("VetoOption%d",iveto));
		dirVeto->cd();
		
		h_elasticity[iveto]->Write();
		h_elasticityConstr[iveto]->Write();
		h_mgg[iveto]->Write();
		h_mggConstr[iveto]->Write();
		h_mggConstr_coh[iveto]->Write();
		h_mm[iveto]->Write();
		h_mm_coh[iveto]->Write();
		h_mm_elas[iveto]->Write();
		h_mm_elas_coh[iveto]->Write();
		h_pt[iveto]->Write();
		h_ptCoh[iveto]->Write();
		
		dirVeto->cd("../");
	}
	h_t_vs_theta->Write();
	
	h_mgg_vs_vertexZ->Write();
	h_mgg_vs_vertexR->Write();
	
	h_scDeltaPhi->Write();
	h_scDeltaPhi_singleHit->Write();
	h_bcalDeltaPhi->Write();
	h_bcalDeltaPhi_singleHit->Write();
	
	h_nSC->Write();
	h_nSC_nobcal->Write();
	h_nSC_onebcal->Write();
	h_nSC_bcalVeto->Write();
	h_nSC_extra->Write();
	
	h_xy1->Write();
	h_xy2->Write();
	
	TDirectory *dirMatrix = new TDirectoryFile("AngularMatrix","AngularMatrix");
	dirMatrix->cd();
	if(h_AngularMatrix_vetos.size()) {
		for(int i=0; i<h_AngularMatrix_vetos.size(); i++) {
			h_AngularMatrix_vetos[i]->Write();
		}
	}
	dirMatrix->cd("../");
	
	printf("  Done.\n");
	return;
}
