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
	
	// Plot SC Timing Distribution for monitoring:
	for(int isc = 0; isc < m_nsc; isc++) {
		double locT  = m_scT[isc] - m_rfTime;
		double locdE = m_scdE[isc];
		if(locdE > 0.0002) {
			h_scRFdt->Fill(locT);
		}
	}
	
	// Plot TOF Timing Distribution for monitoring:
	for(int itof = 0; itof < m_ntof; itof++) {
		TVector3 locPos(m_tofX[itof], m_tofY[itof], m_tofZ[itof]);
		locPos -= m_vertex;
		double locT = m_tofT[itof] - (locPos.Mag()/m_c) - m_rfTime;
		h_tofRFdt->Fill(locT);
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
				
				if((1.0 < locT) && (locT < 9.0) && (locdE > 0.0002)) {
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
				if(isElastic && locVetoOptions[5]) {
					h_beamRFdt_cut->Fill(brfdt);
				}
				
				TLorentzVector k_beam(0.0, 0.0, eb, eb);
				TLorentzVector p_eta(pggx, pggy, pggz, Egg);
				double t    = (k_beam-p_eta)*(k_beam-p_eta);
				double logt = log10(-t);
				
				// acceptance correction:
				double locAcc = 1.0;
				if(h_acceptance != NULL) {
					locAcc = h_acceptance->GetBinContent(h_acceptance->GetXaxis()->FindBin(eb), 
						h_acceptance->GetYaxis()->FindBin(prodTheta));
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
				if(isEta) {
					h_pT_vs_elas->Fill(Egg/etaEnergy, pggt/pTCalc, fillWeight);
					if(locVetoOptions[5]) h_pT_vs_elas_cut->Fill(Egg/etaEnergy, pggt/pTCalc, fillWeight);
				}
				*/
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
				
				if(isEta && isElastic && locVetoOptions[5]) {
					h_e_vs_theta->Fill(acos((px1*pggx + py1*pggy + pz1*pggz)/(e1*sqrt(pow(pggt,2.0)+pow(pggz,2.0)))), e1/Egg);
					h_e_vs_theta->Fill(acos((px2*pggx + py2*pggy + pz2*pggz)/(e2*sqrt(pow(pggt,2.0)+pow(pggz,2.0)))), e2/Egg);
				}
				*/
				if(isEta && isElastic) {
					if(locVetoOptions[5]) {
						h_xy1->Fill(pos1.X(),pos1.Y());
						h_xy2->Fill(pos2.X(),pos2.Y());
						h_t_vs_theta->Fill(prodTheta, logt, fillWeight);
					}
					
					// plot SC-DeltaPhi:
					for(int isc = 0; isc < m_nsc; isc++) {
						double locT  = m_scT[isc] - m_rfTime;
						double locdE = m_scdE[isc];
						if((1.0 < locT) && (locT < 9.0) && (locdE > 0.0002)) {
							double deltaPhi = m_scPhi[isc]-prodPhi;
							h_scDeltaPhi->Fill(prodTheta, deltaPhi);
						}
					}
					
					// plot BCAL-DeltaPhi:
					for(int ibcal = 0; ibcal < locGoodBCALShowers.size(); ibcal++) {
						double deltaPhi = GetBCALPosition(locGoodBCALShowers[ibcal]).Phi() * TMath::RadToDeg() - prodPhi;
						h_bcalDeltaPhi->Fill(prodTheta, deltaPhi);
					}
				}
				
				for(int iveto=0; iveto<m_nVetos; iveto++) {
					
					if(locVetoOptions[iveto]==0) continue;
					if(isEta) {
						h_elasticity[iveto]->Fill(prodTheta, Egg/etaEnergy, fillWeight);
						h_elasticityConstr[iveto]->Fill(prodTheta, elasConstr, fillWeight);
						
						h_pt[iveto]->Fill(prodTheta, pggt - pTCalc, fillWeight);
						h_ptCoh[iveto]->Fill(prodTheta, pggt - pTCalcCoh, fillWeight);
						
						h_mm[iveto]->Fill(prodTheta, mmSq, fillWeight);
						h_mm_coh[iveto]->Fill(prodTheta, (mmSqCoh - pow(ParticleMass(Helium),2.0)), fillWeight);
						if(isElastic) {
							h_mm_elas[iveto]->Fill(prodTheta, mmSq, fillWeight);
							h_mm_elas_coh[iveto]->Fill(prodTheta, (mmSqCoh - pow(ParticleMass(Helium),2.0)), fillWeight);
						}
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
				
				if((m_nmc>0) && locVetoOptions[5] && isElastic) {
					h_mgg_vs_vertex->Fill(m_mcZ[0], invmassConstr, fillWeight);
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

