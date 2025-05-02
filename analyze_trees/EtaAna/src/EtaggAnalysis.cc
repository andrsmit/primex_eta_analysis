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
		m_FCALEnergyCut, m_FCALExtraEnergyCut, 2.0, m_FCALRFCut);
	
	// Apply multiplicity cut on the number of FCAL showers: 
	int locNFCALShowersGood = (int)locGoodFCALShowers.size();
	if((locNFCALShowersEnergyCut!=2) || (locNFCALShowersGood!=2)) return;
	
	//-------------------------------------------//
	// Get list of 'good' BCAL showers:
	
	vector<int> locGoodBCALShowers; locGoodBCALShowers.clear();
	int locNBCALShowers = GetBCALShowerList(locGoodBCALShowers, m_BCALEnergyCut, m_BCALRFCut);
	
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
	// RF Timing distributions:
	
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
	
	for(int ishow=0; ishow<(locNFCALShowersGood-1); ishow++) {
		
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
		
		for(int jshow=(ishow+1); jshow<locNFCALShowersGood; jshow++) {
			
			int show2 = locGoodFCALShowers[jshow];
			TVector3 pos2 = GetFCALPosition(show2);
			
			double  e2 = m_fcalE[show2];
			
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
			
			/*
			I'm moving the TOF veto to later on, because I want to check the events which
			get removed by this veto for a signature of Compton scattering.
			If it is Compton where the photon is converted to e+e- downstream, 
			then one of the showers should follow the Klein-Nishina kinematics.
			
			3/11/24: 
			I saw no such signature of Compton, so I'm reinstating this placement of the TOF veto.
			*/
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
			
			for(int iveto=0; iveto<m_nVetos; iveto++) {
				if(!IsHadronicVeto(iveto, locNBCALShowers, locNBCALShowers_1ns, locNSCHits, locNSCHits_coplanar, 
					prodPhi, locBCALPhi, locBCALTheta, locBCALRFDT)
				) locVetoOptions[iveto] = 1;
			}
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(int igam=0; igam<locNBeamPhotons; igam++) {
				
				int ibeam = locGoodBeamPhotons[igam].first;
				double fillWeight = locGoodBeamPhotons[igam].second;
				
				double eb    = m_beamE[ibeam];
				double brfdt = m_beamT[ibeam] - m_rfTime;
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				double etaEnergyCoh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
				
				// Calculate the energy of the eta meson, assuming production on a free nucleon:
				double etaEnergy = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(Proton));
				
				// Apply a cut on the elasticity
				//  (ratio of measured energy of 2-photons, to the calculated energy above):
				bool isElastic = IsElasticCut(Egg, etaEnergy, prodTheta);
				
				// set a variable to indicate if the two-photon mass is consistent with an eta meson:
				bool isEta = IsEtaCut(invmass);
				
				//----------------------------------------//
				/*
				// Checking for Compton signature:
				
				if(isElastic && isTOFVeto)
				{
					double ecomp1 = eb / (1. + (eb/ParticleMass(Electron))*(1.-cos(pos1.Theta())));
					double ecomp2 = eb / (1. + (eb/ParticleMass(Electron))*(1.-cos(pos2.Theta())));
					
					h_ecomp1->Fill(prodTheta, ecomp1/e1, fillWeight);
					h_ecomp2->Fill(prodTheta, ecomp2/e2, fillWeight);
					if(pos1.Theta() < pos2.Theta()) {
						h_ecomp->Fill(prodTheta, ecomp1/e1, fillWeight);
					}
					else {
						h_ecomp->Fill(prodTheta, ecomp2/e2, fillWeight);
					}
				}
				else if(isElastic)
				{
					// baseline comparison
					
					double ecomp1 = eb / (1. + (eb/ParticleMass(Electron))*(1.-cos(pos1.Theta())));
					double ecomp2 = eb / (1. + (eb/ParticleMass(Electron))*(1.-cos(pos2.Theta())));
					
					h_ecomp1_clean->Fill(prodTheta, ecomp1/e1, fillWeight);
					h_ecomp2_clean->Fill(prodTheta, ecomp2/e2, fillWeight);
					if(pos1.Theta() < pos2.Theta()) {
						h_ecomp_clean->Fill(prodTheta, ecomp1/e1, fillWeight);
					}
					else {
						h_ecomp_clean->Fill(prodTheta, ecomp2/e2, fillWeight);
					}
				}
				if(isTOFVeto) continue;
				*/
				//----------------------------------------//
				
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
				
				// plot opening angle:
				
				if(isElastic) {
					h_openAngle->Fill(prodTheta, cos12, fillWeight);
					if(isEta) {
						h_openAngleCut->Fill(prodTheta, cos12, fillWeight);
					}
				}
				if(isEta) {
					h_openAngleVsElas->Fill(Egg/etaEnergy, cos12, fillWeight);
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
				double e1c_coh = e1/(1.+sigr) + (etaEnergyCoh-e2)/(1.+(1./sigr));
				double e2c_coh = etaEnergyCoh - e1c_coh;
				double invmassConstrCoh = sqrt(2.*e1c_coh*e2c_coh*(1.-cos12)); // energy-constrained invariant mass
				
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
				
				double elasConstr = (e1mc+e2mc)/etaEnergyCoh;
				
				//-----------------------------------------------------//
				// Compare measured transverse momentum of eta to calculated value:
				
				double pTCalc    = sqrt(pow(etaEnergy,2.0)     - pow(ParticleMass(Eta),2.0))*sin(prodTheta*TMath::DegToRad());
				double pTCalcCoh = sqrt(pow(etaEnergyCoh,2.0) - pow(ParticleMass(Eta),2.0))*sin(prodTheta*TMath::DegToRad());
				
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
				// Hybrid Mass
				
				double hmass        = (invmass/ParticleMass(Eta))*cos(TMath::Pi()/4.0) - (Egg/eb)*sin(TMath::Pi()/4.0);
				double hmassCorr    = (invmass/ParticleMass(Eta))*cos(TMath::Pi()/4.0) - (Egg/etaEnergy)*sin(TMath::Pi()/4.0);
				double hmassCorrCoh = (invmass/ParticleMass(Eta))*cos(TMath::Pi()/4.0) - (Egg/etaEnergyCoh)*sin(TMath::Pi()/4.0);
				
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
				
				if(locVetoOptions[6]) {
					if(prodTheta<2.0) {
						h_elasVSmgg_low->Fill(invmass/ParticleMass(Eta), Egg/eb, fillWeight);
						h_elasCorrVSmgg_low->Fill(invmass/ParticleMass(Eta), Egg/etaEnergy, fillWeight);
					} else {
						h_elasVSmgg_high->Fill(invmass/ParticleMass(Eta), Egg/eb, fillWeight);
						h_elasCorrVSmgg_high->Fill(invmass/ParticleMass(Eta), Egg/etaEnergy, fillWeight);
					}
				}
				
				// plot measured 't' divided by calculated 't' as a function of angle:
				/*
				double tCalc = CalculateT(prodTheta, eb);
				
				h_trec_vs_tcalc->Fill(prodTheta, t/tCalc, fillWeight);
				if(isEta) {
					h_trec_vs_tcalc_massCut->Fill(prodTheta, t/tCalc, fillWeight);
					if(isElastic) {
						h_trec_vs_tcalc_massElasCut->Fill(prodTheta, t/tCalc, fillWeight);
					}
				}
				*/
				
				if(isEta && isElastic) {
					
					if(locVetoOptions[6]) {
						h_xy1->Fill(pos1.X(),pos1.Y());
						h_xy2->Fill(pos2.X(),pos2.Y());
						h_t_vs_theta->Fill(prodTheta, logt, fillWeight);
						
						// plot likelihood of TOF match as a function of radius:
						/*
						double rho1 = sqrt(pow(pos1.X(),2.0) + pow(pos1.Y(),2.0));
						double rho2 = sqrt(pow(pos2.X(),2.0) + pow(pos2.Y(),2.0));
						h_tofMatch->Fill(rho1, tof_match1, fillWeight);
						h_tofMatch->Fill(rho2, tof_match2, fillWeight);
						*/
					}
					
					//=======================================//
					// SC Plots:
					
					// Loop over all sc hits and plot timing:
					for(int isc=0; isc<m_nsc; isc++) {
						double locdE = m_scdE[isc];
						double locdT = m_scT[isc] - m_rfTime;
						h_scRFdt->Fill(prodTheta, locdT, fillWeight);
						if(locdE>0.0002) {
							h_scRFdt_cut->Fill(prodTheta, locdT, fillWeight);
						}
					}
					
					// Loop over in-time SC hits and plot energy loss and deltaPhi:
					for(int isc=0; isc<locGoodSCHits.size(); isc++) {
						int hitIndex       = locGoodSCHits[isc];
						double locdE       = 1.e3 * m_scdE[hitIndex]; // [MeV]
						double locDeltaPhi = prodPhi - m_scPhi[hitIndex];
						if(locDeltaPhi>0.0) locDeltaPhi-=180.0;
						else                locDeltaPhi+=180.0;
						
						h_scdE->Fill(prodTheta, locdE, fillWeight);
						h_scDeltaPhi->Fill(prodTheta, locDeltaPhi, fillWeight);
						
						if(locGoodSCHits.size()==1) {
							h_scdE_singleHit->Fill(prodTheta, locdE, fillWeight);
							h_scDeltaPhi_singleHit->Fill(prodTheta, locDeltaPhi, fillWeight);
						}
					}
					
					// Plot number of SC hits:
					h_nSC->Fill(prodTheta, locNSCHits, fillWeight);
					if(locNBCALShowers==0) h_nSC_nobcal->Fill(prodTheta, locNSCHits, fillWeight);
					if(locNBCALShowers==1) h_nSC_onebcal->Fill(prodTheta, locNSCHits, fillWeight);
					if(locVetoOptions[4]) {
						h_nSC_bcalVeto->Fill(prodTheta, locNSCHits, fillWeight);
						h_nSC_extra->Fill(prodTheta, locNSCHits-locNSCHits_coplanar, fillWeight);
					}
					
					//=======================================//
					// BCAL Plots:
					
					// Loop over all BCAL showers and plot timing:
					for(int ishow = 0; ishow < m_nbcal; ishow++) {
						TVector3 locPos = GetBCALPosition(ishow);
						double locT = m_bcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
						h_bcalRFdt->Fill(prodTheta, locT, fillWeight);
					}
					
					// Loop over in-time BCAL showers and plot energy, deltaPhi, theta:
					for(int ibcal = 0; ibcal < locGoodBCALShowers.size(); ibcal++) {
						TVector3 locPos    = GetBCALPosition(locGoodBCALShowers[ibcal]);
						double locTheta    = locPos.Theta() * TMath::RadToDeg();
						double locDeltaPhi = prodPhi - (locPos.Phi() * TMath::RadToDeg());
						if(locDeltaPhi>0.0) locDeltaPhi-=180.0;
						else                locDeltaPhi+=180.0;
						double locE        = m_bcalE[locGoodBCALShowers[ibcal]];
						double locT        = m_bcalT[locGoodBCALShowers[ibcal]] - (locPos.Mag()/m_c) - m_rfTime;
						
						h_bcalEnergy->Fill(prodTheta, locE, fillWeight);
						h_bcalTheta->Fill(prodTheta, locTheta, fillWeight);
						h_bcalDeltaPhi->Fill(prodTheta, locDeltaPhi, fillWeight);
						if(locGoodBCALShowers.size()==1) {
							h_bcalEnergy_singleHit->Fill(prodTheta, locE, fillWeight);
							h_bcalTheta_singleHit->Fill(prodTheta, locTheta, fillWeight);
							h_bcalDeltaPhi_singleHit->Fill(prodTheta, locDeltaPhi, fillWeight);
							
							h_bcalThetaVSTime->Fill(locT, locTheta, fillWeight);
							h_bcalDeltaPhiVSTime->Fill(locT, locDeltaPhi, fillWeight);
							h_bcalThetaVSEnergy->Fill(locE, locTheta, fillWeight);
						}
					}
					
					// Plot number of BCAL showers:
					h_nBCAL->Fill(prodTheta, locNBCALShowers, fillWeight);
					h_nBCAL_1ns->Fill(prodTheta, locNBCALShowers_1ns, fillWeight);
					
					//=======================================//
					// FCAL Plots:
					
					// Plot FCAL-RF dt for accepted events:
					for(int ishow = 0; ishow < m_nfcal; ishow++) {
						TVector3 locPos = GetFCALPosition(ishow);
						double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
						h_fcalRFdt->Fill(locT, fillWeight);
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
					
					h_mgg[iveto]->Fill(prodTheta, invmass, fillWeight);
					h_mggConstr[iveto]->Fill(prodTheta, invmassConstr, fillWeight);
					h_mggConstr_coh[iveto]->Fill(prodTheta, invmassConstrCoh, fillWeight);
					if(isElastic) {
						h_mgg_cut[iveto]->Fill(prodTheta, invmass, fillWeight);
						h_mggConstr_cut[iveto]->Fill(prodTheta, invmassConstr, fillWeight);
						h_mggConstr_coh_cut[iveto]->Fill(prodTheta, invmassConstrCoh, fillWeight);
					}
					
					h_hmass[iveto]->Fill(prodTheta, hmass, fillWeight);
					h_hmassCorr[iveto]->Fill(prodTheta, hmassCorr, fillWeight);
					h_hmassCorrCoh[iveto]->Fill(prodTheta, hmassCorrCoh, fillWeight);
					if(isEta) {
						h_hmass_etaCut[iveto]->Fill(prodTheta, hmass, fillWeight);
						h_hmassCorr_etaCut[iveto]->Fill(prodTheta, hmassCorr, fillWeight);
						h_hmassCorrCoh_etaCut[iveto]->Fill(prodTheta, hmassCorrCoh, fillWeight);
					}
					
					if(m_FillThrown && IsEtaCut(invmassConstr) && isElastic) {
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
	h_tofRFdt      = new TH1F( "tof_rf_dt",      "t_{TOF} - t_{RF}; [ns]", 10000, -100., 100.);
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
		
		h_mgg_cut[iveto] = new TH2F(Form("mgg_cut_veto_%d",iveto), 
			Form("Two-Photon Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mgg_cut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mgg_cut[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma} [GeV/c^{2}]");
		h_mgg_cut[iveto]->Sumw2();
		
		h_mggConstr[iveto] = new TH2F(Form("mgg_const_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr[iveto]->Sumw2();
		
		h_mggConstr_cut[iveto] = new TH2F(Form("mgg_const_cut_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr_cut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_cut[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr_cut[iveto]->Sumw2();
		
		h_mggConstr_coh[iveto] = new TH2F(Form("mgg_const_coh_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr_coh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_coh[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr_coh[iveto]->Sumw2();
		
		h_mggConstr_coh_cut[iveto] = new TH2F(Form("mgg_const_coh_cut_veto_%d",iveto), 
			Form("Energy-Constrained Invariant Mass (Veto Option %d)",iveto), 
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		h_mggConstr_coh_cut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_mggConstr_coh_cut[iveto]->GetYaxis()->SetTitle("m_{#gamma#gamma}^{Constr} [GeV/c^{2}]");
		h_mggConstr_coh_cut[iveto]->Sumw2();
		
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
		
		h_hmass[iveto] = new TH2F(Form("hmass_veto_%d",iveto),
			Form("Hybrid Mass (Veto Option %d)",iveto),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -1.0, 1.0);
		h_hmass[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_hmass[iveto]->GetYaxis()->SetTitle("Hybrid Mass");
		h_hmass[iveto]->Sumw2();
		
		h_hmassCorr[iveto] = new TH2F(Form("hmassCorr_veto_%d",iveto),
			Form("Hybrid Mass (Veto Option %d)",iveto),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -1.0, 1.0);
		h_hmassCorr[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_hmassCorr[iveto]->GetYaxis()->SetTitle("Hybrid Mass");
		h_hmassCorr[iveto]->Sumw2();
		
		h_hmassCorrCoh[iveto] = new TH2F(Form("hmassCorrCoh_veto_%d",iveto),
			Form("Hybrid Mass (Veto Option %d)",iveto),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -1.0, 1.0);
		h_hmassCorrCoh[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_hmassCorrCoh[iveto]->GetYaxis()->SetTitle("Hybrid Mass");
		h_hmassCorrCoh[iveto]->Sumw2();
		
		h_hmass_etaCut[iveto] = new TH2F(Form("hmass_etaCut_veto_%d",iveto),
			Form("Hybrid Mass (Veto Option %d)",iveto),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -1.0, 1.0);
		h_hmass_etaCut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_hmass_etaCut[iveto]->GetYaxis()->SetTitle("Hybrid Mass");
		h_hmass_etaCut[iveto]->Sumw2();
		
		h_hmassCorr_etaCut[iveto] = new TH2F(Form("hmassCorr_etaCut_veto_%d",iveto),
			Form("Hybrid Mass (Veto Option %d)",iveto),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -1.0, 1.0);
		h_hmassCorr_etaCut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_hmassCorr_etaCut[iveto]->GetYaxis()->SetTitle("Hybrid Mass");
		h_hmassCorr_etaCut[iveto]->Sumw2();
		
		h_hmassCorrCoh_etaCut[iveto] = new TH2F(Form("hmassCorrCoh_etaCut_veto_%d",iveto),
			Form("Hybrid Mass (Veto Option %d)",iveto),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, -1.0, 1.0);
		h_hmassCorrCoh_etaCut[iveto]->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
		h_hmassCorrCoh_etaCut[iveto]->GetYaxis()->SetTitle("Hybrid Mass");
		h_hmassCorrCoh_etaCut[iveto]->Sumw2();
		
	}
	
	// Kinematics:
	
	h_t_vs_theta = new TH2F("t_vs_theta", "; #theta_{#gamma#gamma} [#circ]; -t [GeV/c^{2}]", 
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 550, -4.0, 1.0);
	
	h_elasVSmgg_low = new TH2F("elasVSmgg_low", "; m_{#gamma#gamma}/m_{#eta}; E_{#gamma#gamma}/E_{tag}",
		500, 0.0, 2.0, 500, 0.0, 2.0);
	h_elasVSmgg_high = new TH2F("elasVSmgg_high", "; m_{#gamma#gamma}/m_{#eta}; E_{#gamma#gamma}/E_{tag}",
		500, 0.0, 2.0, 500, 0.0, 2.0);
	h_elasCorrVSmgg_low = new TH2F("elasCorrVSmgg_low", "; m_{#gamma#gamma}/m_{#eta}; E_{#gamma#gamma}/E_{#eta}",
		500, 0.0, 2.0, 500, 0.0, 2.0);
	h_elasCorrVSmgg_high = new TH2F("elasCorrVSmgg_high", "; m_{#gamma#gamma}/m_{#eta}; E_{#gamma#gamma}/E_{#eta}",
		500, 0.0, 2.0, 500, 0.0, 2.0);
	
	// in CM frame of eta:
	
	//h_e_vs_theta  = new TH2F("h_e_vs_theta", "; #theta_{cm} [#circ]; E_{CM} [GeV]", 500, 0.0, 1.0, 500, 0.0, 1.0);
	//h_deltaPhi_CM = new TH1F("deltaPhi_CM", "; #Delta#phi_{cm} [#circ]", 3600, -360.0, 360.0);
	
	// Variation of invariant mass vs. thrown z position:
	
	h_mgg_vs_vertexZ = new TH2F("mgg_vs_vertexZ", "; z_{thrown} [cm]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]",
		500, 0.0, 500.0, 300, 0.0, 1.2);
	
	// Variation of invariant mass vs. thrown radial distance:
	
	h_mgg_vs_vertexR = new TH2F("mgg_vs_vertexR", "; r_{thrown} [cm]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]",
		500, 0.0, 5.0, 300, 0.0, 1.2);
	
	h_openAngle       = new TH2F("openAngle", "; #theta_{#gamma#gamma} [#circ]; cos#theta_{12}", 550, 0.0, 5.5, 500, 0.8, 1.0);
	h_openAngleCut    = new TH2F("openAngleCut", "(#eta mass cut); #theta_{#gamma#gamma} [#circ]; cos#theta_{12}", 550, 0.0, 5.5, 500, 0.8, 1.0);
	h_openAngleVsElas = new TH2F("openAngleVsElas", "; E_{#gamma#gamma}/E_{#eta}; cos#theta_{12}", 500, 0.0, 2.0, 500, 0.8, 1.0);
	
	//------------------------------------------------------------------------------------------//
	// SC Plots:
	
	h_scRFdt       = new TH2F("scRFdt", 
		"t_{SC} - t_{RF}; #theta_{#gamma#gamma} [#circ]; [ns]", 100, 0.0, 5.0, 1000, -100., 100.);
	h_scRFdt_cut   = new TH2F("scRFdt_cut", 
		"t_{SC} - t_{RF}; #theta_{#gamma#gamma} [#circ]; [ns]", 100, 0.0, 5.0, 1000, -100., 100.);
	
	h_scdE           = new TH2F("scdE", 
		"SC Energy Loss; #theta_{#gamma#gamma} [#circ]; dE [MeV]", 100, 0.0, 5.0, 2000, 0.0, 10.0);
	h_scdE_singleHit = new TH2F("scdE_singleHit", 
		"SC Energy Loss; #theta_{#gamma#gamma} [#circ]; dE [MeV]", 100, 0.0, 5.0, 2000, 0.0, 10.0);
	
	h_scDeltaPhi           = new TH2F("scDeltaPhi", 
		"#phi_{SC} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		500, 0.0, 5.0, 2000, -360.0, 360.0);
	h_scDeltaPhi_singleHit = new TH2F("scDeltaPhi_singleHit", 
		"#phi_{SC} - #phi_{#gamma#gamma}; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		500, 0.0, 5.0, 2000, -360.0, 360.0);
	
	h_nSC          = new TH2F("nSC",          "Number of 'good' SC hits", 100, 0.0, 5.0, 10, -0.5, 9.5);
	h_nSC_nobcal   = new TH2F("nSC_nobcal",   "Number of 'good' SC hits", 100, 0.0, 5.0, 10, -0.5, 9.5);
	h_nSC_onebcal  = new TH2F("nSC_onebcal",  "Number of 'good' SC hits", 100, 0.0, 5.0, 10, -0.5, 9.5);
	h_nSC_bcalVeto = new TH2F("nSC_bcalVeto", "Number of 'good' SC hits", 100, 0.0, 5.0, 10, -0.5, 9.5);
	h_nSC_extra    = new TH2F("nSC_extra",    "Number of non-coplanar SC hits after applying BCAL veto", 100, 0.0, 5.0, 10, -0.5, 9.5);
	
	//------------------------------------------------------------------------------------------//
	// BCAL Plots:
	
	h_bcalRFdt     = new TH2F("bcalRFdt", 
		"t_{BCAL} - t_{RF}; #theta_{#gamma#gamma} [#circ]; [ns]", 100, 0.0, 5.0, 1000, -100., 100.);
	
	h_bcalEnergy           = new TH2F("bcalEnergy:",
		"BCAL Shower Energy; #theta_{#gamma#gamma} [#circ]; [GeV]", 100, 0.0, 5.0, 2000, 0.0, 2.0);
	h_bcalEnergy_singleHit = new TH2F("bcalEnergy_singleHit:",
		"BCAL Shower Energy; #theta_{#gamma#gamma} [#circ]; [GeV]", 100, 0.0, 5.0, 2000, 0.0, 2.0);
	
	h_bcalTheta           = new TH2F("bcalTheta", 
		"#theta_{BCAL}; #theta_{#gamma#gamma} [#circ]; [#circ]", 100, 0.0, 5.0, 1000, 0.0, 180.0);
	h_bcalTheta_singleHit = new TH2F("bcalTheta_singleHit", 
		"#theta_{BCAL}; #theta_{#gamma#gamma} [#circ]; [#circ]", 100, 0.0, 5.0, 1000, 0.0, 180.0);
	
	h_bcalDeltaPhi           = new TH2F("bcalDeltaPhi", 
		"#phi_{BCAL} - #phi_{#gamma#gamma} [#circ]; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		100, 0.0, 5.0, 2000, -360.0, 360.0);
	h_bcalDeltaPhi_singleHit = new TH2F("bcalDeltaPhi_singleHit", 
		"#phi_{BCAL} - #phi_{#gamma#gamma} [#circ]; #theta_{#gamma#gamma} [#circ]; #Delta#phi [#circ]",
		100, 0.0, 5.0, 2000, -360.0, 360.0);
	
	h_bcalThetaVSTime = new TH2F("bcalThetaVSTime",
		"; t_{BCAL} - t_{RF} [ns]; #theta_{BCAL} [#circ]", 
		1300, -1.0, 12.0, 1000, 0.0, 180.0);
	h_bcalDeltaPhiVSTime = new TH2F("bcalDeltaPhiVSTime",
		"; t_{BCAL} - t_{RF} [ns]; #phi_{BCAL} - #phi_{#gamma#gamma} [#circ]", 
		1300, -1.0, 12.0, 1000, -360.0, 360.0);
	h_bcalThetaVSEnergy = new TH2F("bcalThetaVSEnergy",
		"; t_{BCAL} - t_{RF} [ns]; E_{BCAL} [GeV]", 
		1000, 0.0, 180.0, 1000, 0.0, 2.0);
	
	h_nBCAL     = new TH2F("nBCAL",     "Number of in-time BCAL Showers", 100, 0.0, 5.0, 10, -0.5, 9.5);
	h_nBCAL_1ns = new TH2F("nBCAL_1ns", "Number of <1ns BCAL Showers",    100, 0.0, 5.0, 10, -0.5, 9.5);
	
	//------------------------------------------------------------------------------------------//
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
	
	// Looking for Compton signature in events which get removed by TOF veto:
	/*
	h_ecomp1 = new TH2F("ecomp1", "; #theta_{#gamma#gamma} [#circ]; E_{Comp,1}#left(E_{#gamma},#theta_{1}#right) / E_{1}",
		500, 0.0, 5.0, 500, 0.0, 2.0);
	h_ecomp2 = new TH2F("ecomp2", "; #theta_{#gamma#gamma} [#circ]; E_{Comp,2}#left(E_{#gamma},#theta_{2}#right) / E_{2}",
		500, 0.0, 5.0, 500, 0.0, 2.0);
	h_ecomp  = new TH2F("ecomp",  "; #theta_{#gamma#gamma} [#circ]; E_{Comp}#left(E_{#gamma},#theta#right) / E",
		500, 0.0, 5.0, 500, 0.0, 2.0);
	
	h_ecomp1_clean = new TH2F("ecomp1_clean", "; #theta_{#gamma#gamma} [#circ]; E_{Comp,1}#left(E_{#gamma},#theta_{1}#right) / E_{1}",
		500, 0.0, 5.0, 500, 0.0, 2.0);
	h_ecomp2_clean = new TH2F("ecomp2_clean", "; #theta_{#gamma#gamma} [#circ]; E_{Comp,2}#left(E_{#gamma},#theta_{2}#right) / E_{2}",
		500, 0.0, 5.0, 500, 0.0, 2.0);
	h_ecomp_clean  = new TH2F("ecomp_clean",  "; #theta_{#gamma#gamma} [#circ]; E_{Comp}#left(E_{#gamma},#theta#right) / E",
		500, 0.0, 5.0, 500, 0.0, 2.0);
	*/
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
		h_mgg_cut[iveto]->Reset();
		h_mggConstr[iveto]->Reset();
		h_mggConstr_cut[iveto]->Reset();
		h_mggConstr_coh[iveto]->Reset();
		h_mggConstr_coh_cut[iveto]->Reset();
		
		h_mm[iveto]->Reset();
		h_mm_coh[iveto]->Reset();
		h_mm_elas[iveto]->Reset();
		h_mm_elas_coh[iveto]->Reset();
		
		h_hmass[iveto]->Reset();
		h_hmassCorr[iveto]->Reset();
		h_hmassCorrCoh[iveto]->Reset();
		h_hmass_etaCut[iveto]->Reset();
		h_hmassCorr_etaCut[iveto]->Reset();
		h_hmassCorrCoh_etaCut[iveto]->Reset();
		
		h_pt[iveto]->Reset();
		h_ptCoh[iveto]->Reset();
	}
	h_t_vs_theta->Reset();
	
	h_mgg_vs_vertexZ->Reset();
	h_mgg_vs_vertexR->Reset();
	
	h_elasVSmgg_low->Reset();
	h_elasVSmgg_high->Reset();
	h_elasCorrVSmgg_low->Reset();
	h_elasCorrVSmgg_high->Reset();
	
	h_openAngle->Reset();
	h_openAngleCut->Reset();
	h_openAngleVsElas->Reset();
	
	h_scRFdt->Reset();
	h_scRFdt_cut->Reset();
	h_scdE->Reset();
	h_scdE_singleHit->Reset();
	h_scDeltaPhi->Reset();
	h_scDeltaPhi_singleHit->Reset();
	h_nSC->Reset();
	h_nSC_nobcal->Reset();
	h_nSC_onebcal->Reset();
	h_nSC_bcalVeto->Reset();
	h_nSC_extra->Reset();
	
	h_bcalRFdt->Reset();
	h_bcalEnergy->Reset();
	h_bcalEnergy_singleHit->Reset();
	h_bcalTheta->Reset();
	h_bcalTheta_singleHit->Reset();
	h_bcalDeltaPhi->Reset();
	h_bcalDeltaPhi_singleHit->Reset();
	h_bcalThetaVSTime->Reset();
	h_bcalDeltaPhiVSTime->Reset();
	h_bcalThetaVSEnergy->Reset();
	h_nBCAL->Reset();
	h_nBCAL_1ns->Reset();
	
	h_xy1->Reset();
	h_xy2->Reset();
	
	if(h_AngularMatrix_vetos.size()) {
		for(int i=0; i<h_AngularMatrix_vetos.size(); i++) {
			h_AngularMatrix_vetos[i]->Reset();
		}
	}
	/*
	h_ecomp1->Reset();
	h_ecomp2->Reset();
	h_ecomp->Reset();
	
	h_ecomp1_clean->Reset();
	h_ecomp2_clean->Reset();
	h_ecomp_clean->Reset();
	*/
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
		h_mgg_cut[iveto]->Write();
		h_mggConstr[iveto]->Write();
		h_mggConstr_cut[iveto]->Write();
		h_mggConstr_coh[iveto]->Write();
		h_mggConstr_coh_cut[iveto]->Write();
		
		h_mm[iveto]->Write();
		h_mm_coh[iveto]->Write();
		h_mm_elas[iveto]->Write();
		h_mm_elas_coh[iveto]->Write();
		
		h_hmass[iveto]->Write();
		h_hmassCorr[iveto]->Write();
		h_hmassCorrCoh[iveto]->Write();
		h_hmass_etaCut[iveto]->Write();
		h_hmassCorr_etaCut[iveto]->Write();
		h_hmassCorrCoh_etaCut[iveto]->Write();
		
		h_pt[iveto]->Write();
		h_ptCoh[iveto]->Write();
		
		dirVeto->cd("../");
	}
	h_t_vs_theta->Write();
	
	h_mgg_vs_vertexZ->Write();
	h_mgg_vs_vertexR->Write();
	
	h_elasVSmgg_low->Write();
	h_elasVSmgg_high->Write();
	h_elasCorrVSmgg_low->Write();
	h_elasCorrVSmgg_high->Write();
	
	h_openAngle->Write();
	h_openAngleCut->Write();
	h_openAngleVsElas->Write();
	
	TDirectory *dirSC = new TDirectoryFile("SC","SC");
	dirSC->cd();
	h_scRFdt->Write();
	h_scRFdt_cut->Write();
	h_scdE->Write();
	h_scdE_singleHit->Write();
	h_scDeltaPhi->Write();
	h_scDeltaPhi_singleHit->Write();
	h_nSC->Write();
	h_nSC_nobcal->Write();
	h_nSC_onebcal->Write();
	h_nSC_bcalVeto->Write();
	h_nSC_extra->Write();
	dirSC->cd("../");
	
	TDirectory *dirBCAL = new TDirectoryFile("BCAL","BCAL");
	dirBCAL->cd();
	h_bcalRFdt->Write();
	h_bcalEnergy->Write();
	h_bcalEnergy_singleHit->Write();
	h_bcalTheta->Write();
	h_bcalTheta_singleHit->Write();
	h_bcalDeltaPhi->Write();
	h_bcalDeltaPhi_singleHit->Write();
	h_bcalThetaVSTime->Write();
	h_bcalDeltaPhiVSTime->Write();
	h_bcalThetaVSEnergy->Write();
	h_nBCAL->Write();
	h_nBCAL_1ns->Write();
	dirBCAL->cd("../");
	
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
	
	/*
	TDirectory *dirCompton = new TDirectoryFile("Compton", "Compton");
	dirCompton->cd();
	h_ecomp1->Write();
	h_ecomp2->Write();
	h_ecomp->Write();
	h_ecomp1_clean->Write();
	h_ecomp2_clean->Write();
	h_ecomp_clean->Write();
	dirCompton->cd("../");
	*/
	
	printf("  Done.\n");
	return;
}
