#include "EtaAna.h"

void EtaAna::EtaggAnalysis_special() {
	
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
	int locNFCALShowersTotal, locNFCALShowersEnergyCut, locNFCALShowersGood;
	locNFCALShowersTotal = GetFCALShowerList(locGoodFCALShowers, locNFCALShowersEnergyCut, 
		m_FCALEnergyCut, m_FCALExtraEnergyCut, 2.0, m_FCALRFCut);
	locNFCALShowersGood  = (int)locGoodFCALShowers.size();
	
	// Apply multiplicity cut on the number of FCAL showers: 
	int isExclusive = 0;
	if((locNFCALShowersEnergyCut!=2) || (locNFCALShowersGood!=2)) isExclusive = 1;
	
	//-------------------------------------------//
	// Get list of 'good' BCAL showers:
	
	vector<int> locGoodBCALShowers; locGoodBCALShowers.clear();
	int locNBCALShowers = GetBCALShowerList(locGoodBCALShowers, 0.0, m_BCALRFCut);
	
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
	// Get list of 'good' SC hits:
	
	vector<int> locGoodSCHits; locGoodSCHits.clear();
	int locNSCHits = GetSCHitList(locGoodSCHits);
	
	for(int ishow=0; ishow<(m_nfcal-1); ishow++) {
		
		TVector3 pos1 = GetFCALPosition(ishow);
		
		double  t1 = m_fcalT[ishow] - (pos1.Mag()/m_c);
		double  e1 = m_fcalE[ishow];
		
		double px1 = e1*pos1.X() / pos1.Mag();
		double py1 = e1*pos1.Y() / pos1.Mag();
		double pz1 = e1*pos1.Z() / pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		CheckTOFMatch(pos1, t1, tof_dx1, tof_dy1, tof_dt1, m_TOFRFCut);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		for(int jshow=ishow+1; jshow<m_nfcal; jshow++) {
			
			TVector3 pos2 = GetFCALPosition(jshow);
			
			double  t2 = m_fcalT[jshow] - (pos2.Mag()/m_c);
			double  e2 = m_fcalE[jshow];
			
			double px2 = e2*pos2.X() / pos2.Mag();
			double py2 = e2*pos2.Y() / pos2.Mag();
			double pz2 = e2*pos2.Z() / pos2.Mag();
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			CheckTOFMatch(pos2, t2, tof_dx2, tof_dy2, tof_dt2, m_TOFRFCut);
			double tof_dr2 = sqrt(pow(tof_dx2,2.0)+pow(tof_dy2,2.0));
			
			//-----------------------------------------------------//
			// TOF Veto
			
			int tof_match1 = tof_dr1 < m_FCALTOFCut ? 1.0 : 0.0;
			int tof_match2 = tof_dr2 < m_FCALTOFCut ? 1.0 : 0.0;
			
			// reject combinations of FCAL showers where both showers are near a TOF hit:
			bool isTOFVeto = false;
			if((tof_dr1 < m_FCALTOFCut) && (tof_dr2 < m_FCALTOFCut)) isTOFVeto = true;
			
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
			
			int isHadronicVeto = IsHadronicVeto(m_vetoOption, locNBCALShowers, locNBCALShowers_1ns, locNSCHits, locNSCHits_coplanar, 
				prodPhi, locBCALPhi, locBCALTheta, locBCALRFDT);
			
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
				
				h_mgg_Special[0]->Fill(prodTheta, invmass, fillWeight);
				h_mggConstr_Special[0]->Fill(prodTheta, invmassConstr, fillWeight);
				
				if(e1>0.25 && e2>0.25) {
					h_mgg_Special[1]->Fill(prodTheta, invmass, fillWeight);
					h_mggConstr_Special[1]->Fill(prodTheta, invmassConstr, fillWeight);
					
					if(!FCALFiducialCut(pos1, 2.0) && !FCALFiducialCut(pos2,2.0)) {
						h_mgg_Special[2]->Fill(prodTheta, invmass, fillWeight);
						h_mggConstr_Special[2]->Fill(prodTheta, invmassConstr, fillWeight);
						
						if(locNFCALShowersGood==2) {
							h_mgg_Special[3]->Fill(prodTheta, invmass, fillWeight);
							h_mggConstr_Special[3]->Fill(prodTheta, invmassConstr, fillWeight);
							
							if(locNFCALShowersEnergyCut==2) {
								h_mgg_Special[4]->Fill(prodTheta, invmass, fillWeight);
								h_mggConstr_Special[4]->Fill(prodTheta, invmassConstr, fillWeight);
								
								if(m_nfcal==2) {
									h_mgg_Special[5]->Fill(prodTheta, invmass, fillWeight);
									h_mggConstr_Special[5]->Fill(prodTheta, invmassConstr, fillWeight);
								}
								
								if(!isTOFVeto) {
									h_mgg_Special[6]->Fill(prodTheta, invmass, fillWeight);
									h_mggConstr_Special[6]->Fill(prodTheta, invmassConstr, fillWeight);
									
									if(!isHadronicVeto) {
										h_mgg_Special[7]->Fill(prodTheta, invmass, fillWeight);
										h_mggConstr_Special[7]->Fill(prodTheta, invmassConstr, fillWeight);
										
										if(isElastic) {
											h_mgg_Special[8]->Fill(prodTheta, invmass, fillWeight);
											h_mggConstr_Special[8]->Fill(prodTheta, invmassConstr, fillWeight);
										}
									}
								}
							}
						}
					}
				}
				
			} // end loop over beam photons
		} // end inner loop over fcal showers
	} // end outer loop over fcal showers
}

void EtaAna::InitializeSpecialHists() 
{
	for(int ihist=0; ihist<m_nSpecialHists; ihist++) {
		h_mgg_Special[ihist] = new TH2F(Form("mgg_special_%02d",ihist), 
			"; #theta_{#gamma#gamma} [#circ]; m_{#gamma#gamma} [GeV/c^{2}]",
			500, 0.0, 5.0, 450, 0.3, 1.2);
		h_mggConstr_Special[ihist] = new TH2F(Form("mggConstr_special_%02d",ihist), 
			"; #theta_{#gamma#gamma} [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]",
			500, 0.0, 5.0, 450, 0.3, 1.2);
	}
}

void EtaAna::ResetSpecialHists()
{
	for(int ihist=0; ihist<m_nSpecialHists; ihist++) {
		h_mgg_Special[ihist]->Reset();
		h_mggConstr_Special[ihist]->Reset();
	}
}

void EtaAna::WriteSpecialHists()
{
	for(int ihist=0; ihist<m_nSpecialHists; ihist++) {
		h_mgg_Special[ihist]->Write();
		h_mggConstr_Special[ihist]->Write();
	}
}
