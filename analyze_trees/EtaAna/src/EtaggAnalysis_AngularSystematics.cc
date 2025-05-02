#include "EtaAna.h"

void EtaAna::EtaggAnalysis_AngularSystematics() {
	
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
	
	//=====================================================================================//
	
	// Strict Veto on BCAL:
	if(locNBCALShowers>0) return;
	
	//=====================================================================================//
	
	for(int ishow=0; ishow<(locNFCALShowersGood-1); ishow++) {
		
		int show1 = locGoodFCALShowers[ishow];
		TVector3 pos1 = GetFCALPosition(show1);
		
		double e1 = m_fcalE[show1];
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		CheckTOFMatch(pos1, tof_dx1, tof_dy1, tof_dt1, m_TOFRFCut);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		for(int jshow=(ishow+1); jshow<locNFCALShowersGood; jshow++) {
			
			int show2 = locGoodFCALShowers[jshow];
			TVector3 pos2 = GetFCALPosition(show2);
			
			double e2 = m_fcalE[show2];
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			CheckTOFMatch(pos2, tof_dx2, tof_dy2, tof_dt2, m_TOFRFCut);
			double tof_dr2 = sqrt(pow(tof_dx2,2.0)+pow(tof_dy2,2.0));
			
			//-----------------------------------------------------//
			// TOF Veto
			
			if((tof_dr1 < m_FCALTOFCut) && (tof_dr2 < m_FCALTOFCut)) continue;
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(int igam=0; igam<locNBeamPhotons; igam++) {
				
				int ibeam = locGoodBeamPhotons[igam].first;
				double fillWeight = locGoodBeamPhotons[igam].second;
				
				double eb    = m_beamE[ibeam];
				double brfdt = m_beamT[ibeam] - m_rfTime;
				
				//-----------------------------------------------------//
				// Two-Photon kinematics:
				
				vector<double> locProdTheta, locMgg, locMggConst;
				
				double theta1 = pos1.Theta();
				double phi1   = pos1.Phi();
				double theta2 = pos2.Theta();
				double phi2   = pos2.Phi();
				
				int nAngularShifts = (int)m_AngularShifts.size();
				for(int ishift=0; ishift<nAngularShifts; ishift++)
				{
					double locTheta1 = theta1 * (1.0 + m_AngularShifts[ishift]);
					double locTheta2 = theta2 * (1.0 + m_AngularShifts[ishift]);
					
					TVector3 locPos1 = pos1;
					locPos1.SetTheta(locTheta1);
					
					TVector3 locPos2 = pos2;
					locPos2.SetTheta(locTheta2);
					
					double px1 = e1*locPos1.X() / locPos1.Mag();
					double py1 = e1*locPos1.Y() / locPos1.Mag();
					double pz1 = e1*locPos1.Z() / locPos1.Mag();
					
					double px2 = e2*locPos2.X() / locPos2.Mag();
					double py2 = e2*locPos2.Y() / locPos2.Mag();
					double pz2 = e2*locPos2.Z() / locPos2.Mag();
					
					double Egg  =  e1 +  e2; // energy of 2-photon pair
					double pggx = px1 + px2; // momentum along x-axis
					double pggy = py1 + py2; // momentum along y-axis
					double pggz = pz1 + pz2; // momentum along z-axis
					
					// transverse momentum:
					double pggt = sqrt(pow(pggx,2.0) + pow(pggy,2.0));
					
					// polar angle:
					double prodTheta = atan2(pggt,pggz) * TMath::RadToDeg();
					
					// opening angle:
					double cos12   = (locPos1.X()*locPos2.X() + locPos1.Y()*locPos2.Y() + locPos1.Z()*locPos2.Z()) 
						/ (locPos1.Mag()*locPos2.Mag());
					
					// invariant mass:
					double invmass = sqrt(2.0*e1*e2*(1.-cos12));
					
					// Calculate the energy of the eta meson, assuming a coherent production process:
					double etaEnergyCoh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
					
					// Calculate the energy of the eta meson, assuming production on a free nucleon:
					double etaEnergy = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(Proton));
					
					// Apply a cut on the elasticity
					//  (ratio of measured energy of 2-photons, to the calculated energy above):
					bool isElastic = IsElasticCut(Egg, etaEnergy, prodTheta);
					
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
					
					h_mgg_AngularShift[ishift]->Fill(prodTheta, invmass, fillWeight);
					h_mggConstr_AngularShift[ishift]->Fill(prodTheta, invmassConstr, fillWeight);
					
					if(m_FillThrown && IsEtaCut(invmassConstr) && isElastic) {
						h_AngularMatrix_AngularShift[ishift]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
					}
				}
				
				int nAngularSmears = (int)m_AngularSmears.size();
				for(int ismear=0; ismear<nAngularSmears; ismear++)
				{
					TVector3 locPos1 = pos1;
					locPos1.SetTheta(m_random->Gaus(pos1.Theta(), m_AngularSmears[ismear]*pos1.Theta()));
					
					TVector3 locPos2 = pos2;
					locPos2.SetTheta(m_random->Gaus(pos2.Theta(), m_AngularSmears[ismear]*pos2.Theta()));
					
					double px1 = e1*locPos1.X() / locPos1.Mag();
					double py1 = e1*locPos1.Y() / locPos1.Mag();
					double pz1 = e1*locPos1.Z() / locPos1.Mag();
					
					double px2 = e2*locPos2.X() / locPos2.Mag();
					double py2 = e2*locPos2.Y() / locPos2.Mag();
					double pz2 = e2*locPos2.Z() / locPos2.Mag();
					
					double Egg  =  e1 +  e2; // energy of 2-photon pair
					double pggx = px1 + px2; // momentum along x-axis
					double pggy = py1 + py2; // momentum along y-axis
					double pggz = pz1 + pz2; // momentum along z-axis
					
					// transverse momentum:
					double pggt = sqrt(pow(pggx,2.0) + pow(pggy,2.0));
					
					// polar angle:
					double prodTheta = atan2(pggt,pggz) * TMath::RadToDeg();
					
					// opening angle:
					double cos12   = (locPos1.X()*locPos2.X() + locPos1.Y()*locPos2.Y() + locPos1.Z()*locPos2.Z()) / (locPos1.Mag()*locPos2.Mag());
					
					// invariant mass:
					double invmass = sqrt(2.0*e1*e2*(1.-cos12));
					
					// Calculate the energy of the eta meson, assuming a coherent production process:
					double etaEnergyCoh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
					
					// Calculate the energy of the eta meson, assuming production on a free nucleon:
					double etaEnergy = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(Proton));
					
					// Apply a cut on the elasticity
					//  (ratio of measured energy of 2-photons, to the calculated energy above):
					bool isElastic = IsElasticCut(Egg, etaEnergy, prodTheta);
					
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
					
					h_mgg_AngularSmear[ismear]->Fill(prodTheta, invmass, fillWeight);
					h_mggConstr_AngularSmear[ismear]->Fill(prodTheta, invmassConstr, fillWeight);
					
					if(m_FillThrown && IsEtaCut(invmassConstr) && isElastic) {
						h_AngularMatrix_AngularSmear[ismear]->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
					}
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

void EtaAna::InitializeAngularHists() 
{
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	//----------------------------------------------//
	
	m_AngularShifts.clear();
	m_AngularShifts.push_back(-0.025);
	m_AngularShifts.push_back(-0.020);
	m_AngularShifts.push_back(-0.015);
	m_AngularShifts.push_back(-0.010);
	m_AngularShifts.push_back(-0.005);
	m_AngularShifts.push_back( 0.000);
	m_AngularShifts.push_back( 0.005);
	m_AngularShifts.push_back( 0.010);
	m_AngularShifts.push_back( 0.015);
	m_AngularShifts.push_back( 0.020);
	m_AngularShifts.push_back( 0.025);
	
	int nAngularShifts = (int)m_AngularShifts.size();
	for(int ishift=0; ishift<nAngularShifts; ishift++) {
		
		TH2F *loc_h_mgg = new TH2F(Form("mgg_angularShift_%02d",ishift), 
			Form("#theta_{#gamma} shifted by %.3f%; #theta_{#gamma#gamma} [#circ]; m_{#gamma#gamma} [GeV/c^{2}]", 
				m_AngularShifts[ishift]),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		
		TH2F *loc_h_mggConstr = new TH2F(Form("mgg_const_angularShift_%02d",ishift), 
			Form("#theta_{#gamma} shifted by %.3f%; #theta_{#gamma#gamma} [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]", 
				m_AngularShifts[ishift]),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		
		h_mgg_AngularShift.push_back(loc_h_mgg);
		h_mggConstr_AngularShift.push_back(loc_h_mggConstr);
		
		if(m_FillThrown) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_AngularShift_%02d",ishift), 
				Form("#theta_{#gamma} shifted by %.3f%; #theta(thrown) [#circ]; #theta(rec) [#circ]; E_{#gamma} [GeV]", 
					m_AngularShifts[ishift]),
				nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
				nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
				nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
			hMatrix->SetDirectory(0);
			h_AngularMatrix_AngularShift.push_back(hMatrix);
		}
	}
	
	//----------------------------------------------//
	
	m_AngularSmears.clear();
	m_AngularSmears.push_back( 0.000);
	m_AngularSmears.push_back( 0.005);
	m_AngularSmears.push_back( 0.010);
	m_AngularSmears.push_back( 0.015);
	m_AngularSmears.push_back( 0.020);
	m_AngularSmears.push_back( 0.025);
	m_AngularSmears.push_back( 0.030);
	m_AngularSmears.push_back( 0.035);
	m_AngularSmears.push_back( 0.040);
	m_AngularSmears.push_back( 0.045);
	m_AngularSmears.push_back( 0.050);
	
	int nAngularSmears = (int)m_AngularSmears.size();
	for(int ismear=0; ismear<nAngularSmears; ismear++) {
		
		TH2F *loc_h_mgg = new TH2F(Form("mgg_angularSmear_%02d",ismear), 
			Form("#theta_{#gamma} smeared by %.3f%; #theta_{#gamma#gamma} [#circ]; m_{#gamma#gamma} [GeV/c^{2}]", 
				m_AngularSmears[ismear]),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		
		TH2F *loc_h_mggConstr = new TH2F(Form("mgg_const_angularSmear_%02d",ismear), 
			Form("#theta_{#gamma} smeared by %.3f%; #theta_{#gamma#gamma} [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]", 
				m_AngularSmears[ismear]),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin,
			nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		
		h_mgg_AngularSmear.push_back(loc_h_mgg);
		h_mggConstr_AngularSmear.push_back(loc_h_mggConstr);
		
		if(m_FillThrown) {
			TH3F *hMatrix = new TH3F(Form("AngularMatrix_AngularSmear_%02d",ismear), 
				Form("#theta_{#gamma} smeared by %.3f%; #theta(thrown) [#circ]; #theta(rec) [#circ]; E_{#gamma} [GeV]", 
					m_AngularSmears[ismear]),
				nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
				nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
				nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
			hMatrix->SetDirectory(0);
			h_AngularMatrix_AngularSmear.push_back(hMatrix);
		}
	}
	
	return;
}

void EtaAna::ResetAngularHists() {
	
	int nAngularShifts = (int)m_AngularShifts.size();
	for(int ishift=0; ishift<nAngularShifts; ishift++) {
		h_mgg_AngularShift[ishift]->Reset();
		h_mggConstr_AngularShift[ishift]->Reset();
		if(m_FillThrown) h_AngularMatrix_AngularShift[ishift]->Reset();
	}
	int nAngularSmears = (int)m_AngularSmears.size();
	for(int ismear=0; ismear<nAngularSmears; ismear++) {
		h_mgg_AngularSmear[ismear]->Reset();
		h_mggConstr_AngularSmear[ismear]->Reset();
		if(m_FillThrown) h_AngularMatrix_AngularSmear[ismear]->Reset();
	}
	return;
}

void EtaAna::WriteAngularHists() {
	
	int nAngularShifts = (int)m_AngularShifts.size();
	for(int ishift=0; ishift<nAngularShifts; ishift++) {
		h_mgg_AngularShift[ishift]->Write();
		h_mggConstr_AngularShift[ishift]->Write();
		if(m_FillThrown) h_AngularMatrix_AngularShift[ishift]->Write();
	}
	int nAngularSmears = (int)m_AngularSmears.size();
	for(int ismear=0; ismear<nAngularSmears; ismear++) {
		h_mgg_AngularSmear[ismear]->Write();
		h_mggConstr_AngularSmear[ismear]->Write();
		if(m_FillThrown) h_AngularMatrix_AngularSmear[ismear]->Write();
	}
	return;
}

