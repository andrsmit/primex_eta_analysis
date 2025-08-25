#include "EtaAna.h"

void EtaAna::EtaggAnalysis_matrix() {
	
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
	
	for(int ishow=0; ishow<(locNFCALShowersGood-1); ishow++) {
		
		int show1 = locGoodFCALShowers[ishow];
		TVector3 pos1 = GetFCALPosition(show1);
		
		double  t1 = m_fcalT[show1] - (pos1.Mag()/m_c);
		double  e1 = m_fcalE[show1];
		
		double px1 = e1*pos1.X() / pos1.Mag();
		double py1 = e1*pos1.Y() / pos1.Mag();
		double pz1 = e1*pos1.Z() / pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		CheckTOFMatch(pos1, t1, tof_dx1, tof_dy1, tof_dt1, m_TOFRFCut);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		for(int jshow=(ishow+1); jshow<locNFCALShowersGood; jshow++) {
			
			int show2 = locGoodFCALShowers[jshow];
			TVector3 pos2 = GetFCALPosition(show2);
			
			double t2 = m_fcalT[show2] - (pos2.Mag()/m_c);
			double e2 = m_fcalE[show2];
			
			double px2 = e2*pos2.X() / pos2.Mag();
			double py2 = e2*pos2.Y() / pos2.Mag();
			double pz2 = e2*pos2.Z() / pos2.Mag();
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			CheckTOFMatch(pos2, t2, tof_dx2, tof_dy2, tof_dt2, m_TOFRFCut);
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
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				double etaEnergyCoh = GetEnergyAfterRecoil(eb, prodTheta, ParticleMass(Eta), ParticleMass(m_Target));
				
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
				double constraintEnergy = m_IsCohMC ? etaEnergyCoh : etaEnergy;
				
				double e1c = e1/(1.+sigr) + (constraintEnergy-e2)/(1.+(1./sigr));
				double e2c = constraintEnergy - e1c;
				double invmassConstr = sqrt(2.*e1c*e2c*(1.-cos12));
				
				//-----------------------------------------------------//
				
				if(m_FillThrown && IsEtaCut(invmassConstr) && isElastic) {
					FillAngularMatrix(locThrownBeamEnergy, locThrownAngle, prodTheta, fillWeight);
				}
				
				FillInvmassMatrix(prodTheta, invmassConstr, eb, fillWeight);
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

void EtaAna::InitializeMatrixHists()
{
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	h_invmassMatrix = new TH3F("invmassMatrix", 
		"Invariant Mass Matrix; #theta(rec) [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]; E_{#gamma} [GeV]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin, 
		nBeamEnergyBins, m_minBeamEnergyBin, m_maxBeamEnergyBin);
	h_invmassMatrix->Sumw2();
	h_invmassMatrix->SetDirectory(0);
	
	h_invmassMatrix_prompt = new TH3F("invmassMatrix_prompt",
		"Invariant Mass matrix; #theta(rec) [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]; E_{#gamma} [GeV]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin, 
		nBeamEnergyBins, m_minBeamEnergyBin, m_maxBeamEnergyBin);
	h_invmassMatrix_prompt->SetDirectory(0);
	h_invmassMatrix_prompt->Sumw2();
	
	h_invmassMatrix_acc = new TH3F("invmassMatrix_acc", 
		"Invariant Mass Matrix; #theta(rec) [#circ]; m_{#gamma#gamma}^{Constr} [GeV/c^{2}]; E_{#gamma} [GeV]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 
		nInvmassBins, m_minInvmassBin, m_maxInvmassBin, 
		nBeamEnergyBins, m_minBeamEnergyBin, m_maxBeamEnergyBin);
	h_invmassMatrix_acc->SetDirectory(0);
	h_invmassMatrix_acc->Sumw2();
	
	if(m_FillThrown) {
		h_AngularMatrix = new TH3F("AngularMatrix", "; #theta(thrown) [#circ]; #theta(rec) [#circ]",
			nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
			nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
			nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
	}
	return;
}

void EtaAna::ResetMatrixHists()
{
	if(h_invmassMatrix!=nullptr) h_invmassMatrix->Reset();
	if(h_invmassMatrix_prompt!=nullptr) h_invmassMatrix_prompt->Reset();
	if(h_invmassMatrix_acc!=nullptr) h_invmassMatrix_acc->Reset();
	if(h_AngularMatrix!=nullptr) h_AngularMatrix->Reset();
	
	return;
}

void EtaAna::WriteMatrixHists()
{
	printf("\n  Writing angular and invariant mass matrices...\n");
	
	if(h_invmassMatrix!=nullptr) h_invmassMatrix->Write();
	if(h_invmassMatrix_prompt!=nullptr) h_invmassMatrix_prompt->Write();
	if(h_invmassMatrix_acc!=nullptr) h_invmassMatrix_acc->Write();
	if(h_AngularMatrix!=nullptr) h_AngularMatrix->Write();
	
	printf("  Done.\n");
	return;
}
