#include "EtaAna.h"

void EtaAna::Omega3gAnalysis() {
	
	if(CheckEventMultiplicities()) {
		printf("    Skipping event %d\n", m_event);
		return;
	}
	
	vector<pair<int,double>> locGoodBeamPhotons; locGoodBeamPhotons.clear();
	int locNBeamPhotons = GetBeamPhotonList(locGoodBeamPhotons, m_minBeamEnergyCut, m_maxBeamEnergyCut);
	
	vector<int> locGoodBCALShowers; locGoodBCALShowers.clear();
	int locNBCALShowers = GetBCALShowerList(locGoodBCALShowers, 0.0, m_BCALRFCut);
	
	//=====================================================================================//
	
	// Look for events where there are either 3 photons in FCAL or 2 photons in FCAL + 1 photon in BCAL:
	
	vector<int> locFCALShowers; locFCALShowers.clear();
	int locNFCALShowers     = 0;
	int locNGoodFCALShowers = 0;
	for(int ishow=0; ishow<m_nfcal; ishow++) {
		
		TVector3 locPos = GetFCALPosition(ishow);
		double locT = m_fcalT[ishow] - (locPos.Mag()/m_c) - m_rfTime;
		
		if(fabs(locT) < m_FCALRFCut) {
			locNFCALShowers++;
			locFCALShowers.push_back(ishow);
			if((m_fcalE[ishow] > m_FCALEnergyCut) && !FCALFiducialCut(locPos, 2.0)) {
				locNGoodFCALShowers++;
			}
		}
	}
	
	bool isThreeGammaEvent = false;
	if(locNFCALShowers==3 && locNBCALShowers==0) isThreeGammaEvent = true;
	if(!isThreeGammaEvent) return;
	
	// sort our shower list by energy:
	vector<int> sortedShowerList = {-1, -1, -1};
	
	double e1 = m_fcalE[locFCALShowers[0]];
	double e2 = m_fcalE[locFCALShowers[1]];
	double e3 = m_fcalE[locFCALShowers[2]];
	
	if((e1>e2) && (e1>e3)) {
		sortedShowerList[0] = locFCALShowers[0];
		if(e2 > e3) {
			sortedShowerList[1] = locFCALShowers[1];
			sortedShowerList[2] = locFCALShowers[2];
		} else {
			sortedShowerList[1] = locFCALShowers[2];
			sortedShowerList[2] = locFCALShowers[1];
		}
	} else if((e2>e1) && (e2>e3)) {
		sortedShowerList[0] = locFCALShowers[1];
		if(e1 > e3) {
			sortedShowerList[1] = locFCALShowers[0];
			sortedShowerList[2] = locFCALShowers[2];
		} else {
			sortedShowerList[1] = locFCALShowers[2];
			sortedShowerList[2] = locFCALShowers[0];
		}
	} else {
		sortedShowerList[0] = locFCALShowers[2];
		if(e1 > e2) {
			sortedShowerList[1] = locFCALShowers[0];
			sortedShowerList[2] = locFCALShowers[1];
		} else {
			sortedShowerList[1] = locFCALShowers[1];
			sortedShowerList[2] = locFCALShowers[0];
		}
	}
	
	e1 = m_fcalE[sortedShowerList[0]];
	TVector3 pos1 = GetFCALPosition(sortedShowerList[0]);
	
	e2 = m_fcalE[sortedShowerList[1]];
	TVector3 pos2 = GetFCALPosition(sortedShowerList[1]);
	
	e3 = m_fcalE[sortedShowerList[2]];
	TVector3 pos3 = GetFCALPosition(sortedShowerList[2]);
	
	double t1 = m_fcalT[sortedShowerList[0]] - (pos1.Mag()/m_c);
	double t2 = m_fcalT[sortedShowerList[1]] - (pos2.Mag()/m_c);
	double t3 = m_fcalT[sortedShowerList[2]] - (pos3.Mag()/m_c);
	
	//-------------------------------------------------------//
	// apply Fiducial cut:
	
	if(FCALFiducialCut(pos1, 2.0) || FCALFiducialCut(pos2,2.0) || FCALFiducialCut(pos3,2.0)) return;
	
	//-------------------------------------------------------//
	// apply TOF veto:
	
	double locTOFdx1, locTOFdy1, locTOFdt1;
	double locTOFdx2, locTOFdy2, locTOFdt2;
	double locTOFdx3, locTOFdy3, locTOFdt3;
	
	CheckTOFMatch(pos1, t1, locTOFdx1, locTOFdy1, locTOFdt1, m_TOFRFCut);
	CheckTOFMatch(pos2, t2, locTOFdx2, locTOFdy2, locTOFdt2, m_TOFRFCut);
	CheckTOFMatch(pos3, t3, locTOFdx3, locTOFdy3, locTOFdt3, m_TOFRFCut);
	
	double locTOFdr1 = sqrt(pow(locTOFdx1,2.0)+pow(locTOFdy1,2.0));
	double locTOFdr2 = sqrt(pow(locTOFdx2,2.0)+pow(locTOFdy2,2.0));
	double locTOFdr3 = sqrt(pow(locTOFdx3,2.0)+pow(locTOFdy3,2.0));
	
	int isTOFVeto = false;
	if((locTOFdr1 < m_FCALTOFCut) || (locTOFdr2 < m_FCALTOFCut) || (locTOFdr3 < m_FCALTOFCut)) isTOFVeto = true;
	if(isTOFVeto) return;
	
	//-------------------------------------------------------//
	
	double m12 = CalcInvmass(e1, e2, pos1, pos2);
	double m13 = CalcInvmass(e1, e3, pos1, pos3);
	double m23 = CalcInvmass(e2, e3, pos2, pos3);
	double m3g = CalcInvmass(e1, e2, e3, pos1, pos2, pos3);
	
	h_3gamma_m12->Fill(m12);
	h_3gamma_m13->Fill(m13);
	h_3gamma_m23->Fill(m23);
	h_3gamma_m3g->Fill(m3g);
	
	// require at least one of the pairs to be a pion:
	
	int nPions = 0;
	if((0.11<m12) && (m12<0.16)) nPions++;
	if((0.11<m13) && (m13<0.16)) nPions++;
	if((0.11<m23) && (m23<0.16)) nPions++;
	//if(nPions!=1) return;
	
	// try to find z position that makes this 3-gamma consistent with an omega:
	
	double zVertex = 100.0;
	int nIterations = 0;
	while(1) {
		
		// calculate new 3-gamma mass:
		double epsilon    = 0.1;
		double cost       = CalculateCost(zVertex,         e1, e2, e3, pos1, pos2, pos3);
		double cost_plus  = CalculateCost(zVertex+epsilon, e1, e2, e3, pos1, pos2, pos3);
		double cost_minus = CalculateCost(zVertex-epsilon, e1, e2, e3, pos1, pos2, pos3);
		double gradient   = (cost_plus - cost_minus) / (2.0*epsilon);
		//printf("  iteration %d, gradient = %f\n", nIterations, gradient);
		zVertex = zVertex - 1.e5*gradient;
		
		if(cost<1.e-6) break;
		if(nIterations>100) break;
		nIterations++;
	}
	
	double m12_new = CalcInvmassShifted(zVertex, e1, e2, pos1, pos2);
	double m13_new = CalcInvmassShifted(zVertex, e1, e3, pos1, pos3);
	double m23_new = CalcInvmassShifted(zVertex, e2, e3, pos2, pos3);
	double m3g_new = CalcInvmassShifted(zVertex, e1, e2, e3, pos1, pos2, pos3);
	
	int nPions_new = 0;
	if((0.11<m12_new) && (m12_new<0.16)) nPions_new++;
	if((0.11<m13_new) && (m13_new<0.16)) nPions_new++;
	if((0.11<m23_new) && (m23_new<0.16)) nPions_new++;
	if(nPions_new==1) {
		h_3gamma_vz->Fill(65.0+zVertex);
	}
	
	double fdcz1 = 182.3-65.0;
	double fdcz2 = 240.8-65.0;
	double fdcz3 = 299.4-65.0;
	
	double prodTheta      = CalcProdTheta(e1, e2, e3, pos1, pos2, pos3);
	double prodTheta_fdc1 = CalcProdThetaShifted(fdcz1, e1, e2, e3, pos1, pos2, pos3);
	double prodTheta_fdc2 = CalcProdThetaShifted(fdcz2, e1, e2, e3, pos1, pos2, pos3);
	double prodTheta_fdc3 = CalcProdThetaShifted(fdcz3, e1, e2, e3, pos1, pos2, pos3);
	
	bool isTarget = false;
	bool isFDC1   = false;
	bool isFDC2   = false;
	bool isFDC3   = false;
	if((-15.0<zVertex) && (zVertex<15.0)) isTarget = true;
	else if(((fdcz1-20.0)<zVertex) && (zVertex<(fdcz1+20.0))) isFDC1 = true;
	else if(((fdcz2-20.0)<zVertex) && (zVertex<(fdcz2+20.0))) isFDC2 = true;
	else if(((fdcz3-20.0)<zVertex) && (zVertex<(fdcz3+20.0))) isFDC3 = true;
	
	/*
	printf("Took %d iterations to find zVertex = %f\n", nIterations, zVertex);
	printf("  original mass: %f\n", m3g);
	printf("  new mass: %f\n", m3g_shifted);
	printf("  shifted invariant mass: %f\n", CalcInvmassShifted(zVertex, e1, e2, e3, pos1, pos2, pos3));
	*/
	
	
	
	//-----------------------------------------------------//
	// Loop over Beam photons
	
	for(int igam=0; igam<locNBeamPhotons; igam++) {
		
		int ibeam = locGoodBeamPhotons[igam].first;
		double fillWeight = locGoodBeamPhotons[igam].second;
		
		double eb    = m_beamE[ibeam];
		double brfdt = m_beamT[ibeam] - m_rfTime;
		
		// Apply a cut on the elasticity
		//  (ratio of measured energy of 2-photons, to the calculated energy above):
		bool isElastic = IsElasticCut((e1+e2+e3), eb, 0.0);
		if(!isElastic) continue;
		
		h_3gamma_m12_elas->Fill(m12, fillWeight);
		h_3gamma_m13_elas->Fill(m13, fillWeight);
		h_3gamma_m23_elas->Fill(m23, fillWeight);
		if(nPions_new==1) {
			
			h_3gamma_m3g_elas->Fill(m3g, fillWeight);
			h_3gamma_vz_elas->Fill(65.0+zVertex, fillWeight);
			
			if(isTarget) {
				h_3gamma_theta_targ->Fill(prodTheta, fillWeight);
				h_xy_targ->Fill(pos1.X(), pos1.Y(),  fillWeight);
				h_xy_targ->Fill(pos2.X(), pos2.Y(),  fillWeight);
				h_xy_targ->Fill(pos3.X(), pos3.Y(),  fillWeight);
			} else if(isFDC1) {
				h_3gamma_theta_fdc1->Fill(prodTheta_fdc1, fillWeight);
				h_xy_fdc1->Fill(pos1.X(), pos1.Y(),  fillWeight);
				h_xy_fdc1->Fill(pos2.X(), pos2.Y(),  fillWeight);
				h_xy_fdc1->Fill(pos3.X(), pos3.Y(),  fillWeight);
			} else if(isFDC2) {
				h_3gamma_theta_fdc2->Fill(prodTheta_fdc2, fillWeight);
				h_xy_fdc2->Fill(pos1.X(), pos1.Y(),  fillWeight);
				h_xy_fdc2->Fill(pos2.X(), pos2.Y(),  fillWeight);
				h_xy_fdc2->Fill(pos3.X(), pos3.Y(),  fillWeight);
			} else if(isFDC3) {
				h_3gamma_theta_fdc3->Fill(prodTheta_fdc3, fillWeight);
				h_xy_fdc3->Fill(pos1.X(), pos1.Y(),  fillWeight);
				h_xy_fdc3->Fill(pos2.X(), pos2.Y(),  fillWeight);
				h_xy_fdc3->Fill(pos3.X(), pos3.Y(),  fillWeight);
			}
		}
		
	} // end loop over beam photons
	
	return;
}

double EtaAna::CalcInvmass(double e1, double e2, TVector3 pos1, TVector3 pos2) 
{
	double cos12 = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
	double mgg   = sqrt(2.0*e1*e2*(1.0-cos12));
	return mgg;
}

double EtaAna::CalcInvmass(double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3) 
{
	double cos12 = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
	double cos13 = (pos1.X()*pos3.X() + pos1.Y()*pos3.Y() + pos1.Z()*pos3.Z()) / (pos1.Mag()*pos3.Mag());
	double cos23 = (pos2.X()*pos3.X() + pos2.Y()*pos3.Y() + pos2.Z()*pos3.Z()) / (pos2.Mag()*pos3.Mag());
	
	double m3g = sqrt(2.0*e1*e2*(1.0-cos12) + 2.0*e1*e3*(1.0-cos13) + 2.0*e2*e3*(1.0-cos23));
	return m3g;
}

double EtaAna::CalcInvmassShifted(double deltaZ, double e1, double e2, TVector3 pos1, TVector3 pos2) 
{
	TVector3 pos1_new = pos1;
	TVector3 pos2_new = pos2;
	pos1_new.SetZ(pos1.Z()-deltaZ);
	pos2_new.SetZ(pos2.Z()-deltaZ);
	double mgg = CalcInvmass(e1, e2, pos1_new, pos2_new);
	return mgg;
}

double EtaAna::CalcInvmassShifted(double deltaZ, double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3) 
{
	TVector3 pos1_new = pos1;
	TVector3 pos2_new = pos2;
	TVector3 pos3_new = pos3;
	pos1_new.SetZ(pos1.Z()-deltaZ);
	pos2_new.SetZ(pos2.Z()-deltaZ);
	pos3_new.SetZ(pos3.Z()-deltaZ);
	double m3g = CalcInvmass(e1, e2, e3, pos1_new, pos2_new, pos3_new);
	return m3g;
}

double EtaAna::CalculateCost(double deltaZ, double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3)
{
	double m3g = CalcInvmassShifted(deltaZ, e1, e2, e3, pos1, pos2, pos3);
	double cost = pow(m3g-ParticleMass(omega),2.0);
	return cost;
}

double EtaAna::CalcProdTheta(double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3)
{
	double px = (e1*pos1.X()/pos1.Mag()) + (e2*pos2.X()/pos2.Mag()) + (e3*pos3.X()/pos3.Mag());
	double py = (e1*pos1.Y()/pos1.Mag()) + (e2*pos2.Y()/pos2.Mag()) + (e3*pos3.Y()/pos3.Mag());
	double pz = (e1*pos1.Z()/pos1.Mag()) + (e2*pos2.Z()/pos2.Mag()) + (e3*pos3.Z()/pos3.Mag());
	double pt = sqrt(pow(px,2.0)+pow(py,2.0));
	double prodTheta = atan2(pt,pz) * TMath::RadToDeg();
	return prodTheta;
}

double EtaAna::CalcProdThetaShifted(double deltaZ, double e1, double e2, double e3, TVector3 pos1, TVector3 pos2, TVector3 pos3)
{
	TVector3 pos1_new = pos1;
	TVector3 pos2_new = pos2;
	TVector3 pos3_new = pos3;
	pos1_new.SetZ(pos1.Z()-deltaZ);
	pos2_new.SetZ(pos2.Z()-deltaZ);
	pos3_new.SetZ(pos3.Z()-deltaZ);
	double prodTheta = CalcProdTheta(e1, e2, e3, pos1_new, pos2_new, pos3_new);
	return prodTheta;
}

void EtaAna::InitializeOmegaHists()
{
	// Omega->pi0+gamma:
	
	h_3gamma_m12 = new TH1F("3gamma_m12", ";m_{12} [GeV/c^{2}]", 1200, 0.0, 1.2);
	h_3gamma_m13 = new TH1F("3gamma_m13", ";m_{13} [GeV/c^{2}]", 1200, 0.0, 1.2);
	h_3gamma_m23 = new TH1F("3gamma_m23", ";m_{23} [GeV/c^{2}]", 1200, 0.0, 1.2);
	h_3gamma_m3g = new TH1F("3gamma_m3g", ";m_{3#gamma} [GeV/c^{2}]", 1200, 0.0, 1.2);
	
	h_3gamma_m12_elas = new TH1F("3gamma_m12_elas", ";m_{12} [GeV/c^{2}]", 1200, 0.0, 1.2);
	h_3gamma_m13_elas = new TH1F("3gamma_m13_elas", ";m_{13} [GeV/c^{2}]", 1200, 0.0, 1.2);
	h_3gamma_m23_elas = new TH1F("3gamma_m23_elas", ";m_{23} [GeV/c^{2}]", 1200, 0.0, 1.2);
	h_3gamma_m3g_elas = new TH1F("3gamma_m3g_elas", ";m_{3#gamma} [GeV/c^{2}]", 1200, 0.0, 1.2);
	
	h_3gamma_vz      = new TH1F("3gamma_vz", "Calculated Vertex Z position for #omega#rightarrow#pi^{0}#gamma; z_{vertex} [cm]",
		1000, -500.0, 500.0);
	h_3gamma_vz_elas = new TH1F("3gamma_vz_elas", "Calculated Vertex Z position for #omega#rightarrow#pi^{0}#gamma; z_{vertex} [cm]",
		1000, -500.0, 500.0);
	
	h_3gamma_theta_targ = new TH1F("3gamma_theta_targ", 
		"Angular Distribution of #omega's produced near target; #theta_{#omega} [#circ]",
		1000, 0.0, 6.0);
	h_3gamma_theta_fdc1 = new TH1F("3gamma_theta_fdc1", 
		"Angular Distribution of #omega's produced near FDC Package 1; #theta_{#omega} [#circ]",
		1000, 0.0, 6.0);
	h_3gamma_theta_fdc2 = new TH1F("3gamma_theta_fdc2", 
		"Angular Distribution of #omega's produced near FDC Package 2; #theta_{#omega} [#circ]",
		1000, 0.0, 6.0);
	h_3gamma_theta_fdc3 = new TH1F("3gamma_theta_fdc3", 
		"Angular Distribution of #omega's produced near FDC Package 3; #theta_{#omega} [#circ]",
		1000, 0.0, 6.0);
	
	h_xy_targ = new TH2F("xy_targ", "#omega's from Target Area; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		600, -150.0, 150.0, 600, -150.0, 150.0);
	h_xy_fdc1 = new TH2F("xy_fdc1", "#omega's from 1st FDC Package; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		600, -150.0, 150.0, 600, -150.0, 150.0);
	h_xy_fdc2 = new TH2F("xy_fdc2", "#omega's from 2nd FDC Package; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		600, -150.0, 150.0, 600, -150.0, 150.0);
	h_xy_fdc3 = new TH2F("xy_fdc3", "#omega's from 3rd FDC Package; x_{FCAL} [cm]; y_{FCAL} [cm]", 
		600, -150.0, 150.0, 600, -150.0, 150.0);
	
	return;
}

void EtaAna::ResetOmegaHists()
{
	h_3gamma_m12->Reset();
	h_3gamma_m13->Reset();
	h_3gamma_m23->Reset();
	h_3gamma_m3g->Reset();
	h_3gamma_m12_elas->Reset();
	h_3gamma_m13_elas->Reset();
	h_3gamma_m23_elas->Reset();
	h_3gamma_m3g_elas->Reset();
	h_3gamma_vz->Reset();
	h_3gamma_vz_elas->Reset();
	
	h_3gamma_theta_targ->Reset();
	h_3gamma_theta_fdc1->Reset();
	h_3gamma_theta_fdc2->Reset();
	h_3gamma_theta_fdc3->Reset();
	h_xy_targ->Reset();
	h_xy_fdc1->Reset();
	h_xy_fdc2->Reset();
	h_xy_fdc3->Reset();
	
	return;
}

void EtaAna::WriteOmegaHists()
{
	printf("\n  Writing omega histograms...\n");
	
	h_3gamma_m12->Write();
	h_3gamma_m13->Write();
	h_3gamma_m23->Write();
	h_3gamma_m3g->Write();
	h_3gamma_m12_elas->Write();
	h_3gamma_m13_elas->Write();
	h_3gamma_m23_elas->Write();
	h_3gamma_m3g_elas->Write();
	h_3gamma_vz->Write();
	h_3gamma_vz_elas->Write();
	
	h_3gamma_theta_targ->Write();
	h_3gamma_theta_fdc1->Write();
	h_3gamma_theta_fdc2->Write();
	h_3gamma_theta_fdc3->Write();
	h_xy_targ->Write();
	h_xy_fdc1->Write();
	h_xy_fdc2->Write();
	h_xy_fdc3->Write();
	
	printf("  Done.\n");
	return;
}
