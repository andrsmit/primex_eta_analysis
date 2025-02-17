#include "EtaAna.h"

void EtaAna::EtaggAnalysis_bggen() {
	
	// Skip non-mc events:
	
	if(m_nmc<=0) {
		return;
	}
	
	// Skip events where the initial energy is below cut:
	
	double locMCEnergy = 0.0;
	for(int imc=0; imc<m_nmc; imc++) locMCEnergy += m_mcE[imc];
	if(locMCEnergy < (m_minBeamEnergyCut+ParticleMass(Proton))) return;
	
	// Distribute events realistically inside target:
	
	if(AcceptRejectEvent()) return;
	
	// get reaction index:
	
	int locFinalState = GetFinalState_bggen();
	h_thrown_reactions_bggen->Fill(locFinalState);
	
	double locThrownBeamEnergy = 0.0, locThrownAngle = 0.0;
	GetThrownEnergyAndAngle(locThrownBeamEnergy, locThrownAngle);
	if(locFinalState==0) {
		h_thrown_proton->Fill(locThrownBeamEnergy, locThrownAngle);
	}
	else if(locFinalState==1) {
		h_thrown_neutron->Fill(locThrownBeamEnergy, locThrownAngle);
	}
	
	PlotThrown(locThrownBeamEnergy, locThrownAngle);
	
	
	// Get list of good photons, FCAL showers, and BCAL showers to analyze:
	
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
			
			bool isHadronicVeto = true;
			
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
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a coherently-produced eta meson:
				
				double sig1 = GetFCALEnergyResolution(e1);
				double sig2 = GetFCALEnergyResolution(e2);
				double sigr = sig1/sig2;
				
				// Energy-constrained invariant mass assuming production on free nucleon:
				double e1c = e1/(1.+sigr) + (etaEnergy-e2)/(1.+(1./sigr));
				double e2c = etaEnergy - e1c;
				double invmassConstr = sqrt(2.*e1c*e2c*(1.-cos12));
				
				//-----------------------------------------------------//
				
				if((locFinalState==0) || (locFinalState==1)) {
					// gamma+N -> eta+N
					
					if(isElastic) h_mgg_const_bggen_signal->Fill(prodTheta, invmassConstr, fillWeight);
					h_elas_bggen_signal->Fill(prodTheta, Egg/etaEnergy, fillWeight);
					if(isEta) h_elas_bggen_signal_cut->Fill(prodTheta, Egg/etaEnergy, fillWeight);
				}
				else if((locFinalState>=2) && (locFinalState<=5)) {
					// gamma+N -> eta+pion+N'
					
					if(isElastic) h_mgg_const_bggen_etapion->Fill(prodTheta, invmassConstr, fillWeight);
					h_elas_bggen_etapion->Fill(prodTheta, Egg/etaEnergy, fillWeight);
					if(isEta) h_elas_bggen_etapion_cut->Fill(prodTheta, Egg/etaEnergy, fillWeight);
				}
				else {
					// gamma+N -> eta+pion+N'
					
					if(isElastic) h_mgg_const_bggen_bkgd->Fill(prodTheta, invmassConstr, fillWeight);
					h_elas_bggen_bkgd->Fill(prodTheta, Egg/etaEnergy, fillWeight);
					if(isEta) h_elas_bggen_bkgd_cut->Fill(prodTheta, Egg/etaEnergy, fillWeight);
				}
				
				// make individual plots for every specified channel:
				if(isElastic) {
					h_mgg_const_bggen[locFinalState]->Fill(prodTheta, invmassConstr, fillWeight);
					if(isEta && (locNBCALShowers==1)) {
						h_bcalDeltaPhi_bggen[locFinalState]->Fill(prodTheta, fabs(locBCALPhi-prodPhi), fillWeight);
					}
				}
				
				if(isElastic && isEta) {
					h_rec_reactions_bggen->Fill(locFinalState, fillWeight);
					
					if((m_minBeamEnergyBin<=locThrownBeamEnergy) && (locThrownBeamEnergy<m_maxBeamEnergyBin)) {
						if(locFinalState==0) {
							h_AngularMatrix_proton->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
						}
						else if(locFinalState==1) {
							h_AngularMatrix_neutron->Fill(locThrownAngle, prodTheta, locThrownBeamEnergy, fillWeight);
						}
					}
				}
				
			} // end loop over beam photons
		} // end loop 2 over fcal showers
	} // end loop 1 over fcal showers
	
	return;
}

void EtaAna::InitializeBGGENHists() {
	
	int nInvmassBins     = (int)((m_maxInvmassBin-m_minInvmassBin)/m_invmassBinSize);
	int nRecAngleBins    = (int)((m_maxRecAngleBin-m_minRecAngleBin)/m_recAngleBinSize);
	int nBeamEnergyBins  = (int)((m_maxBeamEnergyBin-m_minBeamEnergyBin)/m_beamEnergyBinSize);
	int nThrownAngleBins = (int)((m_maxThrownAngleBin-m_minThrownAngleBin)/m_thrownAngleBinSize);
	
	int nReactions = (int)m_reaction_types.size();
	
	h_thrown_reactions_bggen   = new TH1F("thrown_reactions_bggen", 
		"Number Thrown; Reaction Type; N_{thrown}", 
		nReactions+1, -0.5, ((double)nReactions)+0.5);
	h_rec_reactions_bggen = new TH1F("rec_reactions_bggn", 
		"Number Accepted; Reaction Type; N_{rec}", 
		nReactions+1, -0.5, ((double)nReactions)+0.5);
	
	for(int irea=0; irea<nReactions; irea++) {
		TString bin_label = "";
		for(int iparticle=0; iparticle<m_reaction_types[irea].size(); iparticle++) {
			bin_label += ParticleName_ROOT(m_reaction_types[irea][iparticle]);
		}
		h_thrown_reactions_bggen->GetXaxis()->SetBinLabel(irea+1, bin_label.Data());
		h_rec_reactions_bggen->GetXaxis()->SetBinLabel(irea+1, bin_label.Data());
	}
	h_thrown_reactions_bggen->GetXaxis()->SetBinLabel(h_thrown_reactions_bggen->GetNbinsX(), "other");
	h_rec_reactions_bggen->GetXaxis()->SetBinLabel(h_rec_reactions_bggen->GetNbinsX(), "other");
	
	//-----------------------------------------------------------------------------//
	
	h_mgg_const_bggen_signal = new TH2F("mgg_const_bggen_signal",
		"Signal (#eta+N#rightarrow#eta+N); #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_const_bggen_signal->Sumw2();
	
	h_mgg_const_bggen_etapion = new TH2F("mgg_const_bggen_etapion",
		"Eta+Pion (#eta+N#rightarrow#eta+#pi+N'); #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_const_bggen_etapion->Sumw2();
	
	h_mgg_const_bggen_bkgd = new TH2F("mgg_const_bggen_bkgd",
		"Other Bkgd; #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
	h_mgg_const_bggen_bkgd->Sumw2();
	
	//-----------------------------------------------------------------------------//
	
	h_elas_bggen_signal = new TH2F("elas_bggen_signal",
		"Signal (#eta+N#rightarrow#eta+N); #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
	
	h_elas_bggen_etapion = new TH2F("elas_bggen_etapion",
		"Eta+Pion (#eta+N#rightarrow#eta+#pi+N'); #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
	
	h_elas_bggen_bkgd = new TH2F("elas_bggen_bkgd",
		"Other Bkgd; #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
	
	//-----------------------------------------------------------------------------//
	
	h_elas_bggen_signal_cut = new TH2F("elas_bggen_signal_cut",
		"Signal (#eta+N#rightarrow#eta+N); #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
	
	h_elas_bggen_etapion_cut = new TH2F("elas_bggen_etapion_cut",
		"Eta+Pion (#eta+N#rightarrow#eta+#pi+N'); #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
	
	h_elas_bggen_bkgd_cut = new TH2F("elas_bggen_bkgd_cut",
		"Other Bkgd; #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]",
		nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1000, 0.0, 2.0);
	
	//-----------------------------------------------------------------------------//
	
	for(int irea=0; irea<=nReactions; irea++) {
		
		TString hist_label = "", hist_title = "";
		if(irea==nReactions) {
			hist_label = "other";
			hist_title = "other";
		}
		else {
			for(int iparticle=0; iparticle<m_reaction_types[irea].size(); iparticle++) {
				hist_label += ShortName(m_reaction_types[irea][iparticle]);
				hist_title += ParticleName_ROOT(m_reaction_types[irea][iparticle]);
			}
		}
		
		TH2F *loc_h2 = new TH2F(Form("mgg_const_bggen_%s",hist_label.Data()),
			Form("%s; #theta_{rec.} [#circ]; m_{#gamma#gamma}^{Constr.} [GeV/c^{2}]", hist_title.Data()),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, nInvmassBins, m_minInvmassBin, m_maxInvmassBin);
		loc_h2->Sumw2();
		h_mgg_const_bggen.push_back(loc_h2);
		
		TH2F *loc_h2_deltaPhi = new TH2F(Form("bcal_deltaPhi_bggen_%s",hist_label.Data()),
			Form("%s; #theta_{rec.} [#circ]; #left|#phi_{#gamma#gamma}-#phi_{BCAL}#right| [#circ]", hist_title.Data()),
			nRecAngleBins, m_minRecAngleBin, m_maxRecAngleBin, 1800, 0.0, 360.0);
		h_bcalDeltaPhi_bggen.push_back(loc_h2_deltaPhi);
	}
	
	//-----------------------------------------------------------------------------//
	
	h_thrown_proton = new TH2F("thrown_proton", 
		"Thrown Angle vs. Thrown Beam Energy (#gamma+p#rightarrow#eta+p); E_{#gamma}(thrown) [GeV]; #theta(thrown) [#circ]",
		nBeamEnergyBins, m_minBeamEnergyBin, m_maxBeamEnergyBin, 
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin);
	
	h_thrown_neutron = new TH2F("thrown_neutron", 
		"Thrown Angle vs. Thrown Beam Energy (#gamma+n#rightarrow#eta+n); E_{#gamma}(thrown) [GeV]; #theta(thrown) [#circ]",
		nBeamEnergyBins, m_minBeamEnergyBin, m_maxBeamEnergyBin, 
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin);
	
	h_AngularMatrix_proton = new TH3F("AngularMatrix_proton", 
		"#gamma+p#rightarrow#eta+p; #theta(thrown) [#circ]; #theta(rec) [#circ]",
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
		nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
	
	h_AngularMatrix_neutron = new TH3F("AngularMatrix_neutron", 
		"#gamma+n#rightarrow#eta+n; #theta(thrown) [#circ]; #theta(rec) [#circ]",
		nThrownAngleBins, m_minThrownAngleBin, m_maxThrownAngleBin, 
		nRecAngleBins,    m_minRecAngleBin,    m_maxRecAngleBin,
		nBeamEnergyBins,  m_minBeamEnergyBin,  m_maxBeamEnergyBin);
	
	return;
}

void EtaAna::ResetBGGENHists() {
	
	h_mgg_const_bggen_signal->Reset();
	h_mgg_const_bggen_etapion->Reset();
	h_mgg_const_bggen_bkgd->Reset();
	
	h_elas_bggen_signal->Reset();
	h_elas_bggen_etapion->Reset();
	h_elas_bggen_bkgd->Reset();
	
	h_elas_bggen_signal_cut->Reset();
	h_elas_bggen_etapion_cut->Reset();
	h_elas_bggen_bkgd_cut->Reset();
	
	h_thrown_reactions_bggen->Reset();
	h_rec_reactions_bggen->Reset();
	
	for(int irea=0; irea<(int)h_mgg_const_bggen.size(); irea++) {
		h_mgg_const_bggen[irea]->Reset();
		h_bcalDeltaPhi_bggen[irea]->Reset();
	}
	
	h_thrown_proton->Reset();
	h_thrown_neutron->Reset();
	h_AngularMatrix_proton->Reset();
	h_AngularMatrix_neutron->Reset();
	return;
}

void EtaAna::WriteBGGENHists() {
	printf("\nWriting BGGEN histograms...\n");
	h_mgg_const_bggen_signal->Write();
	h_mgg_const_bggen_etapion->Write();
	h_mgg_const_bggen_bkgd->Write();
	
	h_elas_bggen_signal->Write();
	h_elas_bggen_etapion->Write();
	h_elas_bggen_bkgd->Write();
	
	h_elas_bggen_signal_cut->Write();
	h_elas_bggen_etapion_cut->Write();
	h_elas_bggen_bkgd_cut->Write();
	
	h_thrown_reactions_bggen->Write();
	h_rec_reactions_bggen->Write();
	
	for(int irea=0; irea<(int)h_mgg_const_bggen.size(); irea++) {
		h_mgg_const_bggen[irea]->Write();
		h_bcalDeltaPhi_bggen[irea]->Write();
	}
	
	h_thrown_proton->Write();
	h_thrown_neutron->Write();
	h_AngularMatrix_proton->Write();
	h_AngularMatrix_neutron->Write();
	printf("Done.\n");
	return;
}

void EtaAna::InitializeReactionTypes() {
	
	m_reaction_types.clear();
	
	m_reaction_types.push_back({Eta, Proton});
	m_reaction_types.push_back({Eta, Neutron});
	m_reaction_types.push_back({Eta, Pi0, Proton});
	m_reaction_types.push_back({Eta, Pi0, Neutron});
	m_reaction_types.push_back({Eta, PiMinus, Proton});
	m_reaction_types.push_back({Eta, PiPlus, Neutron});
	m_reaction_types.push_back({omega, Proton});
	m_reaction_types.push_back({omega, Neutron});
	m_reaction_types.push_back({Rho0, Proton});
	m_reaction_types.push_back({Rho0, Neutron});
	m_reaction_types.push_back({phiMeson, Proton});
	m_reaction_types.push_back({phiMeson, Neutron});
	m_reaction_types.push_back({Eta, Pi0, Pi0, Proton});
	m_reaction_types.push_back({Eta, Pi0, Pi0, Neutron});
	m_reaction_types.push_back({Eta, Pi0, PiPlus, Neutron});
	m_reaction_types.push_back({Eta, Pi0, PiMinus, Proton});
	m_reaction_types.push_back({Eta, PiPlus, PiMinus, Proton});
	m_reaction_types.push_back({Eta, PiPlus, PiMinus, Neutron});
	
	return;
}

int EtaAna::GetFinalState_bggen(int debug) {
	
	vector<Particle_t> thrown_vec;
	thrown_vec.clear();
	for(int imc=0; imc<m_nmc; imc++) {
		thrown_vec.push_back(PDGtoPType((int)m_mcPDGType[imc]));
	}
	
	if(debug) {
		cout << "Thrown list of particles: ";
		for(int imc=0; imc<(m_nmc-1); imc++) {
			cout << ParticleType(thrown_vec[imc]) << ", ";
		}
		cout << ParticleType(thrown_vec[m_nmc-1]) << endl;
	}
	
	/*
	We want to loop over the specified final states stored in the vector 'm_reaction_types'
	and see if any of them match the final state we have in 'thrown_vec'.
	
	This boils down to comparing two unsorted vectors containing Particle_t objects.
	Furthermore, the vectors could have repeated elements.
	
	To do this, we first check that the number of elements in 'thrown_vec' equals the number of elements 
	in the specific vector element of 'm_reaction_types' that we are currently looking at 
	(remember we are looping through these). If they are unequal, we move on to the next element of 'm_reaction_types'. 
	Otherwise, we make a new vector called 'loc_final_state' as a copy of the vector-element of 'm_reaction_types' we're on.
	Then, we loop through each element of 'thrown_vec' and if that element is found in 'loc_final_state', we 
	delete it from 'loc_final_state'.
	If we come across an element of 'thrown_vec' that's not in the current version of 'loc_final_state', the two vectors
	do not match and we keep looping.
	If we loop through all elements of 'thrown_vec' and they all have a match, then we found our final state match and we 
	break the loop.
	*/
	
	int final_state = (int)m_reaction_types.size();
	
	// loop through 'm_reaction_types' vector and see if it matches the thrown final state:
	for(int itype = 0; itype < m_reaction_types.size(); itype++) {
		
		int loc_n_particles = m_reaction_types[itype].size();
		if(m_nmc != loc_n_particles) continue;
		
		vector<Particle_t> loc_final_state;
		loc_final_state.clear();
		for(int i=0; i<loc_n_particles; i++) loc_final_state.push_back(m_reaction_types[itype][i]);
		
		if(debug) {
			cout << "  comparing to: ";
			for(int imc=0; imc<(loc_final_state.size()-1); imc++) {
				cout << ParticleType(loc_final_state[imc]) << ", ";
			}
			cout << ParticleType(loc_final_state[loc_final_state.size()-1]) << "..";
		}
		
		int match_val = 1;
		for(int iparticle = 0; iparticle < thrown_vec.size(); iparticle++) {
			
			// look for this particle in 'loc_final_state':
			int found_val = 0;
			for(int jparticle = 0; jparticle < loc_final_state.size(); jparticle++) {
				if(thrown_vec[iparticle] == loc_final_state[jparticle]) {
					// delete this element from loc_final_state:
					loc_final_state.erase(loc_final_state.begin()+jparticle);
					found_val = 1;
					break;
				}
			}
			// if there's an element of 'thrown_vec' that isn't found in 'loc_final_state', then move on:
			if(found_val==0) {
				match_val = 0;
				break;
			}
		}
		if(match_val==1) {
			final_state = itype;
			if(debug) cout << "match!" << endl;
			break;
		} else {
			if(debug) cout << "not matched" << endl;
		}
	}
	
	return final_state;
}
