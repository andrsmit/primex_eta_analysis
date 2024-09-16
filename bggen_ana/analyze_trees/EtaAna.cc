#include "EtaAna.h"

// Default Constructor:

EtaAna::EtaAna() {
	
	// set default run number to 61434 (first run at 11.2 GeV):
	setRunNumber(61434);
	initializeReactionTypes();
	
	m_random = new TRandom3(0);
	
	// set defaults for cut values:
	m_BEAM_RF_CUT =  2.004;
	m_FCAL_RF_CUT =  2.004;
	m_BCAL_RF_CUT = 12.0;
	m_TOF_RF_CUT  =  1.0;
	
	// default values for the minimum energy cuts:
	m_MIN_FCAL_ENERGY = 0.5; // energy of each FCAL shower used in the analysis
	m_MIN_BEAM_ENERGY = 8.0; // energy of the tagged photon energy
	m_MIN_BCAL_ENERGY = 0.;  // energy sum of BCAL showers within timing cut
	
	// default cut value for selecting a match between the FCAL and TOF:
	m_FCAL_TOF_CUT = 8.0; // distance between fcal shower and closest DTOFPoint
	
	// default value for elasticity cut:
	m_ELAS_CUT_SIGMA = 0.031;
	m_ELAS_CUT_WIDTH = 5.0;
	m_ELAS_CUT_MU_P0 = 1.0; // mu = p0 + p1*E_gamma
	m_ELAS_CUT_MU_P1 = 0.0;
	
	// Set event number to zero on initialization:
	m_event = 0;
}

void EtaAna::initializeReactionTypes() {
	
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

void EtaAna::runAnalysis(TString infname) {
	
	m_infile = new TFile(infname.Data(), "READ");
	int n_events_total = loadTree();
	
	while(m_event < n_events_total) {
		readEvent();
		etaggAnalysis();
		m_event++;
	}
	
	m_infile->Close();
	
	return;
}

void EtaAna::etaggAnalysis() {
	
	// skip events where the initial energy is below 8 GeV:
	double mc_energy = 0.0;
	for(int i=0; i<m_nmc; i++) {
		mc_energy += m_mc_e[i];
	}
	if(mc_energy < (8.0+m_Proton)) return;
	
	int locFinalState = getFinalState_bggen();
	
	h_thrown->Fill(locFinalState);
	
	if(m_nfcal<2) return;
	
	//---------------------------------------------------------------------------//
	// Make a list of good FCAL showers to use in analysis:
	
	vector<int> locGoodFCALShowers;
	locGoodFCALShowers.clear();
	int locNFCALShowers = 0, locNGoodFCALShowers = 0;
	for(int ishow=0; ishow<m_nfcal; ishow++) {
		
		TVector3 loc_pos = getFCALPosition(ishow);
		double loc_t = m_fcal_t[ishow] - (loc_pos.Mag()/m_c) - m_rfTime;
		
		int loc_fid_cut = fcal_fiducial_cut(loc_pos, 2.0);
		if(fabs(loc_t) < m_FCAL_RF_CUT) {
			locNFCALShowers++;
			if(m_fcal_e[ishow] > m_MIN_FCAL_ENERGY) {
				locNGoodFCALShowers++;
				if(!loc_fid_cut) {
					locGoodFCALShowers.push_back(ishow);
				}
			}
		}
	}
	
	//---------------------------------------------------------------------------//
	// Make a list of prompt and selected-sideband beam photons to use in analysis:
	
	vector<int> locGoodBeamPhotons;
	vector<double> locGoodBeamPhotons_weight;
	locGoodBeamPhotons.clear();
	locGoodBeamPhotons_weight.clear();
	
	for(int igam=0; igam<m_nbeam; igam++) {
		
		double loc_dt     = m_beam_t[igam] - m_rfTime;
		double loc_weight = 0.0;
		
		double loc_beam_cut = m_beam_bunches_main*m_BEAM_RF_CUT;
		
		if(fabs(loc_dt) < loc_beam_cut) loc_weight = 1.0;
		else if(
			((m_beam_bunches_main+5.5)*4.008 < fabs(loc_dt)) && 
			(fabs(loc_dt) < (m_beam_bunches_main+5.5+m_beam_bunches_acc)*4.008)
		) loc_weight = -1.0/(2.0*m_beam_bunches_acc);
		else continue;
		
		//if(loc_weight < 0.0) loc_weight *= m_acc_scale_factor[igam];
		if(m_beam_e[igam] > m_MIN_BEAM_ENERGY) {
			locGoodBeamPhotons.push_back(igam);
			locGoodBeamPhotons_weight.push_back(loc_weight);
		}
	}
	
	// BCAL Veto:
	int    locNBCALShowers  = 0, locNBCALShowers_1ns = 0;
	double locBCALEnergySum = 0.;
	double locBCALRFDT = 0., locBCALPhi = 0.; // only useful when there's exactly 1 BCAL shower within timing cut
	for(int ishow=0; ishow<m_nbcal; ishow++) {
		TVector3 loc_pos(m_bcal_x[ishow], m_bcal_y[ishow], m_bcal_z[ishow]);
		loc_pos -= m_vertex;
		double loc_t = m_bcal_t[ishow] - (loc_pos.Mag()/m_c) - m_rfTime;
		if(fabs(loc_t) < 12.0) {
			locBCALEnergySum += m_bcal_e[ishow];
			locNBCALShowers++;
			locBCALRFDT = loc_t;
			locBCALPhi  = loc_pos.Phi() * (180./TMath::Pi());
			if(fabs(loc_t) < 1.0) {
				locNBCALShowers_1ns++;
			}
		}
	}
	
	// FCAL multiplicity cut:
	if(!(locNFCALShowers==2 && locNGoodFCALShowers==2)) return;
	
	//-----------------------------------------------------//
	// Loop over pairs of FCAL showers:
	
	for(int ishow = 0; ishow < locGoodFCALShowers.size(); ishow++) {
		
		TVector3 pos1 = getFCALPosition(ishow);
		double t1 = m_fcal_t[ishow] - (pos1.Mag()/m_c) - m_rfTime;
		double e1 = m_fcal_e[ishow];
		
		// apply minimum energy and RF timing cuts:
		if(e1 < m_MIN_FCAL_ENERGY || fabs(t1) >= m_FCAL_RF_CUT) continue;
		
		// apply fiducial cut to remove the innermost two FCAL layers:
		if(fcal_fiducial_cut(pos1, 2.0)) continue;
		
		double px1 = e1*pos1.X()/pos1.Mag();
		double py1 = e1*pos1.Y()/pos1.Mag();
		double pz1 = e1*pos1.Z()/pos1.Mag();
		
		// check the distance between this shower and the closest (if any) tof hit:
		double tof_dx1, tof_dy1, tof_dt1;
		check_TOF_match(pos1, tof_dx1, tof_dy1, tof_dt1, m_TOF_RF_CUT);
		double tof_dr1 = sqrt(pow(tof_dx1,2.0)+pow(tof_dy1,2.0));
		
		//-----------------------------------------------------//
		
		for(int jshow = ishow+1; jshow < locGoodFCALShowers.size(); jshow++) {
			
			TVector3 pos2 = getFCALPosition(jshow);
			double t2 = m_fcal_t[jshow] - (pos2.Mag()/m_c) - m_rfTime;
			double e2 = m_fcal_e[jshow];
			
			// apply minimum energy and RF timing cuts:
			if(e2 < m_MIN_FCAL_ENERGY || fabs(t2) >= m_FCAL_RF_CUT) continue;
			
			// apply fiducial cut to remove the innermost two FCAL layers:
			if(fcal_fiducial_cut(pos2, 2.0)) continue;
			
			double px2 = e2*pos2.X()/pos2.Mag();
			double py2 = e2*pos2.Y()/pos2.Mag();
			double pz2 = e2*pos2.Z()/pos2.Mag();
			
			// check the distance between this shower and the closest (if any) tof hit:
			double tof_dx2, tof_dy2, tof_dt2;
			check_TOF_match(pos2, tof_dx2, tof_dy2, tof_dt2, m_TOF_RF_CUT);
			double tof_dr2 = sqrt(pow(tof_dx2,2.0)+pow(tof_dy2,2.0));
			
			//-----------------------------------------------------//
			// TOF Veto
			
			if((tof_dr1 < m_FCAL_TOF_CUT) && (tof_dr2 < m_FCAL_TOF_CUT)) continue;
			
			//-----------------------------------------------------//
			// Two-Photon kinematics:
			
			double Egg  =  e1 +  e2; // energy of 2-photon pair
			double pggx = px1 + px2; // momentum along x-axis
			double pggy = py1 + py2; // momentum along y-axis
			double pggz = pz1 + pz2; // momentum along z-axis
			
			// transverse momentum:
			double pggt = sqrt(pow(pggx,2.0) + pow(pggy,2.0));
			
			// polar angle:
			double prod_th = (180./TMath::Pi()) * atan2(pggt,pggz);
			
			// azimuthal angle:
			double prod_phi = (180./TMath::Pi()) * atan2(pggy,pggx);
			
			// opening angle:
			double cos12   = (pos1.X()*pos2.X() + pos1.Y()*pos2.Y() + pos1.Z()*pos2.Z()) / (pos1.Mag()*pos2.Mag());
			
			// invariant mass:
			double invmass = sqrt(2.0*e1*e2*(1.-cos12));
			
			// check different BCAL veto options:
			
			vector<int> loc_bcal_vetos;
			for(int iveto=0; iveto<m_n_bcal_vetos; iveto++) loc_bcal_vetos.push_back(0);
			if(locNBCALShowers==0) 
				loc_bcal_vetos[0] = 1;
			if(locNBCALShowers_1ns==0) 
				loc_bcal_vetos[1] = 1;
			if((locNBCALShowers==0) ||
				(locNBCALShowers==1 && fabs(fabs(locBCALPhi-prod_phi)-180.0) < 30.0))
				loc_bcal_vetos[2] = 1;
			if((locNBCALShowers==0) || 
				(locNBCALShowers==1 && fabs(fabs(locBCALPhi-prod_phi)-180.0) < 30.0 && locBCALRFDT>1.0))
				loc_bcal_vetos[3] = 1;
			
			//-----------------------------------------------------//
			// Loop over Beam photons
			
			for(int igam=0; igam<(int)locGoodBeamPhotons.size(); igam++) {
				
				int loc_beam_index  = locGoodBeamPhotons[igam];
				double eb           = m_beam_e[loc_beam_index];
				double fill_weight  = locGoodBeamPhotons_weight[igam];
				
				int loc_tag_sys     = m_tag_sys[loc_beam_index];
				int loc_tag_counter = m_tag_counter[loc_beam_index];
				
				// Calculate the energy of the eta meson, assuming a coherent production process:
				double eeta = energy_after_recoil(eb, prod_th, m_eta, m_Target);
				
				// Apply a cut on the elasticity
				//  (ratio of measured energy of 2-photons, to the calculated energy above):
				
				bool elas_cut = false;
				double loc_elas_mean  = m_ELAS_CUT_MU_P0 + m_ELAS_CUT_MU_P1*prod_th;
				double loc_elas_width = m_ELAS_CUT_WIDTH * m_ELAS_CUT_SIGMA;
				if(fabs((Egg/eeta)-loc_elas_mean)<loc_elas_width) elas_cut = true;
				
				// set a variable to indicate if the two-photon mass is consistent with an eta meson:
				bool  eta_cut = false;
				if(0.497862<invmass && invmass<0.597862) eta_cut = true;
				
				//-----------------------------------------------------//
				// Energy constraint
				
				// adjust the measured energies of the two-photons to exactly equal the energy of 
				// a coherently-produced eta meson:
				
				double sig1 = fcal_energy_res(e1);
				double sig2 = fcal_energy_res(e2);
				double sigr = pow(sig1/sig2,2.0);
				
				double e1c = e1/(1.+sigr) + (eeta-e2)/(1.+(1./sigr));
				double e2c = eeta - e1c;
				double invmass_const = sqrt(2.*e1c*e2c*(1.-cos12)); // energy-constrained invariant mass
				
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
				double prod_th_const = (180./TMath::Pi()) * atan2(pggtc,pggzc);
				
				//-----------------------------------------------------//
				// Missing Mass 
				//   - This is calculated assuming the reaction took place on a quasifree nucleon
				//     However, we neglect the fermi-momentum of said nucleon
				
				double mmsq = 2.0*m_Proton*eb - 2.0*eb*Egg + m_Proton*m_Proton + m_eta*m_eta 
					- 2.0*m_Proton*Egg + 2.0*eb*cos(prod_th*TMath::Pi()/180.)*sqrt(Egg*Egg - m_eta*m_eta);
				
				//-----------------------------------------------------//
				// Different ways to apply BCAL Veto:
				/*
				 1. Veto events with any BCAL shower within +/-12 ns (strict veto)
				 2. Veto events with any BCAL shower within +/-1 ns
				 3. Allow a single BCAL shower if the deltaPhi from two-photon pair is 180 +/- 30, no cut on timing
				 4. Allow a single BCAL shower if the deltaPhi from two-photon pair is 180 +/- 30 and the timing is >1ns
				*/
				//-----------------------------------------------------//
				
				h_elas_vs_mgg->Fill(invmass/m_eta, Egg/eeta, fill_weight);
				
				if(eta_cut && loc_bcal_vetos[3]) {
					h_elas->Fill(prod_th, Egg/eb, fill_weight);
					h_elas_corr->Fill(prod_th, Egg/eeta, fill_weight);
					if(locFinalState==0 || locFinalState==1) {
						h_elas_signal->Fill(prod_th, Egg/eb, fill_weight);
						h_elas_corr_signal->Fill(prod_th, Egg/eeta, fill_weight);
					} else {
						h_elas_bkgd->Fill(prod_th, Egg/eb, fill_weight);
						h_elas_corr_bkgd->Fill(prod_th, Egg/eeta, fill_weight);
					}
				}
				
				// apply elasticity cut and plot the invariant mass distriubtion:
				if(elas_cut) {
					
					if(loc_bcal_vetos[3]) {
						// invariant mass vs. polar angle:
						h_mgg->Fill(prod_th, invmass, fill_weight);
						
						// energy-constrained invariant mass vs. polar angle:
						h_mgg_const->Fill(prod_th, invmass_const, fill_weight);
						
						if(locFinalState==0 || locFinalState==1) {
							h_mgg_signal->Fill(prod_th, invmass, fill_weight);
							h_mgg_const_signal->Fill(prod_th, invmass_const, fill_weight);
						} else {
							h_mgg_bkgd->Fill(prod_th, invmass, fill_weight);
							h_mgg_const_bkgd->Fill(prod_th, invmass_const, fill_weight);
						}
					}
					if(eta_cut) {
						
						h_theta[locFinalState]->Fill(prod_th, fill_weight);
						if(loc_bcal_vetos[3]) {
							h_accepted->Fill(locFinalState, fill_weight);
						}
						
						for(int iveto=0; iveto<m_n_bcal_vetos; iveto++) {
							if(loc_bcal_vetos[iveto]) h_theta_veto[locFinalState][iveto]->Fill(prod_th, fill_weight);
						}
						
						// plot the distribution of BCAL showers for each reaction type:
						h_nbcal[locFinalState]->Fill(locNBCALShowers, fill_weight);
						h_bcal_energy[locFinalState]->Fill(locBCALEnergySum, fill_weight);
						if(locNBCALShowers==1) {
							h_bcal_energy_single[locFinalState]->Fill(locBCALEnergySum, fill_weight);
							// plot time difference vs angle of eta:
							h_bcal_dt_vs_eta_angle[locFinalState]->Fill(prod_th, locBCALRFDT, fill_weight);
							// plot time difference vs energy of BCAL shower:
							h_bcal_dt_vs_bcal_energy[locFinalState]->Fill(locBCALEnergySum, locBCALRFDT, fill_weight);
							// plot deltaPhi:
							h_bcal_deltaPhi[locFinalState]->Fill(fabs(locBCALPhi-prod_phi), fill_weight);
						}
						
						// plot x-y distribution of showers:
						h_xy_1->Fill(pos1.X(), pos1.Y(), fill_weight);
						h_xy_2->Fill(pos2.X(), pos2.Y(), fill_weight);
					}
				}
				
				// missing mass:
				h_mm_vs_theta->Fill(prod_th, mmsq, fill_weight);
				
				// plot the above distribution with a cut around the eta mass:
				
				if(eta_cut) {
					h_mm_vs_theta_eta_cut->Fill(prod_th, mmsq, fill_weight);
					
					// next, apply elasticity cut:
					if(elas_cut) {
						h_mm_vs_theta_eta_elas_cut->Fill(prod_th, mmsq, fill_weight);
					}
				}
				
			} // end loop over beam photons
		} // end loop over FCAL shower
	} // end loop over FCAL showers
	
	return;
}

int EtaAna::getFinalState_bggen(int debug) {
	
	vector<Particle_t> loc_thrown_vec;
	loc_thrown_vec.clear();
	for(int imc=0; imc<m_nmc; imc++) {
		loc_thrown_vec.push_back(PDGtoPType((int)m_mc_pdgtype[imc]));
	}
	
	if(debug) {
		cout << "Thrown list of particles: ";
		for(int imc=0; imc<(m_nmc-1); imc++) {
			cout << ParticleType(loc_thrown_vec[imc]) << ", ";
		}
		cout << ParticleType(loc_thrown_vec[m_nmc-1]) << endl;
	}
	
	/*
	We want to loop over the specified final states stored in the vector 'm_reaction_types'
	and see if any of them match the final state we have in 'loc_thrown_vec'.
	
	This boils down to comparing two unsorted vectors containing Particle_t objects.
	Furthermore, the vectors could have repeated elements.
	
	To do this, we first check that the number of elements in 'loc_thrown_vec' equals the number of elements 
	in the specific vector element of 'm_reaction_types' that we are currently looking at 
	(remember we are looping through these). If they are unequal, we move on to the next element of 'm_reaction_types'. 
	Otherwise, we make a new vector called 'loc_final_state' as a copy of the vector-element of 'm_reaction_types' we're on.
	Then, we loop through each element of 'loc_thrown_vec' and if that element is found in 'loc_final_state', we 
	delete it from 'loc_final_state'.
	If we come across an element of 'loc_thrown_vec' that's not in the current version of 'loc_final_state', the two vectors
	do not match and we keep looping.
	If we loop through all elements of 'loc_thrown_vec' and they all have a match, then we found our final state match and we 
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
		for(int iparticle = 0; iparticle < loc_thrown_vec.size(); iparticle++) {
			
			// look for this particle in 'loc_final_state':
			int found_val = 0;
			for(int jparticle = 0; jparticle < loc_final_state.size(); jparticle++) {
				if(loc_thrown_vec[iparticle] == loc_final_state[jparticle]) {
					// delete this element from loc_final_state:
					loc_final_state.erase(loc_final_state.begin()+jparticle);
					found_val = 1;
					break;
				}
			}
			// if there's an element of 'loc_thrown_vec' that isn't found in 'loc_final_state', then move on:
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

TVector3 EtaAna::getFCALPosition(int index) {
	
	TVector3 pos(m_fcal_x[index], m_fcal_y[index], m_fcal_z[index]);
	pos = pos + m_fcal_correction - m_vertex;
	
	return pos;
}

int EtaAna::fcal_fiducial_cut(TVector3 pos, double cut_layer) {
	
	int fid_cut = 0;
	
	double fcal_inner_layer_cut = (1.5 + cut_layer) * m_fcal_block_size;
	
	double fcal_face_x = m_vertex.X() + (pos.X() * (m_fcal_face.Z() - m_vertex.Z())/pos.Z()) - m_fcal_face.X();
	double fcal_face_y = m_vertex.Y() + (pos.Y() * (m_fcal_face.Z() - m_vertex.Z())/pos.Z()) - m_fcal_face.Y();
	
	if((fabs(fcal_face_x) < fcal_inner_layer_cut) && (fabs(fcal_face_y) < fcal_inner_layer_cut)) fid_cut = 1;
	
	// only apply the next fiducial cut for runs from phase-I:
	
	if(m_phase_val < 2) {
		if((-32.<fcal_face_y && fcal_face_y<-20.) && (-8.<fcal_face_x && fcal_face_x<4.)) {
			fid_cut = 1;
		}
	}
	
	return fid_cut;
}

void EtaAna::check_TOF_match(TVector3 pos1, double &dx_min, double &dy_min, double &dt_min, double rf_time_cut) {
	
	dx_min = 1000.;
	dy_min = 1000.;
	dt_min = 1000.;
	
	for(int itof=0; itof<m_ntof; itof++) {
		
		double xt = m_tof_x[itof] - m_vertex.X();
		double yt = m_tof_y[itof] - m_vertex.Y();
		double zt = m_tof_z[itof] - m_vertex.Z();
		double rt = sqrt(xt*xt + yt*yt + zt*zt);
		double tt = m_tof_t[itof] - (rt/m_c);
		double dt = tt - m_rfTime;
		xt *= pos1.Z() / zt;
		yt *= pos1.Z() / zt;
		double dx = pos1.X() - xt;
		double dy = pos1.Y() - yt;
		
		if(fabs(dt) < rf_time_cut) {
			if((dx*dx + dy*dy) < (dx_min*dx_min + dy_min*dy_min)) {
				dx_min = dx;
				dy_min = dy;
				dt_min = dt;
			}
		}
	}
	return;
}

double EtaAna::energy_after_recoil(double eb, double theta, double m0, double mp) 
{
	theta *= (TMath::Pi()/180.);
	
	double t1 = eb*cos(theta);
	double t2 = mp+eb;
	double t3 = mp*eb + m0*m0*0.5;
	
	double a = t1*t1-t2*t2;
	double b = 2.*t2*t3;
	double c = -m0*m0*t1*t1-t3*t3;
	double d = b*b - 4.*a*c;
	
	if(d < 0. || a == 0.) {
		cout << "IMAGINARY ETA ENERGY!!!" << endl;
		return 0.;
	}
	
	double energy = (-b-sqrt(d))/2./a;
	return energy;
}

double EtaAna::fcal_energy_res(double e)
{
	// hard-coded values for the FCAL energy resolution (taken from GlueX NIM paper)
	
	double a = 0.062, b = 0.047;
	double sig = (a*a)/e + (b*b);
	sig = sqrt(sig) * e;
	return sig;
}

int EtaAna::getPrimexPhase(int run_number) {
	
	int loc_phase = 0;
	if(run_number<60000) {
		return 0;
	} else if(run_number >=  60000 && run_number <=  69999) {
		return 1;
	} else if(run_number >=  80000 && run_number <=  89999) {
		return 2;
	} else if(run_number >= 110000 && run_number <= 119999) {
		return 3;
	} else {
		return 0;
	}
}

void EtaAna::setRunNumber(int runNum) { 
	
	m_runNumber = runNum;
	m_phase_val = getPrimexPhase(runNum);
	
	setGeometry();
	
	return;
}

void EtaAna::setGeometry() {
	
	// default target is proton:
	m_Target = m_Proton;
	
	switch(m_phase_val) {
		case 1:
		{
			if(m_runNumber<61355) {
				m_Target         = m_Be9;
				
				m_target_length  = 1.7755;
				m_target_density = 1.848;
				m_target_atten   = 0.01172;
				
				m_vertex.SetZ(64.935);
			} else {
				m_Target         = m_He4;
				
				m_target_length  = 29.535;
				m_target_density = 0.1217;
				m_target_atten   = 0.00821;
				
				m_vertex.SetZ(65.0);
			}
			
			// Correction to alignment after Simon updated beam spot:
			
			m_fcal_face.SetXYZ(0.455, -0.032, 624.906);
			m_fcal_correction.SetXYZ(0.455-0.529, -0.032+0.002, 0.0);
			
			if(m_runNumber<61483) {
				m_vertex.SetX( 0.027);
				m_vertex.SetY(-0.128);
			} else if(m_runNumber<61774) {
				m_vertex.SetX( 0.001);
				m_vertex.SetY(-0.077);
			} else {
				m_vertex.SetX( 0.038);
				m_vertex.SetY(-0.095);
			}
			
			break;
		}
		case 2:
		{
			if(m_runNumber<81396) {
				m_Target         = m_Be9;
				
				m_target_length  = 1.7755;
				m_target_density = 1.848;
				m_target_atten   = 0.01172;
				
				m_vertex.SetZ(64.935);
			} else {
				m_Target         = m_He4;
				
				m_target_length  = 29.535;
				m_target_density = 0.1217;
				m_target_atten   = 0.00821;
				
				m_vertex.SetZ(65.0);
			}
			
			m_fcal_face.SetXYZ(0.455, -0.032, 624.906);
			m_fcal_correction.SetXYZ(0.455-0.529, -0.032+0.002, 0.0);
			m_vertex.SetX( 0.129);
			m_vertex.SetY(-0.038);
			break;
		}
		case 3:
		{
			if(m_runNumber<110622) {
				m_Target         = m_Be9;
				
				m_target_length  = 1.7755;
				m_target_density = 1.848;
				m_target_atten   = 0.01172;
				
				m_vertex.SetZ(64.935);
			} else {
				m_Target         = m_He4;
				
				m_target_length  = 29.535;
				m_target_density = 0.1217;
				m_target_atten   = 0.00821;
				
				m_vertex.SetZ(65.0);
			}
			
			m_fcal_face.SetXYZ(0.189, 0.022, 624.32);
			m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
			m_vertex.SetX(-0.00439156);
			m_vertex.SetY(-0.0456728);
			break;
		}
		case 0:
		{
			m_Target = m_Proton;
			// use defaults for Helium target for now:
			m_target_length  = 29.535;
			m_target_density = 0.1217;
			m_target_atten   = 0.00821;
			
			m_fcal_face.SetXYZ(0.189, 0.022, 624.32);
			m_fcal_correction.SetXYZ(0.0, 0.0, 0.0);
			m_vertex.SetXYZ(0.0, 0.0, 65.0);
			break;
		}
	}
	
	return;
}

void EtaAna::setOutputFileName(string name) { m_output_fname = name; return; }

int EtaAna::loadTree() {
	
	m_tree = (TTree*)m_infile->Get("eta_gg");
	if(m_tree==NULL) return 0;
	
	// Reset event count to zero when laoding a new Tree:
	m_event = 0;
	
	return m_tree->GetEntries();
}

void EtaAna::readEvent() {
	
	if(m_event == 0) {
		
		// Set Branch Address on first Event:
		
		m_tree->SetBranchAddress("rfTime",             &m_rfTime);
		m_tree->SetBranchAddress("nbeam",              &m_nbeam);
		m_tree->SetBranchAddress("tag_counter",        &m_tag_counter);
		m_tree->SetBranchAddress("tag_sys",            &m_tag_sys);
		m_tree->SetBranchAddress("beam_e",             &m_beam_e);
		m_tree->SetBranchAddress("beam_t",             &m_beam_t);
		//m_tree->SetBranchAddress("acc_scale_factor",   &m_acc_scale_factor);
		m_tree->SetBranchAddress("nfcal",              &m_nfcal);
		m_tree->SetBranchAddress("fcal_e",             &m_fcal_e);
		m_tree->SetBranchAddress("fcal_x",             &m_fcal_x);
		m_tree->SetBranchAddress("fcal_y",             &m_fcal_y);
		m_tree->SetBranchAddress("fcal_z",             &m_fcal_z);
		m_tree->SetBranchAddress("fcal_t",             &m_fcal_t);
		m_tree->SetBranchAddress("fcal_nblocks",       &m_fcal_nblocks);
		m_tree->SetBranchAddress("nbcal",              &m_nbcal);
		m_tree->SetBranchAddress("bcal_e",             &m_bcal_e);
		m_tree->SetBranchAddress("bcal_x",             &m_bcal_x);
		m_tree->SetBranchAddress("bcal_y",             &m_bcal_y);
		m_tree->SetBranchAddress("bcal_z",             &m_bcal_z);
		m_tree->SetBranchAddress("bcal_t",             &m_bcal_t);
		m_tree->SetBranchAddress("ntof",               &m_ntof);
		m_tree->SetBranchAddress("tof_x",              &m_tof_x);
		m_tree->SetBranchAddress("tof_y",              &m_tof_y);
		m_tree->SetBranchAddress("tof_z",              &m_tof_z);
		m_tree->SetBranchAddress("tof_t",              &m_tof_t);
		m_tree->SetBranchAddress("nmc",                &m_nmc);
		m_tree->SetBranchAddress("mc_pdgtype",         &m_mc_pdgtype);
		m_tree->SetBranchAddress("mc_x",               &m_mc_x);
		m_tree->SetBranchAddress("mc_y",               &m_mc_y);
		m_tree->SetBranchAddress("mc_z",               &m_mc_z);
		m_tree->SetBranchAddress("mc_t",               &m_mc_t);
		m_tree->SetBranchAddress("mc_e",               &m_mc_e);
		m_tree->SetBranchAddress("mc_p",               &m_mc_p);
		m_tree->SetBranchAddress("mc_theta",           &m_mc_theta);
		m_tree->SetBranchAddress("mc_phi",             &m_mc_phi);
	}
	
	m_tree->GetEvent(m_event);
	
	return;
}

void EtaAna::initHistograms() {
	
	int n_reactions = (int)m_reaction_types.size();
	
	h_thrown   = new TH1F("thrown",   "Number Thrown; Reaction Type; N_{thrown}", 
		n_reactions+1, -0.5, ((double)n_reactions)+0.5);
	h_accepted = new TH1F("accepted", "Number Accepted; Reaction Type; N_{rec}", 
		n_reactions+1, -0.5, ((double)n_reactions)+0.5);
	
	for(int itype=0; itype<n_reactions; itype++) {
		TString bin_label = "";
		for(int iparticle=0; iparticle<m_reaction_types[itype].size(); iparticle++) {
			bin_label += ParticleName_ROOT(m_reaction_types[itype][iparticle]);
		}
		h_thrown->GetXaxis()->SetBinLabel(itype+1, bin_label.Data());
		h_accepted->GetXaxis()->SetBinLabel(itype+1, bin_label.Data());
	}
	h_thrown->GetXaxis()->SetBinLabel(h_thrown->GetNbinsX(), "other");
	h_accepted->GetXaxis()->SetBinLabel(h_accepted->GetNbinsX(), "other");
	
	// Plots specific to different reaction types:
	
	for(int itype=0; itype<(n_reactions+1); itype++) {
		
		TString hist_label = "", hist_title = "";
		if(itype==n_reactions) {
			hist_label = "other";
			hist_title = "other";
		}
		else {
			for(int iparticle=0; iparticle<m_reaction_types[itype].size(); iparticle++) {
				hist_label += ShortName(m_reaction_types[itype][iparticle]);
				hist_title += ParticleName_ROOT(m_reaction_types[itype][iparticle]);
			}
		}
		
		// angular yields (with and without BCAL veto):
		
		TH1F *loc_h_theta      = new TH1F(Form("theta_%s",hist_label.Data()), 
			Form("Angular Yield (%s)",hist_title.Data()), 65, 0.0, 6.5);
		h_theta.push_back(loc_h_theta);
		
		vector<TH1F*> loc_h_theta_veto_vec;
		for(int iveto=0; iveto<m_n_bcal_vetos; iveto++) {
			TH1F *loc_h_theta_veto = new TH1F(Form("theta_veto_%d_%s",iveto,hist_label.Data()),
				Form("Angular Yield with BCAL Veto (option %d) (%s)",iveto,hist_title.Data()), 65, 0.0, 6.5);
			loc_h_theta_veto_vec.push_back(loc_h_theta_veto);
		}
		h_theta_veto.push_back(loc_h_theta_veto_vec);
		
		// BCAL shower distributions:
		
		TH1F *loc_h_nbcal = new TH1F(Form("nbcal_%s",hist_label.Data()), 
			Form("Number of BCAL Showers (%s)",hist_title.Data()), 20, -0.5, 19.5);
		
		TH1F *loc_h_bcal_energy = new TH1F(Form("bcal_energy_sum_%s",hist_label.Data()),
			Form("BCAL Energy Sum (%s)",hist_title.Data()), 500, 0.0, 5.0);
		
		TH1F *loc_h_bcal_energy_single = new TH1F(Form("bcal_energy_sum_single_%s",hist_label.Data()),
			Form("BCAL Energy Sum, 1-Shower (%s)",hist_title.Data()), 500, 0.0, 5.0);
		
		TH2F *loc_h_bcal_dt_vs_eta_angle = new TH2F(Form("bcal_dt_vs_eta_angle_%s",hist_label.Data()),
			Form("BCAL-RF Time vs #eta Angle (%s); #theta_{rec} [deg.]; t_{BCAL}-t_{RF} [ns]",hist_title.Data()),
			100, 0.0, 10.0, 400, -20.0, 20.0);
		
		TH2F *loc_h_bcal_dt_vs_bcal_energy = new TH2F(Form("bcal_dt_vs_bcal_energy_%s",hist_label.Data()),
			Form("BCAL-RF Time vs Shower Energy (%s); E_{BCAL} [GeV]; t_{BCAL}-t_{RF} [ns]",hist_title.Data()),
			100, 0.0, 1.0, 400, -20.0, 20.0);
		
		TH1F *loc_h_bcal_deltaPhi = new TH1F(Form("bcal_deltaPhi_%s",hist_label.Data()),
			Form("#left|#phi_{#gamma#gamma} - #phi_{BCAL}#right| (%s); #Delta#phi_{#eta-BCAL} [deg.]",hist_title.Data()),
			360, 0., 360.0);
		
		h_nbcal.push_back(loc_h_nbcal);
		h_bcal_energy.push_back(loc_h_bcal_energy);
		h_bcal_energy_single.push_back(loc_h_bcal_energy_single);
		h_bcal_dt_vs_eta_angle.push_back(loc_h_bcal_dt_vs_eta_angle);
		h_bcal_dt_vs_bcal_energy.push_back(loc_h_bcal_dt_vs_bcal_energy);
		h_bcal_deltaPhi.push_back(loc_h_bcal_deltaPhi);
	}
	
	
	// Elasticity vs. mass ratio:
	h_elas_vs_mgg = new TH2F("elas_vs_mgg", 
		"; M_{#gamma#gamma}/M_{#eta}(PDG); #left(E_{1}+E_{2}#right)/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)", 
		1000, 0., 2., 1000, 0., 2.);
	
	// Elasticity with tagged photon:
	h_elas             = new TH2F("elas", 
		"Elasticity; E_{#gamma#gamma}/E_{#gamma}", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_signal      = new TH2F("elas_signal", 
		"Elasticity; E_{#gamma#gamma}/E_{#gamma}", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_bkgd        = new TH2F("elas_bkgd", 
		"Elasticity; E_{#gamma#gamma}/E_{#gamma}", 650, 0., 6.5, 1000, 0., 2.);
	
	// Elasticity with coherently-produced eta:
	h_elas_corr        = new TH2F("elas_corr", 
		"Elasticity; E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_corr_signal = new TH2F("elas_corr_signal", 
		"Elasticity; E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)", 650, 0., 6.5, 1000, 0., 2.);
	h_elas_corr_bkgd   = new TH2F("elas_corr_bkgd", 
		"Elasticity; E_{#gamma#gamma}/E_{#eta}#left(E_{#gamma},#theta_{#gamma#gamma}#right)", 650, 0., 6.5, 1000, 0., 2.);
	
	// 2-photon invariant mass vs. angle:
	h_mgg = new TH2F("mgg", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg->Sumw2();
	
	h_mgg_signal = new TH2F("mgg_signal", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_signal->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_signal->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg_signal->Sumw2();
	
	h_mgg_bkgd = new TH2F("mgg_bkgd", "Two-Photon Invariant Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_bkgd->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_bkgd->GetYaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
	h_mgg_bkgd->Sumw2();
	
	// Energy-constrained mass vs. angle:
	h_mgg_const = new TH2F("mgg_const", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_const->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const->Sumw2();
	
	h_mgg_const_signal = new TH2F("mgg_const_signal", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_signal->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_const_signal->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const_signal->Sumw2();
	
	h_mgg_const_bkgd = new TH2F("mgg_const_bkgd", "Energy-Constrained Inv Mass", 650, 0., 6.5, 600, 0., 1.2);
	h_mgg_const_bkgd->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mgg_const_bkgd->GetYaxis()->SetTitle("M_{#gamma#gamma}^{constr} [GeV/c^{2}]");
	h_mgg_const_bkgd->Sumw2();
	
	// Missing-mass vs. angle:
	h_mm_vs_theta = new TH2F("mm_vs_theta", "Squared Missing Mass", 
		650, 0., 6.5, 4000, 0., 40.);
	h_mm_vs_theta->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_vs_theta->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_vs_theta->Sumw2();
	
	h_mm_vs_theta_eta_cut = new TH2F("mm_vs_theta_eta_cut", "Squared Missing Mass (m_{#gamma#gamma} cut)", 
		650, 0., 6.5, 4000, 0., 40.);
	h_mm_vs_theta_eta_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_vs_theta_eta_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_vs_theta_eta_cut->Sumw2();
	
	h_mm_vs_theta_eta_elas_cut = new TH2F("mm_vs_theta_eta_elas_cut", 
		"Squared Missing Mass (m_{#gamma#gamma} + elasticity cut)", 650, 0., 6.5, 4000, 0., 40.);
	h_mm_vs_theta_eta_elas_cut->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [#circ]");
	h_mm_vs_theta_eta_elas_cut->GetYaxis()->SetTitle("#Deltam^{2} [GeV^{2}/c^{4}]");
	h_mm_vs_theta_eta_elas_cut->Sumw2();
	
	// x-y position of FCAL showers that survive all cuts:
	h_xy_1 = new TH2F("xy_1", "Postion of Shower 1; x_{1} [cm]; y_{1} [cm]", 500, -100., 100., 500, -100., 100.);
	h_xy_2 = new TH2F("xy_2", "Postion of Shower 2; x_{2} [cm]; y_{2} [cm]", 500, -100., 100., 500, -100., 100.);
	
	return;
}

void EtaAna::resetHistograms() {
	
	h_thrown->Reset();
	h_accepted->Reset();
	
	for(int itype=0; itype<h_theta.size(); itype++) {
		h_theta[itype]->Reset();
		for(int iveto=0; iveto<m_n_bcal_vetos; iveto++) h_theta_veto[itype][iveto]->Reset();
		h_nbcal[itype]->Reset();
		h_bcal_energy[itype]->Reset();
		h_bcal_energy_single[itype]->Reset();
		h_bcal_dt_vs_eta_angle[itype]->Reset();
		h_bcal_dt_vs_bcal_energy[itype]->Reset();
		h_bcal_deltaPhi[itype]->Reset();
	}
	
	h_elas_vs_mgg->Reset();
	
	h_elas->Reset();
	h_elas_corr->Reset();
	h_elas_signal->Reset();
	h_elas_corr_signal->Reset();
	h_elas_bkgd->Reset();
	h_elas_corr_bkgd->Reset();
	
	h_mgg->Reset();
	h_mgg_const->Reset();
	h_mgg_signal->Reset();
	h_mgg_const_signal->Reset();
	h_mgg_bkgd->Reset();
	h_mgg_const_bkgd->Reset();
	
	h_mm_vs_theta->Reset();
	h_mm_vs_theta_eta_cut->Reset();
	h_mm_vs_theta_eta_elas_cut->Reset();
	h_xy_1->Reset();
	h_xy_2->Reset();
	
	return;
}

void EtaAna::writeHistograms() {
	
	cout << "writing histograms to: " << m_output_fname << "..." << endl;
	TFile *fOut = new TFile(m_output_fname.c_str(), "RECREATE");
	fOut->cd();
	
	h_thrown->Write();
	h_accepted->Write();
	
	for(int itype=0; itype<h_theta.size(); itype++) {
		
		TString dir_label = "other";
		if(itype<m_reaction_types.size()) {
			dir_label = "";
			for(int iparticle=0; iparticle<m_reaction_types[itype].size(); iparticle++) {
				dir_label += ShortName(m_reaction_types[itype][iparticle]);
			}
		}
		
		TDirectory *loc_dir = new TDirectoryFile(dir_label.Data(), dir_label.Data());
		loc_dir->cd();
		h_theta[itype]->Write();
		for(int iveto=0; iveto<m_n_bcal_vetos; iveto++) h_theta_veto[itype][iveto]->Write();
		h_nbcal[itype]->Write();
		h_bcal_energy[itype]->Write();
		h_bcal_energy_single[itype]->Write();
		h_bcal_dt_vs_eta_angle[itype]->Write();
		h_bcal_dt_vs_bcal_energy[itype]->Write();
		h_bcal_deltaPhi[itype]->Write();
		loc_dir->cd("../");
	}
	
	h_elas_vs_mgg->Write();
	
	h_elas->Write();
	h_elas_corr->Write();
	h_elas_signal->Write();
	h_elas_corr_signal->Write();
	h_elas_bkgd->Write();
	h_elas_corr_bkgd->Write();
	
	h_mgg->Write();
	h_mgg_const->Write();
	h_mgg_signal->Write();
	h_mgg_const_signal->Write();
	h_mgg_bkgd->Write();
	h_mgg_const_bkgd->Write();
	
	h_mm_vs_theta->Write();
	h_mm_vs_theta_eta_cut->Write();
	h_mm_vs_theta_eta_elas_cut->Write();
	h_xy_1->Write();
	h_xy_2->Write();
	
	fOut->Write();
	
	return;
}
