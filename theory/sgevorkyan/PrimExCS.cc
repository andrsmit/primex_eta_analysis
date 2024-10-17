#include "PrimExCS.h"

PrimExCS::PrimExCS()
{
	// Set default parameters for eta produced off Helium:
	
	setMeson(Eta);
	setTarget(Helium);
	return;
}

PrimExCS::PrimExCS(Particle_t meson, Particle_t target)
{
	setMeson(meson);
	setTarget(target);
	return;
}

void PrimExCS::initialize() {
	
	m_angularBins.clear();
	m_CoulombFF.clear();
	m_StrongFF.clear();
	
	double locTheta = m_minTheta + 0.5*m_thetaBinSize;
	int nThetaBins = 0;
	double locMaxTheta = m_minTheta;
	while(locTheta < m_maxTheta) {
		
		nThetaBins++;
		m_angularBins.push_back(locTheta);
		
		// Coulomb Form Factor:
		m_CoulombFF.push_back({0.,0.});
		
		// Strong Form Factor:
		m_StrongFF.push_back({0.,0.});
		
		locTheta += m_thetaBinSize;
		locMaxTheta = locTheta + 0.5*m_thetaBinSize;
	}
	m_maxTheta = locMaxTheta;
	
	TString m_proc_name[4]  = {"prim", "coh", "inter", "incoh"};
	int loc_colors[4] = {kRed, kBlue, kMagenta, kGreen};
	
	for(int iproc=0; iproc<4; iproc++) {
		h_CrossSection[iproc] = new TH1F(Form("h_%s_xs",m_proc_name[iproc].Data()), 
			"; #theta [#circ]; d#sigma/d#theta [#mub/rad]", nThetaBins, m_minTheta, m_maxTheta);
		StyleHistogram(h_CrossSection[iproc]);
		
		h_CrossSection[iproc]->SetLineColor(loc_colors[iproc]);
		h_CrossSection[iproc]->SetMarkerColor(loc_colors[iproc]);
	}
	
	h_CoulombFF[0] = new TH1F("CoulombFF_real", 
		";#theta [#circ];Re#left(F_{em}#right)", nThetaBins, m_minTheta, m_maxTheta);
	StyleHistogram(h_CoulombFF[0]);
	h_CoulombFF[0]->SetLineColor(kBlue);
	h_CoulombFF[0]->SetMarkerColor(kBlue);
	
	h_CoulombFF[1] = new TH1F("CoulombFF_imag", 
		";#theta [#circ];Im#left(F_{em}#right)", nThetaBins, m_minTheta, m_maxTheta);
	StyleHistogram(h_CoulombFF[1]);
	h_CoulombFF[1]->SetLineColor(kRed);
	h_CoulombFF[1]->SetMarkerColor(kRed);
	
	h_StrongFF[0] = new TH1F("StrongFF_real", 
		";#theta [#circ];Re#left(F_{st}#right)", nThetaBins, m_minTheta, m_maxTheta);
	StyleHistogram(h_StrongFF[0]);
	h_StrongFF[0]->SetLineColor(kBlue);
	h_StrongFF[0]->SetMarkerColor(kBlue);
	
	h_StrongFF[1] = new TH1F("StrongFF_imag", 
		";#theta [#circ];Im#left(F_{st}#right)", nThetaBins, m_minTheta, m_maxTheta);
	StyleHistogram(h_StrongFF[1]);
	h_StrongFF[1]->SetLineColor(kRed);
	h_StrongFF[1]->SetMarkerColor(kRed);
	
	return;
}

void PrimExCS::setTarget(Particle_t particle)
{
	m_target = particle;
	
	// Set the target mass, charge (Z), mass number (A), and nuclear radius (a0):
	m_targetMass   = ParticleMass(particle);
	m_targetZ      = static_cast<double>(ParticleCharge(particle));
	
	if(particle==Helium) {
		m_targetA = 4.0;
	}
	else if(particle==C12) {
		m_targetA = 12.0;
	}
	else {
		// assume A = 2Z:
		m_targetA = 2.0*m_targetZ;
		std::cout << "Unsupported target provided. Setting A = 2Z." << std::endl;
	}
	
	// Print settings to terminal:
	//printf("Target: %s\n  mass   = %f GeV\n  Z      = %d\n  A      = %d\n",
	//	ParticleType(particle), m_targetMass, (int)m_targetZ, (int)m_targetA);
}

//------------------------------------------------------------------------------------//

void PrimExCS::CalculateCrossSection() {
	
	/* 
	Set the bin contents of the supplied histogram equal to the photoproduction cross section.
	Assumes the x-axis gives the production angle in units of degrees.
	*/
	
	int nThetaBins = (int)m_CoulombFF.size();
	for(int iThetaBin=0; iThetaBin<nThetaBins; iThetaBin++) {
		
		double locTheta = h_CrossSection[0]->GetXaxis()->GetBinCenter(iThetaBin+1);
		
		// Primakoff
		
		double  ffCoulomb = norm(m_CoulombFF[iThetaBin]);
		double ampCoulomb = PrimakoffAmplitude(locTheta);
		double  csCoulomb = pow(ampCoulomb,2.0) * ffCoulomb;
		
		// The amplitude-squared has units of GeV^-2. Convert to mb:
		csCoulomb /= m_GeV2mb;
		
		// the above cross section is in mb/sr. Convert to ub/rad:
		csCoulomb *= (1.e3 * sin(locTheta*TMath::DegToRad()) * 2.0 * TMath::Pi());
		
		h_CrossSection[0]->SetBinContent(iThetaBin+1, csCoulomb);
		
		// Strong Coherent:
		
		double  ffStrong = norm(m_StrongFF[iThetaBin]);
		double ampStrong = NuclearCoherentAmplitude(locTheta);
		double  csStrong = pow(ampStrong,2.0) * ffStrong;
		
		// the above cross section is in mb/sr. Convert to ub/rad:
		csStrong *= (1.e3 * sin(locTheta*TMath::DegToRad()) * 2.0 * TMath::Pi());
		
		h_CrossSection[1]->SetBinContent(iThetaBin+1, csStrong);
		
		// Interference between Primakoff and Strong:
		
		double Rem = real(m_CoulombFF[iThetaBin]);
		double Iem = imag(m_CoulombFF[iThetaBin]);
		double Rst = real(m_StrongFF[iThetaBin]);
		double Ist = imag(m_StrongFF[iThetaBin]);
		
		double ff1 = (Rem*Rst + Iem*Ist) * cos(m_cohPhaseAngle*TMath::DegToRad());
		double ff2 = (Rst*Iem - Rem*Ist) * sin(m_cohPhaseAngle*TMath::DegToRad());
		
		double csInterference = ampCoulomb*ampStrong*(ff1 + ff2);
		
		// Convert to ub/rad:
		csInterference /= m_GeV2mb;
		csInterference *= (1.e3 * sin(locTheta*TMath::DegToRad()) * 2.0 * TMath::Pi());
		
		h_CrossSection[2]->SetBinContent(iThetaBin+1, csInterference);
	}
	
	return;
}

double PrimExCS::PrimakoffAmplitude(double angle) {
	
	// Eq.2 from PrimEx Note #85
	
	double mesonMom    = getMesonMomentum(angle);
	double mesonEnergy = sqrt(pow(mesonMom,2.0) + pow(m_mesonMass,2.0));
	double mesonBeta   = mesonMom / mesonEnergy;
	
	/* 
	Mandelstam t:
	double costh = cos(angle*TMath::DegToRad());
	double t = pow(m_mesonMass,2.0) - 2.0*m_beamEnergy*(mesonEnergy - mesonMom*costh);
	double t = 2.0*m_targetMass*(mesonEnergy - m_beamEnergy);
	*/
	double costh = cos(angle*TMath::DegToRad());
	double t = pow(m_mesonMass,2.0) - 2.0*m_beamEnergy*(mesonEnergy - mesonMom*costh);
	
	double sinth = sin(angle*TMath::DegToRad());
	
	double T_prim = m_targetZ * sqrt(8.0*m_alpha*m_RadDecayWidth*m_eV2GeV) * pow(mesonBeta/m_mesonMass,1.5) 
		* (pow(m_beamEnergy,2.0)*sinth / fabs(t));
	
	return T_prim;
}

double PrimExCS::NuclearCoherentAmplitude(double angle) {
	
	// Copied from ds_eta_he4.F in halld_sim/src/programs/Simulation/gen_primex_eta_he4/ds_eta_he4
	
	double T_nc = m_targetA * m_beamEnergy * sqrt(16.0*1.e-3) * sin(angle*TMath::DegToRad());
	return T_nc;
}

double PrimExCS::getMesonMomentum(double angle) {
	
	// Calculates the momentum of the meson from the beam energy and production angle
	// (assumes a coherent production process)
	
	double costh = cos(angle*TMath::Pi()/180.0);
	double f = 0.5*(pow(m_mesonMass,2.0) + 2.0*m_targetMass*m_beamEnergy);
	
	double A = pow(m_beamEnergy*costh,2.0) - pow(m_beamEnergy+m_targetMass,2.0);
	double B = 2.0*f*m_beamEnergy*costh;
	double C = pow(f,2.0) - pow(m_mesonMass*(m_beamEnergy+m_targetMass),2.0);
	
	double D = pow(B,2.0) - 4.0*A*C;
	if(D < 0.0) {
		return -1.0;
	}
	if(A == 0.0) {
		return 0.0;
	}
	double mom = (-B - sqrt(D)) / (2.0*A);
	
	return mom;
}

void PrimExCS::PlotCoulombFF() {
	
	int nThetaBins = (int)m_CoulombFF.size();
	if(nThetaBins==0) return;
	
	for(int itbin=1; itbin<=nThetaBins; itbin++) {
		h_CoulombFF[0]->SetBinContent(itbin, real(m_CoulombFF[itbin-1]));
		h_CoulombFF[1]->SetBinContent(itbin, imag(m_CoulombFF[itbin-1]));
	}
	
	h_CoulombFF[0]->GetYaxis()->SetRangeUser(-0.25,1.25);
	h_CoulombFF[1]->GetYaxis()->SetRangeUser(-0.25,1.25);
	
	TCanvas *c_CoulombFF = new TCanvas("c_CoulombFF","CoulombFF",900,700);
	StyleCanvas(c_CoulombFF);
	h_CoulombFF[0]->Draw("hist");
	h_CoulombFF[1]->Draw("hist same");
	c_CoulombFF->Update();
	c_CoulombFF->Modified();
	
	return;
}

void PrimExCS::PlotStrongFF() {
	
	int nThetaBins = (int)m_StrongFF.size();
	if(nThetaBins==0) return;
	
	for(int itbin=1; itbin<=nThetaBins; itbin++) {
		h_StrongFF[0]->SetBinContent(itbin, real(m_StrongFF[itbin-1]));
		h_StrongFF[1]->SetBinContent(itbin, imag(m_StrongFF[itbin-1]));
	}
	
	h_StrongFF[0]->GetYaxis()->SetRangeUser(-0.25,1.25);
	h_StrongFF[1]->GetYaxis()->SetRangeUser(-0.25,1.25);
	
	TCanvas *c_StrongFF = new TCanvas("c_StrongFF", "StrongFF", 900, 700);
	StyleCanvas(c_StrongFF);
	h_StrongFF[0]->Draw("hist");
	h_StrongFF[1]->Draw("hist same");
	c_StrongFF->Update();
	c_StrongFF->Modified();
	
	return;
}

void PrimExCS::PlotCrossSection() {
	
	TCanvas *cCS = new TCanvas("cCS", "CS", 900, 700);
	StyleCanvas(cCS);
	
	double loc_max = h_CrossSection[0]->GetMaximum();
	for(int ihist=1; ihist<4; ihist++) {
		if(h_CrossSection[ihist]->GetMaximum() > loc_max) {
			loc_max = h_CrossSection[ihist]->GetMaximum();
		}
	}
	for(int ihist=0; ihist<4; ihist++) h_CrossSection[ihist]->GetYaxis()->SetRangeUser(0.0, 1.15*loc_max);
	
	h_CrossSection[0]->Draw("hist");
	for(int ihist=1; ihist<4; ihist++) 
		h_CrossSection[ihist]->Draw("same hist");
	
	cCS->Update();
	cCS->Modified();
	
	return;
}

void PrimExCS::WriteToROOTFile() {
	
	TFile *fOut = new TFile(m_outputFileName.Data(), "RECREATE");
	fOut->cd();
	
	// write out all cross section histograms:
	for(int ihist=0; ihist<4; ihist++) h_CrossSection[ihist]->Write();
	
	// if form factor histograms were created, write them out too:
	
	if(h_CoulombFF[0]!=NULL) h_CoulombFF[0]->Write();
	if(h_CoulombFF[1]!=NULL) h_CoulombFF[1]->Write();
	
	if(h_StrongFF[0]!=NULL) h_StrongFF[0]->Write();
	if(h_StrongFF[1]!=NULL) h_StrongFF[1]->Write();
	
	fOut->Write();
	fOut->Close();
	
	return;
}

void PrimExCS::StyleHistogram(TH1F *h) 
{
	h->SetLineWidth(2);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.0);
	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.0);
	h->GetYaxis()->CenterTitle(true);
	h->SetMarkerStyle(8);
	h->SetMarkerSize(0.5);
	return;
}

void PrimExCS::StyleCanvas(TCanvas *c) 
{
	c->SetTickx();
	c->SetTicky();
	c->SetLeftMargin(0.13); c->SetRightMargin(0.07);
	c->SetBottomMargin(0.13); c->SetTopMargin(0.07);
	c->cd();
	return;
}
