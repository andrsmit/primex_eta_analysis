#ifndef _PRIMEX_ETA_CS_
#define _PRIMEX_ETA_CS_

#include "particleType.h"

#include <math.h>
#include <complex>
#include <iostream>

#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"

using namespace std;

class PrimExCS {
	public:
		PrimExCS();
		PrimExCS(Particle_t meson, Particle_t target);
		
		// Physical Constants:
		static constexpr double m_eV2GeV   = 1.e-9;
		static constexpr double m_alpha    = 7.2973525643e-03; // fine-structure constant
		static constexpr double m_GeV2fm   = 0.1973269603;     // 1 fm = 0.197 GeV-1
		static constexpr double m_GeV2mb   = 2.56819;          // 1 mb = 2.568 GeV-2
		static constexpr double m_rhoMass  = 0.7690;           // mass of rho0 meson
		static constexpr double m_protonRadius = sqrt(2.0/3.0)*0.842;
		
		vector<std::complex<double>> m_CoulombFF, m_StrongFF;
		vector<double> m_angularBins;
		
		void setCohPhaseAngle(double theta) { m_cohPhaseAngle = theta; };
		void setOutputFileName(TString fileName) { m_outputFileName = fileName; };
		
		//----------------------------------------------------------------//
		// Get/Set meson and target properties using Particle_t structure:
		
		void setMeson(Particle_t meson) { m_meson = meson; m_mesonMass = ParticleMass(meson); };
		double getMesonMass() { return m_mesonMass; };
		char*  getMesonName() { return ParticleType(m_meson); };
		
		void setTarget(Particle_t);
		char*  getTargetName() { return ParticleType(m_target); };
		double getTargetZ() { return m_targetZ; };
		double getTargetA() { return m_targetA; };
		double getTargetMass() { return m_targetMass; };
		
		//----------------------------------------------------------------//
		// Get/Set configuration:
		
		void setAngularRange(double minimum, double maximum) {
			m_minTheta = minimum;
			m_maxTheta = maximum;
		};
		void setAngularBinSize(double dtheta) { m_thetaBinSize = dtheta; };
		
		void setBeamEnergy(double energy) { m_beamEnergy = energy; }
		double getBeamEnergy() { return m_beamEnergy; }
		
		//----------------------------------------------------------------//
		
		void initialize();
		
		double getMesonMomentum(double angle);
		
		void CalculateCrossSection();
		
		void PlotCrossSection();
		void PlotCoulombFF();
		void PlotStrongFF();
		void WriteToROOTFile();
	
	private:
		
		double m_RadDecayWidth = 515.0;
		double m_cohPhaseAngle = TMath::RadToDeg();
		
		TH1F *h_CrossSection[4];
		TH1F *h_CoulombFF[2], *h_StrongFF[2];
		
		TString m_outputFileName;
		
		// Parameters Unique to Meson:
		
		Particle_t m_meson;
		double m_mesonMass;
		
		// Parameters Unique to Target:
		
		Particle_t m_target;
		double m_targetMass;
		double m_targetZ;      // atomic number of target
		double m_targetA;      // Number of nucleons in target
		
		// Configuration of cross section calculation:
		
		double m_minTheta = 0.0, m_maxTheta = 10.0;
		double m_thetaBinSize = 0.001;
		
		double m_beamEnergy = 10.0;
		
		double PrimakoffAmplitude(double angle);
		double NuclearCoherentAmplitude(double angle);
		
		void StyleHistogram(TH1F*);
		void StyleCanvas(TCanvas*);
};

#endif
