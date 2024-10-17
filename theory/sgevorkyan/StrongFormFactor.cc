#include "PrimExCS.h"
#include "FormFactor.h"

double StrongFFIntegrationRange_b;
double StrongFFIntegrationRange_z1, StrongFFIntegrationRange_z2;
double StrongFFIntegrationRange_s;

TF2  *f_StrongFFIntegrand_IS_2;
TH2F *h_shadowingMatrix[2];

double StrongFFIntegrand(double *x, double *par) {
	
	double b = x[0];
	double z = x[1];
	double s = x[2];
	
	double isReal = par[0]; // positive for real, negative for imaginary
	double q      = par[1];
	double Delta  = par[2];
	
	double kinFactor = 2.0*TMath::Pi() / q / pow(modelParameters.slopeAp,2.0); // * exp(pow(0.5*targetRadius*q,2.0)/targetA);
	double besJ1     = b*ROOT::Math::cyl_bessel_j(1.0, q*b);
	
	std::complex<double> fAbs = IntegrateAbsorption(b, z);
	
	double besDiff = b*ROOT::Math::cyl_bessel_i(0.0, b*s/modelParameters.slopeAp) 
		- s*ROOT::Math::cyl_bessel_i(1.0, b*s/modelParameters.slopeAp);
	
	double density = GetNuclearDensity(sqrt(pow(s,2.0)+pow(z,2.0)));
	
	double S_b_z_s = s*exp(-(pow(b,2.0)+pow(s,2.0))/(2.0*modelParameters.slopeAp)) * besDiff * density;
	
	double Z_b_z = 0.0;
	if(isReal>=0.0) Z_b_z = real(fAbs)*cos(Delta*z) - imag(fAbs)*sin(Delta*z);
	else            Z_b_z = imag(fAbs)*cos(Delta*z) + real(fAbs)*sin(Delta*z);
	
	double locFF = kinFactor * besJ1 * Z_b_z * S_b_z_s;
	return locFF;
}

double StrongFFIntegrand_IS_1(double *x, double *par) {
	
	double b = x[0];
	double z = x[1];
	double s = x[2];
	
	double isReal    = par[0]; // positive for real, negative for imaginary
	double q         = par[1];
	double Delta     = par[2];
	double Delta_rho = par[3];
	
	double kinFactor = (modelParameters.targetA-1.0) * TMath::Pi() * 
		(2.0*modelParameters.sigmaby2) / q / modelParameters.slopeAs / pow(modelParameters.slopeAp,2.0);
	
	double besJ1     = b*ROOT::Math::cyl_bessel_j(1.0, q*b);
	
	double besDiff = b*ROOT::Math::cyl_bessel_i(0.0, b*s/modelParameters.slopeAp) 
		- s*ROOT::Math::cyl_bessel_i(1.0, b*s/modelParameters.slopeAp);
	
	double density = GetNuclearDensity(sqrt(pow(s,2.0)+pow(z,2.0)));
	double S_b_z_s = s*exp(-(pow(b,2.0)+pow(s,2.0))/(2.0*modelParameters.slopeAp)) * besDiff * density;
	
	std::complex<double> fAbs      = IntegrateAbsorption(b, z, 2.0);
	std::complex<double> trig_func = exp(ii*Delta_rho*z);
	std::complex<double> integral2 = IntegrateInitialState(b, z, q, Delta, Delta_rho);
	
	std::complex<double> locFF = kinFactor * besJ1 * S_b_z_s * fAbs * trig_func * integral2;
	if(isReal>=0.0) return real(locFF);
	else            return imag(locFF);
}

double StrongFFIntegrand_IS_2(double *x, double *par) {
	
	double z2 = x[0];
	double s2 = x[1];
	
	double isReal    = par[0];
	double q         = par[1];
	double Delta     = par[2];
	double Delta_rho = par[3];
	double b         = par[4];
	
	double density   = GetNuclearDensity(sqrt(pow(s2,2.0)+pow(z2,2.0)));
	
	std::complex<double> trig_func = exp(ii*(Delta-Delta_rho)*z2);
	std::complex<double> S_b_z_s   = trig_func *
		s2 * exp(-(pow(b,2.0)+pow(s2,2.0))/(2.0*modelParameters.slopeAs)) * 
		ROOT::Math::cyl_bessel_i(0.0, b*s2/modelParameters.slopeAs) * density;
	
	if(isReal>=0.0) return real(S_b_z_s);
	else            return imag(S_b_z_s);
}

std::complex<double> IntegrateInitialState(double b, double z, double q, double Delta, double Delta_rho) {
	
	if((h_shadowingMatrix[0] == NULL) || (h_shadowingMatrix[1] == NULL)) InitializeShadowingMatrix(q, Delta, Delta_rho);
	
	double realIntegral = Interpolate2D(h_shadowingMatrix[0], b, z);
	double imagIntegral = Interpolate2D(h_shadowingMatrix[1], b, z);
	
	/*
	// Real part of the integration:
	f_StrongFFIntegrand_IS_2->SetParameters(1.0, q, Delta, Delta_rho, b);
	double realIntegral = f_StrongFFIntegrand_IS_2->Integral(z, StrongFFIntegrationRange_z2, 0.0, StrongFFIntegrationRange_s);
	
	// Imaginary part:
	f_StrongFFIntegrand_IS_2->SetParameters(-1.0, q, Delta, Delta_rho, b);
	double imagIntegral = f_StrongFFIntegrand_IS_2->Integral(z, StrongFFIntegrationRange_z2, 0.0, StrongFFIntegrationRange_s);
	*/
	return {realIntegral, imagIntegral};
}

void InitializeShadowingMatrix(double q, double Delta, double Delta_rho) {
	
	for(int ihist=0; ihist<2; ihist++) {
		h_shadowingMatrix[ihist] = new TH2F(Form("shadowingMatrix_%d",ihist), ";b [fm];z [fm]", 500, 0.0, StrongFFIntegrationRange_b, 
			500, StrongFFIntegrationRange_z1, StrongFFIntegrationRange_z2);
	}
	
	f_StrongFFIntegrand_IS_2->SetParameter(1, q);
	f_StrongFFIntegrand_IS_2->SetParameter(2, Delta);
	f_StrongFFIntegrand_IS_2->SetParameter(3, Delta_rho);
	
	std::cout << "Calculating shadowing matrix..." << std::flush;
	
	for(int xbin=1; xbin<=h_shadowingMatrix[0]->GetXaxis()->GetNbins(); xbin++) {
		double b = h_shadowingMatrix[0]->GetXaxis()->GetBinCenter(xbin);
		printf("\nb=%f",b);
		f_StrongFFIntegrand_IS_2->SetParameter(4, b);
		for(int ybin=1; ybin<=h_shadowingMatrix[0]->GetYaxis()->GetNbins(); ybin++) {
			double z = h_shadowingMatrix[0]->GetYaxis()->GetBinCenter(ybin);
			
			f_StrongFFIntegrand_IS_2->SetParameter(0, 1.0);
			double locReIntegral = f_StrongFFIntegrand_IS_2->Integral(z, StrongFFIntegrationRange_z2, 0.0, StrongFFIntegrationRange_s);
			
			f_StrongFFIntegrand_IS_2->SetParameter(0,-1.0);
			double locImIntegral = f_StrongFFIntegrand_IS_2->Integral(z, StrongFFIntegrationRange_z2, 0.0, StrongFFIntegrationRange_s);
			
			h_shadowingMatrix[0]->SetBinContent(xbin, ybin, locReIntegral);
			h_shadowingMatrix[1]->SetBinContent(xbin, ybin, locImIntegral);
		}
	}
	
	std::cout << "done." << std::endl;
	
	return;
}

//-----------------------------------------------------------------------//

void CalculateStrongFF(PrimExCS &csObj) {
	
	StrongFFIntegrationRange_b  =   7.0;
	StrongFFIntegrationRange_z1 =  -5.0;
	StrongFFIntegrationRange_z2 =   5.0;
	StrongFFIntegrationRange_s  =  10.0;
	
	double xmin[] = {0.0,                        StrongFFIntegrationRange_z1, 0.0                       };
	double xmax[] = {StrongFFIntegrationRange_b, StrongFFIntegrationRange_z2, StrongFFIntegrationRange_s};
	ROOT::Math::AdaptiveIntegratorMultiDim locIntegrator;
	
	//----------------------------------------//
	// Initialize the function which will be integrated over the variable b:
	
	TF3 *f_StrongFFIntegrand = new TF3("StrongFFIntegrand", StrongFFIntegrand, 
		0.0, StrongFFIntegrationRange_b, StrongFFIntegrationRange_z1, StrongFFIntegrationRange_z2, 
		0.0, StrongFFIntegrationRange_s, 3);
	
	//----------------------------------------//
	// Initialize the function which will be integrated over the variables b, z1, s1:
	
	TF2 *f_StrongFFIntegrand_IS_1 = new TF3("StrongFFIntegrand_IS_1", StrongFFIntegrand_IS_1, 
		0.0, StrongFFIntegrationRange_b, StrongFFIntegrationRange_z1, StrongFFIntegrationRange_z2, 
		0.0, StrongFFIntegrationRange_s, 4);
	
	//----------------------------------------//
	// Initialize the function which will be integrated over the variables z2, s2:
	
	f_StrongFFIntegrand_IS_2 = new TF2("StrongFFIntegrand_2", StrongFFIntegrand_IS_2, 
		StrongFFIntegrationRange_z1, StrongFFIntegrationRange_z2, 0.0, StrongFFIntegrationRange_s, 5);
	
	double DeltaRho = 0.5 * pow(PrimExCS::m_rhoMass,2.0) / csObj.getBeamEnergy() / PrimExCS::m_GeV2fm;
	
	f_StrongFFIntegrand_IS_1->SetParameter(3, DeltaRho);
	f_StrongFFIntegrand_IS_2->SetParameter(3, DeltaRho);
	
	//----------------------------------------//
	
	int nThetaBins = (int)csObj.m_StrongFF.size();
	
	for(int itbin=0; itbin < nThetaBins; itbin++) {
		
		double locTheta = csObj.m_angularBins[itbin];
		
		// longitudinal component of momentum-transferred:
		double locDelta = 0.5 * pow(csObj.getMesonMass(),2.0) / csObj.getBeamEnergy();
		
		// transverse momentum-transferred:
		double locQ2 = 4.0 * csObj.getBeamEnergy() * csObj.getMesonMomentum(locTheta) * 
			pow(sin(0.5*locTheta*TMath::DegToRad()),2.0);
		double locQ  = sqrt(locQ2);
		
		// express q and Delta in fm^-1:
		locDelta /= PrimExCS::m_GeV2fm;
		locQ     /= PrimExCS::m_GeV2fm;
		
		std::complex<double> FF_S, FF_I;
		
		//==================================================================================//
		// Calculate the form factor without initial state interaction:
		
		f_StrongFFIntegrand->SetParameter(1, locQ);
		f_StrongFFIntegrand->SetParameter(2, locDelta);
		
		//---------------------------------//
		// Real part:
		
		f_StrongFFIntegrand->SetParameter(0, 1.0);
		
		ROOT::Math::WrappedMultiTF1 wf_real(*f_StrongFFIntegrand);
		locIntegrator.SetFunction(wf_real);
		
		double ReFF_S = locIntegrator.Integral(xmin, xmax);
		
		//---------------------------------//
		// Imaginary part:
		
		f_StrongFFIntegrand->SetParameter(0, -1.0);
		
		ROOT::Math::WrappedMultiTF1 wf_imag(*f_StrongFFIntegrand);
		locIntegrator.SetFunction(wf_imag);
		
		double ImFF_S = locIntegrator.Integral(xmin, xmax);
		
		if(ReFF_S != ReFF_S) ReFF_S = 0.0;
		if(ImFF_S != ImFF_S) ImFF_S = 0.0;
		
		FF_S = {ReFF_S, ImFF_S};
		
		//==================================================================================//
		// Next, calculate the form factor from the initial state reaction:
		
		if(modelParameters.shadowingParameter>0.0) {
			
			f_StrongFFIntegrand_IS_1->SetParameter(1, locQ);
			f_StrongFFIntegrand_IS_1->SetParameter(2, locDelta);
			
			//---------------------------------//
			// Real part:
			
			f_StrongFFIntegrand_IS_1->SetParameter(0, 1.0);
			
			ROOT::Math::WrappedMultiTF1 wf_real_IS(*f_StrongFFIntegrand_IS_1);
			locIntegrator.SetFunction(wf_real_IS);
			
			double ReFF_I = locIntegrator.Integral(xmin, xmax);
			
			//---------------------------------//
			// Imaginary part:
			
			f_StrongFFIntegrand_IS_1->SetParameter(0, -1.0);
			
			ROOT::Math::WrappedMultiTF1 wf_imag_IS(*f_StrongFFIntegrand_IS_1);
			locIntegrator.SetFunction(wf_imag_IS);
			
			double ImFF_I = locIntegrator.Integral(xmin, xmax);
			
			if(ReFF_I != ReFF_I) ReFF_I = 0.0;
			if(ImFF_I != ImFF_I) ImFF_I = 0.0;
			
			FF_I = {ReFF_I + modelParameters.alphaIm*ImFF_I, ImFF_I - modelParameters.alphaIm*ReFF_I};
		}
		else {
			FF_I = {0.0, 0.0};
		}	
		//==================================================================================//
		
		printf("theta = %f\n  Re(FF_st) = %f, Im(FF_st) = %f, FF2 = %f\n  Re(FF_st_I) = %f, Im(FF_st_I)  %f\n",
			locTheta, real(FF_S), imag(FF_S), norm(FF_S), real(FF_I), imag(FF_I));
		
		csObj.m_StrongFF[itbin] = FF_S - modelParameters.shadowingParameter*FF_I;
	}
	
	f_StrongFFIntegrand->Delete();
	
	return;
}
