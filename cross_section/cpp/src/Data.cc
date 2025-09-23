#include "EtaAnalyzer.h"
#include "MggFitter.h"
#include "CrossSection.h"

TString    dataDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/analyze_trees/rootFiles";
TString   etamcDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/eta_gg_matrix/analyze_trees/rootFiles";
TString omegamcDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/omega_gammapi0/analyze_trees/rootFiles";
TString   bggenDirectory = "/work/halld/home/andrsmit/primex_eta_analysis/bggen_ana/analyze_trees/rootFiles";

int EtaAnalyzer::LoadDataHistograms()
{
	printf("\nREADING DATA HISTGRAMS...\n");
	
	TString fieldString = "bfield";
	if(m_phase==1) fieldString = "nobfield";
	
	TString anaString = "";
	TString dirString = "";
	switch(m_analysisOption) {
		case 1:
			dirString = "matrix";
			anaString = "_matrix";
			break;
		case 2:
			dirString = "FCAL";
			anaString = "_FCAL";
			break;
		case 3:
			dirString = "BCAL";
			anaString = "_BCAL";
			break;
		case 4:
			dirString = "beam";
			anaString = "_beam";
			break;
		case 5:
			dirString = "TOF";
			anaString = "_TOF";
			break;
		default:
			dirString = "default";
			anaString = "";
			break;
	}
	
	TString vetoStr = "";
	if(m_analysisOption>0) vetoStr = Form("_VetoOption%d", m_vetoOption);
	if(m_analysisOption==8) vetoStr = "";
	
	TString fileTypeString = m_phase==1 ? "EVIO" : "REST";
	
	TString fullTargetFileName  = Form("%s/phase%d/%s/CORRECTED_0.5PERCENT/%s/full_target_%s%s%s.root",  dataDirectory.Data(), m_phase, fileTypeString.Data(), 
		dirString.Data(), fieldString.Data(), anaString.Data(), vetoStr.Data());
	
	TString emptyTargetFileName = Form("%s/phase%d/%s/CORRECTED_0.5PERCENT/%s/empty_target_%s%s%s.root", dataDirectory.Data(), m_phase, fileTypeString.Data(), 
		dirString.Data(), fieldString.Data(), anaString.Data(), vetoStr.Data());
	
	if(gSystem->AccessPathName(fullTargetFileName.Data()) || gSystem->AccessPathName(emptyTargetFileName.Data()))
		return 1;
	
	printf("   Full Target: %s\n",  fullTargetFileName.Data());
	printf("  Empty Target: %s\n", emptyTargetFileName.Data());
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (m_analysisOption=1) is fundamentally different from the rest:
	
	if(m_analysisOption==1) {
		
		//-------------------------------------//
		// Get the full target data:
		
		TFile *fullTargetFile = new TFile(fullTargetFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)fullTargetFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-m_minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-m_maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of data histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_mggVsThetaFull = (TH2F*)h3Full->Project3D("yx")->Clone("mgg_vs_theta_full");
		h_mggVsThetaFull->SetDirectory(0);
		
		fullTargetFile->Close();
		
		//-------------------------------------//
		// Do the same for the empty target data:
		
		TFile *emptyTargetFile = new TFile(emptyTargetFileName.Data(), "READ");
		
		TH3F *h3Empty = (TH3F*)emptyTargetFile->Get("invmassMatrix")->Clone("h3Empty");
		
		locBinSize = h3Empty->GetZaxis()->GetBinWidth(1);
		ebin1 = h3Empty->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		ebin2 = h3Empty->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		h3Empty->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_mggVsThetaEmpty = (TH2F*)h3Empty->Project3D("yx")->Clone("mgg_vs_theta_empty");
		h_mggVsThetaEmpty->SetDirectory(0);
		
		emptyTargetFile->Close();
	}
	else {
		TFile *fullTargetFile = new TFile(fullTargetFileName.Data(), "READ");
		
		h_mggVsThetaFull = (TH2F*)fullTargetFile->Get(Form("%s",m_mggHistName.Data()))->Clone("mgg_vs_theta_full");
		//h_mggVsThetaFull = (TH2F*)fullTargetFile->Get("VetoOption1/mgg_const_cut_veto_1")->Clone("mgg_vs_theta_full");
		h_mggVsThetaFull->SetDirectory(0);
		fullTargetFile->Close();
		
		TFile *emptyTargetFile = new TFile(emptyTargetFileName.Data(), "READ");
		
		h_mggVsThetaEmpty = (TH2F*)emptyTargetFile->Get(Form("%s",m_mggHistName.Data()))->Clone("mgg_vs_theta_empty");
		//h_mggVsThetaEmpty = (TH2F*)emptyTargetFile->Get("VetoOption1/mgg_const_cut_veto_1")->Clone("mgg_vs_theta_empty");
		h_mggVsThetaEmpty->SetDirectory(0);
		emptyTargetFile->Close();
	}
	
	// scale empty target data by photon flux ratio:
	
	h_mggVsThetaEmpty->Scale(m_emptyTargetFluxRatio);
	
	return 0;
}

int EtaAnalyzer::LoadLineshapes()
{
	printf("\nREADING LINESHAPES...\n");
	
	if(LoadEtaLineshape()) return 1;
	if(LoadBGGENLineshape()) return 2;
	if(LoadEtaPionLineshape()) return 3;
	if(LoadOmegaLineshape()) return 4;
	
	// Read the lineshape of omega->pi0+gamma from the first FDC package:
	if(m_fitOption_empty==1 && m_emptyFitOption_fdc>=2)
	{
		if(LoadFDCOmegaLineshape()) return 5;
	}
	return 0;
}

int EtaAnalyzer::LoadEtaLineshape()
{
	TString anaString = "";
	switch(m_analysisOption) {
		case 1:
			anaString = "_matrix";
			break;
		case 2:
			anaString = "_FCAL";
			break;
		case 3:
			anaString = "_BCAL";
			break;
		case 4:
			anaString = "_beam";
			break;
		case 5:
			anaString = "_TOF";
			break;
		case 8:
			anaString = "_angular";
			break;
		default:
			anaString = "";
			break;
	}
	
	TString vetoStr = "";
	if(m_analysisOption>0) vetoStr = Form("_VetoOption%d", m_vetoOption);
	if(m_analysisOption==8) vetoStr = "";
	
	int locPhase = m_phase;
	TString mcFileName = Form("%s/phase%d/phase%d%s%s.root",  etamcDirectory.Data(), locPhase, locPhase, anaString.Data(), vetoStr.Data());
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	printf("  Eta lineshape from %s\n", mcFileName.Data());
	
	//--------------------------------//
	
	// Reading the histograms with the "matrix" option (m_analysisOption=1) is fundamentally different from the rest:
	
	if(m_analysisOption==1) {
		
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		TH3F *h3Full = (TH3F*)mcFile->Get("invmassMatrix")->Clone("h3Full");
		
		double locBinSize = h3Full->GetZaxis()->GetBinWidth(1);
		int ebin1 = h3Full->GetZaxis()->FindBin(m_minBeamEnergy + 0.5*locBinSize);
		int ebin2 = h3Full->GetZaxis()->FindBin(m_maxBeamEnergy - 0.5*locBinSize);
		
		// consistency check to make sure binning of data histogram aligns with specified energy range:
		
		double locMinimumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin1) - 0.5*locBinSize;
		double locMaximumEnergy = h3Full->GetZaxis()->GetBinCenter(ebin2) + 0.5*locBinSize;
		
		if((fabs(locMinimumEnergy-m_minBeamEnergy)>1.e6) || (fabs(locMaximumEnergy-m_maxBeamEnergy)>1.e6)) {
			printf("\n\nWarning! Bin edges of eta mc histogram are inconsistent with specified energy range:\n");
			printf("   Specified energy range: %.5f GeV - %.5f GeV\n", m_minBeamEnergy, m_maxBeamEnergy);
			printf("   Flux Integration Range: %.5f GeV - %.5f GeV\n\n", locMinimumEnergy, locMaximumEnergy);
		}
		
		h3Full->GetZaxis()->SetRange(ebin1, ebin2);
		
		h_etaLineshapeCoh = (TH2F*)h3Full->Project3D("yx")->Clone("etaLineshapeCoh");
		h_etaLineshapeCoh->SetDirectory(0);
		
		mcFile->Close();
	}
	else if(m_analysisOption==8) {
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		TString locHistName = "mgg_const_angularSmear_00";
		if(m_matrixHistName.Contains("_AngularSmear_")) {
			locHistName = Form("mgg_const_angularSmear_%c%c",
				m_matrixHistName[m_matrixHistName.Length()-2], m_matrixHistName[m_matrixHistName.Length()-1]);
		}
		else if(m_matrixHistName.Contains("_AngularShift_")) {
			locHistName = Form("mgg_const_angularShift_%c%c",
				m_matrixHistName[m_matrixHistName.Length()-2], m_matrixHistName[m_matrixHistName.Length()-1]);
		}
		printf("\n\n\n");
		printf("LINESHAPE HISTOGRAM: %s\n", locHistName.Data());
		printf("\n\n\n");
		
		h_etaLineshapeCoh = (TH2F*)mcFile->Get(locHistName.Data())->Clone("etaLineshapeCoh");
		h_etaLineshapeCoh->SetDirectory(0);
		mcFile->Close();
	}
	else {
		TFile *mcFile = new TFile(mcFileName.Data(), "READ");
		
		if(m_analysisOption==0) {
			h_etaLineshapeCoh = (TH2F*)mcFile->Get(Form("VetoOption%d/mgg_const_coh_veto_%d", 
				m_vetoOption, m_vetoOption))->Clone("etaLineshape");
		} else {
			h_etaLineshapeCoh = (TH2F*)mcFile->Get(Form("%s",m_mggHistName.Data()))->Clone("etaLineshapeCoh");
		}
		h_etaLineshapeCoh->SetDirectory(0);
		mcFile->Close();
	}
	
	return 0;
}

int EtaAnalyzer::LoadBGGENLineshape()
{
	TString anaString = "";
	TString locMggHistName = "mgg_const_bggen_signal";
	/*
	if(m_analysisOption==4) {
		int beamCutIndex = 0;
		if(m_mggHistName.Contains("mgg_Elasticity")) {
			int locLength = m_mggHistName.Length();
			TString indexStr(m_mggHistName(locLength-2,locLength));
			beamCutIndex = indexStr.Atoi();
		}
		anaString = "_beam";
		locMggHistName = Form("mgg_Elasticity_signal_%02d",beamCutIndex);
	}
	*/
	TString mcFileName  = Form("%s/phase%d/Helium_VetoOption%d%s.root",  bggenDirectory.Data(), m_phase, m_vetoOption, anaString.Data());
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	printf("  Eta lineshape from %s\n", mcFileName.Data());
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	h_etaLineshapeBGGEN = (TH2F*)mcFile->Get(Form("%s",locMggHistName.Data()))->Clone("etaLineshapeBGGEN");
	h_etaLineshapeBGGEN->SetDirectory(0);
	mcFile->Close();
	return 0;
}

int EtaAnalyzer::LoadEtaPionLineshape()
{
	TString anaString = "";
	TString locMggHistName1 = "mgg_const_bggen_";
	TString locMggHistName2 = "_cut";
	// histName = 'Form(%ssignal%s, locMggHistName1, locMggHistName2)', or 'Form(%setapion%s, locMggHistName1, locMggHistName2)'
	/*
	if(m_analysisOption==4) {
		int beamCutIndex = 0;
		if(m_mggHistName.Contains("mgg_Elasticity")) {
			int locLength = m_mggHistName.Length();
			TString indexStr(m_mggHistName(locLength-2,locLength));
			beamCutIndex = indexStr.Atoi();
		}
		anaString = "_beam";
		locMggHistName1 = "mgg_Elasticity_";
		locMggHistName2 = Form("_%02d",beamCutIndex);
	}
	*/
	TString mcFileName  = Form("%s/phase%d/Helium_VetoOption%d%s.root",  bggenDirectory.Data(), m_phase, m_vetoOption, anaString.Data());
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	printf("  Eta+Pion lineshape from %s\n", mcFileName.Data());
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	// thrown histogram:
	
	TH1F *locThrown            = (TH1F*)mcFile->Get("thrown_reactions_bggen");
	double nThrown             = locThrown->Integral() * 10.0;
	double simulatedLuminosity = nThrown / 1.18673e+08; // 1.18e+08 pb/nucleon is total photoproduction cross section
	
	double scaleFactor = 2.0*(1.e-6)*m_luminosity / simulatedLuminosity;
	
	// eta+pion background channels:
	
	h_eta1PionLineshape = (TH2F*)mcFile->Get(Form("%setapion%s", locMggHistName1.Data(), locMggHistName2.Data()))->Clone("etapionLineshape");
	h_eta1PionLineshape->SetDirectory(0);
	h_eta1PionLineshape->Scale(scaleFactor);
	
	// eta+2pion background channels:
	
	h_eta2PionLineshape = (TH2F*)mcFile->Get(Form("%seta2pion%s", locMggHistName1.Data(), locMggHistName2.Data()))->Clone("eta2pionLineshape");
	h_eta2PionLineshape->SetDirectory(0);
	h_eta2PionLineshape->Scale(scaleFactor);
	
	// eta+3pion background channels:
	
	h_eta3PionLineshape = (TH2F*)mcFile->Get(Form("%seta3pion%s", locMggHistName1.Data(), locMggHistName2.Data()))->Clone("eta3pionLineshape");
	h_eta3PionLineshape->SetDirectory(0);
	h_eta3PionLineshape->Scale(scaleFactor);
	
	// other hadronic background channels:
	
	h_hadronicBkgdLineshape = (TH2F*)mcFile->Get(Form("%sbkgd%s", locMggHistName1.Data(), locMggHistName2.Data()))->Clone("bkgdLineshape");
	h_hadronicBkgdLineshape->SetDirectory(0);
	h_hadronicBkgdLineshape->Scale(scaleFactor);
	/*
	TH2F *h_rhoBkgdLineshape = (TH2F*)mcFile->Get(Form("%srho%s", locMggHistName1.Data(), locMggHistName2.Data()))->Clone("rhoBkgdLineshape");
	h_rhoBkgdLineshape->Scale(scaleFactor);
	h_hadronicBkgdLineshape->Add(h_rhoBkgdLineshape);
	*/
	//-------------------------------------------------------------------//
	// Store fractions of each background:
	
	TH2F *h2_exclusive = (TH2F*)mcFile->Get(Form("%ssignal%s", locMggHistName1.Data(), locMggHistName2.Data()));
	h2_exclusive->Scale(scaleFactor);
	
	//----------------------//
	// Without cut on invariant mass:
	
	int minMggCutBin = h2_exclusive->GetYaxis()->FindBin(m_minFitRange);
	int maxMggCutBin = h2_exclusive->GetYaxis()->FindBin(m_maxFitRange);
	
	TH1F *h1_exclusive = (TH1F*)h2_exclusive->ProjectionX("h1_exclusive",       minMggCutBin, maxMggCutBin);
	TH1F *h1_eta1pion  = (TH1F*)h_eta1PionLineshape->ProjectionX("h1_eta1pion", minMggCutBin, maxMggCutBin);
	TH1F *h1_eta2pion  = (TH1F*)h_eta2PionLineshape->ProjectionX("h1_eta2pion", minMggCutBin, maxMggCutBin);
	TH1F *h1_eta3pion  = (TH1F*)h_eta3PionLineshape->ProjectionX("h1_eta3pion", minMggCutBin, maxMggCutBin);
	TH1F *h1_bkgd      = (TH1F*)h_hadronicBkgdLineshape->ProjectionX("h1_bkgd", minMggCutBin, maxMggCutBin);
	
	h1_exclusive->Rebin(m_rebinsTheta);
	h1_eta1pion->Rebin(m_rebinsTheta);
	h1_eta2pion->Rebin(m_rebinsTheta);
	h1_eta3pion->Rebin(m_rebinsTheta);
	h1_bkgd->Rebin(m_rebinsTheta);
	
	//----------------------//
	// With cut on invariant mass:
	
	minMggCutBin = h2_exclusive->GetYaxis()->FindBin(0.5);
	maxMggCutBin = h2_exclusive->GetYaxis()->FindBin(0.6);
	
	TH1F *h1_exclusive_cut = (TH1F*)h2_exclusive->ProjectionX("h1_exclusive_cut",       minMggCutBin, maxMggCutBin);
	TH1F *h1_eta1pion_cut  = (TH1F*)h_eta1PionLineshape->ProjectionX("h1_eta1pion_cut", minMggCutBin, maxMggCutBin);
	TH1F *h1_eta2pion_cut  = (TH1F*)h_eta2PionLineshape->ProjectionX("h1_eta2pion_cut", minMggCutBin, maxMggCutBin);
	TH1F *h1_eta3pion_cut  = (TH1F*)h_eta3PionLineshape->ProjectionX("h1_eta3pion_cut", minMggCutBin, maxMggCutBin);
	TH1F *h1_bkgd_cut      = (TH1F*)h_hadronicBkgdLineshape->ProjectionX("h1_bkgd_cut", minMggCutBin, maxMggCutBin);
	
	h1_exclusive_cut->Rebin(m_rebinsTheta);
	h1_eta1pion_cut->Rebin(m_rebinsTheta);
	h1_eta2pion_cut->Rebin(m_rebinsTheta);
	h1_eta3pion_cut->Rebin(m_rebinsTheta);
	h1_bkgd_cut->Rebin(m_rebinsTheta);
	
	if(m_fitOption_signal<9) {
		// sum all hadronic background channels into a single histogram:
		
		h_hadronicBkgdLineshape->Add(h_eta1PionLineshape);
		h_hadronicBkgdLineshape->Add(h_eta2PionLineshape);
		h_hadronicBkgdLineshape->Add(h_eta3PionLineshape);
		
		h1_bkgd->Add(h1_eta1pion);
		h1_bkgd->Add(h1_eta2pion);
		h1_bkgd->Add(h1_eta3pion);
		h1_bkgd->Divide(h1_exclusive);
		
		h_HadronicBkgdFraction_bggen = (TH1F*)h1_bkgd->Clone("hadronicBkgdFraction_bggen");
		h_HadronicBkgdFraction_bggen->SetDirectory(0);
		
		h1_bkgd_cut->Add(h1_eta1pion_cut);
		h1_bkgd_cut->Add(h1_eta2pion_cut);
		h1_bkgd_cut->Add(h1_eta3pion_cut);
		h1_bkgd_cut->Divide(h1_exclusive_cut);
		
		h_HadronicBkgdFraction_bggen_cut = (TH1F*)h1_bkgd_cut->Clone("hadronicBkgdFraction_bggen_cut");
		h_HadronicBkgdFraction_bggen_cut->SetDirectory(0);
	}
	else {
		// separate histograms for eta+pion and other backgrounds:
		
		h_hadronicBkgdLineshape->Add(h_eta2PionLineshape);
		h_hadronicBkgdLineshape->Add(h_eta3PionLineshape);
		
		h1_bkgd->Add(h1_eta2pion);
		h1_bkgd->Add(h1_eta3pion);
		
		h1_bkgd->Divide(h1_exclusive);
		h_HadronicBkgdFraction_bggen = (TH1F*)h1_bkgd->Clone("hadronicBkgdFraction_bggen");
		h_HadronicBkgdFraction_bggen->SetDirectory(0);
		
		h1_eta1pion->Divide(h1_exclusive);
		h_EtaPionBkgdFraction_bggen = (TH1F*)h1_eta1pion->Clone("etaPionFraction_bggen");
		h_EtaPionBkgdFraction_bggen->SetDirectory(0);
		
		h1_bkgd_cut->Add(h1_eta2pion_cut);
		h1_bkgd_cut->Add(h1_eta3pion_cut);
		
		h1_bkgd_cut->Divide(h1_exclusive_cut);
		h_HadronicBkgdFraction_bggen_cut = (TH1F*)h1_bkgd_cut->Clone("hadronicBkgdFraction_bggen_cut");
		h_HadronicBkgdFraction_bggen_cut->SetDirectory(0);
		
		h1_eta1pion_cut->Divide(h1_exclusive_cut);
		h_EtaPionBkgdFraction_bggen_cut = (TH1F*)h1_eta1pion_cut->Clone("etaPionFraction_bggen_cut");
		h_EtaPionBkgdFraction_bggen_cut->SetDirectory(0);
	}
	
	mcFile->Close();
	
	return 0;
}

int EtaAnalyzer::LoadOmegaLineshape()
{
	// Use the lineshape from bggen for the omega+other backgrounds:
	
	int locVetoOption = 4;
	if((m_vetoOption==1) || (m_vetoOption==6) || (m_vetoOption==7)) locVetoOption = m_vetoOption;
	
	TString subDirString = "";
	//if(m_phase==3) subDirString="/phase3_cuts";
	TString mcFileName = Form("%s/phase%d%s/Helium_VetoOption%d.root",  bggenDirectory.Data(), 1, subDirString.Data(), locVetoOption);
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	printf("  Omega lineshape from %s\n", mcFileName.Data());
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	h_omegaLineshape = (TH2F*)mcFile->Get("mgg_const_bggen_omega_cut")->Clone("omegaLineshape");
	h_omegaLineshape->SetDirectory(0);
	
	h_rhoLineshape = (TH2F*)mcFile->Get("mgg_const_bggen_rho_cut")->Clone("rhoLineshape");
	h_rhoLineshape->SetDirectory(0);
	
	// thrown histogram:
	
	TH1F *locThrown            = (TH1F*)mcFile->Get("thrown_reactions_bggen");
	double nThrown             = locThrown->Integral() * 10.0;
	double simulatedLuminosity = nThrown / 1.18673e+08; // 1.18e+08 pb/nucleon is total photoproduction cross section
	
	double scaleFactor = 2.0*(1.e-6)*m_luminosity / simulatedLuminosity;
	
	h_omegaLineshape->Scale(scaleFactor);
	h_rhoLineshape->Scale(scaleFactor);
	
	if(m_fitOption_omega!=4) {
		h_omegaLineshape->Add(h_rhoLineshape);
		
		TH2F *loc_h2_bkgd = (TH2F*)mcFile->Get("mgg_const_bggen_bkgd_cut")->Clone("bkgdLineshape");
		loc_h2_bkgd->Scale(scaleFactor);
		h_omegaLineshape->Add(loc_h2_bkgd);
		
	}
	
	mcFile->Close();
	return 0;
}

int EtaAnalyzer::LoadFDCOmegaLineshape()
{
	// Use the lineshape from bggen for the omega+other backgrounds:
	
	TString mcFileName  = Form("%s/phase3_1st-FDC-pack.root",  omegamcDirectory.Data());
	
	if(gSystem->AccessPathName(mcFileName.Data())) return 1;
	
	TFile *mcFile = new TFile(mcFileName.Data(), "READ");
	
	h_fdcOmegaLineshape = (TH2F*)mcFile->Get("mgg_const_veto_5")->Clone("fdcOmegaLineshape");
	h_fdcOmegaLineshape->SetDirectory(0);
	mcFile->Close();
	return 0;
}
