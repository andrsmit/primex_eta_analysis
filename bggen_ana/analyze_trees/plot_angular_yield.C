
#include "particleType.h"

vector<pair<int,vector<Particle_t>>> m_reaction_types;

void initializeReactionTypes();

void plot_angular_yield() {
	
	gStyle->SetOptStat(0);
	
	//-------------------------------------------------------------------//
	
	TFile *fIn_proton  = new TFile("rootFiles/phase1/Helium-proton.root",  "READ");
	TFile *fIn_neutron = new TFile("rootFiles/phase1/Helium-neutron.root", "READ");
	
	TH1F *h_thrown[2];
	h_thrown[0] = (TH1F*)fIn_proton->Get("thrown")->Clone("thrown_proton");
	h_thrown[1] = (TH1F*)fIn_neutron->Get("thrown")->Clone("thrown_neutron");
	
	h_thrown[0]->SetDirectory(0);
	h_thrown[1]->SetDirectory(0);
	
	fIn_proton->Close();
	fIn_neutron->Close();
	
	//-------------------------------------------------------------------//
	
	TFile *fIn = new TFile("rootFiles/phase1/Helium.root", "READ");
	
	initializeReactionTypes();
	
	// list of reactions we want to plot:
	vector<int> reaction_list = {0, 1, 4, 5, 16, 17};
	vector<pair<TString,TH1F*>> h_theta;
	
	TCanvas *c = new TCanvas("c", "c", 950, 700);
	c->SetTickx(); c->SetTicky();
	
	double hist_max = 0.;
	for(int ireaction=0; ireaction<reaction_list.size(); ireaction++) {
		
		int loc_reaction = reaction_list[ireaction];
		double n_thrown = h_thrown[m_reaction_types[loc_reaction].first]->Integral();
		
		TString dir_label = "";
		for(int iparticle=0; iparticle<m_reaction_types[loc_reaction].second.size(); iparticle++) {
			dir_label += ShortName(m_reaction_types[loc_reaction].second[iparticle]);
		}
		TH1F *loc_h1 = (TH1F*)fIn->Get(Form("%s/theta_veto_%s",dir_label.Data(),dir_label.Data()));
		loc_h1->Scale(1.0/n_thrown);
		if(loc_h1->GetMaximum() > hist_max) hist_max = loc_h1->GetMaximum();
		
		TString loc_reaction_name = "";
		for(int iparticle=0; iparticle<m_reaction_types[loc_reaction].second.size(); iparticle++) {
			loc_reaction_name += ParticleName_ROOT(m_reaction_types[loc_reaction].second[iparticle]);
		}
		h_theta.push_back({loc_reaction_name, loc_h1});
	}
	
	int first = 1;
	
	TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
	
	for(int ihist=0; ihist<h_theta.size(); ihist++) {
		h_theta[ihist].second->SetLineColor(1+ihist);
		h_theta[ihist].second->SetTitle("Angular Yield");
		h_theta[ihist].second->GetXaxis()->SetTitle("#theta_{#gamma#gamma} [deg.]");
		h_theta[ihist].second->GetXaxis()->SetTitleSize(0.05);
		h_theta[ihist].second->GetXaxis()->CenterTitle(true);
		
		h_theta[ihist].second->GetYaxis()->SetRangeUser(0., 1.2*hist_max);
		c->cd();
		if(first) {
			h_theta[ihist].second->Draw("PE");
			first = 0;
		}
		else h_theta[ihist].second->Draw("P same");
		
		leg->AddEntry(h_theta[ihist].second, h_theta[ihist].first.Data(), "PE");
	}
	leg->Draw();
	
	return;
}

void initializeReactionTypes() {
	
	m_reaction_types.clear();
	
	m_reaction_types.push_back({0, {Eta, Proton}});                   //  0
	m_reaction_types.push_back({1, {Eta, Neutron}});                  //  1
	m_reaction_types.push_back({0, {Eta, Pi0, Proton}});              //  2
	m_reaction_types.push_back({1, {Eta, Pi0, Neutron}});             //  3
	m_reaction_types.push_back({1, {Eta, PiMinus, Proton}});          //  4
	m_reaction_types.push_back({0, {Eta, PiPlus, Neutron}});          //  5
	m_reaction_types.push_back({0, {omega, Proton}});                 //  6
	m_reaction_types.push_back({1, {omega, Neutron}});                //  7
	m_reaction_types.push_back({0, {Rho0, Proton}});                  //  8
	m_reaction_types.push_back({1, {Rho0, Neutron}});                 //  9
	m_reaction_types.push_back({0, {phiMeson, Proton}});              // 10
	m_reaction_types.push_back({1, {phiMeson, Neutron}});             // 11
	m_reaction_types.push_back({0, {Eta, Pi0, Pi0, Proton}});         // 12
	m_reaction_types.push_back({1, {Eta, Pi0, Pi0, Neutron}});        // 13 
 	m_reaction_types.push_back({0, {Eta, Pi0, PiPlus, Neutron}});     // 14
	m_reaction_types.push_back({1, {Eta, Pi0, PiMinus, Proton}});     // 15
	m_reaction_types.push_back({0, {Eta, PiPlus, PiMinus, Proton}});  // 16
	m_reaction_types.push_back({1, {Eta, PiPlus, PiMinus, Neutron}}); // 17
	
	return;
}
