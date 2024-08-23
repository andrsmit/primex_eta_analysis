#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/eta_inc.h"
#include "/work/halld/home/andrsmit/primex_eta_analysis/cross_section/include/style.cc"

void getEtaMean(TF1 *f1, double &mean, double &mean_err);

void plotFitResults() {
	
	double *angle     = new double[m_n_bins_fit];
	double *angle_err = new double[m_n_bins_fit];
	
	double *eta_mean     = new double[m_n_bins_fit];
	double *eta_mean_err = new double[m_n_bins_fit];
	
	int mean_par_index = 1;
	if(m_signal_fit_option==2) mean_par_index = 2;
	
	for(int ibin=0; ibin<m_n_bins_fit; ibin++) {
		angle[ibin]     = m_angular_bin[ibin].first;
		angle_err[ibin] = m_angular_bin[ibin].second;
		
		getEtaMean(f_mgg_fit_functions[ibin].second, eta_mean[ibin], eta_mean_err[ibin]);
	}
	
	TGraphErrors *gEtaMean = new TGraphErrors(m_n_bins_fit, angle, eta_mean, angle_err, eta_mean_err);
	gEtaMean->GetXaxis()->SetTitle("Polar Angle, #theta_{#gamma#gamma} [#circ]");
	gEtaMean->GetXaxis()->SetTitleSize(0.05);
	gEtaMean->GetXaxis()->SetTitleOffset(1.0);
	gEtaMean->GetXaxis()->CenterTitle("");
	gEtaMean->GetYaxis()->SetTitle("#mu_{#eta} (fit result) [GeV/c^{2}]");
	gEtaMean->GetYaxis()->SetTitleSize(0.05);
	gEtaMean->GetYaxis()->SetTitleOffset(1.0);
	gEtaMean->GetYaxis()->CenterTitle("");
	gEtaMean->SetTitle("");
	gEtaMean->SetMarkerStyle(4);
	gEtaMean->SetMarkerSize(0.7);
	gEtaMean->SetMarkerColor(kBlue+2);
	gEtaMean->SetLineColor(kBlue+2);
	gEtaMean->SetLineWidth(1);
	
	TCanvas *cEtaMean = new TCanvas("cEtaMean", "Eta Mean", 700, 500);
	styleCanvas(cEtaMean);
	
	cEtaMean->cd();
	gEtaMean->Draw("APE");
	
	return;
}

void getEtaMean(TF1 *f1, double &mean, double &mean_err) {
	
	switch(m_signal_fit_option) {
		case 1:
			mean      = f1->GetParameter(1);
			mean_err  = f1->GetParError(1);
			break;
		case 2:
		{
			double n1 = f1->GetParameter(0);
			double n2 = f1->GetParameter(1);
			mean      = (n1*f1->GetParameter(2) + n2*(f1->GetParameter(3)+f1->GetParameter(2))) / (n1+n2);
			mean_err  = sqrt(pow(f1->GetParError(2),2.0) + pow(f1->GetParError(3),2.0));
			break;
		}
		case 3:
			mean      = f1->GetParameter(1);
			mean_err  = f1->GetParError(1);
			break;
		case 4:
			mean      = f1->GetParameter(1);
			mean_err  = f1->GetParError(1);
			break;
	}
	
	return;
}
