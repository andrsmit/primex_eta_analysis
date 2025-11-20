double R = 14.0;
double pf = 0.26;
double Gp(double *x, double *par) {
	double q = x[0];
	double G = (1.0 + pow(q*R/sqrt(15.0),4.0))*exp(-2.0*pow(q*R,2.0)/15.0);
	return (1.0-G);
}
double G2(double *x, double *par) {
	
	double q = x[0];
	if(q>2.0*pf) return 1.0;
	double G = (3.0/4.0)*(q/pf) - (1.0/16.0)*pow(q/pf,3.0);
	return G;
}
void plot() {
	
	TF1 *f1 = new TF1("f1",Gp,0.,1.0,0);
	f1->Draw();
	TF1 *f2 = new TF1("f2",G2,0.,1.0,0);
	f2->SetLineColor(kBlue);
	f2->SetLineStyle(2);
	f2->Draw("same");
	
}
