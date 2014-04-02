using std::cout;

#define FAST_TEST 0

enum {MAX_SAMPLE=1000000};
double G_eta [MAX_SAMPLE];
double G_reco[MAX_SAMPLE];
double G_new [MAX_SAMPLE];
double G_old [MAX_SAMPLE];
double G_pu  [MAX_SAMPLE];
int    G_ieta[MAX_SAMPLE];
int    G_num;
// fit flags:
int    G_use_old; //C-Style Bool
double G_eta_bin = -1; // no requirement

double G_eta_min = 0.0;
double G_eta_max = 5.0;

// number and type of eta bin fits:
#define NUM_ETA_BINS 22
#define MIN_ETA_BIN 0

//Calculate jet area by Eta Code
//see RCT eta index https://twiki.cern.ch/twiki/pub/CMS/RCTMap/trigtowermap.png

get_jet_area(int l1gEtaCode){
	jetarea = 0.0;
	if (l1gEtaCode==0) jetarea = 0.5*0.348*6.0;
	if (l1gEtaCode==1) jetarea = 0.5*0.348*9.0;
	if (l1gEtaCode==2) jetarea = 0.5*0.348*9.0;
	if (l1gEtaCode==3) jetarea = 0.5*0.348*6.0+0.828*0.348*3.0;
	if (l1gEtaCode==4) jetarea = 0.828*0.348*3.0+0.432*0.348*3.0+0.5*0.348*3.0;
	if (l1gEtaCode==5) jetarea = 0.432*0.348*3.0+3.0*0.432*0.348+3.0*0.828*0.348;
	if (l1gEtaCode==6) jetarea = 0.348*0.348*6.0+3.0*0.432*0.348;
	if (l1gEtaCode==7) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==8) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==9) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==10) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==21) jetarea = 0.5*0.348*6.0;
	if (l1gEtaCode==20) jetarea = 0.5*0.348*9.0;
	if (l1gEtaCode==19) jetarea = 0.5*0.348*9.0;
	if (l1gEtaCode==18) jetarea = 0.5*0.348*6.0+0.828*0.348*3.0;
	if (l1gEtaCode==17) jetarea = 0.828*0.348*3.0+0.432*0.348*3.0+0.5*0.348*3.0;
	if (l1gEtaCode==16) jetarea = 0.432*0.348*3.0+3.0*0.432*0.348+3.0*0.828*0.348;
	if (l1gEtaCode==15) jetarea = 0.348*0.348*6.0+3.0*0.432*0.348;
	if (l1gEtaCode==14) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==13) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==12) jetarea = 0.348*0.348*9.0;
	if (l1gEtaCode==11) jetarea = 0.348*0.348*9.0;
	return jetarea;
}  


//results of full calibration:

double R_ogamma, R_dogamma;
double R_alpha[NUM_ETA_BINS],  R_beta[NUM_ETA_BINS],  R_gamma[NUM_ETA_BINS];
double R_dalpha[NUM_ETA_BINS], R_dbeta[NUM_ETA_BINS], R_dgamma[NUM_ETA_BINS];
double R_old[MAX_SAMPLE];
double R_new[MAX_SAMPLE];

//Event with an eta range
int event_passes(int i){
	double abseta = fabs(G_eta[i]);

	if (abseta < G_eta_min) return 0;
	if (abseta > G_eta_max) return 0;

	if (G_eta_bin >= 0){
		if (G_ieta[i] != G_eta_bin) return 0;
	}

	return 1;
}



// choose pt function old or new
double choose_pt(int i){
	if (G_use_old) return G_old[i];
	return G_new[i];
}


//What is this?
//defining a graph
//with a y_min value for y target?


double interpolate_increasing(double ytarget, TGraph * g){  
	int n = g->GetN();
	double * x = g->GetX();
	double * y = g->GetY();
	for (int i=0; i<n; i++){
		//cout << x[i] << " " << y[i] << "\n";

		if (y[i] > ytarget){
			if (i==0) return x[0];
			double ya  = y[i-1]; //previous y value
			double yb  = y[i]; //current y value
			double dya = ytarget - ya; //(delta y previous) y target- previous y value
			double dyb = yb - ytarget; //(delta y current) y target- current y value
			double tot = yb - ya; //change from current y to previous y  
			if (tot == 0.0) return x[i];  //if current= previous return the x value    
			return (x[i-1] * dya + x[i] * dyb) / tot; // why??? if else? return previous x*(target y-previousy) + current x*(target y - previous y)
		} 
	} 
	return x[n-1];//return number of x bins?
}




double cross_calibrate(double thresh, double start, double delta, const char * plotname){
	const int ntobins = 40;
	double thresho       = thresh;  
	double threshn       = thresh;  

	TGraphAsymmErrors * tobef = NULL;
	TGraphAsymmErrors * toaft = NULL;

	double best_diff   = 1000.0;
	double best_offset = 0.0;
	double offset = start;

	//MODES:  INIT, SCAN, RECALC
	for (int mode=-1; mode <=1; ){
		if (delta == 0) mode = 1;

		if (mode == -1) mode = 0;
		else if (mode == 0){
			offset += delta;
			if (fabs(offset) > fabs(start + 20*delta)) mode = 1;
		}
		if (mode == 1){
			offset = best_offset;
		}

		TH1F * hturnondenom  = new TH1F("hfdifb","",ntobins, 0.0,200.0);
		TH1F * hturnonbef    = new TH1F("hturnon1bef","",ntobins, 0.0,200.0);
		TH1F * hturnonaft    = new TH1F("hturnon1aft","",ntobins, 0.0,200.0);

		for (int i=0; i<G_num; i++){
			if (! event_passes(i)) continue;

			double ptrec = G_reco[i];
			if (ptrec <=0.0) continue;
			double ptold  = R_old[i];
			double ptnew  = R_new[i];

			hturnondenom->Fill(ptrec);
			if (ptold >= thresho) hturnonbef->Fill(ptrec);
			if (ptnew >= threshn-offset) hturnonaft->Fill(ptrec);

		}    

		if (tobef) delete tobef;
		if (toaft) delete toaft;
		TGraphAsymmErrors * tobef = new TGraphAsymmErrors(hturnonbef, hturnondenom);
		TGraphAsymmErrors * toaft = new TGraphAsymmErrors(hturnonaft, hturnondenom);    
		double bef50 = interpolate_increasing(0.50, tobef);
		double aft50 = interpolate_increasing(0.50, toaft);
		//cout << "bef50:  " << bef50 << "\n";
		//cout << "aft50:  " << aft50 << "\n";
		//cout << "offset:  " << offset << "  delta:  " << bef50 - aft50 << "\n";    

		delete hturnondenom;
		delete hturnonbef;
		delete hturnonaft;

		if (mode == 1) mode++;
		else {
			double diff = fabs(bef50-aft50); 
			if (diff < best_diff) {
				best_diff  = diff;
				best_offset = offset; 
			}
		}
	}
	cout << "FOUND best offset value:  " << best_offset << " with diff " << best_diff << "\n";

	double bef95 = interpolate_increasing(0.95, tobef);
	double aft95 = interpolate_increasing(0.95, toaft);

	cout << "bef95:  " << bef95 << "\n";
	cout << "aft95:  " << aft95 << "\n";

	tobef->SetLineColor(2);
	toaft->SetLineColor(4);
	TCanvas * c3 = new TCanvas();
	toaft->Draw("APE");
	tobef->Draw("PESAME");
	toaft->Draw("PESAME");
	c3->SaveAs(plotname);
}

void control_plots(double thresh, double offset){
	TH1F * hietab = new TH1F("hietab","",23, 0.0, 23.0);
	TH1F * hietaa = new TH1F("hietaa","",23, 0.0, 23.0);

	cout << "running control plot!\n";

	for (int i=0; i<G_num; i++){

		if (! event_passes(i)) continue;    
		//cout << G_ieta[i] << "\n";

		double ptrec = G_reco[i];
		if (ptrec <=0.0) continue;
		double ptold  = R_old[i];
		double ptnew  = R_new[i];

		if (ptold >= thresh){
			hietab->Fill(G_ieta[i]);
		} 
		if (ptnew >= thresh-offset){
			hietaa->Fill(G_ieta[i]);
		}
	}

	TCanvas * c = new TCanvas();
	hietab->Draw("H");
	hietaa->Draw("HSAME");
}

void write_calibrated(double alpha, double gamma){
	for (int i=0; i<G_num; i++){
		if (! event_passes(i)) continue;
		if (G_use_old){
			if (G_old[i] > 0.0) R_old[i] = G_old[i] + gamma;
			//cout << R_old[i] << "\n";
		} else {
			if (G_new[i] > 0.0) R_new[i] = alpha*(G_new[i]) + gamma;
		}
	}
}




/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
//  Change Code below this line                  //
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/ 


void control_plots(){
	//	pass1 settings
	//	double offset20 = cross_calibrate(20, 0.,  0.0, "jet20.png");
	//	double offset36 = cross_calibrate(36, 1.,  0.1, "jet36.png");
	//	double offset48 = cross_calibrate(48, 5.0, 0.1, "jet48.png");  
	//

	G_eta_bin = -1; // no requirement, one factor!
	//	double offset20 = cross_calibrate(20, 0.0,  0.0, "jet20.png");
	//	double offset36 = cross_calibrate(36, 0.0,  0.0, "jet36.png");
	//	double offset48 = cross_calibrate(48, 0.0, 0.0, "jet48.png");  
	//	double offset80 = cross_calibrate(80, 0.0, 0.0, "jet80.png");

	//	control_plots(20, offset20);
	// 	control_plots(36, offset36);
	//	control_plots(48, offset48);
	//	control_plots(80, offset80);


	ofstream output("calib_npvsall_pum0_SF.txt");
	for (int etabin=0; etabin < NUM_ETA_BINS; etabin++){
		output << R_alpha[etabin] << ", " << R_beta[etabin]  << ", " << R_gamma[etabin] << ",\n"; 
	}
	output.close();
}



//where funciton to be minimized is defined
void fcn_chisq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *x, Int_t iflag){
	double alpha  = x[0];
	double gamma  = x[1];

	double chisq = 0.0;
	for (int i=0; i<G_num; i++){
		if (! event_passes(i)) continue;
		double pt   = choose_pt(i);
		double reco = G_reco[i];
		double pu   = G_pu[i];
		if (reco <= 0.0) continue;
		if (pt <= 0.0) continue;
		double delta = reco - ((alpha * pt)+gamma);
		chisq += delta*delta;
	}
	//cout << chisq << "\n";
	f = chisq;
}

void testcalib(){

	fill_calib_sample();

//	G_use_old = 1;
//	G_eta_bin = -1; // no requirement 
//	write_calibrated(1.0, 0.0, 0.0); //?

	int ierflg;
	TMinuit *gMinuit = new TMinuit(10);  //initialize TMinuit with a maximum of 10  params
	gMinuit->mninit(5,6,7);
	gMinuit->SetFCN(fcn_chisq);

	for (int etapack=0; etapack < NUM_ETA_BINS; etapack++){
		G_eta_bin = MIN_ETA_BIN + etapack;
		G_use_old = 0;


		double gamma = 0.0; // default Subtraction

		//double beta = -get_jet_area(G_eta_bin); // PU UIC Uncomment line when minimizing pile up. Provides a starting variable for beta
		gMinuit->mnparm(0, "alpha",    1.0, 0.1, 0,0,ierflg);
		gMinuit->mnparm(1, "gamma",    0.0, 0.1, 0,0,ierflg);
		gMinuit->FixParameter(0);
		gMinuit->FixParameter(1);

		gMinuit->Release(0);
	//	gMinuit->Release(1);

		gMinuit->Migrad();

		gMinuit->GetParameter(0, R_alpha [etapack], R_dalpha [etapack]);
		gMinuit->GetParameter(1, R_gamma  [etapack], R_dgamma  [etapack]);

		write_calibrated(R_alpha[etapack], R_gamma[etapack]);

	}  

	control_plots();
}

void fill_calib_sample(){

	//file used to get calibrations
	TFile f("/nfs_scratch/laura/March4_BUGFIX_SFGamma_allet.root");
	//change to corrjetEfficiency when running calibrations
	TTree * tree = (TTree *) f.Get("jetEfficiency/Ntuple");

	float eta, phi, rpt, opt, npt, lumi, npvs, pu;
	int ieta;

	tree->SetBranchAddress("recoEta", &eta);
	tree->SetBranchAddress("recoPhi",    &phi);
	tree->SetBranchAddress("recoPt",     &rpt);
	tree->SetBranchAddress("l1Pt",       &opt);
	tree->SetBranchAddress("l1gPt",      &npt);
	tree->SetBranchAddress("lumi",       &lumi);
	tree->SetBranchAddress("nPVs",       &npvs);
	//tree->SetBranchAddress("l1gPUUIC",   &pu);
	tree->SetBranchAddress("l1gPU",      &pu);
	//tree->SetBranchAddress("l1gPUM",      &pum);
	tree->SetBranchAddress("l1gEtaCode", &ieta);

	int n = tree->GetEntries();

	if (FAST_TEST){
		if (n > 1000000) n = 1000000;
	}

	G_num = 0;
	int update = n / 20;
	if (update < 1) update=1;  
	cout << "Progress:  ";
	for (int i=0; i<n; i++){
		if ((i+1)%update==0) { cout << "."; cout.flush(); }    
		tree->GetEntry(i);

		if (rpt < 20.0) continue;
		//if (npt <= 0.0) continue;
		//if (npt <= 0.0) continue;
		//if (reco > 100.0) continue;
		//if (fabs(eta) > 2.0) continue;
		//if (ghd <= 0.0) continue;

		//initial cuts can be made here.


		//cout << "lumi:  " << lumi << "\n";
	//	if (npvs > 10) continue;

		G_eta [G_num] = eta;
		G_reco[G_num] = rpt;
		G_new [G_num] = npt;
		G_old [G_num] = opt;
		G_pu  [G_num] = pu;
		G_ieta[G_num] = ieta;
		G_num++;

		if (G_num >= MAX_SAMPLE){
			cout << "\n";
			cout << "ERROR:  calibration sample reached maximum size.\n";
			exit(0);
		}    
	}
	cout <<"\n";
	cout << "Passed:  " << G_num << "\n";
}

