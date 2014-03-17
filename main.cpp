/**
  * Some kinematic tests for primarily eta prime decays performed, create plots from the generated pluto files
  * 
  * data generated using my simulation scripts, arguments are: number of total files, total events, number of channels (begin counting with 0), number of array entries for the simulation process (array contains 4 entries per channel, i. e. channel name, number of files, events per file, number of existing files (needed for number at the end of the file name to prevent overwriting))
  * ./sim.sh 6 60000000 5 24 etap_e+e-g 1 10000000 1 etap_pi0pi0eta 1 10000000 0 etap_pi0pi0pi0 1 10000000 0 etap_pi+pi-pi0 1 10000000 0 etap_omegag 1 10000000 0 omega_etag 1 10000000 0 > sim_log_10M & disown
  */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <vector>
#include <map>
#include <algorithm>  // for_each
//#include <initializer_list>  // C++11, usage of -std=gnu++11 or -std=c++11 required

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <THStack.h>
#include <TList.h>

typedef std::map<int, const char*> IntCharMap;
typedef std::pair<int, const char*> ICPair;
typedef std::map<int, const char*>::iterator ICIter;

static const double MASS_PROTON = 938.272;
static int count = 0;

void prepare_hist(TH1 *h, const char* x_name, const char* y_name = "#Events", Int_t color = 3)
{
	h->GetXaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelSize(.045);
	h->GetXaxis()->SetTitleSize(.048);
	h->GetXaxis()->SetTitleFont(42);
	h->GetXaxis()->SetTitle(x_name);
	//h->GetYaxis()->SetRange(10,320);
	h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelSize(.045);
	h->GetYaxis()->SetTitleSize(.048);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitle(y_name);
	h->GetZaxis()->SetLabelFont(42);
	h->GetZaxis()->SetLabelSize(.042);
	h->GetZaxis()->SetTitleSize(.035);
	h->GetZaxis()->SetTitleFont(42);
	h->SetLineColor(color);
	//h->SetTitle("");
	h->GetYaxis()->SetTitleOffset(1.3);
	h->GetYaxis()->SetLabelOffset(.008);  // per cent of pad width; standard is .005
	h->GetYaxis()->SetDecimals();  //show e. g. 1.0 instead of just 1 (same decimals for every label)
	/*h->GetYaxis()->SetNoExponent(false);  //show exponent for very large or small values (seems only to work for values smaller than e-5)
	h->GetXaxis()->SetTitleOffset(1.05);
	h->GetXaxis()->SetLabelOffset(.01);  // per cent of pad width; standard is .005
	*/
}

TH1F* proton_energy(const char* file, Int_t color)
{
	// declare the histogram BEFORE opening the root file, otherwise it will cause a segfault for some strange reason
	char name[4];
	sprintf(name, "h%d", count++);  // prevent "Potential memory leak"
	TH1F* h = new TH1F(name, "", 1000, 0, 1000);

	TFile f(file, "READ");
	if (!f.IsOpen()) {
		fprintf(stderr, "Error opening file %s: %s\n", file, strerror(errno));
		exit(1);
	}

	TTree* MCTree = (TTree*)f.Get("data");
	if (!MCTree) {
		perror("Error opening TTree 'data'");
		exit(1);
	}
	Long64_t treeSize = MCTree->GetEntries();
	printf("%d events in file %s\n", treeSize, file);

	Int_t part;
	Int_t pid[20];
	Double_t fE[20];
	Float_t pE;

	MCTree->SetMakeClass(1);
    MCTree->SetBranchAddress("Particles", &part);
	MCTree->SetBranchAddress("Particles.pid", pid);
	MCTree->SetBranchAddress("Particles.fE", fE);

	for (UInt_t i = 0; i < treeSize; i++) {
		MCTree->GetEntry(i);
		pE = 1000*fE[1] - MASS_PROTON;  // proton has always the index 1
		(*h).Fill(pE);
	}

	prepare_hist(h, "E_{p} [GeV]", "#Events", color);
	f.Close();

	return h;
}

TH1F* proton_theta(const char* file, Int_t color)
{
	char name[4];
	sprintf(name, "h%d", count++);
	TH1F* h = new TH1F(name, "", 120, 0, 60);

	TFile f(file, "READ");
	if (!f.IsOpen()) {
		fprintf(stderr, "Error opening file %s: %s\n", file, strerror(errno));
		exit(1);
	}

	TTree* MCTree = (TTree*)f.Get("data");
	if (!MCTree) {
		perror("Error opening TTree 'data'");
		exit(1);
	}
	Long64_t treeSize = MCTree->GetEntries();
	printf("%d events in file %s\n", treeSize, file);

	Int_t part;
	Int_t pid[20];
	Double_t fE[20];
	Double_t fPx[20];
	Double_t fPy[20];
	Double_t fPz[20];
	TLorentzVector p4;

	MCTree->SetMakeClass(1);
    MCTree->SetBranchAddress("Particles", &part);
	MCTree->SetBranchAddress("Particles.pid", pid);
	MCTree->SetBranchAddress("Particles.fP.fX", fPx);
	MCTree->SetBranchAddress("Particles.fP.fY", fPy);
	MCTree->SetBranchAddress("Particles.fP.fZ", fPz);
	MCTree->SetBranchAddress("Particles.fE", fE);

	for (UInt_t i = 0; i < treeSize; i++) {
		MCTree->GetEntry(i);
		p4.SetXYZT(1000*fPx[1], 1000*fPy[1], 1000*fPz[1], 1000*fE[1]);  // proton has always the index 1
		(*h).Fill(p4.Theta()*TMath::RadToDeg());
	}

	prepare_hist(h, "#vartheta_{p} [#circ]", "#Events", color);
	f.Close();

	return h;
}

THStack* energies_stack(const char* file, Int_t* color, std::vector<int> partIdx)//std::initializer_list<int> partIdx)
{
	const int nParticles = partIdx.size();
	char name[5];
	sprintf(name, "h%d", count++);

	THStack *hs = new THStack(name, "");

	TH1F* h[nParticles];
	for (int i = 0; i < nParticles; i++) {
		sprintf(name, "h%d.%d", count-1, i);
		h[i] = new TH1F(name, "", 1000, 0, 1000);
	}

	TFile f(file, "READ");
	if (!f.IsOpen()) {
		fprintf(stderr, "Error opening file %s: %s\n", file, strerror(errno));
		exit(1);
	}

	TTree* MCTree = (TTree*)f.Get("data");
	if (!MCTree) {
		perror("Error opening TTree 'data'");
		exit(1);
	}
	Long64_t treeSize = MCTree->GetEntries();
	printf("%d events in file %s\n", treeSize, file);

	Int_t part;
	Int_t pid[20];
	Double_t fE[20];
	Double_t fPx[20];
	Double_t fPy[20];
	Double_t fPz[20];
	TLorentzVector p4[nParticles];

	MCTree->SetMakeClass(1);
    MCTree->SetBranchAddress("Particles", &part);
	MCTree->SetBranchAddress("Particles.pid", pid);
	MCTree->SetBranchAddress("Particles.fP.fX", fPx);
	MCTree->SetBranchAddress("Particles.fP.fY", fPy);
	MCTree->SetBranchAddress("Particles.fP.fZ", fPz);
	MCTree->SetBranchAddress("Particles.fE", fE);

	for (UInt_t i = 0; i < treeSize; i++) {
		MCTree->GetEntry(i);
		for (int i = 0; i < nParticles; i++) {
			p4[i].SetXYZT(1000*fPx[partIdx[i]], 1000*fPy[partIdx[i]], 1000*fPz[partIdx[i]], 1000*fE[partIdx[i]]);
			(*h[i]).Fill(p4[i].T());
		}
	}

	for (int i = 0; i < nParticles; i++) {
		prepare_hist(h[i], "E [GeV]", "#Events", color[i]);
		h[i]->SetFillColor(color[i]);
		hs->Add(h[i]);
	}
	f.Close();

	// the following calls let root crash, Paint method has to be calles first -.-
	//hs->GetXaxis()->SetTitle("E [#GeV]");
	//hs->GetYaxis()->SetTitle("#Events");

	return hs;
}

TList* theta_vs_energy(const char* file, std::vector<int> partIdx)
{
	const int nParticles = partIdx.size();
	char name[5];
	//sprintf(name, "h%d", count++);
	count++;

	TList *l = new TList();

	TH2F* h[nParticles];
	for (int i = 0; i < nParticles; i++) {
		sprintf(name, "h%d.%d", count-1, i);
		h[i] = new TH2F(name, "", 200, 0, 1000, 180, 0, 180);
	}

	TFile f(file, "READ");
	if (!f.IsOpen()) {
		fprintf(stderr, "Error opening file %s: %s\n", file, strerror(errno));
		exit(1);
	}

	TTree* MCTree = (TTree*)f.Get("data");
	if (!MCTree) {
		perror("Error opening TTree 'data'");
		exit(1);
	}
	Long64_t treeSize = MCTree->GetEntries();
	printf("%d events in file %s\n", treeSize, file);

	Int_t part;
	Int_t pid[20];
	Double_t fE[20];
	Double_t fPx[20];
	Double_t fPy[20];
	Double_t fPz[20];
	TLorentzVector p4[nParticles];

	MCTree->SetMakeClass(1);
    MCTree->SetBranchAddress("Particles", &part);
	MCTree->SetBranchAddress("Particles.pid", pid);
	MCTree->SetBranchAddress("Particles.fP.fX", fPx);
	MCTree->SetBranchAddress("Particles.fP.fY", fPy);
	MCTree->SetBranchAddress("Particles.fP.fZ", fPz);
	MCTree->SetBranchAddress("Particles.fE", fE);

	for (UInt_t i = 0; i < treeSize; i++) {
		MCTree->GetEntry(i);
		for (int i = 0; i < nParticles; i++) {
			p4[i].SetXYZT(1000*fPx[partIdx[i]], 1000*fPy[partIdx[i]], 1000*fPz[partIdx[i]], 1000*fE[partIdx[i]]);
			(*h[i]).Fill(p4[i].T(), p4[i].Theta()*TMath::RadToDeg());
		}
	}

	for (int i = 0; i < nParticles; i++) {
		prepare_hist(h[i], "E [GeV]", "#vartheta [#circ]");
		l->Add(h[i]);
	}
	f.Close();

	return l;
}

TH1F* etapEnergy_etap_eeg(const char* file)
{
	TFile f(file, "READ");
	if (!f.IsOpen()) {
		fprintf(stderr, "Error opening file %s: %s\n", file, strerror(errno));
		exit(1);
	}

	TTree* MCTree = (TTree*)f.Get("data");
	if (!MCTree) {
		perror("Error opening TTree 'data'");
		exit(1);
	}
	Long64_t treeSize = MCTree->GetEntries();
	printf("%d events in file %s\n", treeSize, file);

	TH1F* h = new TH1F("h", "E_{#eta'} #eta'#rightarrowe^{+}e^{-}#gamma", 1000, 0, 1000);
	TLorentzVector virtGP4, elecP4, posiP4;

	Int_t part;
	Int_t pid[20];
	Double_t fE[20];
	Double_t fPx[20];
	Double_t fPy[20];
	Double_t fPz[20];

	MCTree->SetMakeClass(1);
	MCTree->SetBranchAddress("Particles", &part);
	MCTree->SetBranchAddress("Particles.pid", pid);
	MCTree->SetBranchAddress("Particles.fP.fX", fPx);
	MCTree->SetBranchAddress("Particles.fP.fY", fPy);
	MCTree->SetBranchAddress("Particles.fP.fZ", fPz);
	MCTree->SetBranchAddress("Particles.fE", fE);

	for (UInt_t i = 0; i < treeSize; i++) {
		MCTree->GetEntry(i);
		virtGP4.SetXYZT(-1000*fPx[3], -1000*fPy[3], -1000*fPz[3], 1000*fE[3]);
		elecP4.SetXYZT(1000*fPx[5], 1000*fPy[5], 1000*fPz[5], 1000*fE[5]);
		posiP4.SetXYZT(1000*fPx[6], 1000*fPy[6], 1000*fPz[6], 1000*fE[6]);
		(*h).Fill((elecP4+posiP4+virtGP4).T());
	}

	prepare_hist(h, "E_{#eta'} [GeV]");
	f.Close();

	return h;
}

int main(int argc, char **argv)
{
	char buffer[50];  // buffer for temporary operations
	double hMax;  // temporary variable to store maximum value of the histogram
	int iMax;  // index of histogram with maximum
	TH1 *h_tmp;  // for temporary histogram usage
	int j;  // counter used for several plots
	TIter *iter;  // Iterator for TList, used to iterate through THStack

	// vectors used for dynamically changes for histogram stacking and file naming
	std::vector<int> indices;
	std::vector<const char*> particles;
	std::vector<const char*> names;

	enum chan {
		etap_pi0pi0eta,
		etap_pi0pi0pi0,
		etap_pipipi0,
		etap_omegag,
		etap_eeg,
		omega_etag,
		unknown
	};

	IntCharMap channel;
	channel.insert(ICPair(etap_pi0pi0eta, "etap_pi0pi0eta"));
	channel.insert(ICPair(etap_pi0pi0pi0, "etap_pi0pi0pi0"));
	channel.insert(ICPair(etap_pipipi0, "etap_pi+pi-pi0"));
	channel.insert(ICPair(etap_omegag, "etap_omegag"));
	channel.insert(ICPair(etap_eeg, "etap_e+e-g"));
	channel.insert(ICPair(omega_etag, "omega_etag"));

	IntCharMap identifier;
	identifier.insert(ICPair(etap_pi0pi0eta, "etap_pi0pi0eta"));
	identifier.insert(ICPair(etap_pi0pi0pi0, "etap_pi0pi0pi0"));
	identifier.insert(ICPair(etap_pipipi0, "etap_pipipi0"));
	identifier.insert(ICPair(etap_omegag, "etap_omegag"));
	identifier.insert(ICPair(etap_eeg, "etap_eeg"));
	identifier.insert(ICPair(omega_etag, "omega_etag"));

	IntCharMap legend;
	legend.insert(ICPair(etap_pi0pi0eta, "#eta'#rightarrow#pi^{0}#pi^{0}#eta"));
	legend.insert(ICPair(etap_pi0pi0pi0, "#eta'#rightarrow#pi^{0}#pi^{0}#pi^{0}"));
	legend.insert(ICPair(etap_pipipi0, "#eta'#rightarrow#pi^{+}#pi^{-}#pi^{0}"));
	legend.insert(ICPair(etap_omegag, "#eta'#rightarrow#omega#gamma"));
	legend.insert(ICPair(etap_eeg, "#eta'#rightarrowe^{+}e^{-}#gamma"));
	legend.insert(ICPair(omega_etag, "#omega#rightarrow#eta#gamma"));

	// general settings: set canvas background to white and hide stat box
	gStyle->SetCanvasColor(0);
	gStyle->SetOptStat(0);

	// change contour for 2d plots
	const Int_t nRGBs = 5;
	const Int_t nCont = 255;

	Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
	gStyle->SetNumberContours(nCont);

	// create canvas for all plots
	TCanvas *c = new TCanvas("c", "1D Plots", 10, 10, 700, 600);
	c->SetBorderSize(2);
	c->SetFrameFillColor(0);
	c->SetFrameBorderMode(0);
	c->SetLeftMargin(.12);
	c->SetRightMargin(.05);
	c->SetBottomMargin(.12);
	c->SetTopMargin(.05);
	
	TCanvas *c2 = new TCanvas("c2", "2D Plots", 10, 10, 700, 600);
	c2->SetBorderSize(2);
	c2->SetFrameFillColor(0);
	c2->SetFrameBorderMode(0);
	c2->SetLeftMargin(.12);
	c2->SetRightMargin(.13);
	c2->SetBottomMargin(.12);
	c2->SetTopMargin(.1);
	// char array containing the file names
	const int nChannels = channel.size();
	const int nFiles = 1;//5;
	char sim_files[nChannels*nFiles][100];
	const char* path = "files";//"/data/simulation/background/channels/new_triggerTesting";
	int n = 0;
	for (ICIter it = channel.begin(); it != channel.end(); ++it)
		for (UInt_t i = 1; i <= nFiles; i++)
			sprintf(sim_files[n++], "%s/sim_%s_%02d.root", path, it->second, i);

	Int_t color[6] = {kRed+1, kAzure, kGreen+2, kOrange-3, kSpring-8, kCyan-3};

	// proton energies
	c->cd();
	TH1F* h_Ep[nChannels];
	for (UInt_t i = 0; i < nChannels; i++)
		h_Ep[i] = proton_energy(sim_files[i], color[i]);
	TLegend *leg = new TLegend(.62, .6, .92, .92);
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	// determine histogram with maximum y value for drawing (a bit messy, but for the moment it works)
	hMax = 0;
	for (ICIter it = legend.begin(); it != legend.end(); ++it)
		if (hMax < h_Ep[it->first]->GetBinContent(h_Ep[it->first]->GetMaximumBin())) {
			hMax = h_Ep[it->first]->GetBinContent(h_Ep[it->first]->GetMaximumBin());
			iMax = it->first;
		}
	h_Ep[iMax]->Draw();
	for (ICIter it = legend.begin(); it != legend.end(); ++it) {
		h_Ep[it->first]->Draw("SAME");
		sprintf(buffer, "E_{p} %s", it->second);
		leg->AddEntry(h_Ep[it->first], buffer, "l");
	}
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/proton_energy.pdf");

	// proton theta
	c->Clear();
	TH1F* h_thetap[nChannels];
	for (UInt_t i = 0; i < nChannels; i++)
		h_thetap[i] = proton_theta(sim_files[i], color[i]);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	hMax = 0;
	for (ICIter it = legend.begin(); it != legend.end(); ++it)
		if (hMax < h_thetap[it->first]->GetBinContent(h_thetap[it->first]->GetMaximumBin())) {
			hMax = h_thetap[it->first]->GetBinContent(h_thetap[it->first]->GetMaximumBin());
			iMax = it->first;
		}
	h_thetap[iMax]->Draw();
	for (ICIter it = legend.begin(); it != legend.end(); ++it) {
		h_thetap[it->first]->Draw("SAME");
		sprintf(buffer, "#vartheta_{p} %s", it->second);
		leg->AddEntry(h_thetap[it->first], buffer, "l");
	}
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/proton_theta.pdf");

	// energies decay particles
	// channel eta' --> pi0 pi0 eta
	c->Clear();
	//int indices[] = {6, 7, 8, 9, 10, 11};
	indices.resize(6);
	// shorter way of assigning values, possible since C++11
	indices = {6, 7, 8, 9, 10, 11};
	particles.resize(6);
	particles = {"E_{#gamma,1}(#pi^{0}_{1})", "E_{#gamma,2}(#pi^{0}_{1})", "E_{#gamma,3}(#pi^{0}_{2})", "E_{#gamma,4}(#pi^{0}_{2})", "E_{#gamma,5}(#eta)", "E_{#gamma,6}(#eta)"};
	THStack *hs_E_etap_pi0pi0eta = energies_stack(sim_files[0], color, indices);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	iter = new TIter(hs_E_etap_pi0pi0eta->GetHists());
	j = 0;
	// no idea how to access the histogram from the iterator, maybe I have to define some struct for this...
	//std::for_each (iter.Begin(), TIter::End(), leg->AddEntry((TH1*)iter, particles[j++], "l"));
	// I somehow prefer the following way
	//TObject *o;
	//while (o = iter.Next())
	//	leg->AddEntry((TH1*)o, particles[j++], "l");
	//TH1 *h_tmp; 
	while (h_tmp = (TH1*)iter->Next())
		leg->AddEntry(h_tmp, particles[j++], "l");
	hs_E_etap_pi0pi0eta->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	hs_E_etap_pi0pi0eta->GetXaxis()->SetTitle("E [GeV]");
	hs_E_etap_pi0pi0eta->GetYaxis()->SetTitle("#Events");
	hs_E_etap_pi0pi0eta->Draw();
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/energies_etap_pi0pi0eta.pdf");
	// this fuckin' root object ownership :-@
	//delete[] indices;
	//delete[] particles;

	// channel eta' --> pi0 pi0 pi0
	c->Clear();
	indices = {6, 7, 8, 9, 10, 11};
	particles = {"E_{#gamma,1}(#pi^{0}_{1})", "E_{#gamma,2}(#pi^{0}_{1})", "E_{#gamma,3}(#pi^{0}_{2})", "E_{#gamma,4}(#pi^{0}_{2})", "E_{#gamma,5}(#pi^{0}_{3})", "E_{#gamma,6}(#pi^{0}_{3})"};
	THStack *hs_E_etap_3pi0 = energies_stack(sim_files[1], color, indices);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	delete iter;
	iter = new TIter(hs_E_etap_3pi0->GetHists());
	j = 0;
	while (h_tmp = (TH1*)iter->Next())
		leg->AddEntry(h_tmp, particles[j++], "l");
	hs_E_etap_3pi0->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	hs_E_etap_3pi0->GetXaxis()->SetTitle("E [GeV]");
	hs_E_etap_3pi0->GetYaxis()->SetTitle("#Events");
	hs_E_etap_3pi0->Draw();
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/energies_etap_pi0pi0pi0.pdf");

	// channel eta' --> pi+ pi- pi0
	c->Clear();
	indices.resize(4);
	particles.resize(4);
	indices = {3, 4, 6, 7};
	particles = {"E_{#pi^{+}}", "E_{#pi^{-}}", "E_{#gamma,1}(#pi^{0})", "E_{#gamma,2}(#pi^{0})"};
	THStack *hs_E_etap_pipipi0 = energies_stack(sim_files[2], color, indices);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	delete iter;
	iter = new TIter(hs_E_etap_pipipi0->GetHists());
	j = 0;
	while (h_tmp = (TH1*)iter->Next())
		leg->AddEntry(h_tmp, particles[j++], "l");
	hs_E_etap_pipipi0->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	hs_E_etap_pipipi0->GetXaxis()->SetTitle("E [GeV]");
	hs_E_etap_pipipi0->GetYaxis()->SetTitle("#Events");
	hs_E_etap_pipipi0->Draw();
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/energies_etap_pipipi0.pdf");

	// channel eta' --> omega gamma
	c->Clear();
	indices = {4, 6, 7, 8};
	particles = {"E_{#gamma,1}", "E_{#gamma,2}(#omega)", "E_{#gamma,3}(#eta)", "E_{#gamma,4}(#eta)"};
	THStack *hs_E_etap_omegag = energies_stack(sim_files[3], color, indices);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	delete iter;
	iter = new TIter(hs_E_etap_omegag->GetHists());
	j = 0;
	while (h_tmp = (TH1*)iter->Next())
		leg->AddEntry(h_tmp, particles[j++], "l");
	hs_E_etap_omegag->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	hs_E_etap_omegag->GetXaxis()->SetTitle("E [GeV]");
	hs_E_etap_omegag->GetYaxis()->SetTitle("#Events");
	hs_E_etap_omegag->Draw();
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/energies_etap_omegag.pdf");

	// channel eta' --> e+ e- gamma
	c->Clear();
	indices.resize(3);
	particles.resize(3);
	indices = {5, 6, 4};
	particles = {"E_{e^{+}}", "E_{e^{-}}", "E_{#gamma}"};
	THStack *hs_E_etap_eeg = energies_stack(sim_files[4], color, indices);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	delete iter;
	iter = new TIter(hs_E_etap_eeg->GetHists());
	j = 0;
	while (h_tmp = (TH1*)iter->Next())
		leg->AddEntry(h_tmp, particles[j++], "l");
	hs_E_etap_eeg->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	hs_E_etap_eeg->GetXaxis()->SetTitle("E [GeV]");
	hs_E_etap_eeg->GetYaxis()->SetTitle("#Events");
	hs_E_etap_eeg->Draw();
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/energies_etap_eeg.pdf");

	// channel omega --> eta gamma
	c->Clear();
	indices = {4, 5, 6};
	particles = {"E_{#gamma,1}(#omega)", "E_{#gamma,2}(#eta)", "E_{#gamma,3}(#eta)"};
	THStack *hs_E_omega_etag = energies_stack(sim_files[5], color, indices);
	leg->Clear();
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);
	delete iter;
	iter = new TIter(hs_E_omega_etag->GetHists());
	j = 0;
	while (h_tmp = (TH1*)iter->Next())
		leg->AddEntry(h_tmp, particles[j++], "l");
	hs_E_omega_etag->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	hs_E_omega_etag->GetXaxis()->SetTitle("E [GeV]");
	hs_E_omega_etag->GetYaxis()->SetTitle("#Events");
	hs_E_omega_etag->Draw();
	leg->Draw("SAME");
	c->Update();
	c->Print("plots/energies_omega_etag.pdf");

	// theta vs energy angle 2d plots for single final state particles
	// channel eta' --> pi0 pi0 eta
	c2->cd();
	indices.resize(6);
	indices = {6, 7, 8, 9, 10, 11};
	particles.resize(6);
	particles = {"#gamma_{1}(#pi^{0}_{1})", "#gamma_{2}(#pi^{0}_{1})", "#gamma_{3}(#pi^{0}_{2})", "#gamma_{4}(#pi^{0}_{2})", "#gamma_{5}(#eta)", "#gamma_{6}(#eta)"};
	names.resize(6);
	names = {"gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6"};
	//TList *l_etap_pi0pi0eta = theta_vs_energy(sim_files[etap_pi0pi0eta], indices);
	//iter = new TIter(l_etap_pi0pi0eta);
	delete iter;
	iter = new TIter(theta_vs_energy(sim_files[etap_pi0pi0eta], indices));
	j = 0;
	while (h_tmp = (TH2F*)iter->Next()) {
		c2->Clear();
		h_tmp->SetTitle(particles[j]);
		h_tmp->Draw("COLZ");
		sprintf(buffer, "plots/theta_vs_energy_%s_%s.pdf", identifier.find(etap_pi0pi0eta)->second, names[j++]);
		c2->Update();
		c2->Print(buffer);
	}

	// channel eta' --> pi0 pi0 pi0
	c2->cd();
	indices = {6, 7, 8, 9, 10, 11};
	particles = {"#gamma_{1}(#pi^{0}_{1})", "#gamma_{2}(#pi^{0}_{1})", "#gamma_{3}(#pi^{0}_{2})", "#gamma_{4}(#pi^{0}_{2})", "#gamma_{5}(#pi^{0}_{3})", "#gamma_{6}(#pi^{0}_{3})"};
	names = {"gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6"};
	delete iter;
	iter = new TIter(theta_vs_energy(sim_files[etap_pi0pi0pi0], indices));
	j = 0;
	while (h_tmp = (TH2F*)iter->Next()) {
		c2->Clear();
		h_tmp->SetTitle(particles[j]);
		h_tmp->Draw("COLZ");
		sprintf(buffer, "plots/theta_vs_energy_%s_%s.pdf", identifier.find(etap_pi0pi0pi0)->second, names[j++]);
		c2->Update();
		c2->Print(buffer);
	}

	// channel eta' --> pi+ pi- pi0
	c2->cd();
	indices.resize(4);
	particles.resize(4);
	names.resize(4);
	indices = {3, 4, 6, 7};
	particles = {"#pi^{+}", "#pi^{-}", "#gamma_{1}", "#gamma_{2}"};
	names = {"pi1", "pi2", "gamma1", "gamma2"};
	delete iter;
	iter = new TIter(theta_vs_energy(sim_files[etap_pipipi0], indices));
	j = 0;
	while (h_tmp = (TH2F*)iter->Next()) {
		c2->Clear();
		h_tmp->SetTitle(particles[j]);
		h_tmp->Draw("COLZ");
		sprintf(buffer, "plots/theta_vs_energy_%s_%s.pdf", identifier.find(etap_pipipi0)->second, names[j++]);
		c2->Update();
		c2->Print(buffer);
	}

	// channel eta' --> omega gamma
	c2->cd();
	indices = {4, 6, 7, 8};
	particles = {"#gamma_{1}", "#gamma_{2}(#omega)", "#gamma_{3}(#eta)", "#gamma_{4}(#eta)"};
	names = {"gamma1", "gamma2", "gamma3", "gamma4"};
	delete iter;
	iter = new TIter(theta_vs_energy(sim_files[etap_omegag], indices));
	j = 0;
	while (h_tmp = (TH2F*)iter->Next()) {
		c2->Clear();
		h_tmp->SetTitle(particles[j]);
		h_tmp->Draw("COLZ");
		sprintf(buffer, "plots/theta_vs_energy_%s_%s.pdf", identifier.find(etap_omegag)->second, names[j++]);
		c2->Update();
		c2->Print(buffer);
	}

	// channel eta' --> e+ e- gamma
	c2->cd();
	indices.resize(3);
	particles.resize(3);
	names.resize(3);
	indices = {5, 6, 4};
	particles = {"e^{+}", "e^{-}", "#gamma"};
	names = {"e1", "e2", "gamma"};
	delete iter;
	iter = new TIter(theta_vs_energy(sim_files[etap_eeg], indices));
	j = 0;
	while (h_tmp = (TH2F*)iter->Next()) {
		c2->Clear();
		h_tmp->SetTitle(particles[j]);
		h_tmp->Draw("COLZ");
		sprintf(buffer, "plots/theta_vs_energy_%s_%s.pdf", identifier.find(etap_eeg)->second, names[j++]);
		c2->Update();
		c2->Print(buffer);
	}

	// channel omega --> eta gamma
	c2->cd();
	indices = {4, 5, 6};
	particles = {"#gamma_{1}(#omega)", "#gamma_{2}(#eta)", "#gamma_{3}(#eta)"};
	names = {"gamma1", "gamma2", "gamma3"};
	delete iter;
	iter = new TIter(theta_vs_energy(sim_files[omega_etag], indices));
	j = 0;
	while (h_tmp = (TH2F*)iter->Next()) {
		c2->Clear();
		h_tmp->SetTitle(particles[j]);
		h_tmp->Draw("COLZ");
		sprintf(buffer, "plots/theta_vs_energy_%s_%s.pdf", identifier.find(omega_etag)->second, names[j++]);
		c2->Update();
		c2->Print(buffer);
	}

	return 0;
}

