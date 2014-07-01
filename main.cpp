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
// accessing files and directories
#include <sys/types.h>
#include <sys/stat.h>

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
// for handling of different final state particles
typedef std::map<int, std::vector<int>> IntVecintMap;
typedef std::map<int, std::vector<const char*>> IntVecharMap;
typedef std::pair<int, std::vector<int>> IViPair;
typedef std::pair<int, std::vector<const char*>> IVcPair;
typedef std::map<int, std::vector<int>>::iterator IViIter;
typedef std::map<int, std::vector<int>>::const_iterator constIViIter;  // const_iterator needed while iterating through const map
typedef std::map<int, std::vector<const char*>>::iterator IVcIter;
// one map containing all 4-vectors (reading all information from MC Tree file only once required)
typedef std::map<int, std::vector<std::vector<TLorentzVector>>> IntP4Map;
typedef std::pair<int, std::vector<std::vector<TLorentzVector>>> IP4Pair;
typedef std::map<int, std::vector<std::vector<TLorentzVector>>>::iterator IP4Iter;
typedef std::vector<std::vector<TLorentzVector>> VVP4;
typedef VVP4::const_iterator VVP4Iter;

//static const double MASS_PROTON = 938.272;
static int count = 0;  // counter used for individual histogram naming
static const int READ_LIMIT = 1000000;  // limit to which number events are read per file; number smaller than zero for all events, e. g. -1

int collect_particles(IntP4Map& p4, const IntVecintMap& idx, const char files[][100], const int nFiles = 1);  // structure of two-dimensional char array has to be char a[][n] or, equivalent, char (a*)[n]
void prepare_hist(TH1 *h, const char* x_name, const char* y_name = "#Events", Int_t color = 3);
void prepare_hist(THStack *h, const char* x_name, const char* y_name = "#Events");
TList* energies(const VVP4& p4, const std::vector<int> partIdx);
TList* thetas(const VVP4& p4, const std::vector<int> partIdx);
TList* theta_vs_energy(const VVP4& p4, const std::vector<int> partIdx);
TH1F* etapEnergy_etap_eeg(const char* file);

int main(int argc, char **argv)
{
	char buffer[50];  // buffer for temporary operations
	double hMax;  // temporary variable to store maximum value of the histogram
	int iMax;  // index of histogram with maximum
	TH1 *h_tmp;  // for temporary histogram usage
	int j, p;  // counter used for several plots etc.
	TIter *iter;  // Iterator for TList, used to iterate through THStack and TList

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
		omega_eepi0,
		unknown
	};

	IntCharMap channel;
	channel.insert(ICPair(etap_pi0pi0eta, "etap_pi0pi0eta"));
	channel.insert(ICPair(etap_pi0pi0pi0, "etap_pi0pi0pi0"));
	channel.insert(ICPair(etap_pipipi0, "etap_pi+pi-pi0"));
	channel.insert(ICPair(etap_omegag, "etap_omegag"));
	channel.insert(ICPair(etap_eeg, "etap_e+e-g"));
	channel.insert(ICPair(omega_etag, "omega_etag"));
	channel.insert(ICPair(omega_eepi0, "omega_e+e-pi0"));

	IntCharMap identifier;
	identifier.insert(ICPair(etap_pi0pi0eta, "etap_pi0pi0eta"));
	identifier.insert(ICPair(etap_pi0pi0pi0, "etap_pi0pi0pi0"));
	identifier.insert(ICPair(etap_pipipi0, "etap_pipipi0"));
	identifier.insert(ICPair(etap_omegag, "etap_omegag"));
	identifier.insert(ICPair(etap_eeg, "etap_eeg"));
	identifier.insert(ICPair(omega_etag, "omega_etag"));
	identifier.insert(ICPair(omega_eepi0, "omega_eepi0"));

	// legend entry when only one value of a channel is used in a plot
	IntCharMap legend;
	legend.insert(ICPair(etap_pi0pi0eta, "#eta'#rightarrow#pi^{0}#pi^{0}#eta"));
	legend.insert(ICPair(etap_pi0pi0pi0, "#eta'#rightarrow#pi^{0}#pi^{0}#pi^{0}"));
	legend.insert(ICPair(etap_pipipi0, "#eta'#rightarrow#pi^{+}#pi^{-}#pi^{0}"));
	legend.insert(ICPair(etap_omegag, "#eta'#rightarrow#omega#gamma"));
	legend.insert(ICPair(etap_eeg, "#eta'#rightarrowe^{+}e^{-}#gamma"));
	legend.insert(ICPair(omega_etag, "#omega#rightarrow#eta#gamma"));
	legend.insert(ICPair(omega_eepi0, "#omega#rightarrowe^{+}e^{-}#pi^{0}"));

	// maps containing final state particle information for more automated behaviour
	IntVecintMap indicesFS;
	// recoil proton has always the index 1
	indicesFS.insert(IViPair(etap_pi0pi0eta, {6, 7, 8, 9, 10, 11, 1}));
	indicesFS.insert(IViPair(etap_pi0pi0pi0, {6, 7, 8, 9, 10, 11, 1}));
	indicesFS.insert(IViPair(etap_pipipi0, {3, 4, 6, 7, 1}));
	indicesFS.insert(IViPair(etap_omegag, {4, 6, 7, 8, 1}));
	indicesFS.insert(IViPair(etap_eeg, {5, 6, 4, 1}));
	indicesFS.insert(IViPair(omega_etag, {4, 5, 6, 1}));
	indicesFS.insert(IViPair(omega_eepi0, {7, 8, 5, 6, 1}));
	IntVecharMap particlesFS;
	particlesFS.insert(IVcPair(etap_pi0pi0eta, {"#gamma_{1}(#pi^{0}_{1})", "#gamma_{2}(#pi^{0}_{1})", "#gamma_{3}(#pi^{0}_{2})", "#gamma_{4}(#pi^{0}_{2})", "#gamma_{5}(#eta)", "#gamma_{6}(#eta)", "p"}));
	particlesFS.insert(IVcPair(etap_pi0pi0pi0, {"#gamma_{1}(#pi^{0}_{1})", "#gamma_{2}(#pi^{0}_{1})", "#gamma_{3}(#pi^{0}_{2})", "#gamma_{4}(#pi^{0}_{2})", "#gamma_{5}(#pi^{0}_{3})", "#gamma_{6}(#pi^{0}_{3})", "p"}));
	particlesFS.insert(IVcPair(etap_pipipi0, {"#pi^{+}", "#pi^{-}", "#gamma_{1}", "#gamma_{2}", "p"}));
	// for omega --> eta g
	//particlesFS.insert(IVcPair(etap_omegag, {"#gamma_{1}", "#gamma_{2}(#omega)", "#gamma_{3}(#eta)", "#gamma_{4}(#eta)", "p"}));
	// for omega --> pi0 g
	particlesFS.insert(IVcPair(etap_omegag, {"#gamma_{1}", "#gamma_{2}(#omega)", "#gamma_{3}(#pi^{0})", "#gamma_{4}(#pi^{0})", "p"}));
	particlesFS.insert(IVcPair(etap_eeg, {"e^{+}", "e^{-}", "#gamma", "p"}));
	particlesFS.insert(IVcPair(omega_etag, {"#gamma_{1}(#omega)", "#gamma_{2}(#eta)", "#gamma_{3}(#eta)", "p"}));
	particlesFS.insert(IVcPair(omega_eepi0, {"e^{+}", "e^{-}", "#gamma_{1}(#pi^{0})", "#gamma_{2}(#pi^{0})", "p"}));
	IntVecharMap namesFS;
	namesFS.insert(IVcPair(etap_pi0pi0eta, {"gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6", "proton"}));
	namesFS.insert(IVcPair(etap_pi0pi0pi0, {"gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6", "proton"}));
	namesFS.insert(IVcPair(etap_pipipi0, {"pi1", "pi2", "gamma1", "gamma2", "proton"}));
	namesFS.insert(IVcPair(etap_omegag, {"gamma1", "gamma2", "gamma3", "gamma4", "proton"}));
	namesFS.insert(IVcPair(etap_eeg, {"e1", "e2", "gamma", "proton"}));
	namesFS.insert(IVcPair(omega_etag, {"gamma1", "gamma2", "gamma3", "proton"}));
	namesFS.insert(IVcPair(omega_eepi0, {"e1", "e2", "gamma1", "gamma2", "proton"}));

	std::cout << "[INFO] Channel initialisation done!" << std::endl
	<< "The following channels will be analysed:" << std::endl;
	for (ICIter it = channel.begin(); it != channel.end(); ++it)
		printf( "  %s, %d final state particles\n", it->second, indicesFS.find(it->first)->second.size());
	std::cout << std::endl;

	// char array containing the file names
	const int nChannels = channel.size();
	const int nFiles = 1;//5;
	char sim_files[nChannels*nFiles][100];
	const char* path = "/data/simulation/background/channels/new_triggerTesting_5M";
	const char* ext = "png";
	const char* save = "plots";
	// READ_LIMIT is set above before method declarations
	std::cout << "The following data path will be used: " << path << std::endl
	<< "Number of files per channel: " << nFiles << std::endl
	<< "Maximum number of events read per file: ";
	if (READ_LIMIT < 0)
		std::cout << "all events" << std::endl;
	else
		std::cout << READ_LIMIT << std::endl;
	for (ICIter it = channel.begin(); it != channel.end(); ++it)
		for (UInt_t i = 1; i <= nFiles; i++)
			sprintf(sim_files[nFiles*it->first+i-1], "%s/sim_%s_%02d.root", path, it->second, i);
	std::cout << "Plots will be saved as " << ext << std::endl;
	//std::cout << "Save directory is " << save << std::endl << std::endl;

	// check if specified save directory exists
	struct stat s;
	if (stat(save, &s)) {  // directory cannot be accessed, seems to not exist, create it
		std::cout << "Create save directory " << save << std::endl;
		if (!mkdir(save, S_IRWXU))
			std::cout << "Directory created successfully." << std::endl << std::endl;
		else {
			perror("Creating directory failed: ");
			exit(1);
		}
	} else if (s.st_mode & S_IFDIR)
		printf("Save directory is \"%s\"\n\n", save);
	else {
		printf("%s is not a directory!\n", save);
		exit(1);
	}

	// colors which will be used for the 1D histograms
	Int_t color[7] = {kRed+1, kAzure, kGreen+2, kOrange-3, kSpring-8, kCyan-3, kRed+2};

	// general settings: set canvas background to white and hide stat box
	gStyle->SetCanvasColor(0);
	gStyle->SetOptStat(0);

	// change contour for 2D plots
	const Int_t nRGBs = 5;
	const Int_t nCont = 255;

	Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
	gStyle->SetNumberContours(nCont);

	// create canvas for all plots 1D plots
	TCanvas *c = new TCanvas("c", "1D Plots", 10, 10, 700, 600);
	c->SetBorderSize(2);
	c->SetFrameFillColor(0);
	c->SetFrameBorderMode(0);
	c->SetLeftMargin(.12);
	c->SetRightMargin(.05);
	c->SetBottomMargin(.12);
	c->SetTopMargin(.05);
	// and one for the 2D ones
	TCanvas *c2 = new TCanvas("c2", "2D Plots", 10, 10, 700, 600);
	c2->SetBorderSize(2);
	c2->SetFrameFillColor(0);
	c2->SetFrameBorderMode(0);
	c2->SetLeftMargin(.12);
	c2->SetRightMargin(.13);
	c2->SetBottomMargin(.12);
	c2->SetTopMargin(.1);

	// gather all needed particle information (4-vectors) from the generated files
	IntP4Map p4FS;
	for (IViIter it = indicesFS.begin(); it != indicesFS.end(); ++it)
		p4FS.insert(IP4Pair(it->first, std::vector<std::vector<TLorentzVector>>(it->second.size())));
	if (!collect_particles(p4FS, indicesFS, sim_files, nFiles))
		printf("\n[INFO] All particles collected!\n\n");
	else
		printf("\nSome error occurred...\n\n");

	// legend used in some of the histograms
	TLegend *leg = new TLegend(.64, .6, .94, .94);
	leg->SetFillColor(0);
	leg->SetBorderSize(1);
	leg->SetTextFont(42);
	leg->SetTextSize(.035);

	// use THStack for plotting several histograms into one canvas
	THStack *hs;

	// energies final state
	std::cout << "[INFO] Create plots for energies in the final state" << std::endl;
	c->cd();
	THStack *hs_E = new THStack("hs_E", "");
	TList *l_p = new TList();  // list containing proton energy histograms
	TList *l_cc = new TList();  // list containing 2D histograms of #patricles in CB vs. theta constrained ESum
	TList *l_ct = new TList();  // list containing 2D histograms of #patricles in TAPS vs. theta constrained ESum
	for (IViIter it = indicesFS.begin(); it != indicesFS.end(); ++it) {
		c->Clear();
		leg->Clear();
		// change legend height according to number of entries that it fits better
		if (it->second.size() > 5)
			leg->SetY1NDC(.45);
		else
			leg->SetY1NDC(.58);
		leg->SetHeader("Energies FS");
		if (!hs)
			delete hs;
		sprintf(buffer, "h%d%d", count-1, it->first);
		hs = new THStack(buffer, "");
		if (!iter)
			delete iter;
		iter = new TIter(energies(p4FS.find(it->first)->second, it->second));
		j = 0;
		while (h_tmp = (TH1*)iter->Next()) {
			if (strstr(h_tmp->GetTitle(), "p")) {  // proton
				h_tmp->SetTitle(legend.find(it->first)->second);  // store current decay channel in histogram title to use it later for the legend entry
				l_p->Add(h_tmp);
			} else if (strstr(h_tmp->GetTitle(), "Energy Sum")) {
				c->Clear();
				h_tmp->SetLineColor(color[1]);
				h_tmp->SetTitle("");
				h_tmp->Draw();
				c->Update();
				sprintf(buffer, "%s/energy_sum_%s.%s", save, identifier.find(it->first)->second, ext);
				c->Print(buffer);
			} else if (strstr(h_tmp->GetTitle(), "ESum thetaConstr")) {
				h_tmp->SetLineColor(color[2]);  // as the full energy sum is saved before the constrained one in the list, first draw both combined before deleting the ESum of the FS from the canvas
				h_tmp->SetTitle("");
				h_tmp->Draw("SAME");
				c->Update();
				sprintf(buffer, "%s/energy_sum_combined_%s.%s", save, identifier.find(it->first)->second, ext);
				c->Print(buffer);
				c->Clear();
				h_tmp->SetLineColor(color[1]);
				h_tmp->SetTitle("");
				h_tmp->Draw();
				c->Update();
				sprintf(buffer, "%s/energy_sum_thetaConstr_%s.%s", save, identifier.find(it->first)->second, ext);
				c->Print(buffer);
				// now change the color and fill style and add the histogram to a stack
				h_tmp->SetLineColor(color[it->first]);
				h_tmp->SetFillColor(color[it->first]);
				hs_E->Add(h_tmp);
			} else if (strstr(h_tmp->GetTitle(), "ESum_nPart")) {
				c2->cd();
				c2->Clear();
				char det[5] = "test";
				if (strstr(h_tmp->GetTitle(), "CB")) {
					strcpy(det, "CB");
					h_tmp->SetTitle("");
					l_cc->Add(h_tmp);
				}
				else if (strstr(h_tmp->GetTitle(), "TAPS")) {
					strcpy(det, "TAPS");
					h_tmp->SetTitle("");
					l_ct->Add(h_tmp);
				}
				h_tmp->Draw("COLZ");
				c2->Update();
				sprintf(buffer, "%s/nPart_vs_ESumConstr_%s_%s.%s", save, det, identifier.find(it->first)->second, ext);
				c2->Print(buffer);
				c->cd();
			} else {  // histograms of decay particles don't have a histogram title
				h_tmp->SetLineColor(color[j]);
				h_tmp->SetFillColor(color[j]);
				hs->Add(h_tmp);
				leg->AddEntry(h_tmp, particlesFS.find(it->first)->second[j++], "l");
			}
		}
		// draw the energies of the decay particles
		c->Clear();
		hs->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
		prepare_hist(hs, "E [MeV]");
		hs->Draw();
		leg->Draw("SAME");
		c->Update();
		sprintf(buffer, "%s/energies_%s.%s", save, identifier.find(it->first)->second, ext);
		c->Print(buffer);
	}
	// after iterating over all channels draw the energies of the protons from the different channels now
	c->Clear();
	leg->Clear();
	leg->SetY1NDC(.6);
	if (!iter)
		delete iter;
	iter = new TIter(l_p);
	// first find the histogram with the maximum bin (used for proper drawing)
	j = 0;
	hMax = 0;
	while (h_tmp = (TH1*)iter->Next()) {
		h_tmp->SetLineColor(color[j++]);
		leg->AddEntry(h_tmp, h_tmp->GetTitle(), "l");
		h_tmp->SetTitle("");
		if (hMax < h_tmp->GetBinContent(h_tmp->GetMaximumBin())) {
			hMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
			iMax = j-1;}
		//h_tmp->Draw("SAME");
	}
	l_p->At(iMax)->Draw();  // draw max hist
	delete iter;
	iter = new TIter(l_p);
	// now draw all the histograms
	while (h_tmp = (TH1*)iter->Next())
		h_tmp->Draw("SAME");
	leg->Draw("SAME");
	c->Update();
	sprintf(buffer, "%s/proton_energies.%s", save, ext);
	c->Print(buffer);
	// draw the stack with the theta constrained energy sum
	c->Clear();
	hs_E->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
	prepare_hist(hs_E, "E_{sum} CB [MeV]");
	hs_E->Draw();
	// move legend to the left
	leg->SetX1NDC(.17);
	leg->SetY1NDC(.6);
	leg->SetX2NDC(.47);
	leg->SetY2NDC(.94);
	leg->Draw("SAME");  // legend should be the same as in the above case
	c->Update();
	sprintf(buffer, "%s/energy_sums_theta_constraint.%s", save, ext);
	c->Print(buffer);
	// change legend position back
	leg->SetX1NDC(.64);
	leg->SetY1NDC(.6);
	leg->SetX2NDC(.94);
	leg->SetY2NDC(.94);
	// lastly draw a summed up 2D histogram of all particle count vs. ESum hists
	// CB
	c2->cd();
	c2->Clear();
	h_tmp = (TH2F*)l_cc->First()->Clone("h_tmp");
	h_tmp->Reset();
	h_tmp->Merge(l_cc);
	h_tmp->Draw("COLZ");
	c2->Update();
	sprintf(buffer, "%s/nPart_vs_ESumConstr_sum_CB.%s", save, ext);
	c2->Print(buffer);
	// TAPS
	c2->Clear();
	h_tmp = (TH2F*)l_ct->First()->Clone("h_tmp");
	h_tmp->Reset();
	h_tmp->Merge(l_ct);
	h_tmp->Draw("COLZ");
	c2->Update();
	sprintf(buffer, "%s/nPart_vs_ESumConstr_sum_TAPS.%s", save, ext);
	c2->Print(buffer);
	c->cd();

	delete hs_E;
	delete l_p;
	delete l_cc;
	delete l_ct;

	// theta angles final state
	std::cout << "[INFO] Create plots for theta angles in the final state" << std::endl;
	l_p = new TList();
	for (IViIter it = indicesFS.begin(); it != indicesFS.end(); ++it) {
		c->Clear();
		leg->Clear();
		// change legend height according to number of entries that it fits better
		if (it->second.size() > 5)
			leg->SetY1NDC(.45);
		else
			leg->SetY1NDC(.58);
		leg->SetHeader("#vartheta FS");
		if (!hs)
			delete hs;
		sprintf(buffer, "h%d%d", count-1, it->first);
		hs = new THStack(buffer, "");
		if (!iter)
			delete iter;
		iter = new TIter(thetas(p4FS.find(it->first)->second, it->second));
		j = 0;
		while (h_tmp = (TH1*)iter->Next()) {
			if (strstr(h_tmp->GetTitle(), "p")) {  // proton
				h_tmp->SetTitle(legend.find(it->first)->second);  // store current decay channel in histogram title to use it later for the legend entry
				l_p->Add(h_tmp);
			} else {  // histograms of decay particles don't have a histogram title
				h_tmp->SetLineColor(color[j]);
				h_tmp->SetFillColor(color[j]);
				hs->Add(h_tmp);
				leg->AddEntry(h_tmp, particlesFS.find(it->first)->second[j++], "l");
			}
		}
		// draw the angles of the decay particles
		c->Clear();
		hs->Paint();  // TAxis objects of THStack are only created when the Paint function is called, otherwise a segfault occurs
		prepare_hist(hs, "#vartheta [#circ]");
		hs->Draw();
		leg->Draw("SAME");
		c->Update();
		sprintf(buffer, "%s/thetas_%s.%s", save, identifier.find(it->first)->second, ext);
		c->Print(buffer);
	}
	// after iterating over all channels draw the energies of the protons from the different channels now
	c->Clear();
	leg->Clear();
	leg->SetY1NDC(.6);
	if (!iter)
		delete iter;
	iter = new TIter(l_p);
	j = 0;
	hMax = 0;
	while (h_tmp = (TH1*)iter->Next()) {
		h_tmp->SetLineColor(color[j++]);
		leg->AddEntry(h_tmp, h_tmp->GetTitle(), "l");
		h_tmp->SetTitle("");
		if (hMax < h_tmp->GetBinContent(h_tmp->GetMaximumBin())) {
			hMax = h_tmp->GetBinContent(h_tmp->GetMaximumBin());
			iMax = j-1;}
		//h_tmp->Draw("SAME");
	}
	l_p->At(iMax)->Draw();
	delete iter;
	iter = new TIter(l_p);
	while (h_tmp = (TH1*)iter->Next())
		h_tmp->Draw("SAME");
	leg->Draw("SAME");
	c->Update();
	sprintf(buffer, "%s/proton_thetas.%s", save, ext);
	c->Print(buffer);

	delete l_p;

	// theta vs energy angle 2d plots for single final state particles
	std::cout << "[INFO] Create 2D plots theta angle vs. energy" << std::endl;
	c2->cd();
	for (IViIter it = indicesFS.begin(); it != indicesFS.end(); ++it) {
		if (!iter)
			delete iter;
		//iter = new TIter(theta_vs_energy(sim_files[it->first], it->second));
		iter = new TIter(theta_vs_energy(p4FS.find(it->first)->second, it->second));
		j = 0;
		while (h_tmp = (TH2F*)iter->Next()) {
			c2->Clear();
			h_tmp->SetTitle(particlesFS.find(it->first)->second[j]);
			if (particlesFS.find(it->first)->second[j] == "p") {
				h_tmp->GetYaxis()->SetRangeUser(0, 50);
				h_tmp->SetTitle("");
			}
			h_tmp->Draw("COLZ");
			sprintf(buffer, "%s/theta_vs_energy_%s_%s.%s", save, identifier.find(it->first)->second, namesFS.find(it->first)->second[j++], ext);
			c2->Update();
			c2->Print(buffer);
		}
	}

	return 0;
}


int collect_particles(IntP4Map& p4, const IntVecintMap& idx, const char files[][100], const int nFiles)
{
	printf("[INFO] Start collecting final state particles for %d channels . . .\n\n", p4.size());

	TTree* MCTree;
	Long64_t treeSize;
	Int_t part;
	Int_t pid[20];
	Double_t fE[20];
	Double_t fPx[20];
	Double_t fPy[20];
	Double_t fPz[20];

	for (int n = 0; n < nFiles; n++) {
		for (constIViIter it = idx.begin(); it != idx.end(); ++it) {
			TFile f(files[nFiles*it->first+n], "READ");
			if (!f.IsOpen()) {
				fprintf(stderr, "Error opening file %s: %s\n", files[it->first], strerror(errno));
				exit(1);
			}

			MCTree = (TTree*)f.Get("data");
			if (!MCTree) {
				perror("Error opening TTree 'data'");
				exit(1);
			}
			treeSize = MCTree->GetEntries();
			printf("%d events in file %s\n", treeSize, files[it->first]);
			/* i will be an iterator over an two-dimensional vector of objects of the type TLorentzVector, as p4.find(it->first)->second delivers this vector inside the p4 map regarding the actually processed channel. The outer vector should have the size of the number of final state particles as described in the second part of the indices map, like the vector is initialized earlier in the main method. So the for loop with the iterator i should have e. g. 6 steps for a decay with 6 particles in the final state. The size of the inner vector, which can be accessed via the iterator i, a pointer to this vector, should be zero in the first run because it isn't initialized in the main method. What we want now is to reserve the expected size of this vector that all the events from the tree file can be stored in it without reallocating memory after every few push_backs. Speeds up the system execution time. Therefore we use the initial size (0 for the first run or the amount of stored events in later loop runs) plus the new size. */
			for (std::vector<std::vector<TLorentzVector>>::iterator i = p4.find(it->first)->second.begin(); i != p4.find(it->first)->second.end(); ++i)
				i->reserve(i->size()+treeSize);
				//std::cout << i->size() << std::endl;  // should be zero for all entries as the vector inside the vector isn't initialized yet

			MCTree->SetMakeClass(1);
			MCTree->SetBranchAddress("Particles", &part);
			MCTree->SetBranchAddress("Particles.pid", pid);
			MCTree->SetBranchAddress("Particles.fE", fE);
			MCTree->SetBranchAddress("Particles.fP.fX", fPx);
			MCTree->SetBranchAddress("Particles.fP.fY", fPy);
			MCTree->SetBranchAddress("Particles.fP.fZ", fPz);

			for (int i = 0; i < treeSize; i++) {
				if (i == READ_LIMIT) break;  // limiting events read per file to READ_LIMIT due to heavy ram usage...
				MCTree->GetEntry(i);
				for (int j = 0; j < it->second.size(); j++)
					p4.find(it->first)->second[j].push_back(TLorentzVector(1000*fPx[it->second[j]], 1000*fPy[it->second[j]], 1000*fPz[it->second[j]], 1000*fE[it->second[j]]));
			}

			f.Close();
		}
	}

	std::cout << "Finished processing all files." << std::endl;

	return 0;
}

void prepare_hist(TH1 *h, const char* x_name, const char* y_name, Int_t color)
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

void prepare_hist(THStack *h, const char* x_name, const char* y_name)
{
	h->GetXaxis()->SetLabelFont(42);
	h->GetXaxis()->SetLabelSize(.045);
	h->GetXaxis()->SetTitleSize(.048);
	h->GetXaxis()->SetTitleFont(42);
	h->GetXaxis()->SetTitle(x_name);
	h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->SetLabelSize(.045);
	h->GetYaxis()->SetTitleSize(.048);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitle(y_name);
	h->SetTitle("");
	h->GetYaxis()->SetTitleOffset(1.3);
	h->GetYaxis()->SetLabelOffset(.008);  // per cent of pad width; standard is .005
	h->GetYaxis()->SetDecimals();  //show e. g. 1.0 instead of just 1 (same decimals for every label)
}

TList* energies(const VVP4& p4, const std::vector<int> partIdx)
{
	const int nParticles = partIdx.size();
	char name[5];
	count++;

	TList *l = new TList();

	TH1F* h[nParticles+2];  // number of particles plus two additional energy sum histograms
	for (int i = 0; i < nParticles; i++) {
		sprintf(name, "h%d.%d", count-1, i);
		h[i] = new TH1F(name, "", 1000, 0, 1000);
		prepare_hist(h[i], "E [MeV]", "#Events");
		if (partIdx[i] == 1) {  // mark histogram with proton energy for later usage and set x-axis title to E_p
			h[i]->SetTitle("p");
			prepare_hist(h[i], "E_{p} [MeV]", "#Events");
		}
	}
	sprintf(name, "hes%d", count-1);
	h[nParticles] = new TH1F(name, "Energy Sum", 950, 650, 1600);
	prepare_hist(h[nParticles], "E_{sum} FS [MeV]", "#Events");
	sprintf(name, "hec%d", count-1);
	h[nParticles+1] = new TH1F(name, "ESum thetaConstr", 1600, 0, 1600);
	prepare_hist(h[nParticles+1], "E_{sum} CB [MeV]", "#Events");
	// at last two histograms that count the number of particles in the CB and TAPS range
	sprintf(name, "h2c%d", count-1);
	TH2F *h2 = new TH2F(name, "ESum_nPart_CB", 400, 0, 1600, partIdx.size(), 0, partIdx.size());
	prepare_hist(h2, "E_{sum} CB [MeV]", "#particles CB");
	sprintf(name, "h2t%d", count);
	TH2F *h3 = new TH2F(name, "ESum_nPart_TAPS", 400, 0, 1600, partIdx.size(), 0, partIdx.size());
	prepare_hist(h3, "E_{sum} TAPS [MeV]", "#particles TAPS");

	double esum, esum_constrCB, esum_constrTAPS, e;
	int c, t;  // counter for CB/TAPS particles
	for (int j = 0; j < p4[0].size(); j++) {
		esum = esum_constrCB = esum_constrTAPS = c = t = 0;
		for (int i = 0; i < nParticles; i++) {
			e = p4[i][j].E()-p4[i][j].M();  // indices of particles (which are passed to this method) are coupled to the four-momenta due to collection process, therefore the usage of a simple for-loop with accessing the momenta via i is possible
			h[i]->Fill(e);
			if (partIdx[i] != 1) {  // exclude proton (has always id 1) from energy sum
				esum += e;
				if (p4[i][j].Theta()*TMath::RadToDeg() > 20. && p4[i][j].Theta()*TMath::RadToDeg() < 160.) {  // only particles in CB range
					esum_constrCB += e;
					c++;
				} else if (p4[i][j].Theta()*TMath::RadToDeg() <= 20.) {  // only particles in TAPS
					esum_constrTAPS += e;
					t++;
				}
			}
		}
		h[nParticles]->Fill(esum);
		/* only fill spectra when esum != 0, i. e. c or t counter greater than zero */
		if (c) {
			h[nParticles+1]->Fill(esum_constrCB);
			h2->Fill(esum_constrCB, c);
		}
		if (t)
			h3->Fill(esum_constrTAPS, t);
	}

	for (int i = 0; i < nParticles+2; i++)
		l->Add(h[i]);
	l->Add(h2);
	l->Add(h3);

	return l;
}

TList* thetas(const VVP4& p4, const std::vector<int> partIdx)
{
	const int nParticles = partIdx.size();
	char name[5];
	count++;

	TList *l = new TList();

	TH1F* h[nParticles];
	for (int i = 0; i < nParticles; i++) {
		sprintf(name, "h%d.%d", count-1, i);
		if (partIdx[i] == 1) {  // other dimensions needed for proton theta; mark histogram for later usage
			h[i] = new TH1F(name, "p", 120, 0, 60);
			prepare_hist(h[i], "#vartheta_{p} [#circ]", "#Events");
		} else {
			h[i] = new TH1F(name, "", 360, 0, 180);
			prepare_hist(h[i], "#vartheta [#circ]", "#Events");
		}
	}

	for (int j = 0; j < p4[0].size(); j++)
		for (int i = 0; i < nParticles; i++)
			h[i]->Fill(p4[i][j].Theta()*TMath::RadToDeg());

	for (int i = 0; i < nParticles; i++)
		l->Add(h[i]);

	return l;
}

TList* theta_vs_energy(const VVP4& p4, const std::vector<int> partIdx)
{
	const int nParticles = partIdx.size();
	char name[5];
	count++;

	TList *l = new TList();

	TH2F* h[nParticles];
	for (int i = 0; i < nParticles; i++) {
		sprintf(name, "h%d.%d", count-1, i);
		h[i] = new TH2F(name, "", 200, 0, 1000, 180, 0, 180);
		prepare_hist(h[i], "E [MeV]", "#vartheta [#circ]");
	}

	for (int j = 0; j < p4[0].size(); j++)
		for (int i = 0; i < nParticles; i++)
			h[i]->Fill(p4[i][j].E()-p4[i][j].M(), p4[i][j].Theta()*TMath::RadToDeg());

	for (int i = 0; i < nParticles; i++)
		l->Add(h[i]);

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

	prepare_hist(h, "E_{#eta'} [MeV]");
	f.Close();

	return h;
}

