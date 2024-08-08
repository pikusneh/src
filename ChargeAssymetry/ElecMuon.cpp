
#ifdef _WIN32
#include <direct.h>
// MSDN recommends against using getcwd & chdir names
#define cwd _getcwd
// #define cd _chdir
#else
#include "unistd.h"
#define cwd getcwd
// #define cd chdir
#endif
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <signal.h> //  our new library
#include <algorithm>
// #include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <new>
#include <climits>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TH3F.h"
#include "TH3.h"
#include "TStyle.h"
#include "TAttFill.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TPostScript.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TPaletteAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TVectorD.h"
#include "TBox.h"
#include "TArrow.h"
#include <functional>
#include <TApplication.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <THttpServer.h>
#include <TMemFile.h>
#include <TNtuple.h>
#include "TROOT.h"
#include <cstring>
#include <TLorentzVector.h>
#include <TChain.h>

#include "Event.hh"
#include "TopVar.hh"
#include "Analysis.hh"
using namespace TopAnalysis;
using namespace std;
double ww = 0;

int main(int argc, char **argv)
{

	// gROOT->SetBatch(true);
	TApplication app("app", &argc, argv);
	TChain *tr = new TChain("AnalysisTree");

	// tr->Add("../data/TTGamma_SingleLept_2016_AnalysisNtuple.root");
	tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
	tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
	tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");
	double CAssym = 0;
	// Global variables to select the mass and PT Bin

	// Variables For Electron Channel

	const int NTopMassBin_Ele_ch = 8;
	const int NTopPtBin_Ele_ch = 7;
	const int NPhotonPtBin_Ele_ch = 5;
	std::vector<double> TopMassBin_Ele_ch{300, 400, 500, 600, 700, 800, 900, 1000};
	std::vector<double> TopPtBin_Ele_ch{-200, 0, 200, 400, 600, 800, 1000};
	std::vector<double> PhotonPtBin_Ele_ch{-100.0, 100.0, 300.0, 500.0, 600.0};

	std::vector<double> TopMassBin_sel_event_eta_plus_Ele_ch(NTopMassBin_Ele_ch, 0);
	std::vector<double> TopMassBin_sel_event_eta_minus_Ele_ch(NTopMassBin_Ele_ch, 0);

	std::vector<double> TopPtBin_sel_event_eta_plus_Ele_ch(NTopPtBin_Ele_ch, 0);
	std::vector<double> TopPtBin_sel_event_eta_minus_Ele_ch(NTopPtBin_Ele_ch, 0);

	std::vector<double> PhotonPtBin_sel_event_eta_plus_Ele_ch(NPhotonPtBin_Ele_ch, 0);
	std::vector<double> PhotonPtBin_sel_event_eta_minus_Ele_ch(NPhotonPtBin_Ele_ch, 0);
	std::vector<double> AvgTopMass_Ele_ch(NPhotonPtBin_Ele_ch, 0);

	// Variables for the Muon Channel
	const int NTopMassBin_Mu_ch = 8;
	const int NTopPtBin_Mu_ch = 7;
	const int NPhotonPtBin_Mu_ch = 5;

	std::vector<double> TopMassBin_Mu_ch{300, 400, 500, 600, 700, 800, 900, 1000};
	std::vector<double> TopPtBin_Mu_ch{-200, 0, 200, 400, 600, 800, 1000};
	std::vector<double> PhotonPtBin_Mu_ch{-100.0, 100.0, 300.0, 500.0, 600.0};

	std::vector<double> TopMassBin_sel_event_eta_plus_Mu_ch(NTopMassBin_Mu_ch, 0);
	std::vector<double> TopMassBin_sel_event_eta_minus_Mu_ch(NTopMassBin_Mu_ch, 0);

	std::vector<double> TopPtBin_sel_event_eta_plus_Mu_ch(NTopPtBin_Mu_ch, 0);
	std::vector<double> TopPtBin_sel_event_eta_minus_Mu_ch(NTopPtBin_Mu_ch, 0);

	std::vector<double> PhotonPtBin_sel_event_eta_plus_Mu_ch(NPhotonPtBin_Mu_ch, 0);
	std::vector<double> PhotonPtBin_sel_event_eta_minus_Mu_ch(NPhotonPtBin_Mu_ch, 0);

	// Mix/Either channel
	const int NTopMassBin_Ele_Mu_ch = 8;
	const int NTopPtBin_Ele_Mu_ch = 7;
	const int NPhotonPtBin_Ele_Mu_ch = 5;
	std::vector<double> TopMassBin_Ele_Mu_ch{300, 400, 500, 600, 700, 800, 900, 1000};
	std::vector<double> TopPtBin_Ele_Mu_ch{-200, 0, 200, 400, 600, 800, 1000};
	std::vector<double> PhotonPtBin_Ele_Mu_ch{-100.0, 100.0, 300.0, 500.0, 600.0};

	std::vector<double> TopMassBin_sel_event_eta_plus_Ele_Mu_ch(NTopMassBin_Ele_Mu_ch, 0);
	std::vector<double> TopMassBin_sel_event_eta_minus_Ele_Mu_ch(NTopMassBin_Ele_Mu_ch, 0);

	std::vector<double> TopPtBin_sel_event_eta_plus_Ele_Mu_ch(NTopPtBin_Ele_Mu_ch, 0);
	std::vector<double> TopPtBin_sel_event_eta_minus_Ele_Mu_ch(NTopPtBin_Ele_Mu_ch, 0);

	std::vector<double> PhotonPtBin_sel_event_eta_plus_Ele_Mu_ch(NPhotonPtBin_Ele_Mu_ch, 0);
	std::vector<double> PhotonPtBin_sel_event_eta_minus_Ele_Mu_ch(NPhotonPtBin_Ele_Mu_ch, 0);

	vector<std::string> Bintype{"TopMass", "TopPt", "PhoPt"};
	vector<std::string> Ch_type{"Ele", "Mu", "Ele_Mu"};

	Event *events = new Event(tr);

	long int numevents = events->tr->GetEntries();
	cout << "numevents:: " << numevents << std::endl;

	// events->SetDebugMode();
	for (long int i = 0; i < 50000; i++) // Event Loop Start
	{
		events->tr->GetEntry(i);

		// Electron channel calculations
		events->SetChannelAndCuts(Ch_type[0], 4, 1, 1);
		bool masstrue = events->FindSelectedEvents(TopMassBin_Ele_ch, TopMassBin_sel_event_eta_plus_Ele_ch, TopMassBin_sel_event_eta_minus_Ele_ch, Bintype[0]);
		bool pttrue = events->FindSelectedEvents(TopPtBin_Ele_ch, TopPtBin_sel_event_eta_plus_Ele_ch, TopPtBin_sel_event_eta_minus_Ele_ch, Bintype[1]);
		bool pttruePh = events->FindSelectedEvents(PhotonPtBin_Ele_ch, PhotonPtBin_sel_event_eta_plus_Ele_ch, PhotonPtBin_sel_event_eta_minus_Ele_ch, Bintype[2]);

		// Muon channel calculations
		events->SetChannelAndCuts(Ch_type[1], 4, 1, 1);
		bool masstrueMu = events->FindSelectedEvents(TopMassBin_Mu_ch, TopMassBin_sel_event_eta_plus_Mu_ch, TopMassBin_sel_event_eta_minus_Mu_ch, Bintype[0]);
		bool pttrueMu = events->FindSelectedEvents(TopPtBin_Mu_ch, TopPtBin_sel_event_eta_plus_Mu_ch, TopPtBin_sel_event_eta_minus_Mu_ch, Bintype[1]);
		bool pttruePhMu = events->FindSelectedEvents(PhotonPtBin_Mu_ch, PhotonPtBin_sel_event_eta_plus_Mu_ch, PhotonPtBin_sel_event_eta_minus_Mu_ch, Bintype[2]);

		// Electron or Muon channel calculations
		events->SetChannelAndCuts(Ch_type[2], 4, 1, 1);
		bool masstrueEleMu = events->FindSelectedEvents(TopMassBin_Ele_Mu_ch, TopMassBin_sel_event_eta_plus_Ele_Mu_ch, TopMassBin_sel_event_eta_minus_Ele_Mu_ch, Bintype[0]);
		bool pttrueEleMu = events->FindSelectedEvents(TopPtBin_Ele_Mu_ch, TopPtBin_sel_event_eta_plus_Ele_Mu_ch, TopPtBin_sel_event_eta_minus_Ele_Mu_ch, Bintype[1]);
		bool pttruePhEleMu = events->FindSelectedEvents(PhotonPtBin_Ele_Mu_ch, PhotonPtBin_sel_event_eta_plus_Ele_Mu_ch, PhotonPtBin_sel_event_eta_minus_Ele_Mu_ch, Bintype[2]);

	} // End of Event loop

	// Plotting
	Analysis *ana = new Analysis();
	/*  int bin = h14->FindBin(NetMass4);
	h14->SetBinContent(bin, Ac_ele4);
	h14->SetBinError(bin, Delta_Ac_ele4);*/

	// ################################################################ Visual OUTPUTS ######################################################
	cout << "####################################################################################################################\n\n";
	cout << "                                              Channel:: Ele"
		 << "\n\n";
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
	vector<double> Cassym_top_mass_ele;
	vector<double> ErrorCassym_top_mass_ele;
	events->ChargeAssymetry(TopMassBin_sel_event_eta_plus_Ele_ch, TopMassBin_sel_event_eta_minus_Ele_ch, Cassym_top_mass_ele, ErrorCassym_top_mass_ele);

	for (int i = 0; i < TopMassBin_sel_event_eta_plus_Ele_ch.size(); i++)
	{
		cout << std::fixed << setprecision(7) << "MASS BIN::" << i << "\tFor +ve Rapidity diff::" << TopMassBin_sel_event_eta_plus_Ele_ch[i];
		cout << std::fixed << setprecision(7) << "\tFor -ve Rapidity diff::" << TopMassBin_sel_event_eta_minus_Ele_ch[i] << "\t Charge Assymetry:: " << Cassym_top_mass_ele[i] << endl;
		if (i > 0 && (Cassym_top_mass_ele[i] != -99999.0 && ErrorCassym_top_mass_ele[i] != -99999.0))
		{
			ana->ChargeAssy_vs_TopMassBin->SetBinContent(i, Cassym_top_mass_ele[i]);
			// ana->ChargeAssy_vs_TopMassBin->SetBinContent(Cassym_top_mass_ele[i]);
			ana->ChargeAssy_vs_TopMassBin->SetBinError(i, ErrorCassym_top_mass_ele[i]);
		}
	}
	cout << "-----------------------------------------------------------------------------------------------------------------------\n";
	// vector<double>Cassym_top_pt_ele = events->ChargeAssymetry(TopPtBin_sel_event_eta_plus_Ele_ch,TopPtBin_sel_event_eta_minus_Ele_ch);
	vector<double> Cassym_top_pt_ele;
	vector<double> ErrorCassym_top_pt_ele;
	events->ChargeAssymetry(TopPtBin_sel_event_eta_plus_Ele_ch, TopPtBin_sel_event_eta_minus_Ele_ch, Cassym_top_pt_ele, ErrorCassym_top_pt_ele);

	for (int i = 0; i < TopPtBin_sel_event_eta_plus_Ele_ch.size(); i++)
	{
		cout << std::fixed << setprecision(7) << "Top Pt BIN::" << i << "\tFor +ve value::" << TopPtBin_sel_event_eta_plus_Ele_ch[i];
		cout << std::fixed << setprecision(7) << "\t\tFor -ve Rapidity diff::" << TopPtBin_sel_event_eta_minus_Ele_ch[i] << "\t Charge Assymetry:: " << Cassym_top_pt_ele[i] << endl;
	}
	cout << "-------------------------------------------------------------------------------------------------------------------------\n";
	// vector<double>Cassym_pho_pt_ele = events->ChargeAssymetry(PhotonPtBin_sel_event_eta_plus_Ele_ch,PhotonPtBin_sel_event_eta_minus_Ele_ch);
	/*
	for(int i=0;i<PhotonPtBin_sel_event_eta_plus_Ele_ch.size();i++)
	{
		  cout<<std::fixed<<setprecision(6)<<std::setfill('0')<<"Pho Pt BIN::"<<i<<"\tFor +ve Rapidity diff::"<<PhotonPtBin_sel_event_eta_plus_Ele_ch[i];
		  cout<<std::fixed<<setprecision(6)<<std::setfill('0')<<"\tFor -ve Rapidity diff::"<<PhotonPtBin_sel_event_eta_minus_Ele_ch[i]<<"\t Charge Assymetry:: "<<Cassym_pho_pt_ele[i]<<endl;

	}
   cout<<"####################################################################################################################\n\n";

   cout<<"                                                 Channel:: Mu"<<"\n\n";

   cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
   // vector<double>Cassym_top_mass_mu = events->ChargeAssymetry(TopMassBin_sel_event_eta_plus_Mu_ch,TopMassBin_sel_event_eta_minus_Mu_ch);
	for(int i=0;i<TopMassBin_sel_event_eta_plus_Mu_ch.size();i++)
	{
		  cout<<std::fixed<<setprecision(7)<<std::setfill('0')<<"MASS BIN::"<<i<<"\tFor +ve Rapidity diff::"<<TopMassBin_sel_event_eta_plus_Mu_ch[i];
		  cout<<std::fixed<<setprecision(7)<<std::setfill('0')<<"\tFor -ve Rapidity diff::"<<TopMassBin_sel_event_eta_minus_Mu_ch[i]<<"\t Charge Assymetry:: "<<Cassym_top_mass_mu[i]<<endl;

	}
	cout<<"-------------------------------------------------------------------------------------------------------------------\n";
   // vector<double>Cassym_top_pt_mu = events->ChargeAssymetry(TopPtBin_sel_event_eta_plus_Mu_ch,TopPtBin_sel_event_eta_minus_Mu_ch);
	for(int i=0;i<TopPtBin_sel_event_eta_plus_Mu_ch.size();i++)
	{
		  cout<<std::fixed<<setprecision(7)<<std::setfill('0')<<"Top Pt BIN::"<<i<<"\tFor +ve Rapidity diff::"<<TopPtBin_sel_event_eta_plus_Mu_ch[i];
		  cout<<std::fixed<<setprecision(7)<<std::setfill('0')<<"\tFor -ve Rapidity diff::"<<TopPtBin_sel_event_eta_minus_Mu_ch[i]<<"\t Charge Assymetry:: "<<Cassym_top_pt_mu[i]<<endl;

	}
	cout<<"--------------------------------------------------\n";
	//vector<double>Cassym_pho_pt_mu = events->ChargeAssymetry(PhotonPtBin_sel_event_eta_plus_Mu_ch,PhotonPtBin_sel_event_eta_minus_Mu_ch);
	for(int i=0;i<PhotonPtBin_sel_event_eta_plus_Mu_ch.size();i++)
	{
		  cout<<std::fixed<<setprecision(6)<<std::setfill('0')<<"Pho Pt BIN::"<<i<<"\tFor +ve Rapidity diff::"<<PhotonPtBin_sel_event_eta_plus_Mu_ch[i];
		  cout<<std::fixed<<setprecision(6)<<std::setfill('0')<<"\tFor -ve Rapidity diff::"<<PhotonPtBin_sel_event_eta_minus_Mu_ch[i]<<"\t Charge Assymetry:: "<<Cassym_pho_pt_mu[i]<<endl;

	}
	 cout<<"####################################################################################################################\n\n";
	 cout<<"                                                 Channel:: Ele_Mu"<<"\n\n";
	 cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
	// vector<double>Cassym_top_mass_elmu = events->ChargeAssymetry(TopMassBin_sel_event_eta_plus_Ele_Mu_ch,TopMassBin_sel_event_eta_minus_Ele_Mu_ch);
	 for(int i=0;i<TopMassBin_sel_event_eta_plus_Ele_Mu_ch.size();i++)
	 {
		  cout<<std::fixed<<setprecision(7)<<"MASS BIN::"<<i<<"\tFor +ve Rapidity diff::"<<TopMassBin_sel_event_eta_plus_Ele_Mu_ch[i];
		  cout<<std::fixed<<setprecision(7)<<"\tFor -ve Rapidity diff::"<<TopMassBin_sel_event_eta_minus_Ele_Mu_ch[i]<<"\t Charge Assymetry:: "<<Cassym_top_mass_elmu[i]<<endl;

	 }
	 cout<<"------------------------------------------------------------------------------------------------------------------\n";
	 //vector<double>Cassym_top_pt_elmu = events->ChargeAssymetry(TopPtBin_sel_event_eta_plus_Ele_Mu_ch,TopPtBin_sel_event_eta_minus_Ele_Mu_ch);
	 for(int i=0;i<TopPtBin_sel_event_eta_plus_Ele_Mu_ch.size();i++)
	 {
		  cout<<std::fixed<<setprecision(7)<<"Top Pt BIN::"<<i<<"\tFor +ve Rapidity diff::"<<TopPtBin_sel_event_eta_plus_Ele_Mu_ch[i];
		  cout<<std::fixed<<setprecision(7)<<"\tFor -ve Rapidity diff::"<<TopPtBin_sel_event_eta_minus_Ele_Mu_ch[i]<<"\t Charge Assymetry:: "<<Cassym_top_pt_elmu[i]<<endl;

	 }
	 cout<<"----------------------------------------------------------------------------------------------------------------\n";
	 //vector<double>Cassym_pho_pt_elmu = events->ChargeAssymetry(PhotonPtBin_sel_event_eta_plus_Ele_Mu_ch,PhotonPtBin_sel_event_eta_minus_Ele_Mu_ch);
	 for(int i=0;i<PhotonPtBin_sel_event_eta_plus_Ele_Mu_ch.size();i++)
	 {
		  cout<<std::fixed<<setprecision(7)<<"Pho Pt BIN::"<<i<<"\tFor +ve Rapidity diff::"<<PhotonPtBin_sel_event_eta_plus_Ele_Mu_ch[i];
		  cout<<std::fixed<<setprecision(7)<<"\tFor -ve Rapidity diff::"<<PhotonPtBin_sel_event_eta_minus_Ele_Mu_ch[i]<<"\t Charge Assymetry:: "<<Cassym_pho_pt_elmu[i]<<endl;

	 } */
	auto *cnv = new TCanvas();
	// ana->ChargeAssy_vs_TopMassBin->Draw("E1");
	gStyle->SetEndErrorSize(1);
	gStyle->SetErrorX(.3);
	ana->ChargeAssy_vs_TopMassBin->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
	ana->ChargeAssy_vs_TopMassBin->GetXaxis()->SetTitleSize(0.04);
	ana->ChargeAssy_vs_TopMassBin->GetXaxis()->SetLabelSize(0.03);
	ana->ChargeAssy_vs_TopMassBin->GetYaxis()->SetTitle("A_{c_{ele}}");
	ana->ChargeAssy_vs_TopMassBin->GetYaxis()->SetTitleOffset(1.1);
	ana->ChargeAssy_vs_TopMassBin->GetYaxis()->SetTitleSize(0.04);
	ana->ChargeAssy_vs_TopMassBin->GetYaxis()->SetLabelSize(0.03);

	ana->ChargeAssy_vs_TopMassBin->SetMarkerSize(1);
	ana->ChargeAssy_vs_TopMassBin->SetMarkerColor(kRed);
	ana->ChargeAssy_vs_TopMassBin->SetLineWidth(1);
	ana->ChargeAssy_vs_TopMassBin->Draw("E1");
	app.Run(true);
}
