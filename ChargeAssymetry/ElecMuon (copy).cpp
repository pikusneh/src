
#ifdef _WIN32
#include <direct.h>
// MSDN recommends against using getcwd & chdir names
#define cwd _getcwd
//#define cd _chdir
#else
#include "unistd.h"
#define cwd getcwd
//#define cd chdir
#endif
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <signal.h> //  our new library 
#include <algorithm>
//#include <iostream>
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
#include<cstring>
#include <TLorentzVector.h>
#include<TChain.h>

#include "Event.hh"
#include "TopVar.hh"
using namespace TopAnalysis;
using namespace std;




void ChargeAssymetrySingleChan(Event *events, double &CAssym,string SelChParticle="Electron",int nJetCut=4, int nBJetCut=1, int nPhoCut=1)
{
	


	//double TotalWeight = 0;
	
	//long int numevents=events->tr->GetEntries();
	//cout << "numevents" << numevents << std::endl;
	
	
	
	ParticleType Ch_particle;
	Ch_particle.type=SelChParticle;
	
	//for(long int i=0; i<numevents; i++) // Event Loop Start
	//{
	//	events->tr->GetEntry(i);
		
		
			
		//TotalWeight = events->evtWeight * events->PUweight * events->muEffWeight * events->eleEffWeight * events->btagWeight_1a * events->prefireSF;
		
		//cout<<"TotalWeight:: "<<TotalWeight<<endl;
		bool evtpass=0;
		
		
		
		if(Ch_particle.type=="Electron")
		{
		 	evtpass=events->passPresel_Ele;	
		
		}
		
		else if(Ch_particle.type=="Muon")
		{
			evtpass=events->passPresel_Mu;
		
		}
		
		//decide the observed channel based on the lepton charge at the final state
		
		
		
		
       	if ( evtpass && events->nJet >= nJetCut && events->nBJet >= nBJetCut && events->nPho == nPhoCut)
        	{
        		cout<<"TotalWeight After Cut:: "<<TotalWeight<<"\n"<<endl;
        		
        	
        	}
		
	
	} // End of Event loop
	
	
	

}


int main()
{

	gROOT->SetBatch(true);

	TChain *tr = new TChain("AnalysisTree");

	tr->Add("../data/TTGamma_SingleLept_2016_AnalysisNtuple.root");
	//tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
	//tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
	//tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");
        double CAssym=0;
	Event *events=new Event(tr);
	double TotalWeight = 0;
	
	long int numevents=events->tr->GetEntries();
	cout << "numevents" << numevents << std::endl;
	
	
	

	
	for(long int i=0; i<numevents; i++) // Event Loop Start
	{
		events->tr->GetEntry(i);
		TotalWeight = events->evtWeight * events->PUweight * events->muEffWeight * events->eleEffWeight * events->btagWeight_1a * events->prefireSF;
		
		ChargeAssymetrySingleChan(events, CAssym);//,std::string="Electron",int nJetCut=4, int nBJetCut=1, int nPhoCut=1);



	}









