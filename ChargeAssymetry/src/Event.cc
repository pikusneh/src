
#include "Event.hh"



using namespace std;
namespace TopAnalysis{

	Event::Event(TChain *trev){
	
	///////////////////////////////////////////////////////////////Set Branch Status//////////////////////////////////////////////////////////////////
	tr=dynamic_cast<TChain*>(trev);

	tr->SetBranchStatus("*", 0);
	tr->SetBranchStatus("elePt", 1);
	tr->SetBranchStatus("eleEta", 1);
	tr->SetBranchStatus("elePhi", 1);
	tr->SetBranchStatus("muPt", 1);
	tr->SetBranchStatus("muEta", 1);
	tr->SetBranchStatus("muPhi", 1);
	tr->SetBranchStatus("jetPt", 1);
	tr->SetBranchStatus("phoEt", 1);
	tr->SetBranchStatus("phoEta", 1);
	tr->SetBranchStatus("phoPhi", 1);
	tr->SetBranchStatus("nEle", 1);
	tr->SetBranchStatus("nMu", 1);
	tr->SetBranchStatus("nJet", 1);
	tr->SetBranchStatus("nBJet", 1);
	tr->SetBranchStatus("nPho", 1);
	tr->SetBranchStatus("TopHad_pt", 1);
	tr->SetBranchStatus("TopLep_pt", 1);
	tr->SetBranchStatus("TopHad_eta", 1);
	tr->SetBranchStatus("TopLep_eta", 1);
	tr->SetBranchStatus("TopHad_phi", 1);
	tr->SetBranchStatus("TopLep_phi", 1);
	tr->SetBranchStatus("M_bjj", 1);
	tr->SetBranchStatus("Mt_blgammaMET", 1);
	tr->SetBranchStatus("TopLep_charge", 1);
	tr->SetBranchStatus("passPresel_Ele", 1);
	tr->SetBranchStatus("passPresel_Mu", 1);
	tr->SetBranchStatus("lumis", 1);
	tr->SetBranchStatus("evtWeight", 1);
	tr->SetBranchStatus("btagWeight_1a", 1);
	tr->SetBranchStatus("prefireSF", 1);
	tr->SetBranchStatus("muEffWeight", 1);
	tr->SetBranchStatus("eleEffWeight", 1);
	tr->SetBranchStatus("PUweight", 1);


	/////////////////////////////////////////////////////////////Set Branch Address/////////////////////////////////////////////////////////////////
	tr->SetBranchAddress("elePt", &elePt);
	tr->SetBranchAddress("elePt", &eleEta);
	tr->SetBranchAddress("elePt", &elePhi);
	tr->SetBranchAddress("muPt", &muPt);
	tr->SetBranchAddress("muPt", &muEta);
	tr->SetBranchAddress("muPt", &muPhi);
	tr->SetBranchAddress("jetPt", &jetPt);
	tr->SetBranchAddress("phoEt", &phoEt);
	tr->SetBranchAddress("phoEta", &phoEta);
	tr->SetBranchAddress("phoPhi", &phoPhi);
	tr->SetBranchAddress("nEle", &nEle);
	tr->SetBranchAddress("nMu", &nMu);
	tr->SetBranchAddress("nJet", &nJet);
	tr->SetBranchAddress("nBJet", &nBJet);
	tr->SetBranchAddress("nPho", &nPho);
	tr->SetBranchAddress("Mt_blgammaMET", &Mt_blgammaMET); // toplepmass
	tr->SetBranchAddress("M_bjj", &M_bjj);                 // tophadmass
	tr->SetBranchAddress("TopHad_pt", &TopHad_pt);
	tr->SetBranchAddress("TopLep_pt", &TopLep_pt);
	tr->SetBranchAddress("TopHad_eta", &TopHad_eta);
	tr->SetBranchAddress("TopLep_eta", &TopLep_eta);
	tr->SetBranchAddress("TopHad_phi", &TopHad_phi);
	tr->SetBranchAddress("TopLep_phi", &TopLep_phi);
	tr->SetBranchAddress("TopLep_charge", &TopLep_charge);
	tr->SetBranchAddress("passPresel_Ele", &passPresel_Ele);
	tr->SetBranchAddress("passPresel_Mu", &passPresel_Mu);
	tr->SetBranchAddress("lumis", &lumis);
	tr->SetBranchAddress("btagWeight_1a", &btagWeight_1a);
	tr->SetBranchAddress("evtWeight", &evtWeight);
	tr->SetBranchAddress("prefireSF", &prefireSF);
	tr->SetBranchAddress("muEffWeight", &muEffWeight);
	tr->SetBranchAddress("eleEffWeight", &eleEffWeight);
	tr->SetBranchAddress("PUweight", &PUweight);
	}
	
	Event::~Event(){}
	

	void Event::TopAntiTopPair(TLorentzVector &TopAntiTop,Double_t &RapiDiff)
	{
		
      		RapiDiff=TMath::Abs(GetTop().Rapidity())-TMath::Abs(GetAntiTop().Rapidity());
      		TopAntiTop=GetTop()+GetAntiTop();
      		
	}
	float Event::SetChannelAndCuts(std::string SelChParticle_temp,int nJetCut_temp, int nBJetCut_temp, int nPhoCut_temp){
	
		 SelChParticle=SelChParticle_temp;
		 nJetCut=nJetCut_temp;
		 nBJetCut=nBJetCut_temp;
		 nPhoCut=nPhoCut_temp;  	
		
	}
	bool Event::FindSelectedEvents(std::vector<double> &BinArray,std::vector<double> &sel_event_eta_plus,std::vector<double> &sel_event_eta_minus,std::string BinType)
	{

		ParticleType Ch_particle;
		Ch_particle.type=SelChParticle;

		bool channel=0;
		double type=0.0;

		
		double TotalWeight = this->evtWeight * this->PUweight * this->muEffWeight * this->eleEffWeight * this->btagWeight_1a * this->prefireSF;	
		
		//decide the observed channel based on the lepton charge at the final state
		TLorentzVector TopAntiTop;
		Double_t RapiDiff=0;
		this->TopAntiTopPair(TopAntiTop,RapiDiff);
		
		if(BinType=="TopMass")
			type=TopAntiTop.M();
		else if(BinType=="TopPt")
			type=TopAntiTop.Pt();
		else if(BinType=="PhoPt" && this->nPho>=1)
			{
			type=this->GetPhoton().Pt();

			}
		else return 0;	

		if(SelChParticle=="Ele" && this->passPresel_Ele)
			{
			channel=1;

			}
		else if (SelChParticle=="Mu" && this->passPresel_Mu)
			channel=1;
		else if (SelChParticle=="Ele_Mu" && (this->passPresel_Mu || this->passPresel_Ele))
			channel=1;
		else return 0;		

		if (channel && this->nJet >= nJetCut && this->nBJet >= nBJetCut && this->nPho == nPhoCut)
		{

			if(debug)
			cout<<"RapiDiff :: "<<RapiDiff<< " channel:: "<<Ch_particle.type<<" TotalWeight :: "<<TotalWeight<<endl;
			//cout<<"TotalWeight After Cut:: "<<TotalWeight<<"\n"<<endl;

			for(int i=0; i<BinArray.size();i++)	
			{
			 
			  if(RapiDiff>0)
			  {
			      if(i==0)sel_event_eta_plus[i]=sel_event_eta_plus[i]+TotalWeight;
			      else if(type>=BinArray[i-1] && type<BinArray[i]) 
			      {
			      	      sel_event_eta_plus[i]=sel_event_eta_plus[i]+TotalWeight;
			      	      //NetMass
			      }
			  }
			  else if(RapiDiff<0)
			  {
			      if(i==0)sel_event_eta_minus[i]=sel_event_eta_minus[i]+TotalWeight;
			      else if(type>=BinArray[i-1] && type<BinArray[i]) 
			      {
			      	     sel_event_eta_minus[i]=sel_event_eta_minus[i]+TotalWeight;
			      }
			  }
			  
			
			
			     /* if(i==0)	{
			      		if(RapiDiff>0)
			  		{
			      		  sel_event_eta_plus[i]=sel_event_eta_plus[i]+TotalWeight;
			      		}
			      		else if(RapiDiff<0)
			  		{
			  		  sel_event_eta_minus[i]=sel_event_eta_minus[i]+TotalWeight;
			  		}
			      		
			      }
			      else if(type>=BinArray[i-1] && type<BinArray[i]) {
			      		if (RapiDiff>0)
			  		{
			      		  sel_event_eta_plus[i]=sel_event_eta_plus[i]+TotalWeight;
			      		}
			      		else if (RapiDiff<0)
			  		{
			  		  sel_event_eta_minus[i]=sel_event_eta_minus[i]+TotalWeight;
			  		}
			      		
			      }*/
			      
			
			}

		return 1;  
		}

	   return 0;     	
	}
	void Event::ChargeAssymetry(std::vector<double> &sel_event_eta_plus,std::vector<double> &sel_event_eta_minus,std::vector<double>&CAssym,std::vector<double> &BinError)
	{
		 // std::cout << "Asymmetry for electron channel vs phopt1 = " << Ac_ele1phopt << std::endl;
 
	     int BinRange=sel_event_eta_plus.size(); 
	    // std::vector<double>CAssym(BinRange,-99999);
	     for(int i=0;i<BinRange;i++)
	     {
	        double sum=sel_event_eta_plus[i]+sel_event_eta_minus[i];
	        double diff=sel_event_eta_plus[i]-sel_event_eta_minus[i];
		
		Float_t CAssym_tmp=1.0*(diff/sum);
		Float_t Delta_diff = sqrt(sum);
    		Float_t Delta_sum = sqrt(sum);
    		Float_t Delta_N_plus = sqrt(sel_event_eta_plus[i]);
    		Float_t Delta_N_minus = sqrt(sel_event_eta_minus[i]);
    		Float_t Delta_Ac = CAssym_tmp * sqrt((Delta_diff/ diff) * (Delta_diff/ diff) + (Delta_sum / sum) * (Delta_sum/ sum));
    
    	        
	       if(TMath::Abs(sum)>0)
	       {
	       	CAssym.push_back(CAssym_tmp);
	       	BinError.push_back(Delta_Ac);
	       }
	       	
	       else
	       {
	           CAssym.push_back(-99999.0);
	           BinError.push_back(-99999.0);	           
	       }	
	       	
	     }
	
	}
}
