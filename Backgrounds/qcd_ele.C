#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TGraphErrors.h>
// use namespace std;
// gROOT.SetBatch(true);
void qcd_ele()
{

  gROOT->SetBatch(true);
  ////////////////////////////////// defining histograms//////////////////////////////////////////

  TH1F *h1 = new TH1F("h1", "elept", 100, 0, 400);
  TH1F *h2 = new TH1F("h2", "muonpt", 100, 0, 400);
  TH1F *h3 = new TH1F("h3", "jetspt", 100, 0, 500);
  // TH1F *h4 = new TH1F("h4","tophadpt",100,-11000,2000);
  // TH1F *h5 = new TH1F("h5","topleppt",100,-11000,2000);
  TH1F *h6 = new TH1F("h6", "tophadeta", 100, -11000, 2000);
  TH1F *h7 = new TH1F("h7", "toplepeta", 100, -11000, 2000);
  TH1F *h8 = new TH1F("h8", "toplepcharge", 100, -11000, 2000);
  TH1F *h9 = new TH1F("h9", "top", 100, -3, 3);
  TH1F *h10 = new TH1F("h10", "antitop", 100, -3, 3);
  TH1F *h11phi = new TH1F("h11", "topphi", 50, -3, 3);
  TH1F *h12phi = new TH1F("h12", "antitopphi", 50, -3, 3);
  TH1F *h16 = new TH1F("h16", "Acele", 50, 0, 1);
  TH1F *h15 = new TH1F("h15", "Mass_top", 50, 0, 500);
  TH1F *h20 = new TH1F("h20", "PhoPt", 5, -150, 800);
  TH1F *h23 = new TH1F("h23", "YPho", 5, 0, 800);
  //////////////////////////////////////////////////////////
  TH1F *h9ele = new TH1F("h9ele", "top", 50, -3, 3);
  TH1F *h10ele = new TH1F("h10ele", "antitop", 50, -3, 3);
  TH1F *h9mu = new TH1F("h9mu", "top", 50, -3, 3);
  TH1F *h10mu = new TH1F("h10mu", "antitop", 50, -3, 3);
  /////////////////////////Necessary Plots/////////////////////////////////////
  TH1F *h14 = new TH1F("h14", "Charge Asymmetry vs M_{t#bar{t}} Electronchannel ", 7, 300, 1000);
  TH1F *h17 = new TH1F("h17", "Charge Asymmetry vs M_{t#bar{t}} Muonchannel ", 7, 300, 1000);
  TH1F *h18 = new TH1F("h18", "Charge Asymmetry vs t#bar{t} Transverse Momentum Electronchannel ", 4, -200, 1000);
  TH1F *h19 = new TH1F("h19", "Charge Asymmetry vs t#bar{t} Transverse Momentum Muonchannel ", 4, -200, 1000);
  TH1F *h21 = new TH1F("h21", "Charge Asymmetry vs #gamma Transverse Momentum Electronchannel ", 5, -150, 800);
  TH1F *h22 = new TH1F("h22", "Charge Asymmetry vs #gamma Transverse Momentum Muonchannel ", 5, -150, 800);
  TH1F *h25 = new TH1F("h25", "Rapidity_ele", 4, -3, 3);
  TH1F *h26 = new TH1F("h26", "Rapidity_mu", 4, -3, 3);
  ///////////////////////////////////////////////////////////////////////////////////
  
  h14->SetLineWidth(2);
  h17->SetLineWidth(2);
  h18->SetLineWidth(2);
  h19->SetLineWidth(2);
  h21->SetLineWidth(2);
  h22->SetLineWidth(2);

  h9ele->SetOption("hist");
  h9mu->SetOption("hist");
  h10ele->SetOption("hist");
  h10mu->SetOption("hist");
  h11phi->SetOption("hist");
  h12phi->SetOption("hist");
  h16->SetOption("hist");
  h25->SetOption("hist");
  h26->SetOption("hist");
  // Set line color
  h9ele->SetLineColor(kBlack);
  h9mu->SetLineColor(kBlack);
  h10ele->SetLineColor(kBlack);
  h10mu->SetLineColor(kBlack);

  // h26->SetLineColor(kRed);
  // Set fill color
  h9ele->SetFillColor(kRed);
  h9mu->SetFillColor(kCyan);
  h10ele->SetFillColor(kBlue);
  h10mu->SetFillColor(kYellow);

  // saving histogram in root files
  TFile *fl = new TFile("qcd_ele.root", "RECREATE");
  TChain *tr = new TChain("AnalysisTree");
  // saving N+ and N- value in root
  TFile *file1 = new TFile("output_qcd_ele.root", "RECREATE");

  TTree *tree1 = new TTree("tree1", "tree1");

  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt20to30_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt30to50_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt50to80_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt80to120_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt120to170_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt170to300_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt300toInf_Ele_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt20to30_bcToE_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt30to80_bcToE_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt80to170_bcToE_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt170to250_bcToE_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/QCD_Pt250toInf_bcToE_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/GJets_HT100To200_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/GJets_HT200To400_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/GJets_HT400To600_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/GJets_HT40To100_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/GJets_HT600ToInf_2016_AnalysisNtuple.root");
  
  // TTree *tr = (TTree*)fl-> Get("AnalysisTree");
  Int_t numevents = tr->GetEntries();
  //   numevents = 100;
  //  Define the variables
  std::vector<Float_t> *elePt = 0;
  std::vector<Float_t> *muPt = 0;
  std::vector<Float_t> *jetPt = 0;
  std::vector<Float_t> *phoEt = 0;
  std::vector<Float_t> *phoEta = 0;
  std::vector<Float_t> *phoPhi = 0;

  Int_t nEle = 0;
  Int_t nMu = 0;
  Int_t nJet = 0;
  Int_t nBJet = 0;
  Int_t nPho = 0;
  Int_t lumis = 0;
  Float_t evtWeight = 0;
  Float_t btagWeight_1a = 0;
  Float_t prefireSF = 0;
  Float_t muEffWeight = 0;
  Float_t eleEffWeight = 0;
  Float_t TopHad_pt = 0;
  Float_t TopLep_pt = 0;
  Float_t TopHad_eta = 0;
  Float_t TopLep_eta = 0;
  Float_t TopLep_phi = 0;
  Float_t TopHad_phi = 0;
  Float_t Mt_blgammaMET = 0; // toplepmass
  Float_t M_bjj = 0;         // tophadmass
  Float_t TopLep_charge = 0;
  Float_t PUweight = 0;
  bool passPresel_Ele;
  bool passPresel_Mu;
  Float_t Rapidity_T_ele; // taking T for top quark
  Float_t Rapidity_t_ele; // taking t for antitop quark
  Float_t Rapidity_T_mu;
  Float_t Rapidity_t_mu;
  Float_t Rapidity_T1; // selection 1
  Float_t Rapidity_t1; // selection 1
  Float_t Rapidity_T2; // selection 2
  Float_t Rapidity_t2; // selection 2
  Float_t PhoRapidity;
  Float_t Mass_top;
  Float_t Mass_antitop;
  Float_t ttbar_mass;

  ///////////////////////////////////////////////////////////////Set Branch Status//////////////////////////////////////////////////////////////////

  tr->SetBranchStatus("*", 0);
  tr->SetBranchStatus("elePt", 1);
  tr->SetBranchStatus("muPt", 1);
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
  cout << "numevents" << numevents << std::endl;

  /////////////////////////////////////////////////////////////Set Branch Address/////////////////////////////////////////////////////////////////
  tr->SetBranchAddress("elePt", &elePt);
  tr->SetBranchAddress("muPt", &muPt);
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

  double weight = 0;

  /////////////////////////////////////////////////////Initialize//////////////////////////////////////////////////////////////////////
  float N_plus_ele1 = 0;
  float N_plus_ele2 = 0;
  float N_plus_ele3 = 0;
  float N_plus_ele4 = 0;
  float N_plus_ele5 = 0;
  float N_plus_ele6 = 0;
  float N_plus_ele7 = 0;
  float N_plus_ele8 = 0;
  float N_minus_ele1 = 0;
  float N_minus_ele2 = 0;
  float N_minus_ele3 = 0;
  float N_minus_ele4 = 0;
  float N_minus_ele5 = 0;
  float N_minus_ele6 = 0;
  float N_minus_ele7 = 0;
  float N_minus_ele8 = 0;
  float N_plus_ele = 0;
  float N_minus_ele = 0;
  float N_plus_mu1 = 0;
  float N_plus_mu2 = 0;
  float N_plus_mu3 = 0;
  float N_plus_mu4 = 0;
  float N_plus_mu5 = 0;
  float N_plus_mu6 = 0;
  float N_plus_mu7 = 0;
  float N_plus_mu8 = 0;
  float N_minus_mu1 = 0;
  float N_minus_mu2 = 0;
  float N_minus_mu3 = 0;
  float N_minus_mu4 = 0;
  float N_minus_mu5 = 0;
  float N_minus_mu6 = 0;
  float N_minus_mu7 = 0;
  float N_minus_mu8 = 0;
  float N_plus_mu = 0;
  float N_minus_mu = 0;
  float N_plus1 = 0;
  float N_minus1 = 0;
  float N_plus2 = 0;
  float N_minus2 = 0;
  ///////////////elept(tt~)//////////////////
  float N_plus_ele1pt = 0;
  float N_plus_ele2pt = 0;
  float N_plus_ele3pt = 0;
  float N_plus_ele4pt = 0;
  float N_minus_ele1pt = 0;
  float N_minus_ele2pt = 0;
  float N_minus_ele3pt = 0;
  float N_minus_ele4pt = 0;
  ///////////////mupt(tt~)///////////////////
  float N_plus_mu1pt = 0;
  float N_plus_mu2pt = 0;
  float N_plus_mu3pt = 0;
  float N_plus_mu4pt = 0;
  float N_minus_mu1pt = 0;
  float N_minus_mu2pt = 0;
  float N_minus_mu3pt = 0;
  float N_minus_mu4pt = 0;
  ////////////////elept(pho)////////////////
  float N_plus_ele1phopt = 0;
  float N_plus_ele2phopt = 0;
  float N_plus_ele3phopt = 0;
  float N_plus_ele4phopt = 0;
  float N_minus_ele1phopt = 0;
  float N_minus_ele2phopt = 0;
  float N_minus_ele3phopt = 0;
  float N_minus_ele4phopt = 0;
  ////////////////mupt(pho)///////////////////
  float N_plus_mu1phopt = 0;
  float N_plus_mu2phopt = 0;
  float N_plus_mu3phopt = 0;
  float N_plus_mu4phopt = 0;
  float N_minus_mu1phopt = 0;
  float N_minus_mu2phopt = 0;
  float N_minus_mu3phopt = 0;
  float N_minus_mu4phopt = 0;

  tree1->Branch("N_plus1", &N_plus1, "N_plus1/F");
  tree1->Branch("N_minus1", &N_minus1, "N_minus1/F");
  tree1->Branch("N_plus2", &N_plus2, "N_plus2/F");
  tree1->Branch("N_minus2", &N_minus2, "N_minus2/F");

  TLorentzVector TopLep;
  TLorentzVector TopHad;
  TLorentzVector AntiTopLep;
  TLorentzVector AntiTopHad;
  TLorentzVector Top;
  TLorentzVector AntiTop;

  // TLorentzVector TopAntiTop = Top + AntiTop;
  Float_t YT_ele;
  Float_t Yt_ele;
  Float_t rapidity_diff_ele;
  Float_t YT_mu;
  Float_t Yt_mu;
  Float_t rapidity_diff_mu;
  Float_t YT1, YT2;
  Float_t Yt1, Yt2;
  Float_t YPho;
  Float_t rapidity_diff1, rapidity_diff2;
  Float_t Mass, NetMass, NetMass1, NetMass2, NetMass3, NetMass4, NetMass5, NetMass6, NetMass7, NetMass8;
  Float_t Pt_ele, Pt_ele1, Pt_ele2, Pt_ele3, Pt_ele4;
  Float_t Pt_mu, Pt_mu1, Pt_mu2, Pt_mu3, Pt_mu4;
  Float_t PhoPt_ele, PhoPt_ele1, PhoPt_ele2, PhoPt_ele3, PhoPt_ele4;
  Float_t PhoPt_mu, PhoPt_mu1, PhoPt_mu2, PhoPt_mu3, PhoPt_mu4;
  Float_t PhoPt;
  float Ac1, Ac2;

  //////////////////////////////////////////////Set the loop//////////////////////////////////////////////////////////////////////////

  // std::cout << "Mass of top" << Mass << std::endl;
  for (Int_t i = 0; i < numevents; i++)
  {
    tr->GetEntry(i);
    /////////////////////////////////////define weight/////////////////////////////////////////////////////////////////////////
    weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;

    /////////////////////////////////////Electron Channel////////////////////////////////////////////////////////////////////

    if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1) //////////////////////////selection for electron channel////////////////////////////////////////

    {
      //////////////Tlorentz for photon///////////////////////
      TLorentzVector Photon;
      Photon.SetPtEtaPhiM(phoEt->at(0), phoEta->at(0), phoPhi->at(0), 0.0);
      PhoRapidity = Photon.Rapidity();
      YPho = TMath::Abs(PhoRapidity);
      PhoPt_ele = Photon.Pt();
      if (YPho >= 1)
      {
        if (TopLep_charge > 0)
        {
          Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
          AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
        }
        else if (TopLep_charge < 0)
        {
          AntiTop.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
          Top.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
        }
        TLorentzVector TopAntiTop = Top + AntiTop; // this is function of ttbar
        Mass = TopAntiTop.M();
        Pt_ele = TopAntiTop.Pt();
        NetMass = Mass;
        ///////////////////////////////rapidity calculation/////////////////////////////////////////////////////////////
        Rapidity_T_ele = Top.Rapidity();
        Rapidity_t_ele = AntiTop.Rapidity();
        // absolute value of rapidity
        YT_ele = TMath::Abs(Rapidity_T_ele); // top
        Yt_ele = TMath::Abs(Rapidity_t_ele); // antitop
        rapidity_diff_ele = YT_ele - Yt_ele;           // difference of rapidity; delta y
        // std::cout << "Rapidity = " << rapidity_diff_ele << std::endl;
        h25->Fill(rapidity_diff_ele, weight);

        // Fill N+ and N- with weight for overall electron channel/////////////////
        // we will get only one value of N+ and N-

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele = N_plus_ele + weight;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele = N_minus_ele + weight;
        }

        // Fill N+ and N- for different mass range of Top+antitop

        if (NetMass >= 300 && NetMass < 400)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele2 = N_plus_ele2 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele2 = N_minus_ele2 + weight;
          }
          NetMass2 = NetMass;
        }

        if (NetMass >= 400 && NetMass < 500)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele3 = N_plus_ele3 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele3 = N_minus_ele3 + weight;
          }
          NetMass3 = NetMass;
        }

        if (NetMass >= 500 && NetMass < 600)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele4 = N_plus_ele4 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele4 = N_minus_ele4 + weight;
          }
          NetMass4 = NetMass;
        }
        if (NetMass >= 600 && NetMass < 700)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele5 = N_plus_ele5 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele5 = N_minus_ele5 + weight;
          }
          NetMass5 = NetMass;
        }
        if (NetMass >= 700 && NetMass < 800)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele6 = N_plus_ele6 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele6 = N_minus_ele6 + weight;
          }
          NetMass6 = NetMass;
        }
        if (NetMass >= 800 && NetMass < 900)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele7 = N_plus_ele7 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele7 = N_minus_ele7 + weight;
          }
          NetMass7 = NetMass;
        }
        if (NetMass >= 900 && NetMass < 1000)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele8 = N_plus_ele8 + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele8 = N_minus_ele8 + weight;
          }
          NetMass8 = NetMass;
        }
        ///////////////////////////////////////////////////////////////End///////////////////////////////////////////////////////////////////////////////////////

        /////////////////////////////////////////////////////fill N+ and N- value in different range of Pt for electron channel///////////////////////////////////
        if (Pt_ele >= -200 && Pt_ele < 0)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele1pt = N_plus_ele1pt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele1pt = N_minus_ele1pt + weight;
          }
          Pt_ele1 = Pt_ele;
        }

        if (Pt_ele >= 0 && Pt_ele < 200)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele2pt = N_plus_ele2pt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele2pt = N_minus_ele2pt + weight;
          }
          Pt_ele2 = Pt_ele;
        }

        if (Pt_ele >= 200 && Pt_ele < 400)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele3pt = N_plus_ele3pt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele3pt = N_minus_ele3pt + weight;
          }
          Pt_ele3 = Pt_ele;
        }

        if (Pt_ele >= 400 && Pt_ele < 1000)
        {

          if (rapidity_diff_ele > 0)
          {
            N_plus_ele4pt = N_plus_ele4pt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele4pt = N_minus_ele4pt + weight;
          }
          Pt_ele4 = Pt_ele;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////

        /////////////////////////////////////////////////////fill N+ and N- value in different range of PhoPt for electron channel///////////////////////////////////
        if (PhoPt_ele >= -100 && PhoPt_ele < 100)
        {
          if (rapidity_diff_ele > 0)
          {
            N_plus_ele1phopt = N_plus_ele1phopt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele1phopt = N_minus_ele1phopt + weight;
          }
          PhoPt_ele1 = PhoPt_ele;
        }

        if (PhoPt_ele >= 100 && PhoPt_ele < 300)
        {
          if (rapidity_diff_ele > 0)
          {
            N_plus_ele2phopt = N_plus_ele2phopt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele2phopt = N_minus_ele2phopt + weight;
          }
          PhoPt_ele2 = PhoPt_ele;
        }

        if (PhoPt_ele >= 300 && PhoPt_ele < 500)
        {
          if (rapidity_diff_ele > 0)
          {
            N_plus_ele3phopt = N_plus_ele3phopt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele3phopt = N_minus_ele3phopt + weight;
          }
          PhoPt_ele3 = PhoPt_ele;
        }

        if (PhoPt_ele >= 500 && PhoPt_ele < 600)
        {
          if (rapidity_diff_ele > 0)
          {
            N_plus_ele4phopt = N_plus_ele4phopt + weight;
          }
          if (rapidity_diff_ele < 0)
          {
            N_minus_ele4phopt = N_minus_ele4phopt + weight;
          }
          PhoPt_ele4 = PhoPt_ele;
        }
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // h20->Fill(PhoPt_ele);
      h15->Fill(Mt_blgammaMET);
      h15->Fill(M_bjj);
      h1->Fill(elePt->at(0), weight);

      // TopLep.SetPtEtaPhiM(TopLep_pt, TopHad_eta, TopLep_phi, TopLep_mass);
      h9ele->Fill(TopLep_eta, weight); // top quark
      // AntiTopHad.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, TopHad_mass);
      h10ele->Fill(TopHad_eta, weight); // anti top quark
      h11phi->Fill(TopLep_phi, weight);
    }

    //////////////////////////////////////Muon channel//////////////////////////////////////////////////////////////////////
    if (passPresel_Mu && nJet >= 4 && nBJet >= 1 && nPho == 1) // Muon channel selection

    {
      //////////////Tlorentz for photon///////////////////////
      TLorentzVector Photon;
      Photon.SetPtEtaPhiM(phoEt->at(0), phoEta->at(0), phoPhi->at(0), 0.0);
      PhoRapidity = Photon.Rapidity();
      YPho = TMath::Abs(PhoRapidity);
      PhoPt_mu = Photon.Pt();

      if (YPho >= 1)
      {

        if (TopLep_charge > 0)
        {
          Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
          AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
        }
        else if (TopLep_charge < 0)
        {
          AntiTop.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
          Top.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
        }
        TLorentzVector TopAntiTop = Top + AntiTop;
        Mass = TopAntiTop.M();
        Pt_mu = TopAntiTop.Pt();
        NetMass = Mass;
        /////rapidity calculation////////
        Rapidity_T_mu = Top.Rapidity();
        Rapidity_t_mu = AntiTop.Rapidity();
        // absolute value of rapidity
        YT_mu = TMath::Abs(Rapidity_T_mu); // top
        Yt_mu = TMath::Abs(Rapidity_t_mu); // antitop
        rapidity_diff_mu = YT_mu - Yt_mu;
        h26->Fill(rapidity_diff_mu, weight);

        // Fill N+ and N- with weight for overall muon channel/////////////////

        if (rapidity_diff_mu > 0)
        {
          N_plus_mu = N_plus_mu + weight;
        }
        if (rapidity_diff_mu < 0)
        {
          N_minus_mu = N_minus_mu + weight;
        }

        if (NetMass >= 300 && NetMass < 400)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu2 = N_plus_mu2 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu2 = N_minus_mu2 + weight;
          }
          NetMass2 = NetMass;
        }

        if (NetMass >= 400 && NetMass < 500)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu3 = N_plus_mu3 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu3 = N_minus_mu3 + weight;
          }
          NetMass3 = NetMass;
        }

        if (NetMass >= 500 && NetMass < 600)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu4 = N_plus_mu4 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu4 = N_minus_mu4 + weight;
          }
          NetMass4 = NetMass;
        }
        if (NetMass >= 600 && NetMass < 700)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu5 = N_plus_mu5 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu5 = N_minus_mu5 + weight;
          }
          NetMass5 = NetMass;
        }
        if (NetMass >= 700 && NetMass < 800)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu6 = N_plus_mu6 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu6 = N_minus_mu6 + weight;
          }
          NetMass6 = NetMass;
        }
        if (NetMass >= 800 && NetMass < 900)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu7 = N_plus_mu7 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu7 = N_minus_mu7 + weight;
          }
          NetMass7 = NetMass;
        }
        if (NetMass >= 900 && NetMass < 1000)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu8 = N_plus_mu8 + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu8 = N_minus_mu8 + weight;
          }
          NetMass8 = NetMass;
        }

        /////////////////////////////////////////////////////fill N+ and N- value in difeerent range of Pt for muon channel///////////////////////////////////
        if (Pt_mu >= -200 && Pt_mu < 0)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu1pt = N_plus_mu1pt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu1pt = N_minus_mu1pt + weight;
          }
          Pt_mu1 = Pt_mu;
        }

        if (Pt_mu >= 0 && Pt_mu < 200)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu2pt = N_plus_mu2pt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu2pt = N_minus_mu2pt + weight;
          }
          Pt_mu2 = Pt_mu;
        }

        if (Pt_mu >= 200 && Pt_mu < 400)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu3pt = N_plus_mu3pt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu3pt = N_minus_mu3pt + weight;
          }
          Pt_mu3 = Pt_mu;
        }

        if (Pt_mu >= 400 && Pt_mu < 1000)
        {

          if (rapidity_diff_mu > 0)
          {
            N_plus_mu4pt = N_plus_mu4pt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu4pt = N_minus_mu4pt + weight;
          }
          Pt_mu4 = Pt_mu;
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////fill N+ and N- value in different range of PhoPt for muon channel/
        if (PhoPt_mu >= -100 && PhoPt_mu < 100)
        {
          if (rapidity_diff_mu > 0)
          {
            N_plus_mu1phopt = N_plus_mu1phopt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu1phopt = N_minus_mu1phopt + weight;
          }
          PhoPt_mu1 = PhoPt_mu;
        }

        if (PhoPt_mu >= 100 && PhoPt_mu < 300)
        {
          if (rapidity_diff_mu > 0)
          {
            N_plus_mu2phopt = N_plus_mu2phopt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu2phopt = N_minus_mu2phopt + weight;
          }
          PhoPt_mu2 = PhoPt_mu;
        }

        if (PhoPt_mu >= 300 && PhoPt_mu < 500)
        {
          if (rapidity_diff_mu > 0)
          {
            N_plus_mu3phopt = N_plus_mu3phopt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu3phopt = N_minus_mu3phopt + weight;
          }
          PhoPt_mu3 = PhoPt_mu;
        }

        if (PhoPt_mu >= 500 && PhoPt_mu < 600)
        {
          if (rapidity_diff_mu > 0)
          {
            N_plus_mu4phopt = N_plus_mu4phopt + weight;
          }
          if (rapidity_diff_mu < 0)
          {
            N_minus_mu4phopt = N_minus_mu4phopt + weight;
          }
          PhoPt_mu4 = PhoPt_mu;
        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      h2->Fill(muPt->at(0), weight);
      h9mu->Fill(TopHad_eta, weight);  // top quark
      h10mu->Fill(TopLep_eta, weight); // anti top quark
      h12phi->Fill(TopHad_phi, weight);
    }
    //////////////////////////////////////////Selection 1 (without photon)////////////////////////////////////
    if ((passPresel_Mu || passPresel_Ele) && nJet >= 4 && nBJet >= 1 && nPho == 0) // selection for both electron and muon channel

    {
      if (TopLep_charge > 0)
      {
        Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
        AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
      }
      else if (TopLep_charge < 0)
      {
        AntiTop.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
        Top.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
      }
      TLorentzVector TopAntiTop = Top + AntiTop;
      Mass = TopAntiTop.M();
      NetMass = Mass;
      Rapidity_T1 = Top.Rapidity();
      Rapidity_t1 = AntiTop.Rapidity();
      YT1 = TMath::Abs(Rapidity_T1); // top
      Yt1 = TMath::Abs(Rapidity_t1); // antitop
      rapidity_diff1 = YT1 - Yt1;
      if (rapidity_diff1 > 0)
      {
        N_plus1 = N_plus1 + weight; // single one value of N+ including both channel
      }
      if (rapidity_diff1 < 0)
      {
        N_minus1 = N_minus1 + weight; // single one value of N- including both channel
      }

      tree1->Fill();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////Selection 2 : with photon//////////////////////////////////////////////

    if ((passPresel_Mu || passPresel_Ele) && nJet >= 4 && nBJet >= 1 && nPho == 1)

    {
      //////////////Tlorentz for photon///////////////////////
      TLorentzVector Photon;
      Photon.SetPtEtaPhiM(phoEt->at(0), phoEta->at(0), phoPhi->at(0), 0.0);
      PhoRapidity = Photon.Rapidity();
      YPho = TMath::Abs(PhoRapidity);
      PhoPt = Photon.Pt();
      if (YPho >= 1)
      {

        // cout<<"testing 1"<<endl;
        if (TopLep_charge > 0)
        {
          Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
          AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
        }
        else if (TopLep_charge < 0)
        {
          AntiTop.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
          Top.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
        }
        // cout<<"testing 2"<<endl;
        TLorentzVector TopAntiTop = Top + AntiTop;
        Mass = TopAntiTop.M();
        NetMass = Mass;
        Rapidity_T2 = Top.Rapidity();
        Rapidity_t2 = AntiTop.Rapidity();
        YT2 = TMath::Abs(Rapidity_T2); // top
        Yt2 = TMath::Abs(Rapidity_t2); // antitop
        rapidity_diff2 = YT2 - Yt2;

        if (rapidity_diff2 > 0)
        {
          N_plus2 = N_plus2 + weight; // single one value of N+ including both channel
        }
        if (rapidity_diff2 < 0)
        {
          N_minus2 = N_minus2 + weight; // single one value of N- including both channel
        }
        //  cout<<"testing 3"<<endl;
        // h23->Fill(YPho);
        tree1->Fill();
      }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (nJet > 0)
    {
      h3->Fill(jetPt->at(0), weight);
    }

    if (TopLep_charge > 0)
    {
      h9->Fill(TopLep_eta, weight);  // top quark
      h10->Fill(TopHad_eta, weight); // anti top quark
    }

    if (TopLep_charge < 0)
    {
      h9->Fill(TopHad_eta, weight);  // top quark
      h10->Fill(TopLep_eta, weight); // antitop quark
    }

    //   h4->Fill(TopHad_pt);
    //   h5->Fill(TopLep_pt);
    //   h6->Fill(TopHad_eta);
    //   h7->Fill(TopLep_eta);
    h8->Fill(TopLep_charge, weight);
  }
  tree1->Fill();
  std::cout << "Rapidity = " << rapidity_diff_ele << std::endl;
  std::cout << "N_ele+ = " << N_plus_ele << std::endl;
  std::cout << "N_ele- = " << N_minus_ele << std::endl;
  std::cout << "N_mu+ = " << N_plus_mu << std::endl;
  std::cout << "N_mu- = " << N_minus_mu << std::endl;
  std::cout << "N+1 = " << N_plus1 << std::endl;
  std::cout << "N-1 = " << N_minus1 << std::endl;
  std::cout << "N+2 = " << N_plus2 << std::endl;
  std::cout << "N-2 = " << N_minus2 << std::endl;
  std::cout << "N_phoele1+ = " << N_plus_ele1phopt << std::endl;
  std::cout << "N_phoele2+ = " << N_plus_ele2phopt << std::endl;
  std::cout << "N_phoele3+ = " << N_plus_ele3phopt << std::endl;
  std::cout << "N_phoele4+ = " << N_plus_ele4phopt << std::endl;
  std::cout << "N_elepho1- = " << N_minus_ele1phopt << std::endl;
  std::cout << "N_elepho2- = " << N_minus_ele2phopt << std::endl;
  std::cout << "N_elepho3- = " << N_minus_ele3phopt << std::endl;
  std::cout << "N_elepho4- = " << N_minus_ele4phopt << std::endl;
  // calculate charge asymmetry/
  /////////////////////////////////////////////////// Electron channel////////////////////////////////////////////
  Float_t Ac, Ac_ele, Ac_ele1, Ac_ele2, Ac_ele3, Ac_ele4, Ac_ele5, Ac_ele6, Ac_ele7, Ac_ele8; //////////for mass
  Float_t sum_ele1 = N_plus_ele1 + N_minus_ele1;
  Float_t diff_ele1 = N_plus_ele1 - N_minus_ele1;
  Float_t sum_ele2 = N_plus_ele2 + N_minus_ele2;
  Float_t diff_ele2 = N_plus_ele2 - N_minus_ele2;
  Float_t sum_ele3 = N_plus_ele3 + N_minus_ele3;
  Float_t diff_ele3 = N_plus_ele3 - N_minus_ele3;
  Float_t sum_ele4 = N_plus_ele4 + N_minus_ele4;
  Float_t diff_ele4 = N_plus_ele4 - N_minus_ele4;
  Float_t sum_ele5 = N_plus_ele5 + N_minus_ele5;
  Float_t diff_ele5 = N_plus_ele5 - N_minus_ele5;
  Float_t sum_ele6 = N_plus_ele6 + N_minus_ele6;
  Float_t diff_ele6 = N_plus_ele6 - N_minus_ele6;
  Float_t sum_ele7 = N_plus_ele7 + N_minus_ele7;
  Float_t diff_ele7 = N_plus_ele7 - N_minus_ele7;
  Float_t sum_ele8 = N_plus_ele8 + N_minus_ele8;
  Float_t diff_ele8 = N_plus_ele8 - N_minus_ele8;
  Float_t sum_ele = N_plus_ele + N_minus_ele;
  Float_t diff_ele = N_plus_ele - N_minus_ele;
  ////////////////////////////////////////////
  Float_t Ac_ele1pt, Ac_ele2pt, Ac_ele3pt, Ac_ele4pt; //////////for pt_ele
  Float_t sum_ele1pt = N_plus_ele1pt + N_minus_ele1pt;
  Float_t diff_ele1pt = N_plus_ele1pt - N_minus_ele1pt;
  Float_t sum_ele2pt = N_plus_ele2pt + N_minus_ele2pt;
  Float_t diff_ele2pt = N_plus_ele2pt - N_minus_ele2pt;
  Float_t sum_ele3pt = N_plus_ele3pt + N_minus_ele3pt;
  Float_t diff_ele3pt = N_plus_ele3pt - N_minus_ele3pt;
  Float_t sum_ele4pt = N_plus_ele4pt + N_minus_ele4pt;
  Float_t diff_ele4pt = N_plus_ele4pt - N_minus_ele4pt;
  ////////////////////////////////////////////
  ////////////////////////////////////for PhoPt_ele//////////////////////////////////
  Float_t Ac_ele1phopt, Ac_ele2phopt, Ac_ele3phopt, Ac_ele4phopt; //////////for phopt_ele
  Float_t sum_ele1phopt = N_plus_ele1phopt + N_minus_ele1phopt;
  Float_t diff_ele1phopt = N_plus_ele1phopt - N_minus_ele1phopt;
  Float_t sum_ele2phopt = N_plus_ele2phopt + N_minus_ele2phopt;
  Float_t diff_ele2phopt = N_plus_ele2phopt - N_minus_ele2phopt;
  Float_t sum_ele3phopt = N_plus_ele3phopt + N_minus_ele3phopt;
  Float_t diff_ele3phopt = N_plus_ele3phopt - N_minus_ele3phopt;
  Float_t sum_ele4phopt = N_plus_ele4phopt + N_minus_ele4phopt;
  Float_t diff_ele4phopt = N_plus_ele4phopt - N_minus_ele4phopt;
  /////////////////////////////////////////////////

  Float_t sum = N_plus2 + N_minus2;
  Float_t diff = N_plus2 - N_minus2;

  if (sum != 0)
  {
    Ac = diff / sum;
    Float_t Delta_diff = sqrt(sum);
    Float_t Delta_sum = sqrt(sum);
    Float_t Delta_N_plus = sqrt(N_plus2);
    Float_t Delta_N_minus = sqrt(N_minus2);
    Float_t Delta_Ac = Ac * sqrt((Delta_diff / diff) * (Delta_diff / diff) + (Delta_sum / sum) * (Delta_sum / sum));
    std::cout << "Asymmetry for TTGamma = " << Ac << " ± " << Delta_Ac << std::endl;
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  //////////////////////////////electron channel//////////////////////////////////////////

  if (sum_ele != 0)
  {
    Ac_ele = diff_ele / sum_ele;
    Float_t Delta_diff_ele = sqrt(sum_ele);
    Float_t Delta_sum_ele = sqrt(sum_ele);
    Float_t Delta_N_plus1 = sqrt(N_plus_ele);
    Float_t Delta_N_minus1 = sqrt(N_minus_ele);
    Float_t Delta_Ac_ele = Ac_ele * sqrt((Delta_diff_ele / diff_ele) * (Delta_diff_ele / diff_ele) + (Delta_sum_ele / sum_ele) * (Delta_sum_ele / sum_ele));
    std::cout << "Asymmetry for electron channel = " << Ac_ele << " ± " << Delta_Ac_ele << std::endl;
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  //////////////////////////////////Fill h14/////////////////////////////////
  Float_t Delta_Ac_ele2 = 0.0;
  if (sum_ele2 != 0)
  {
    Ac_ele2 = diff_ele2 / sum_ele2;
    Float_t Delta_diff_ele2 = sqrt(sum_ele2);
    Float_t Delta_sum_ele2 = sqrt(sum_ele2);
    Float_t Delta_N_plus_ele2 = sqrt(N_plus_ele2);
    Float_t Delta_N_minus_ele2 = sqrt(N_minus_ele2);
    Float_t Delta_Ac_ele2 = Ac_ele2 * sqrt((Delta_diff_ele2 / diff_ele2) * (Delta_diff_ele2 / diff_ele2) + (Delta_sum_ele2 / sum_ele2) * (Delta_sum_ele2 / sum_ele2));
    // std::cout << "Asymmetry for electron channel = " << Ac_ele2 << std::endl;
    // std::cout << "NetMass2 =" << NetMass2 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass2), Ac_ele2);
    int bin = h14->FindBin(NetMass2);
    h14->SetBinContent(bin, Ac_ele2);
    h14->SetBinError(bin, Delta_Ac_ele2);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass2, Ac_ele2);
  Float_t Delta_Ac_ele3 = 0.0;
  if (sum_ele3 != 0)
  {
    Ac_ele3 = diff_ele3 / sum_ele3;
    Float_t Delta_diff_ele3 = sqrt(sum_ele3);
    Float_t Delta_sum_ele3 = sqrt(sum_ele3);
    Float_t Delta_N_plus_ele3 = sqrt(N_plus_ele3);
    Float_t Delta_N_minus_ele3 = sqrt(N_minus_ele3);
    Float_t Delta_Ac_ele3 = Ac_ele3 * sqrt((Delta_diff_ele3 / diff_ele3) * (Delta_diff_ele3 / diff_ele3) + (Delta_sum_ele3 / sum_ele3) * (Delta_sum_ele3 / sum_ele3));
    // std::cout << "Asymmetry for electron channel = " << Ac_ele3 << std::endl;
    // std::cout << "NetMass3 =" << NetMass3 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass3), Ac_ele3);
    int bin = h14->FindBin(NetMass3);
    h14->SetBinContent(bin, Ac_ele3);
    h14->SetBinError(bin, Delta_Ac_ele3);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass3, Ac_ele3);
  Float_t Delta_Ac_ele4 = 0.0;
  if (sum_ele4 != 0)
  {
    Ac_ele4 = diff_ele4 / sum_ele4;
    Float_t Delta_diff_ele4 = sqrt(sum_ele4);
    Float_t Delta_sum_ele4 = sqrt(sum_ele4);
    Float_t Delta_N_plus_ele4 = sqrt(N_plus_ele4);
    Float_t Delta_N_minus_ele4 = sqrt(N_minus_ele4);
    Float_t Delta_Ac_ele4 = Ac_ele4 * sqrt((Delta_diff_ele4 / diff_ele4) * (Delta_diff_ele4 / diff_ele4) + (Delta_sum_ele4 / sum_ele4) * (Delta_sum_ele4 / sum_ele4));
    // std::cout << "Asymmetry for electron channel = " << Ac_ele4 << std::endl;
    // std::cout << "NetMass4 =" << NetMass4 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass4), Ac_ele4);
    int bin = h14->FindBin(NetMass4);
    h14->SetBinContent(bin, Ac_ele4);
    h14->SetBinError(bin, Delta_Ac_ele4);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  Float_t Delta_Ac_ele5 = 0.0;
  if (sum_ele5 != 0)
  {
    Ac_ele5 = diff_ele5 / sum_ele5;
    Float_t Delta_diff_ele5 = sqrt(sum_ele5);
    Float_t Delta_sum_ele5 = sqrt(sum_ele5);
    Float_t Delta_N_plus_ele5 = sqrt(N_plus_ele5);
    Float_t Delta_N_minus_ele5 = sqrt(N_minus_ele5);
    Float_t Delta_Ac_ele5 = Ac_ele5 * sqrt((Delta_diff_ele5 / diff_ele5) * (Delta_diff_ele5 / diff_ele5) + (Delta_sum_ele5 / sum_ele5) * (Delta_sum_ele5 / sum_ele5));
    //  std::cout << "Asymmetry for electron channel = " << Ac_ele5 << std::endl;
    // std::cout << "NetMass5 =" << NetMass5 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass5), Ac_ele5);
    int bin = h14->FindBin(NetMass5);
    h14->SetBinContent(bin, Ac_ele5);
    h14->SetBinError(bin, Delta_Ac_ele5);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  // h14->Fill(NetMass1, Ac_ele1);
  Float_t Delta_Ac_ele6 = 0.0;
  if (sum_ele6 != 0)
  {
    Ac_ele6 = diff_ele6 / sum_ele6;
    Float_t Delta_diff_ele6 = sqrt(sum_ele6);
    Float_t Delta_sum_ele6 = sqrt(sum_ele6);
    Float_t Delta_N_plus_ele6 = sqrt(N_plus_ele6);
    Float_t Delta_N_minus_ele6 = sqrt(N_minus_ele6);
    Float_t Delta_Ac_ele6 = Ac_ele6 * sqrt((Delta_diff_ele6 / diff_ele6) * (Delta_diff_ele6 / diff_ele6) + (Delta_sum_ele6 / sum_ele6) * (Delta_sum_ele6 / sum_ele6));
    //  std::cout << "Asymmetry for electron channel = " << Ac_ele6 << std::endl;
    // std::cout << "NetMass6 =" << NetMass6 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass6), Ac_ele6);
    int bin = h14->FindBin(NetMass6);
    h14->SetBinContent(bin, Ac_ele6);
    h14->SetBinError(bin, Delta_Ac_ele6);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass2, Ac_ele2);
  Float_t Delta_Ac_ele7 = 0.0;
  if (sum_ele7 != 0)
  {
    Ac_ele7 = diff_ele7 / sum_ele7;
    Float_t Delta_diff_ele7 = sqrt(sum_ele7);
    Float_t Delta_sum_ele7 = sqrt(sum_ele7);
    Float_t Delta_N_plus_ele7 = sqrt(N_plus_ele7);
    Float_t Delta_N_minus_ele7 = sqrt(N_minus_ele7);
    Float_t Delta_Ac_ele7 = Ac_ele7 * sqrt((Delta_diff_ele7 / diff_ele7) * (Delta_diff_ele7 / diff_ele7) + (Delta_sum_ele7 / sum_ele7) * (Delta_sum_ele7 / sum_ele7));
    //  std::cout << "Asymmetry for electron channel = " << Ac_ele7 << std::endl;
    //  std::cout << "NetMass7 =" << NetMass7 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass7), Ac_ele7);
    int bin = h14->FindBin(NetMass7);
    h14->SetBinContent(bin, Ac_ele7);
    h14->SetBinError(bin, Delta_Ac_ele7);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass3, Ac_ele3);
  Float_t Delta_Ac_ele8 = 0.0;
  if (sum_ele8 != 0)
  {
    Ac_ele8 = diff_ele8 / sum_ele8;
    Float_t Delta_diff_ele8 = sqrt(sum_ele8);
    Float_t Delta_sum_ele8 = sqrt(sum_ele8);
    Float_t Delta_N_plus_ele8 = sqrt(N_plus_ele8);
    Float_t Delta_N_minus_ele8 = sqrt(N_minus_ele8);
    Float_t Delta_Ac_ele8 = Ac_ele8 * sqrt((Delta_diff_ele8 / diff_ele8) * (Delta_diff_ele8 / diff_ele8) + (Delta_sum_ele8 / sum_ele8) * (Delta_sum_ele8 / sum_ele8));
    //  std::cout << "Asymmetry for electron channel = " << Ac_ele8 << std::endl;
    //  std::cout << "NetMass8 =" << NetMass8 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass8), Ac_ele8);
    int bin = h14->FindBin(NetMass8);
    h14->SetBinContent(bin, Ac_ele8);
    h14->SetBinError(bin, Delta_Ac_ele8);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  //////////////////////////////////////for Ac_ele vs Pt(tt~)//////////////////////////////////////////////////
  if (sum_ele1pt != 0)
  {
    Ac_ele1pt = diff_ele1pt / sum_ele1pt;
    Float_t Delta_diff_ele1pt = sqrt(sum_ele1pt);
    Float_t Delta_sum_ele1pt = sqrt(sum_ele1pt);
    Float_t Delta_N_plus_ele1pt = sqrt(N_plus_ele1pt);
    Float_t Delta_N_minus_ele1pt = sqrt(N_minus_ele1pt);
    Float_t Delta_Ac_ele1pt = Ac_ele1pt * sqrt((Delta_diff_ele1pt / diff_ele1pt) * (Delta_diff_ele1pt / diff_ele1pt) + (Delta_sum_ele1pt / sum_ele1pt) * (Delta_sum_ele1pt / sum_ele1pt));
    // std::cout << "Asymmetry for electron channel vs pt1 = " << Ac_ele1pt << std::endl;
    // std::cout << "Pt1 =" << Pt_ele1 << std::endl;
    // h18->SetBinContent(h18->FindBin(Pt_ele1), Ac_ele1pt);
    int bin = h18->FindBin(Pt_ele1);
    h18->SetBinContent(bin, Ac_ele1pt);
    h18->SetBinError(bin, Delta_Ac_ele1pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele2pt != 0)
  {
    Ac_ele2pt = diff_ele2pt / sum_ele2pt;
    Float_t Delta_diff_ele2pt = sqrt(sum_ele2pt);
    Float_t Delta_sum_ele2pt = sqrt(sum_ele2pt);
    Float_t Delta_N_plus_ele2pt = sqrt(N_plus_ele2pt);
    Float_t Delta_N_minus_ele2pt = sqrt(N_minus_ele2pt);
    Float_t Delta_Ac_ele2pt = Ac_ele2pt * sqrt((Delta_diff_ele2pt / diff_ele2pt) * (Delta_diff_ele2pt / diff_ele2pt) + (Delta_sum_ele2pt / sum_ele2pt) * (Delta_sum_ele2pt / sum_ele2pt));
    // std::cout << "Asymmetry for electron channel vs pt2 = " << Ac_ele2pt << std::endl;
    // std::cout << "Pt2 =" << Pt_ele2 << std::endl;
    // h18->SetBinContent(h18->FindBin(Pt_ele2), Ac_ele2pt);
    int bin = h18->FindBin(Pt_ele2);
    h18->SetBinContent(bin, Ac_ele2pt);
    h18->SetBinError(bin, Delta_Ac_ele2pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele3pt != 0)
  {
    Ac_ele3pt = diff_ele3pt / sum_ele3pt;
    Float_t Delta_diff_ele3pt = sqrt(sum_ele3pt);
    Float_t Delta_sum_ele3pt = sqrt(sum_ele3pt);
    Float_t Delta_N_plus_ele3pt = sqrt(N_plus_ele3pt);
    Float_t Delta_N_minus_ele3pt = sqrt(N_minus_ele3pt);
    Float_t Delta_Ac_ele3pt = Ac_ele3pt * sqrt((Delta_diff_ele3pt / diff_ele3pt) * (Delta_diff_ele3pt / diff_ele3pt) + (Delta_sum_ele3pt / sum_ele3pt) * (Delta_sum_ele3pt / sum_ele3pt));
    // std::cout << "Asymmetry for electron channel vs pt3 = " << Ac_ele3pt << std::endl;
    // std::cout << "Pt3 =" << Pt_ele3 << std::endl;
    // h18->SetBinContent(h18->FindBin(Pt_ele3), Ac_ele3pt);
    int bin = h18->FindBin(Pt_ele3);
    h18->SetBinContent(bin, Ac_ele3pt);
    h18->SetBinError(bin, Delta_Ac_ele3pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele4pt != 0)
  {
    Ac_ele4pt = diff_ele4pt / sum_ele4pt;
    Float_t Delta_diff_ele4pt = sqrt(sum_ele4pt);
    Float_t Delta_sum_ele4pt = sqrt(sum_ele4pt);
    Float_t Delta_N_plus_ele4pt = sqrt(N_plus_ele4pt);
    Float_t Delta_N_minus_ele4pt = sqrt(N_minus_ele4pt);
    Float_t Delta_Ac_ele4pt = Ac_ele4pt * sqrt((Delta_diff_ele4pt / diff_ele4pt) * (Delta_diff_ele4pt / diff_ele4pt) + (Delta_sum_ele4pt / sum_ele4pt) * (Delta_sum_ele4pt / sum_ele4pt));
    // std::cout << "Asymmetry for electron channel vs pt4 = " << Ac_ele4pt << std::endl;
    // std::cout << "Pt4 =" << Pt_ele4 << std::endl;
    // h18->SetBinContent(h18->FindBin(Pt_ele4), Ac_ele4pt);
    int bin = h18->FindBin(Pt_ele4);
    h18->SetBinContent(bin, Ac_ele4pt);
    h18->SetBinError(bin, Delta_Ac_ele4pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  ///////////////////////////////////////for Ac_ele vs PhoPt_ele///////////////////////////////////

  if (sum_ele1phopt != 0)
  {
    Ac_ele1phopt = diff_ele1phopt / sum_ele1phopt;
    Float_t Delta_diff_ele1phopt = sqrt(sum_ele1phopt);
    Float_t Delta_sum_ele1phopt = sqrt(sum_ele1phopt);
    Float_t Delta_N_plus_ele1phopt = sqrt(N_plus_ele1phopt);
    Float_t Delta_N_minus_ele1phopt = sqrt(N_minus_ele1phopt);
    Float_t Delta_Ac_ele1phopt = Ac_ele1phopt * sqrt((Delta_diff_ele1phopt / diff_ele1phopt) * (Delta_diff_ele1phopt / diff_ele1phopt) + (Delta_sum_ele1phopt / sum_ele1phopt) * (Delta_sum_ele1phopt / sum_ele1phopt));
    // std::cout << "Asymmetry for electron channel vs phopt1 = " << Ac_ele1phopt << std::endl;
    // std::cout << "phoPtele1 =" << PhoPt_ele1 << std::endl;
    // h21->SetBinContent(h21->FindBin(PhoPt_ele1), Ac_ele1phopt);
    int bin = h21->FindBin(PhoPt_ele1);
    h21->SetBinContent(bin, Ac_ele1phopt);
    h21->SetBinError(bin, Delta_Ac_ele1phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele2phopt != 0)
  {
    Ac_ele2phopt = diff_ele2phopt / sum_ele2phopt;
    Float_t Delta_diff_ele2phopt = sqrt(sum_ele2phopt);
    Float_t Delta_sum_ele2phopt = sqrt(sum_ele2phopt);
    Float_t Delta_N_plus_ele2phopt = sqrt(N_plus_ele2phopt);
    Float_t Delta_N_minus_ele2phopt = sqrt(N_minus_ele2phopt);
    Float_t Delta_Ac_ele2phopt = Ac_ele2phopt * sqrt((Delta_diff_ele2phopt / diff_ele2phopt) * (Delta_diff_ele2phopt / diff_ele2phopt) + (Delta_sum_ele2phopt / sum_ele2phopt) * (Delta_sum_ele2phopt / sum_ele2phopt));
    // std::cout << "Asymmetry for electron channel vs phopt2 = " << Ac_ele2phopt << std::endl;
    //  std::cout << "phoPtele2 =" << PhoPt_ele2 << std::endl;
    // h21->SetBinContent(h21->FindBin(PhoPt_ele2), Ac_ele2phopt);
    int bin = h21->FindBin(PhoPt_ele2);
    h21->SetBinContent(bin, Ac_ele2phopt);
    h21->SetBinError(bin, Delta_Ac_ele2phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele3phopt != 0)
  {
    Ac_ele3phopt = diff_ele3phopt / sum_ele3phopt;
    Float_t Delta_diff_ele3phopt = sqrt(sum_ele3phopt);
    Float_t Delta_sum_ele3phopt = sqrt(sum_ele3phopt);
    Float_t Delta_N_plus_ele3phopt = sqrt(N_plus_ele3phopt);
    Float_t Delta_N_minus_ele3phopt = sqrt(N_minus_ele3phopt);
    Float_t Delta_Ac_ele3phopt = Ac_ele3phopt * sqrt((Delta_diff_ele3phopt / diff_ele3phopt) * (Delta_diff_ele3phopt / diff_ele3phopt) + (Delta_sum_ele3phopt / sum_ele3phopt) * (Delta_sum_ele3phopt / sum_ele3phopt));
    // std::cout << "Asymmetry for electron channel vs phopt1 = " << Ac_ele3phopt << std::endl;
    // std::cout << "phoPtele3 =" << PhoPt_ele3 << std::endl;
    // h21->SetBinContent(h21->FindBin(PhoPt_ele3), Ac_ele3phopt);
    int bin = h21->FindBin(PhoPt_ele3);
    h21->SetBinContent(bin, Ac_ele3phopt);
    h21->SetBinError(bin, Delta_Ac_ele3phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele4phopt != 0)
  {
    Ac_ele4phopt = diff_ele4phopt / sum_ele4phopt;
    Float_t Delta_diff_ele4phopt = sqrt(sum_ele4phopt);
    Float_t Delta_sum_ele4phopt = sqrt(sum_ele4phopt);
    Float_t Delta_N_plus_ele4phopt = sqrt(N_plus_ele4phopt);
    Float_t Delta_N_minus_ele4phopt = sqrt(N_minus_ele4phopt);
    Float_t Delta_Ac_ele4phopt = Ac_ele4phopt * sqrt((Delta_diff_ele4phopt / diff_ele4phopt) * (Delta_diff_ele4phopt / diff_ele4phopt) + (Delta_sum_ele4phopt / sum_ele4phopt) * (Delta_sum_ele4phopt / sum_ele4phopt));
    // std::cout << "Asymmetry for electron channel vs phopt4 = " << Ac_ele4phopt << std::endl;
    // std::cout << "phoPtele4 =" << PhoPt_ele4 << std::endl;
    // h21->SetBinContent(h21->FindBin(PhoPt_ele4), Ac_ele4phopt);
    int bin = h21->FindBin(PhoPt_ele4);
    h21->SetBinContent(bin, Ac_ele4phopt);
    h21->SetBinError(bin, Delta_Ac_ele4phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////

  /* for (int i = 1; i <= h14->GetNbinsX(); i++)
   {
     std::cout << "Bin " << i << " content: " << h14->GetBinContent(i) << std::endl;
   }*/

  /////////////////////////////////////////Muon channel////////////////////////////////////////////////////////////
  Float_t Ac_mu1, Ac_mu2, Ac_mu3, Ac_mu4, Ac_mu5, Ac_mu6, Ac_mu7, Ac_mu8;
  Float_t sum_mu1 = N_plus_mu1 + N_minus_mu1;
  Float_t diff_mu1 = N_plus_mu1 - N_minus_mu1;
  Float_t sum_mu2 = N_plus_mu2 + N_minus_mu2;
  Float_t diff_mu2 = N_plus_mu2 - N_minus_mu2;
  Float_t sum_mu3 = N_plus_mu3 + N_minus_mu3;
  Float_t diff_mu3 = N_plus_mu3 - N_minus_mu3;
  Float_t sum_mu4 = N_plus_mu4 + N_minus_mu4;
  Float_t diff_mu4 = N_plus_mu4 - N_minus_mu4;
  Float_t sum_mu5 = N_plus_mu5 + N_minus_mu5;
  Float_t diff_mu5 = N_plus_mu5 - N_minus_mu5;
  Float_t sum_mu6 = N_plus_mu6 + N_minus_mu6;
  Float_t diff_mu6 = N_plus_mu6 - N_minus_mu6;
  Float_t sum_mu7 = N_plus_mu7 + N_minus_mu7;
  Float_t diff_mu7 = N_plus_mu7 - N_minus_mu7;
  Float_t sum_mu8 = N_plus_mu8 + N_minus_mu8;
  Float_t diff_mu8 = N_plus_mu8 - N_minus_mu8;
  Float_t sum_mu = N_plus_mu + N_minus_mu;
  Float_t diff_mu = N_plus_mu - N_minus_mu;
  ///////////////////////////////////////////////////////////////
  Float_t Ac_mu1pt, Ac_mu2pt, Ac_mu3pt, Ac_mu4pt; //////////for pt_mu
  Float_t sum_mu1pt = N_plus_mu1pt + N_minus_mu1pt;
  Float_t diff_mu1pt = N_plus_mu1pt - N_minus_mu1pt;
  Float_t sum_mu2pt = N_plus_mu2pt + N_minus_mu2pt;
  Float_t diff_mu2pt = N_plus_mu2pt - N_minus_mu2pt;
  Float_t sum_mu3pt = N_plus_mu3pt + N_minus_mu3pt;
  Float_t diff_mu3pt = N_plus_mu3pt - N_minus_mu3pt;
  Float_t sum_mu4pt = N_plus_mu4pt + N_minus_mu4pt;
  Float_t diff_mu4pt = N_plus_mu4pt - N_minus_mu4pt;
  ///////////////////////////////////////////////////
  ////////////////////////////////////for PhoPt_mu//////////////////////////////////
  Float_t Ac_mu1phopt, Ac_mu2phopt, Ac_mu3phopt, Ac_mu4phopt; //////////for phopt_mu
  Float_t sum_mu1phopt = N_plus_mu1phopt + N_minus_mu1phopt;
  Float_t diff_mu1phopt = N_plus_mu1phopt - N_minus_mu1phopt;
  Float_t sum_mu2phopt = N_plus_mu2phopt + N_minus_mu2phopt;
  Float_t diff_mu2phopt = N_plus_mu2phopt - N_minus_mu2phopt;
  Float_t sum_mu3phopt = N_plus_mu3phopt + N_minus_mu3phopt;
  Float_t diff_mu3phopt = N_plus_mu3phopt - N_minus_mu3phopt;
  Float_t sum_mu4phopt = N_plus_mu4phopt + N_minus_mu4phopt;
  Float_t diff_mu4phopt = N_plus_mu4phopt - N_minus_mu4phopt;
  /////////////////////////////////////////////////

  if (sum_mu != 0)
  {
    Float_t Ac_mu = diff_mu / sum_mu;
    Float_t Delta_diff_mu = sqrt(sum_mu);
    Float_t Delta_sum_mu = sqrt(sum_mu);
    Float_t Delta_N_plus = sqrt(N_plus2);
    Float_t Delta_N_minus = sqrt(N_minus2);
    Float_t Delta_Ac_mu = Ac_mu * sqrt((Delta_diff_mu / diff_mu) * (Delta_diff_mu / diff_mu) + (Delta_sum_mu / sum_mu) * (Delta_sum_mu / sum_mu));
    std::cout << "Asymmetry for muon channel = " << Ac_mu << " ± " << Delta_Ac_mu << std::endl;
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  /////////////////////////////////Fill h17/////////////////////////////////////////////////////
  if (sum_mu2 != 0)
  {
    Ac_mu2 = diff_mu2 / sum_mu2;
    Float_t Delta_diff_mu2 = sqrt(sum_mu2);
    Float_t Delta_sum_mu2 = sqrt(sum_mu2);
    Float_t Delta_N_plus_mu2 = sqrt(N_plus_mu2);
    Float_t Delta_N_minus_mu2 = sqrt(N_minus_mu2);
    Float_t Delta_Ac_mu2 = Ac_mu2 * sqrt((Delta_diff_mu2 / diff_mu2) * (Delta_diff_mu2 / diff_mu2) + (Delta_sum_mu2 / sum_mu2) * (Delta_sum_mu2 / sum_mu2));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu2 << std::endl;
    // std::cout << "NetMass2 =" << NetMass2 << std::endl;
    // h17->SetBinContent(h17->FindBin(NetMass2), Ac_mu2);
    int bin = h17->FindBin(NetMass2);
    h17->SetBinContent(bin, Ac_mu2);
    h17->SetBinError(bin, Delta_Ac_mu2);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h17->Fill(NetMass2, Ac_mu2);

  if (sum_mu3 != 0)
  {
    Ac_mu3 = diff_mu3 / sum_mu3;
    Float_t Delta_diff_mu3 = sqrt(sum_mu3);
    Float_t Delta_sum_mu3 = sqrt(sum_mu3);
    Float_t Delta_N_plus_mu3 = sqrt(N_plus_mu3);
    Float_t Delta_N_minus_mu3 = sqrt(N_minus_mu3);
    Float_t Delta_Ac_mu3 = Ac_mu3 * sqrt((Delta_diff_mu3 / diff_mu3) * (Delta_diff_mu3 / diff_mu3) + (Delta_sum_mu3 / sum_mu3) * (Delta_sum_mu3 / sum_mu3));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu3 << std::endl;
    // std::cout << "NetMass3 =" << NetMass3 << std::endl;
    // h17->SetBinContent(h17->FindBin(NetMass3), Ac_mu3);
    int bin = h17->FindBin(NetMass3);
    h17->SetBinContent(bin, Ac_mu3);
    h17->SetBinError(bin, Delta_Ac_mu3);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h17->Fill(NetMass3, Ac_mu3);

  if (sum_mu4 != 0)
  {
    Ac_mu4 = diff_mu4 / sum_mu4;
    Float_t Delta_diff_mu4 = sqrt(sum_mu4);
    Float_t Delta_sum_mu4 = sqrt(sum_mu4);
    Float_t Delta_N_plus_mu4 = sqrt(N_plus_mu4);
    Float_t Delta_N_minus_mu4 = sqrt(N_minus_mu4);
    Float_t Delta_Ac_mu4 = Ac_mu4 * sqrt((Delta_diff_mu4 / diff_mu4) * (Delta_diff_mu4 / diff_mu4) + (Delta_sum_mu4 / sum_mu4) * (Delta_sum_mu4 / sum_mu4));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu4 << std::endl;
    // std::cout << "NetMass4 =" << NetMass4 << std::endl;
    // h17->SetBinContent(h17->FindBin(NetMass4), Ac_mu4);
    int bin = h17->FindBin(NetMass4);
    h17->SetBinContent(bin, Ac_mu4);
    h17->SetBinError(bin, Delta_Ac_mu4);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_mu5 != 0)
  {
    Ac_mu5 = diff_mu5 / sum_mu5;
    Float_t Delta_diff_mu5 = sqrt(sum_mu5);
    Float_t Delta_sum_mu5 = sqrt(sum_mu5);
    Float_t Delta_N_plus_mu5 = sqrt(N_plus_mu5);
    Float_t Delta_N_minus_mu5 = sqrt(N_minus_mu5);
    Float_t Delta_Ac_mu5 = Ac_mu5 * sqrt((Delta_diff_mu5 / diff_mu5) * (Delta_diff_mu5 / diff_mu5) + (Delta_sum_mu5 / sum_mu5) * (Delta_sum_mu5 / sum_mu5));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu5 << std::endl;
    // std::cout << "NetMass5 =" << NetMass5 << std::endl;
    // h17->SetBinContent(h17->FindBin(NetMass5), Ac_mu5);
    int bin = h17->FindBin(NetMass5);
    h17->SetBinContent(bin, Ac_mu5);
    h17->SetBinError(bin, Delta_Ac_mu5);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  // h17->Fill(NetMass1, Ac_mu1);

  if (sum_mu6 != 0)
  {
    Ac_mu6 = diff_mu6 / sum_mu6;
    Float_t Delta_diff_mu6 = sqrt(sum_mu6);
    Float_t Delta_sum_mu6 = sqrt(sum_mu6);
    Float_t Delta_N_plus_mu6 = sqrt(N_plus_mu6);
    Float_t Delta_N_minus_mu6 = sqrt(N_minus_mu6);
    Float_t Delta_Ac_mu6 = Ac_mu6 * sqrt((Delta_diff_mu6 / diff_mu6) * (Delta_diff_mu6 / diff_mu6) + (Delta_sum_mu6 / sum_mu6) * (Delta_sum_mu6 / sum_mu6));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu6 << std::endl;
    // std::cout << "NetMass6 =" << NetMass6 << std::endl;
    // h17->SetBinContent(h17->FindBin(NetMass6), Ac_mu6);
    int bin = h17->FindBin(NetMass6);
    h17->SetBinContent(bin, Ac_mu6);
    h17->SetBinError(bin, Delta_Ac_mu6);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h17->Fill(NetMass2, Ac_mu2);

  if (sum_mu7 != 0)
  {
    Ac_mu7 = diff_mu7 / sum_mu7;
    Float_t Delta_diff_mu7 = sqrt(sum_mu7);
    Float_t Delta_sum_mu7 = sqrt(sum_mu7);
    Float_t Delta_N_plus_mu7 = sqrt(N_plus_mu7);
    Float_t Delta_N_minus_mu7 = sqrt(N_minus_mu7);
    Float_t Delta_Ac_mu7 = Ac_mu7 * sqrt((Delta_diff_mu7 / diff_mu7) * (Delta_diff_mu7 / diff_mu7) + (Delta_sum_mu7 / sum_mu7) * (Delta_sum_mu7 / sum_mu7));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu7 << std::endl;
    // std::cout << "NetMass7 =" << NetMass7 << std::endl;
    // h17->SetBinContent(h17->FindBin(NetMass7), Ac_mu7);
    int bin = h17->FindBin(NetMass7);
    h17->SetBinContent(bin, Ac_mu7);
    h17->SetBinError(bin, Delta_Ac_mu7);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h17->Fill(NetMass3, Ac_mu3);

  if (sum_mu8 != 0)
  {
    Ac_mu8 = diff_mu8 / sum_mu8;
    Float_t Delta_diff_mu8 = sqrt(sum_mu8);
    Float_t Delta_sum_mu8 = sqrt(sum_mu8);
    Float_t Delta_N_plus_mu8 = sqrt(N_plus_mu8);
    Float_t Delta_N_minus_mu8 = sqrt(N_minus_mu8);
    Float_t Delta_Ac_mu8 = Ac_mu8 * sqrt((Delta_diff_mu8 / diff_mu8) * (Delta_diff_mu8 / diff_mu8) + (Delta_sum_mu8 / sum_mu8) * (Delta_sum_mu8 / sum_mu8));
    // std::cout << "Asymmetry for muon channel = " << Ac_mu8 << std::endl;
    // std::cout << "NetMass8 =" << NetMass8 << std::endl;
    // h17->SetBinContent(h18->FindBin(NetMass8), Ac_mu8);
    int bin = h17->FindBin(NetMass8);
    h17->SetBinContent(bin, Ac_mu8);
    h17->SetBinError(bin, Delta_Ac_mu8);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  ////////////////////////////////////asymmetry vs pt(tt~) for muon channel////////////////////////////////////////////////////
  if (sum_mu1pt != 0)
  {
    Ac_mu1pt = diff_mu1pt / sum_mu1pt;
    Float_t Delta_diff_mu1pt = sqrt(sum_mu1pt);
    Float_t Delta_sum_mu1pt = sqrt(sum_mu1pt);
    Float_t Delta_N_plus_mu1pt = sqrt(N_plus_mu1pt);
    Float_t Delta_N_minus_mu1pt = sqrt(N_minus_mu1pt);
    Float_t Delta_Ac_mu1pt = Ac_mu1pt * sqrt((Delta_diff_mu1pt / diff_mu1pt) * (Delta_diff_mu1pt / diff_mu1pt) + (Delta_sum_mu1pt / sum_mu1pt) * (Delta_sum_mu1pt / sum_mu1pt));
    // std::cout << "Asymmetry for muon channel vs pt1 = " << Ac_mu1pt << std::endl;
    // std::cout << "Pt1 =" << Pt_mu1 << std::endl;
    // h19->SetBinContent(h19->FindBin(Pt_mu1), Ac_mu1pt);
    int bin = h19->FindBin(Pt_mu1);
    h19->SetBinContent(bin, Ac_mu1pt);
    h19->SetBinError(bin, Delta_Ac_mu1pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry for pt1_muon: sum is zero" << std::endl;
  }

  if (sum_mu2pt != 0)
  {
    Ac_mu2pt = diff_mu2pt / sum_mu2pt;
    Float_t Delta_diff_mu2pt = sqrt(sum_mu2pt);
    Float_t Delta_sum_mu2pt = sqrt(sum_mu2pt);
    Float_t Delta_N_plus_mu2pt = sqrt(N_plus_mu2pt);
    Float_t Delta_N_minus_mu2pt = sqrt(N_minus_mu2pt);
    Float_t Delta_Ac_mu2pt = Ac_mu2pt * sqrt((Delta_diff_mu2pt / diff_mu2pt) * (Delta_diff_mu2pt / diff_mu2pt) + (Delta_sum_mu2pt / sum_mu2pt) * (Delta_sum_mu2pt / sum_mu2pt));
    // std::cout << "Asymmetry for muon channel vs pt2 = " << Ac_mu2pt << std::endl;
    // std::cout << "Pt2 =" << Pt_mu2 << std::endl;
    // h19->SetBinContent(h19->FindBin(Pt_mu2), Ac_mu2pt);
    int bin = h19->FindBin(Pt_mu2);
    h19->SetBinContent(bin, Ac_mu2pt);
    h19->SetBinError(bin, Delta_Ac_mu2pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry for pt2_muon: sum is zero" << std::endl;
  }

  if (sum_mu3pt != 0)
  {
    Ac_mu3pt = diff_mu3pt / sum_mu3pt;
    Float_t Delta_diff_mu3pt = sqrt(sum_mu3pt);
    Float_t Delta_sum_mu3pt = sqrt(sum_mu3pt);
    Float_t Delta_N_plus_mu3pt = sqrt(N_plus_mu3pt);
    Float_t Delta_N_minus_mu3pt = sqrt(N_minus_mu3pt);
    Float_t Delta_Ac_mu3pt = Ac_mu3pt * sqrt((Delta_diff_mu3pt / diff_mu3pt) * (Delta_diff_mu3pt / diff_mu3pt) + (Delta_sum_mu3pt / sum_mu3pt) * (Delta_sum_mu3pt / sum_mu3pt));
    // std::cout << "Asymmetry for muon channel vs pt3 = " << Ac_mu3pt << std::endl;
    // std::cout << "Pt3 =" << Pt_mu3 << std::endl;
    // h19->SetBinContent(h19->FindBin(Pt_mu3), Ac_mu3pt);
    int bin = h19->FindBin(Pt_mu3);
    h19->SetBinContent(bin, Ac_mu3pt);
    h19->SetBinError(bin, Delta_Ac_mu3pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry for pt3_muon: sum is zero" << std::endl;
  }

  if (sum_mu4pt != 0)
  {
    Ac_mu4pt = diff_mu4pt / sum_mu4pt;
    Float_t Delta_diff_mu4pt = sqrt(sum_mu4pt);
    Float_t Delta_sum_mu4pt = sqrt(sum_mu4pt);
    Float_t Delta_N_plus_mu4pt = sqrt(N_plus_mu4pt);
    Float_t Delta_N_minus_mu4pt = sqrt(N_minus_mu4pt);
    Float_t Delta_Ac_mu4pt = Ac_mu4pt * sqrt((Delta_diff_mu4pt / diff_mu4pt) * (Delta_diff_mu4pt / diff_mu4pt) + (Delta_sum_mu4pt / sum_mu4pt) * (Delta_sum_mu4pt / sum_mu4pt));
    // std::cout << "Asymmetry for muon channel vs pt4 = " << Ac_mu4pt << std::endl;
    // std::cout << "Pt4 =" << Pt_mu4 << std::endl;
    // h19->SetBinContent(h19->FindBin(Pt_mu4), Ac_mu4pt);
    int bin = h19->FindBin(Pt_mu4);
    h19->SetBinContent(bin, Ac_mu4pt);
    h19->SetBinError(bin, Delta_Ac_mu4pt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry for pt4_muon: sum is zero" << std::endl;
  }

  ///////////////////////////////////////for Ac vs PhoPt_mu///////////////////////////////////

  if (sum_mu1phopt != 0)
  {
    Ac_mu1phopt = diff_mu1phopt / sum_mu1phopt;
    Float_t Delta_diff_mu1phopt = sqrt(sum_mu1phopt);
    Float_t Delta_sum_mu1phopt = sqrt(sum_mu1phopt);
    Float_t Delta_N_plus_mu1phopt = sqrt(N_plus_mu1phopt);
    Float_t Delta_N_minus_mu1phopt = sqrt(N_minus_mu1phopt);
    Float_t Delta_Ac_mu1phopt = Ac_mu1phopt * sqrt((Delta_diff_mu1phopt / diff_mu1phopt) * (Delta_diff_mu1phopt / diff_mu1phopt) + (Delta_sum_mu1phopt / sum_mu1phopt) * (Delta_sum_mu1phopt / sum_mu1phopt));
    // std::cout << "Asymmetry for muon channel vs phopt1 = " << Ac_mu1phopt << std::endl;
    // std::cout << "phoPtmu2 =" << PhoPt_mu2 << std::endl;
    // h22->SetBinContent(h22->FindBin(PhoPt_mu2), Ac_mu2phopt);
    int bin = h22->FindBin(PhoPt_mu1);
    h22->SetBinContent(bin, Ac_mu1phopt);
    h22->SetBinError(bin, Delta_Ac_mu1phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_mu2phopt != 0)
  {
    Ac_mu2phopt = diff_mu2phopt / sum_mu2phopt;
    Float_t Delta_diff_mu2phopt = sqrt(sum_mu2phopt);
    Float_t Delta_sum_mu2phopt = sqrt(sum_mu2phopt);
    Float_t Delta_N_plus_mu2phopt = sqrt(N_plus_mu2phopt);
    Float_t Delta_N_minus_mu2phopt = sqrt(N_minus_mu2phopt);
    Float_t Delta_Ac_mu2phopt = Ac_mu2phopt * sqrt((Delta_diff_mu2phopt / diff_mu2phopt) * (Delta_diff_mu2phopt / diff_mu2phopt) + (Delta_sum_mu2phopt / sum_mu2phopt) * (Delta_sum_mu2phopt / sum_mu2phopt));
    // std::cout << "Asymmetry for muon channel vs phopt1 = " << Ac_mu2phopt << std::endl;
    // std::cout << "phoPtmu2 =" << PhoPt_mu2 << std::endl;
    // h22->SetBinContent(h22->FindBin(PhoPt_mu2), Ac_mu2phopt);
    int bin = h22->FindBin(PhoPt_mu2);
    h22->SetBinContent(bin, Ac_mu2phopt);
    h22->SetBinError(bin, Delta_Ac_mu2phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_mu3phopt != 0)
  {
    Ac_mu3phopt = diff_mu3phopt / sum_mu3phopt;
    Float_t Delta_diff_mu3phopt = sqrt(sum_mu3phopt);
    Float_t Delta_sum_mu3phopt = sqrt(sum_mu3phopt);
    Float_t Delta_N_plus_mu3phopt = sqrt(N_plus_mu3phopt);
    Float_t Delta_N_minus_mu3phopt = sqrt(N_minus_mu3phopt);
    Float_t Delta_Ac_mu3phopt = Ac_mu3phopt * sqrt((Delta_diff_mu3phopt / diff_mu3phopt) * (Delta_diff_mu3phopt / diff_mu3phopt) + (Delta_sum_mu3phopt / sum_mu3phopt) * (Delta_sum_mu3phopt / sum_mu3phopt));
    // std::cout << "Asymmetry for muon channel vs phopt1 = " << Ac_mu3phopt << std::endl;
    // std::cout << "phoPtmu3 =" << PhoPt_mu3 << std::endl;
    // h22->SetBinContent(h22->FindBin(PhoPt_mu3), Ac_mu3phopt);
    int bin = h22->FindBin(PhoPt_mu3);
    h22->SetBinContent(bin, Ac_mu3phopt);
    h22->SetBinError(bin, Delta_Ac_mu3phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_mu4phopt != 0)
  {
    Ac_mu4phopt = diff_mu4phopt / sum_mu4phopt;
    Float_t Delta_diff_mu4phopt = sqrt(sum_mu4phopt);
    Float_t Delta_sum_mu4phopt = sqrt(sum_mu4phopt);
    Float_t Delta_N_plus_mu4phopt = sqrt(N_plus_mu4phopt);
    Float_t Delta_N_minus_mu4phopt = sqrt(N_minus_mu4phopt);
    Float_t Delta_Ac_mu4phopt = Ac_mu4phopt * sqrt((Delta_diff_mu4phopt / diff_mu4phopt) * (Delta_diff_mu4phopt / diff_mu4phopt) + (Delta_sum_mu4phopt / sum_mu4phopt) * (Delta_sum_mu4phopt / sum_mu4phopt));
    // std::cout << "Asymmetry for muon channel vs phopt1 = " << Ac_mu4phopt << std::endl;
    // std::cout << "phoPtmu4 =" << PhoPt_mu4 << std::endl;
    // h22->SetBinContent(h22->FindBin(PhoPt_mu4), Ac_mu4phopt);
    int bin = h22->FindBin(PhoPt_mu4);
    h22->SetBinContent(bin, Ac_mu4phopt);
    h22->SetBinError(bin, Delta_Ac_mu4phopt);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  //////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////

  // h9->Write();
  // h10->Write();

  h9ele->Write();
  h10ele->Write();
  h9mu->Write();
  h10mu->Write();
  h14->Write();
  h17->Write();
  h18->Write();
  h19->Write();
  h21->Write();
  h22->Write();
  h25->Write();
  h26->Write();

  tree1->Write();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1->Draw();
  c1->Update();
  c1->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h1.png");
  TCanvas *c2 = new TCanvas();
  c2->cd();
  h2->Draw();
  c2->Update();
  c2->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h2.png");
  TCanvas *c3 = new TCanvas();
  c3->cd();
  h3->Draw();
  c3->Update();
  c3->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h3.png");
  /* TCanvas *c4 = new TCanvas();
 c4->cd();
 h4->Draw();
 c4->Update();
 c4->SaveAs("h4.png");
 TCanvas *c5 = new TCanvas();
 c5->cd();
 h5->Draw();
 c5->Update();
 c5->SaveAs("h5.png");
 TCanvas *c6 = new TCanvas();
 c6->cd();
 h6->Draw();
 c6->Update();
 c6->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h6.png");
 TCanvas *c7 = new TCanvas();
 c7->cd();
 h7->Draw();
 c7->Update();
 c7->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h7.png");*/
  TCanvas *c8 = new TCanvas();
  c8->cd();
  h8->Draw();
  c8->Update();
  c8->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h8.png");
  TCanvas *c9 = new TCanvas();
  c9->cd();
  h9->Draw();
  c9->Update();
  c9->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h9.png");
  TCanvas *c10 = new TCanvas();
  c10->cd();
  h10->Draw();
  c10->Update();
  c10->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h10.png");

  TCanvas *c11 = new TCanvas();
  c11->cd();
  h9ele->Draw();
  c11->Update();
  c11->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h9electrontop.png");

  TCanvas *c12 = new TCanvas();
  c12->cd();
  h10ele->Draw();
  c12->Update();
  c12->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h10electronantitop.png");

  TCanvas *c13 = new TCanvas();
  c13->cd();
  h9mu->Draw();
  c13->Update();
  c13->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h9mutop.png");

  TCanvas *c14 = new TCanvas();
  c14->cd();
  h10mu->Draw();
  c14->Update();
  c14->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h10muantitop.png");

  TCanvas *c15 = new TCanvas();
  c15->cd();
  h11phi->Draw();
  c15->Update();
  c15->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h11phi.png");

  TCanvas *c16 = new TCanvas();
  c16->cd();
  h12phi->Draw();
  c16->Update();
  c16->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/h12phi.png");

  TCanvas *c17 = new TCanvas();
  c17->cd();
  h14->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
  h14->GetXaxis()->SetTitleSize(0.04);
  h14->GetXaxis()->SetLabelSize(0.03);
  h14->GetYaxis()->SetTitle("A_{c_{ele}}");
  h14->GetYaxis()->SetTitleOffset(1.1);
  h14->GetYaxis()->SetTitleSize(0.04);
  h14->GetYaxis()->SetLabelSize(0.03);

  h14->SetMarkerSize(1);
  h14->SetMarkerColor(kRed);
  h14->SetLineWidth(1);
  gStyle->SetEndErrorSize(1);
  gStyle->SetErrorX(.3);

  // h14->SetOption("E3");

  h14->Draw("E1");
  // create a legend object and set its properties
  TLegend *leg1 = new TLegend(0.6, 0.7, 0.8, 0.6); // position in NDC coordinates
  leg1->SetBorderSize(0);                          // remove border
  leg1->SetFillColor(0);                           // transparent background

  // add legend entries with the selection criteria
  leg1->AddEntry("", "Selection Criteria:", ""); // empty string for plot object
  leg1->AddEntry("", "nJets #geq 4", "");
  leg1->AddEntry("", "nBjet #geq 1", "");
  leg1->AddEntry("", "nPho = 1", "");
  leg1->Draw();
  c17->Update();
  c17->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Ac_ele vs ttbarmass.png");

  TCanvas *c18 = new TCanvas();
  c18->cd();
  h15->Draw();
  c18->Update();
  c18->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/masstop.png");

  TCanvas *c19 = new TCanvas();
  c19->cd();
  h16->Draw();
  c19->Update();
  c19->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Acele.png");

  TCanvas *c20 = new TCanvas();
  c20->cd();
  h17->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
  h17->GetXaxis()->SetTitleSize(0.04);
  h17->GetXaxis()->SetLabelSize(0.03);
  h17->GetYaxis()->SetTitle("A_{c_{mu}}");
  h17->GetYaxis()->SetTitleSize(0.04);
  h17->GetYaxis()->SetTitleOffset(1.1);
  h17->GetYaxis()->SetLabelSize(0.03);
  h17->SetMarkerSize(1);
  h17->SetMarkerColor(kRed);
  h17->SetLineWidth(1);
  gStyle->SetEndErrorSize(1);
  gStyle->SetErrorX(.3);
  // h17->Draw("E1");
  h17->SetOption("E3");
  h17->Draw();
  // create a legend object and set its properties
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.8, 0.6); // position in NDC coordinates
  leg2->SetBorderSize(0);                          // remove border
  leg2->SetFillColor(0);                           // transparent background

  // add legend entries with the selection criteria
  leg2->AddEntry("", "Selection Criteria:", ""); // empty string for plot object
  leg2->AddEntry("", "nJets #geq 4", "");
  leg2->AddEntry("", "nBjet #geq 1", "");
  leg2->AddEntry("", "nPho = 1", "");
  leg2->Draw();
  c20->Update();
  c20->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Ac_mu vs ttbarmass.png");

  TCanvas *c21 = new TCanvas();
  c21->cd();
  h18->GetXaxis()->SetTitle("P_{t}(GeV/c)");
  h18->GetXaxis()->SetTitleSize(0.04);
  h18->GetXaxis()->SetLabelSize(0.03);
  h18->GetYaxis()->SetTitle("A_{c_{ele}}");
  h18->GetYaxis()->SetTitleSize(0.04);
  h18->GetYaxis()->SetTitleOffset(1.1);
  h18->GetYaxis()->SetLabelSize(0.03);
  h18->SetMarkerSize(1);
  h18->SetMarkerColor(kRed);
  h18->SetLineWidth(1);
  // gStyle->SetEndErrorSize(1);
  // gStyle->SetErrorX(.3);
  // h18->Draw("E1");
  h18->SetOption("E");
  c21->Update();
  c21->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Ac_ele vs Pt_ele.png");

  TCanvas *c22 = new TCanvas();
  c22->cd();
  h19->GetXaxis()->SetTitle("P_{t}(GeV/c)");
  h19->GetXaxis()->SetTitleSize(0.04);
  h19->GetXaxis()->SetLabelSize(0.03);
  h19->GetYaxis()->SetTitle("A_{c_{mu}}");
  h19->GetYaxis()->SetTitleSize(0.04);
  h19->GetYaxis()->SetLabelSize(0.03);
  h19->GetYaxis()->SetTitleOffset(1.1);
  h19->SetMarkerSize(1);
  h19->SetMarkerColor(kRed);
  h19->SetLineWidth(1);
  gStyle->SetEndErrorSize(1);
  gStyle->SetErrorX(.3);
  // h19->Draw("E1");
  h19->SetOption("E1");
  c22->Update();
  c22->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Ac_mu vs Pt_mu.png");

  TCanvas *c23 = new TCanvas();
  c23->cd();
  h20->Draw();
  c23->Update();
  c23->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Pho_pt_ele.png");

  TCanvas *c24 = new TCanvas();
  c24->cd();
  h21->GetXaxis()->SetTitle("#gamma_{p_{t}} (GeV/c)");
  h21->GetXaxis()->SetTitleSize(0.04);
  h21->GetXaxis()->SetLabelSize(0.03);
  h21->GetYaxis()->SetTitle("A_{c_{ele}}");
  h21->GetYaxis()->SetTitleSize(0.04);
  h21->GetYaxis()->SetTitleOffset(1.1);
  h21->GetYaxis()->SetLabelSize(0.03);
  h21->SetMarkerSize(1);
  h21->SetMarkerColor(kRed);
  h21->SetLineWidth(1);
  gStyle->SetEndErrorSize(1);
  gStyle->SetErrorX(.3);
  // h21->Draw("E1");
  h21->SetOption("E1");
  c24->Update();
  c24->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Ac_ele vs Pho_pt_ele.png");

  TCanvas *c25 = new TCanvas();
  c25->cd();
  h22->GetXaxis()->SetTitle("#gamma_{p_{t}} (GeV/c)");
  h22->GetXaxis()->SetTitleSize(0.04);
  h22->GetXaxis()->SetLabelSize(0.03);
  h22->GetYaxis()->SetTitle("A_{c_{mu}}");
  h22->GetYaxis()->SetTitleSize(0.04);
  h22->GetYaxis()->SetLabelSize(0.03);
  h22->GetYaxis()->SetTitleOffset(1.1);
  h22->SetMarkerSize(1);
  h22->SetMarkerColor(kRed);
  h22->SetLineWidth(1);
  gStyle->SetEndErrorSize(1);
  gStyle->SetErrorX(.3);
  // h22->Draw("E1");
  h22->SetOption("E1");
  c25->Update();
  c25->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Ac_mu vs Pho_pt_mu.png");

  TCanvas *c26 = new TCanvas();
  c26->cd();
  h23->Draw();
  c26->Update();
  c26->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/YPho.png");

  TCanvas *c27 = new TCanvas();
  c27->cd();
  h25->Draw();
  c27->Update();
  c27->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Rapidity_ele.png");

  TCanvas *c28 = new TCanvas();
  c28->cd();
  h26->Draw();
  c28->Update();
  c28->SaveAs("/eos/user/s/ssnehshu/plots/QCD_ele/Rapidity_mu.png");

  fl->Close();
  file1->Close();
  std::cout << "numevents" << numevents << std::endl;

  // std::cout << "nEle" << nEle << std::endl;
}
