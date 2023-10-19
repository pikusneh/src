#include "TROOT.h"
#include "TH1.h"

// use namespace std;
// gROOT.SetBatch(true);
void DataElectron()
{

  gROOT->SetBatch(true);

  TH1F *h1 = new TH1F("h1", "elept", 100, 0, 400);
  TH1F *h2 = new TH1F("h2", "muonpt", 100, 0, 400);
  TH1F *h3 = new TH1F("h3", "jetspt", 100, 0, 500);
  // TH1F *h4 = new TH1F("h4","tophadpt",100,-11000,2000);
  // TH1F *h5 = new TH1F("h5","topleppt",100,-11000,2000);
  TH1F *h6 = new TH1F("h6", "tophadeta", 100, -11000, 2000);
  TH1F *h7 = new TH1F("h7", "tophlepeta", 100, -11000, 2000);
  TH1F *h8 = new TH1F("h8", "toplepcharge", 100, -11000, 2000);
  TH1F *h9 = new TH1F("h9", "top", 100, -3, 3);
  TH1F *h10 = new TH1F("h10", "antitop", 100, -3, 3);
  TH1F *h9ele = new TH1F("h9ele", "top", 100, -3, 3);
  TH1F *h10ele = new TH1F("h10ele", "antitop", 100, -3, 3);
  TH1F *h9mu = new TH1F("h9mu", "top", 100, -3, 3);
  TH1F *h10mu = new TH1F("h10mu", "antitop", 100, -3, 3);
  TH1F *h11 = new TH1F("h11", "Rapidity", 4, -3,3);

  TH1F *h14 = new TH1F("h14", "Charge Asymmetry vs M_{t#bar{t}} Electronchannel ", 7, 300, 1000);

  // h9->SetLineColor(kYellow+10);
  //  h10->SetLineColor(kCyan+3);

  h9ele->SetMarkerSize(0.5);
  h10ele->SetMarkerSize(0.5);
  h9mu->SetMarkerSize(0.5);
  h10mu->SetMarkerSize(0.5);

  h9ele->SetMarkerStyle(20);
  h10ele->SetMarkerStyle(20);
  h9mu->SetMarkerStyle(20);
  h10mu->SetMarkerStyle(20);

  TFile *fl2 = new TFile("Ele.root", "RECREATE");

  // TFile *fl = new TFile("TTGamma_Dilepton_2016_AnalysisNtuple.root", "read");
  // TFile *fl1 = new TFile("TTGamma_Dilep_Hist",
  TChain *tr = new TChain("AnalysisTree");
  //  tr->Add("TTGamma_Dilepton_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_b_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_b_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_b_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_b_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_b_2016_AnalysisNtuple_5of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_c_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_c_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_c_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_c_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_c_2016_AnalysisNtuple_5of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_d_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_d_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_d_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_d_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_d_2016_AnalysisNtuple_5of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_e_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_e_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_e_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_e_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_e_2016_AnalysisNtuple_5of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_f_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_f_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_f_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_f_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_f_2016_AnalysisNtuple_5of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_g_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_g_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_g_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_g_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_g_2016_AnalysisNtuple_5of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_h_2016_AnalysisNtuple_1of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_h_2016_AnalysisNtuple_2of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_h_2016_AnalysisNtuple_3of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_h_2016_AnalysisNtuple_4of5.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/Data_SingleEle_h_2016_AnalysisNtuple_5of5.root");
  // TTree *tr = (TTree*)fl-> Get("AnalysisTree");
  Int_t numevents = tr->GetEntries();
  // numevents = 1000;
  // Define the variables
  std::vector<Float_t> *elePt = 0;
  std::vector<Float_t> *muPt = 0;
  std::vector<Float_t> *jetPt = 0;
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
  Float_t PUweight = 0;
  Float_t TopHad_pt = 0;
  Float_t TopLep_pt = 0;
  Float_t TopHad_eta = 0;
  Float_t TopLep_eta = 0;
  Float_t TopLep_phi = 0;
  Float_t TopHad_phi = 0;
  Float_t Mt_blgammaMET = 0; // toplepmass
  Float_t M_bjj = 0;         // tophadmass
  Float_t TopLep_charge = 0;
  bool passPresel_Ele;
  bool passPresel_Mu;
  Float_t Rapidity_T_ele; // taking T for top quark
  Float_t Rapidity_t_ele; // taking t for antitop quark
  Float_t Rapidity_T_mu;
  Float_t Rapidity_t_mu;
  Float_t Rapidity_T;
  Float_t Rapidity_t;
  Float_t Mass_top;
  Float_t Mass_antitop;
  Float_t ttbar_mass;

  ////////////////////////////////Set Branch Status////////////////////////////////////////////////

  tr->SetBranchStatus("*", 0);
  tr->SetBranchStatus("elePt", 1);
  tr->SetBranchStatus("muPt", 1);
  tr->SetBranchStatus("jetPt", 1);
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

  // Set Branch Address
  tr->SetBranchAddress("elePt", &elePt);
  tr->SetBranchAddress("muPt", &muPt);
  tr->SetBranchAddress("jetPt", &jetPt);
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

  float N_plus_ele = 0;
  float N_minus_ele = 0;
  float N_plus_mu = 0;
  float N_minus_mu = 0;
  float N_plus = 0;
  float N_minus = 0;
  TLorentzVector TopLep;
  TLorentzVector TopHad;
  TLorentzVector AntiTopLep;
  TLorentzVector AntiTopHad;
  TLorentzVector Top;
  TLorentzVector AntiTop;
  Float_t YT_ele;
  Float_t Yt_ele;
  Float_t rapidity_diff_ele;
  Float_t YT_mu;
  Float_t Yt_mu;
  Float_t rapidity_diff_mu;
  Float_t YT;
  Float_t Yt;
  Float_t rapidity_diff;
  Float_t Mass, NetMass,NetMass1,NetMass2,NetMass3,NetMass4,NetMass5,NetMass6,NetMass7;

  float N_plus_ele1 = 0;
  float N_plus_ele2 = 0;
  float N_plus_ele3 = 0;
  float N_plus_ele4 = 0;
  float N_plus_ele5 = 0;
  float N_plus_ele6 = 0;
  float N_plus_ele7 = 0;
  
  float N_minus_ele1 = 0;
  float N_minus_ele2 = 0;
  float N_minus_ele3 = 0;
  float N_minus_ele4 = 0;
  float N_minus_ele5 = 0;
  float N_minus_ele6 = 0;
  float N_minus_ele7 = 0;



  // Set the loop
  for (Int_t i = 0; i < numevents; i++)
  {
    tr->GetEntry(i);
    // define weight
   // weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;

    /////////////////////////////////Electron channel////////////////////////////////////////
    if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1)

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

      // rapidity calculation
      Rapidity_T_ele = Top.Y();
      Rapidity_t_ele = AntiTop.Y();
      // absolute value of rapidity
      YT_ele = TMath::Abs(Rapidity_T_ele); // top
      Yt_ele = TMath::Abs(Rapidity_t_ele); // antitop
      rapidity_diff_ele = YT_ele - Yt_ele;
      h11->Fill(rapidity_diff_ele);
      if (rapidity_diff_ele > 0)

      {
        N_plus_ele++;
      }
      if (rapidity_diff_ele < 0)
      {
        N_minus_ele++;
      }
      // Fill N+ and N- for different mass range of Top+antitop

      if (NetMass >= 300 && NetMass < 400)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele1++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele1++;
        }
        NetMass1 = NetMass;
      }

      if (NetMass >= 400 && NetMass < 500)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele2++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele2++;
        }
        NetMass2 = NetMass;
      }

      if (NetMass >= 500 && NetMass < 600)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele3++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele3++;
        }
        NetMass3 = NetMass;
      }
      if (NetMass >= 600 && NetMass < 700)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele4++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele4++;
        }
        NetMass4 = NetMass;
      }
      if (NetMass >= 700 && NetMass < 800)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele5++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele5++;
        }
        NetMass5 = NetMass;
      }
      if (NetMass >= 800 && NetMass < 900)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele6 ++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele6++;
        }
        NetMass6 = NetMass;
      }
      if (NetMass >= 900 && NetMass < 1000)
      {

        if (rapidity_diff_ele > 0)
        {
          N_plus_ele7++;
        }
        if (rapidity_diff_ele < 0)
        {
          N_minus_ele7++;
        }
        NetMass7 = NetMass;
      }
      ///////////////////////////////////////////////////////////////End/////////////////////////

      h1->Fill(elePt->at(0));
      h9ele->Fill(TopLep_eta);  // top quark
      h10ele->Fill(TopHad_eta); // anti top quark

      //////////////////////////////////////Muon Channel////////////////////////////////////////
    }
    if (passPresel_Mu && nJet >= 4 && nBJet >= 1 && nPho == 1)
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
      Rapidity_T_mu = Top.Y();
      Rapidity_t_mu = AntiTop.Y();
      YT_mu = TMath::Abs(Rapidity_T_mu); // top
      Yt_mu = TMath::Abs(Rapidity_t_mu); // antitop
      rapidity_diff_mu = YT_mu - Yt_mu;
      if (rapidity_diff_mu > 0)
      {
        N_plus_mu++;
      }
      if (rapidity_diff_mu < 0)
      {
        N_minus_mu++;
      }

      h2->Fill(muPt->at(0));
      h9mu->Fill(TopHad_eta);  // top quark
      h10mu->Fill(TopLep_eta); // anti top quark
    }
    //////////////////////////////////////////Both channel////////////////////////////////////
    if (passPresel_Mu || (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1))

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
      Rapidity_T = Top.Y();
      Rapidity_t = AntiTop.Y();
      YT = TMath::Abs(Rapidity_T); // top
      Yt = TMath::Abs(Rapidity_t); // antitop
      rapidity_diff = YT - Yt;
      if (rapidity_diff > 0)
      {
        N_plus++;
      }
      if (rapidity_diff < 0)
      {
        N_minus++;
      }
    }

    if (nJet > 0)
    {
      h3->Fill(jetPt->at(0));
    }

    if (TopLep_charge > 0)
    {
      h9->Fill(TopLep_eta);  // top quark
      h10->Fill(TopHad_eta); // anti top quark
    }

    if (TopLep_charge < 0)
    {
      h9->Fill(TopHad_eta);  // top quark
      h10->Fill(TopLep_eta); // antitop quark
    }

    //   h4->Fill(TopHad_pt);
    //  h5->Fill(TopLep_pt);
    //   h6->Fill(TopHad_eta);
    //   h7->Fill(TopLep_eta);
    h8->Fill(TopLep_charge);
  }

  //////////////////////////Value of N+ and N-//////////////////////////////////
  
  // std::cout << "N_ele+ = " << N_plus_ele << std::endl;
  // std::cout << "N_ele- = " << N_minus_ele << std::endl;
  // std::cout << "N+ = " << N_plus << std::endl;
  // std::cout << "N- = " << N_minus << std::endl;

  ////////////////////////Calculate charge asymmetry/////////////////////////////
  ///////////////////////////////////Electron Channel///////////////////////////////

  Float_t Ac, Ac_ele;
  Float_t sum_ele = N_plus_ele + N_minus_ele;
  Float_t diff_ele = N_plus_ele - N_minus_ele;
  Float_t sum = N_plus + N_minus;
  Float_t diff = N_plus - N_minus;

  if (sum_ele != 0)
  {
    Ac_ele = diff_ele / sum_ele;
    std::cout << "Asymmetry for electron channel = " << Ac_ele << std::endl;
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  ////////////////////////////////////////Muon Channel////////////////////////
  Float_t sum_mu = N_plus_mu + N_minus_mu;
  Float_t diff_mu = N_plus_mu - N_minus_mu;

  if (sum_mu != 0)
  {
    Float_t Ac_mu = diff_mu / sum_mu;
    std::cout << "Asymmetry for muon channel = " << Ac_mu << std::endl;
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  /////////////////////////////////////Both Channel////////////////////////////////////////////
  if (sum != 0)
  {
    Ac = diff / sum;
    std::cout << "Asymmetry for DataEle = " << Ac << std::endl;
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
/////////////////////////////////////////////////// Electron channel////////////////////////////////////////////
  Float_t  Ac_ele1, Ac_ele2, Ac_ele3, Ac_ele4, Ac_ele5, Ac_ele6, Ac_ele7;  //////////for mass
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

  //////////////////////////////////Fill h14/////////////////////////////////
  Float_t Delta_Ac_ele1 = 0.0;
  if (sum_ele1 != 0)
  {
    Ac_ele1 = diff_ele1 / sum_ele1;
    Float_t Delta_diff_ele1 = sqrt(sum_ele1);
    Float_t Delta_sum_ele1 = sqrt(sum_ele1);
    Float_t Delta_N_plus_ele1 = sqrt(N_plus_ele1);
    Float_t Delta_N_minus_ele1 = sqrt(N_minus_ele1);
    Float_t Delta_Ac_ele1 = Ac_ele1 * sqrt((Delta_diff_ele1 / diff_ele1) * (Delta_diff_ele1 / diff_ele1) + (Delta_sum_ele1 / sum_ele1) * (Delta_sum_ele1 / sum_ele1));
    // std::cout << "Asymmetry for electron channel = " << Ac_ele1 << std::endl;
    // std::cout << "NetMass1 =" << NetMass1 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass1), Ac_ele1);
    int bin = h14->FindBin(NetMass1);
    h14->SetBinContent(bin, Ac_ele1);
    h14->SetBinError(bin, Delta_Ac_ele1);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  
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
    // h14->SetBinContent(h13->FindBin(NetMass3), Ac_ele3);
    int bin = h14->FindBin(NetMass3);
    h14->SetBinContent(bin, Ac_ele3);
    h14->SetBinError(bin, Delta_Ac_ele3);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  Float_t Delta_Ac_ele4 = 0.0;
  if (sum_ele4 != 0)
  {
    Ac_ele4 = diff_ele4 / sum_ele4;
    Float_t Delta_diff_ele4 = sqrt(sum_ele4);
    Float_t Delta_sum_ele4 = sqrt(sum_ele4);
    Float_t Delta_N_plus_ele4 = sqrt(N_plus_ele4);
    Float_t Delta_N_minus_ele4 = sqrt(N_minus_ele4);
    Float_t Delta_Ac_ele4 = Ac_ele4 * sqrt((Delta_diff_ele4 / diff_ele4) * (Delta_diff_ele4 / diff_ele4) + (Delta_sum_ele4 / sum_ele4) * (Delta_sum_ele4 / sum_ele4));
    //  std::cout << "Asymmetry for electron channel = " << Ac_ele4 << std::endl;
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
    //  std::cout << "NetMass6 =" << NetMass6 << std::endl;
    // h14->SetBinContent(h14->FindBin(NetMass6), Ac_ele6);
    int bin = h14->FindBin(NetMass6);
    h14->SetBinContent(bin, Ac_ele6);
    h14->SetBinError(bin, Delta_Ac_ele6);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

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
 


  h9->Write();
  h10->Write();
  h9ele->Write();
  h10ele->Write();
  h9mu->Write();
  h10mu->Write();
  h11->Write();
  h14->Write();
  // Canvas
  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1->Draw();
  c1->Update();
  c1->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h1e.png");
  TCanvas *c2 = new TCanvas();
  c2->cd();
  h2->Draw();
  c2->Update();
  c2->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h2e.png");
  TCanvas *c3 = new TCanvas();
  c3->cd();
  h3->Draw();
  c3->Update();
  c3->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h3e.png");
  /* TCanvas *c4 = new TCanvas();
   c4->cd();
   h4->Draw();
   c4->Update();
   c4->SaveAs("h4e.png");
   TCanvas *c5 = new TCanvas();
   c5->cd();
   h5->Draw();
   c5->Update();
   c5->SaveAs("h5e.png");*/
  TCanvas *c6 = new TCanvas();
  c6->cd();
  h6->Draw();
  c6->Update();
  c6->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h6e.png");
  TCanvas *c7 = new TCanvas();
  c7->cd();
  h7->Draw();
  c7->Update();
  c7->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h7e.png");
  TCanvas *c8 = new TCanvas();
  c8->cd();
  h8->Draw();
  c8->Update();
  c8->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h8e.png");
  TCanvas *c9 = new TCanvas();
  c9->cd();
  h9->Draw();
  c9->Update();
  c9->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h9e.png");
  TCanvas *c10 = new TCanvas();
  c10->cd();
  h10->Draw();
  c10->Update();
  c10->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h10e.png");

  TCanvas *c11 = new TCanvas();
  c11->cd();
  h9ele->Draw("P");
  c11->Update();
  c11->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h9electrontop.png");

  TCanvas *c12 = new TCanvas();
  c12->cd();
  h10ele->Draw("P");
  c12->Update();
  c12->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h10electronantitop.png");

  TCanvas *c13 = new TCanvas();
  c13->cd();
  h9mu->Draw("P");
  c13->Update();
  c13->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h9mutop.png");

  TCanvas *c14 = new TCanvas();
  c14->cd();
  h10mu->Draw("P");
  c14->Update();
  c14->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/h10muantitop.png");

  TCanvas *c15 = new TCanvas();
  c15->cd();
  h11->Draw();
  c15->Update();
  c15->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/Rapidity.png");


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
  c17->Update();
  c17->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/Ac_ele vs ttbarmass.png");

  fl2->Close();

  std::cout << "numevents" << numevents << std::endl;
  std::cout << "nEle" << nEle << std::endl;
}
