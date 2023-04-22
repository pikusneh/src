

#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
// use namespace std;
// gROOT.SetBatch(true);
void TTGamma()
{

  gROOT->SetBatch(true);
  // defining histograms

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
  TH1F *h9ele = new TH1F("h9ele", "top", 50, -3, 3);
  TH1F *h10ele = new TH1F("h10ele", "antitop", 50, -3, 3);
  TH1F *h9mu = new TH1F("h9mu", "top", 50, -3, 3);
  TH1F *h10mu = new TH1F("h10mu", "antitop", 50, -3, 3);
  TH1F *h11phi = new TH1F("h11", "topphi", 50, -3, 3);
  TH1F *h12phi = new TH1F("h12", "antitopphi", 50, -3, 3);
  // TH1F *h13 = new TH1F("h13", "Charge Asymmetry Electron channel", 50, -1, 1);
  //TH2F *h14 = new TH2F("h14", "Charge Asymmetry Electronchannel", 7, 300, 1000, 50, 0.08, 0.1);
 TH1F *h14 = new TH1F("h14", "Charge Asymmetry Electronchannel", 7,300,1000);

  TH1F *h15 = new TH1F("h15", "Mass_top", 50, 0, 500);

  TH1F *h16 = new TH1F("h16", "Acele", 50, 0, 1);

  h14->SetMarkerStyle(20);
  // h14->SetMarkerSize(10);
  h9ele->SetOption("hist");
  h9mu->SetOption("hist");
  h10ele->SetOption("hist");
  h10mu->SetOption("hist");
  h11phi->SetOption("hist");
  h12phi->SetOption("hist");
  h16->SetOption("hist");
  // Set line color
  h9ele->SetLineColor(kBlack);
  h9mu->SetLineColor(kBlack);
  h10ele->SetLineColor(kBlack);
  h10mu->SetLineColor(kBlack);

  // Set fill color
  h9ele->SetFillColor(kRed);
  h9mu->SetFillColor(kCyan);
  h10ele->SetFillColor(kBlue);
  h10mu->SetFillColor(kYellow);

  // saving histogram in root files
  TFile *fl = new TFile("TTGamma.root", "RECREATE");

  // TFile *fl = new TFile("TTGamma_Dilepton_2016_AnalysisNtuple.root", "read");
  // TFile *fl1 = new TFile("TTGamma_Dilep_Hist",
  TChain *tr = new TChain("AnalysisTree");
  //  tr->Add("TTGamma_Dilepton_2016_AnalysisNtuple.root");

  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");
  // TTree *tr = (TTree*)fl-> Get("AnalysisTree");
  Int_t numevents = tr->GetEntries();
  //   numevents = 100;
  //  Define the variables
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
  Float_t Mass_top;
  Float_t Mass_antitop;
  Float_t ttbar_mass;

  // Set Branch Status

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
  cout << "numevents" << numevents << std::endl;

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

  // Set the loop

  double weight = 0;

  // Initialize
  float N_plus_ele1=0;
  float N_plus_ele2=0;
  float N_plus_ele3=0;
  float N_plus_ele4=0;
  float N_plus_ele5=0;
  float N_plus_ele6=0;
  float N_plus_ele7=0;
  float N_plus_ele8=0;
  float N_minus_ele1=0;
  float N_minus_ele2=0;
  float N_minus_ele3=0;
  float N_minus_ele4=0;
  float N_minus_ele5=0;
  float N_minus_ele6=0;
  float N_minus_ele7=0;
  float N_minus_ele8=0;
  float N_plus_mu=0;
  float N_minus_mu=0;
  TLorentzVector TopLep;
  TLorentzVector TopHad;
  TLorentzVector AntiTopLep;
  TLorentzVector AntiTopHad;
  TLorentzVector Top;
  TLorentzVector AntiTop;

  // TLorentzVector TopAntiTop = Top + AntiTop;
  Float_t YT_ele;
  Float_t Yt_ele;
  Float_t sub_ele;
  Float_t YT_mu;
  Float_t Yt_mu;
  Float_t sub_mu;
  Float_t Mass, NetMass, NetMass1, NetMass2, NetMass3, NetMass4, NetMass5, NetMass6, NetMass7, NetMass8;
  float Ac1, Ac2;

  // std::cout << "Mass of top" << Mass << std::endl;
  for (Int_t i = 0; i < numevents; i++)
  {
    tr->GetEntry(i);
    // define weight
    weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;

    /////////////////////////////////////Electron Channel////////////////////////////////////////////////////////////////////

    if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1) // selection

    {
      if (TopLep_charge > 0)
      {
        Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
        AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
      }
      else
        (TopLep_charge < 0);
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
      sub_ele = YT_ele - Yt_ele;

     /* if (NetMass < 50)
      {

        if (sub_ele > 0)
        {
          N_plus_ele1++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele1++;
        }
        NetMass1 = NetMass;
      }*/

      if (NetMass >= 300 && NetMass < 400)
      {

        if (sub_ele > 0)
        {
          N_plus_ele2++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele2++;
        }
        NetMass2 = NetMass;
      }

      if (NetMass >= 400 && NetMass < 500)
      {

        if (sub_ele > 0)
        {
          N_plus_ele3++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele3++;
        }
        NetMass3 = NetMass;
      }

      if (NetMass >= 500 && NetMass < 600)
      {

        if (sub_ele > 0)
        {
          N_plus_ele4++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele4++;
        }
        NetMass4 = NetMass;
      }
      if (NetMass >= 600 && NetMass < 700)
      {

        if (sub_ele > 0)
        {
          N_plus_ele5++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele5++;
        }
        NetMass5 = NetMass;
      }
      if (NetMass >= 700 && NetMass < 800)
      {

        if (sub_ele > 0)
        {
          N_plus_ele6++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele6++;
        }
        NetMass6 = NetMass;
      }
      if (NetMass >= 800 && NetMass < 900)
      {

        if (sub_ele > 0)
        {
          N_plus_ele7++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele7++;
        }
        NetMass7 = NetMass;
      }
      if (NetMass >= 900 && NetMass < 1000)
      {

        if (sub_ele > 0)
        {
          N_plus_ele8++;
        }
        if (sub_ele < 0)
        {
          N_minus_ele8++;
        }
        NetMass8 = NetMass;
      }

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
    if (passPresel_Mu && nJet >= 4 && nBJet >= 1 && nPho == 1)

    {
      if (TopLep_charge > 0)
      {
        Top.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
        AntiTop.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
      }

      else
        (TopLep_charge < 0);
      {
        AntiTop.SetPtEtaPhiM(TopLep_pt, TopLep_eta, TopLep_phi, Mt_blgammaMET);
        Top.SetPtEtaPhiM(TopHad_pt, TopHad_eta, TopHad_phi, M_bjj);
      }
      Rapidity_T_mu = Top.Y();
      Rapidity_t_mu = AntiTop.Y();
      YT_mu = TMath::Abs(Rapidity_T_mu); // top
      Yt_mu = TMath::Abs(Rapidity_t_mu); // antitop
      sub_mu = YT_mu - Yt_mu;
      if (sub_mu > 0)
      {
        N_plus_mu++;
      }
      if (sub_mu < 0)
      {
        N_minus_mu++;
      }
      h2->Fill(muPt->at(0), weight);
      h9mu->Fill(TopHad_eta, weight);  // top quark
      h10mu->Fill(TopLep_eta, weight); // anti top quark
      h12phi->Fill(TopHad_phi, weight);
    }
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

  // calculate charge asymmetry/
  ////////////////////////////// Electron channel////////////////////////////////////////////
  Float_t Ac_ele1, Ac_ele2, Ac_ele3, Ac_ele4, Ac_ele5, Ac_ele6, Ac_ele7, Ac_ele8;
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

 /* if (sum_ele1 != 0)
  {
    Ac_ele1 = diff_ele1 / sum_ele1;
    std::cout << "Asymmetry for electron channel = " << Ac_ele1 << std::endl;
    std::cout << "NetMass1 =" << NetMass1 << std::endl;
    //h14->Fill(NetMass1, Ac_ele1);
  
  //  h14->SetBinContent(1, Ac_ele1);
    h14->SetBinContent(h14->FindBin(NetMass1,Ac_ele1), Ac_ele1);

  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }*/

  // h14->Fill(NetMass1, Ac_ele1);

  if (sum_ele2 != 0)
  {
    Ac_ele2 = diff_ele2 / sum_ele2;
    std::cout << "Asymmetry for electron channel = " << Ac_ele2 << std::endl;
    std::cout << "NetMass2 =" << NetMass2 << std::endl;
    //h14->Fill(NetMass2, Ac_ele2);
    //h14->SetBinContent(2, Ac_ele2);
     h14->SetBinContent(h14->FindBin(NetMass2), Ac_ele2);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass2, Ac_ele2);

  if (sum_ele3 != 0)
  {
    Ac_ele3 = diff_ele3 / sum_ele3;
    std::cout << "Asymmetry for electron channel = " << Ac_ele3 << std::endl;
    std::cout << "NetMass3 =" << NetMass3 << std::endl;
    //h14->Fill(NetMass3, Ac_ele3);
   // h14->SetBinContent(3, Ac_ele3);
     h14->SetBinContent(h14->FindBin(NetMass3), Ac_ele3);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass3, Ac_ele3);

  if (sum_ele4 != 0)
  {
    Ac_ele4 = diff_ele4 / sum_ele4;
    std::cout << "Asymmetry for electron channel = " << Ac_ele4 << std::endl;
    std::cout << "NetMass4 =" << NetMass4 << std::endl;
   // h14->Fill(NetMass4, Ac_ele4);
   // h14->SetBinContent(4, Ac_ele4);
    h14->SetBinContent(h14->FindBin(NetMass4), Ac_ele4);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  if (sum_ele5 != 0)
  {
    Ac_ele5 = diff_ele5 / sum_ele5;
    std::cout << "Asymmetry for electron channel = " << Ac_ele5 << std::endl;
    std::cout << "NetMass5 =" << NetMass5 << std::endl;
    //h14->Fill(NetMass5, Ac_ele5);
    //h14->SetBinContent(5, Ac_ele5);
     h14->SetBinContent(h14->FindBin(NetMass5), Ac_ele5);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  // h14->Fill(NetMass1, Ac_ele1);

  if (sum_ele6 != 0)
  {
    Ac_ele6 = diff_ele6 / sum_ele6;
    std::cout << "Asymmetry for electron channel = " << Ac_ele6 << std::endl;
    std::cout << "NetMass6 =" << NetMass6 << std::endl;
   // h14->Fill(NetMass6, Ac_ele6);
    h14->SetBinContent(6, Ac_ele6);
     h14->SetBinContent(h14->FindBin(NetMass6), Ac_ele6);
  }

  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass2, Ac_ele2);

  if (sum_ele7 != 0)
  {
    Ac_ele7 = diff_ele7 / sum_ele7;
    std::cout << "Asymmetry for electron channel = " << Ac_ele7 << std::endl;
    std::cout << "NetMass7 =" << NetMass7 << std::endl;
   // h14->Fill(NetMass7, Ac_ele3);
     h14->SetBinContent(7, Ac_ele7);
      h14->SetBinContent(h14->FindBin(NetMass7), Ac_ele7);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  // h14->Fill(NetMass3, Ac_ele3);

  if (sum_ele8 != 0)
  {
    Ac_ele8 = diff_ele8 / sum_ele8;
    std::cout << "Asymmetry for electron channel = " << Ac_ele8 << std::endl;
    std::cout << "NetMass8 =" << NetMass8 << std::endl;
    h14->Fill(NetMass8, Ac_ele8);
     h14->SetBinContent(h14->FindBin(NetMass8), Ac_ele8);
  }
  else
  {
    std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  for (int i = 1; i <= h14->GetNbinsX(); i++)
  {
    std::cout << "Bin " << i << " content: " << h14->GetBinContent(i) << std::endl;
  }

  /////////////////////////////////////////Muon channel////////////////////////////////////////////////////////////
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

  // h9->Write();
  // h10->Write();

  h9ele->Write();
  h10ele->Write();
  h9mu->Write();
  h10mu->Write();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1->Draw();
  c1->Update();
  c1->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h1.png");
  TCanvas *c2 = new TCanvas();
  c2->cd();
  h2->Draw();
  c2->Update();
  c2->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h2.png");
  TCanvas *c3 = new TCanvas();
  c3->cd();
  h3->Draw();
  c3->Update();
  c3->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h3.png");
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
 c6->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h6.png");
 TCanvas *c7 = new TCanvas();
 c7->cd();
 h7->Draw();
 c7->Update();
 c7->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h7.png");*/
  TCanvas *c8 = new TCanvas();
  c8->cd();
  h8->Draw();
  c8->Update();
  c8->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h8.png");
  TCanvas *c9 = new TCanvas();
  c9->cd();
  h9->Draw();
  c9->Update();
  c9->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h9.png");
  TCanvas *c10 = new TCanvas();
  c10->cd();
  h10->Draw();
  c10->Update();
  c10->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h10.png");

  TCanvas *c11 = new TCanvas();
  c11->cd();
  h9ele->Draw();
  c11->Update();
  c11->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h9electrontop.png");

  TCanvas *c12 = new TCanvas();
  c12->cd();
  h10ele->Draw();
  c12->Update();
  c12->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h10electronantitop.png");

  TCanvas *c13 = new TCanvas();
  c13->cd();
  h9mu->Draw();
  c13->Update();
  c13->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h9mutop.png");

  TCanvas *c14 = new TCanvas();
  c14->cd();
  h10mu->Draw();
  c14->Update();
  c14->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h10muantitop.png");

  TCanvas *c15 = new TCanvas();
  c15->cd();
  h11phi->Draw();
  c15->Update();
  c15->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h11phi.png");

  TCanvas *c16 = new TCanvas();
  c16->cd();
  h12phi->Draw();
  c16->Update();
  c16->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/h12phi.png");

  TCanvas *c17 = new TCanvas();
  c17->cd();
  h14->Draw();
  c17->Update();
  c17->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/chargeasymm_ele.png");

  TCanvas *c18 = new TCanvas();
  c18->cd();
  h15->Draw();
  c18->Update();
  c18->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/masstop.png");

  TCanvas *c19 = new TCanvas();
  c19->cd();
  h16->Draw();
  c19->Update();
  c19->SaveAs("/eos/user/s/ssnehshu/plots/TTGamma/Acele.png");

  fl->Close();
  std::cout << "numevents" << numevents << std::endl;

  // std::cout << "nEle" << nEle << std::endl;
}
