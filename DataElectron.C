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
  TH1F *h11 = new TH1F("h14", "N+ and N-", 7, 300, 1000);

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
  Float_t sub_ele;
  Float_t YT_mu;
  Float_t Yt_mu;
  Float_t sub_mu;
  Float_t YT;
  Float_t Yt;
  Float_t sub;
  Float_t Mass, NetMass;

  // Set the loop
  for (Int_t i = 0; i < numevents; i++)
  {
    tr->GetEntry(i);
    // define weight
    weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;

    /////////////////////////////////Electron channel////////////////////////////////////////
    if (passPresel_Ele && nJet >= 4 && nBJet >= 1 && nPho == 1)

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
      if (sub_ele > 0)
      {
        N_plus_ele++;
      }
      if (sub_ele < 0)
      {
        N_minus_ele++;
      }

      h1->Fill(elePt->at(0), weight);
      h9ele->Fill(TopLep_eta, weight);  // top quark
      h10ele->Fill(TopHad_eta, weight); // anti top quark

      //////////////////////////////////////Muon Channel////////////////////////////////////////
    }
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
      TLorentzVector TopAntiTop = Top + AntiTop;
      Mass = TopAntiTop.M();
      NetMass = Mass;
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
      else
        (TopLep_charge < 0);
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
      sub = YT - Yt;
      if (sub > 0)
      {
        N_plus++;
      }
      if (sub < 0)
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
  h11->Fill(N_plus_ele);
  h11->Fill(N_minus_ele);
  std::cout << "N_ele+ = " << N_plus_ele << std::endl;
  std::cout << "N_ele- = " << N_minus_ele << std::endl;
  std::cout << "N+ = " << N_plus << std::endl;
  std::cout << "N- = " << N_minus << std::endl;

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

  h9->Write();
  h10->Write();
  h9ele->Write();
  h10ele->Write();
  h9mu->Write();
  h10mu->Write();
  h11->Write();
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
  c15->SaveAs("/eos/user/s/ssnehshu/plots/DataEle/N+and N-.png");

  fl2->Close();

  std::cout << "numevents" << numevents << std::endl;
  std::cout << "nEle" << nEle << std::endl;
}
