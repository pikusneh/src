#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
//#include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"

void GenLevel()
{
  
  TH1F *h1 = new TH1F("h1", "toppt", 100, 0, 400);
  TH1F *h2 = new TH1F("h2", "topantitop_pt", 100, 0, 400);
  TH1F *h3 = new TH1F("h3", "topantitop_mass", 100, 200, 1200);
  TH1F *h4 = new TH1F("h4", "Ac vs topantitop_mass", 4, 300, 1100);
  TH1F *h5 = new TH1F("h5", "Top Rapidity", 6, -3, 3);
  TH1F *h6 = new TH1F("h6", "AntiTop Rapidity", 6, -3, 3);
  TH1F *Rapidity_diff_Gen = new TH1F("Rapidity_diff_Gen", "Rapidity_Diff_Gen", 20, -3,3);

  Rapidity_diff_Gen->SetOption("hist");
//  Rapidity_diff_Gen->SetLineColor(kRed);

  TFile *fl = new TFile("GenLevel.root", "RECREATE");
  TChain *tr = new TChain("AnalysisTree");

  TFile *file1 = new TFile("output_GenLevel.root", "RECREATE");
  TTree *tree1 = new TTree("tree1", "tree1");

  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_SingleLept_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Hadronic_2016_AnalysisNtuple.root");
  tr->Add("/eos/user/s/ssnehshu/topquarksample/2016/TTGamma_Dilepton_2016_AnalysisNtuple.root");

  int numevents = tr->GetEntries();

  int nGenPart;
  std::vector<int> *genPDGID = 0;
  std::vector<int> *genMomIdx = 0;
  std::vector<Float_t> *genMass = 0;
  std::vector<Float_t> *genPhi = 0;
  std::vector<Float_t> *genEta = 0;
  std::vector<Float_t> *genPt = 0;
  float evtWeight = 0;
  double weight = 0;

  tr->SetBranchStatus("*", 0);

  tr->SetBranchStatus("nGenPart", 1);
  tr->SetBranchStatus("genPDGID", 1);
  tr->SetBranchStatus("genMomIdx", 1);
  tr->SetBranchStatus("genMass", 1);
  tr->SetBranchStatus("genPhi", 1);
  tr->SetBranchStatus("genEta", 1);
  tr->SetBranchStatus("genPt", 1);
  tr->SetBranchStatus("evtWeight", 1);

  tr->SetBranchAddress("nGenPart", &nGenPart);
  tr->SetBranchAddress("genPDGID", &genPDGID);
  tr->SetBranchAddress("genMomIdx", &genMomIdx);
  tr->SetBranchAddress("genMass", &genMass);
  tr->SetBranchAddress("genPhi", &genPhi);
  tr->SetBranchAddress("genEta", &genEta);
  tr->SetBranchAddress("genPt", &genPt);
  tr->SetBranchAddress("evtWeight", &evtWeight);

  float N_plus, N_minus; //for single value calculation

  float N_plus1, N_minus1;
  float N_plus2, N_minus2;
  float N_plus3, N_minus3;
  float N_plus4, N_minus4;
  float N_plus5, N_minus5;
  float N_plus6, N_minus6;
  float N_plus7, N_minus7;

  //T is for top
  //t is for anti top

  float Mass_Top, Mass_AntiTop, Rapidity_T, Rapidity_t, Pt_top, Pt_antitop, Pt_topantitop;
  float NetMass, NetMass1, NetMass2, NetMass3, NetMass4, NetMass5, NetMass6, NetMass7;
  float Ac1, Ac2, Ac3, Ac4, Ac5, Ac6, Ac7, Ac;

  for (int j = 0; j < numevents; j++)
  {
    int topindex ;
    int antitopindex ;
    tr->GetEntry(j);
    weight = evtWeight;
    TLorentzVector Top;
    TLorentzVector AntiTop;

    // std::cout << "nGenPart: genpdgidSize" << nGenPart<<"     " <<genPDGID->size() <<std::endl;
    for (int i = 0; i < nGenPart; i++)
    {
      //        std::cout << "PDGID[" << i << "]: " << (*genPDGID)[i] << std::endl;
      if (genPDGID->at(i) == 6)
      {
        topindex = i;
      }
      if (genPDGID->at(i) == -6)
      {
        antitopindex = i;
      }
    }

    Top.SetPtEtaPhiM((*genPt)[topindex], (*genEta)[topindex], (*genPhi)[topindex], (*genMass)[topindex]);
    AntiTop.SetPtEtaPhiM((*genPt)[antitopindex], (*genEta)[antitopindex], (*genPhi)[antitopindex], (*genMass)[antitopindex]);
    TLorentzVector TopAntiTop = Top + AntiTop;

    Mass_Top = Top.M();
    Mass_AntiTop = AntiTop.M();
    Pt_top = Top.Pt();
    Pt_antitop = AntiTop.Pt();
    Pt_topantitop = TopAntiTop.Pt();
    NetMass = TopAntiTop.M();
   
    ///////////////////////////////rapidity calculation/////////////////////////////////////////////////////////////
    Rapidity_T = Top.Rapidity();
    Rapidity_t = AntiTop.Rapidity();
    float YT = TMath::Abs(Rapidity_T); // top rapidity
    float Yt = TMath::Abs(Rapidity_t); // antitop rapidity
    h5->Fill(YT,weight);
    h6->Fill(Yt,weight);
    float rapidity_diff = YT - Yt;  
    Rapidity_diff_Gen->Fill(rapidity_diff,weight);

    Rapidity_diff_Gen->Scale(1.0 / Rapidity_diff_Gen->Integral());

    // std::cout<<"rapidity_diff = "<< rapidity_diff <<std::endl;

    if (NetMass >= 300 && NetMass < 400)
    {
      if (rapidity_diff > 0)
      {
        N_plus1 = N_plus1 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus1 = N_minus1 + weight;
      }
      NetMass1 = NetMass;
    }

    if (NetMass >= 400 && NetMass < 500)
    {
      if (rapidity_diff > 0)
      {
        N_plus2 = N_plus2 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus2 = N_minus2 + weight;
      }
      NetMass2 = NetMass;
    }
    if (NetMass >= 500 && NetMass < 600)
    {
      if (rapidity_diff > 0)
      {
        N_plus3 = N_plus3 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus3 = N_minus3 + weight;
      }
      NetMass3 = NetMass;
    }

    if (NetMass >= 600 && NetMass < 700)
    {
      if (rapidity_diff > 0)
      {
        N_plus4 = N_plus4 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus4 = N_minus4 + weight;
      }
      NetMass4 = NetMass;
    }

    if (NetMass >= 700 && NetMass < 800)
    {
      if (rapidity_diff > 0)
      {
        N_plus5 = N_plus5 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus5 = N_minus5 + weight;
      }
      NetMass5 = NetMass;
    }

    if (NetMass >= 800 && NetMass < 900)
    {
      if (rapidity_diff > 0)
      {
        N_plus6 = N_plus6 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus6 = N_minus6 + weight;
      }
      NetMass6 = NetMass;
    }

    if (NetMass >= 900 && NetMass < 1000)
    {
      if (rapidity_diff > 0)
      {
        N_plus7 = N_plus7 + weight;
      }
      if (rapidity_diff < 0)
      {
        N_minus7 = N_minus7 + weight;
      }
      NetMass7 = NetMass;
    }

    //  std::cout<<"N_plus1 = "<< N_plus1 <<std::endl;
    //  std::cout<<"N_minus1 = "<< N_minus1 <<std::endl;
  }
  // std::cout << "Rapidity diff = " << rapidity_diff << std::endl;
  // std::cout << "N+1 = " << N_plus1 << std::endl;
  // std::cout << "N-1 = " << N_minus1 << std::endl;
  // std::cout << "N+2 = " << N_plus2 << std::endl;
  // std::cout << "N-2 = " << N_minus2 << std::endl;
  // std::cout << "N+3 = " << N_plus3 << std::endl;
  // std::cout << "N-3 = " << N_minus3 << std::endl;
  // std::cout << "N+4 = " << N_plus4 << std::endl;
  // std::cout << "N-4 = " << N_minus4 << std::endl;
  // std::cout << "N+5 = " << N_plus5 << std::endl;
  // std::cout << "N-5 = " << N_minus5 << std::endl;
  // std::cout << "N+6 = " << N_plus6 << std::endl;
  // std::cout << "N-6 = " << N_minus6 << std::endl;
  // std::cout << "N+7 = " << N_plus7 << std::endl;
  // std::cout << "N-7 = " << N_minus7 << std::endl;

  float sum = N_plus + N_minus;
  float diff = N_plus - N_minus;
  ///////////////////////////////////////////
  float sum1 = N_plus1 + N_minus1;
  float diff1 = N_plus1 - N_minus1;
  float sum2 = N_plus2 + N_minus2;
  float diff2 = N_plus2 - N_minus2;
  float sum3 = N_plus3 + N_minus3;
  float diff3 = N_plus3 - N_minus3;
  float sum4 = N_plus4 + N_minus4;
  float diff4 = N_plus4 - N_minus4;
  float sum5 = N_plus5 + N_minus5;
  float diff5 = N_plus5 - N_minus5;
  float sum6 = N_plus6 + N_minus6;
  float diff6 = N_plus6 - N_minus6;
  float sum7 = N_plus7 + N_minus7;
  float diff7 = N_plus7 - N_minus7;

  // std::cout << "diff1 = " << diff1 << std::endl;
  // std::cout << "sum1 = " << sum1 << std::endl;
  //  std::cout<<"diff2 = "<< diff2 <<std::endl;
  //  std::cout<<"sum2 = "<< sum2 <<std::endl;
  //  std::cout<<"diff3 = "<< diff3 <<std::endl;
  //  std::cout<<"sum3 = "<< sum3 <<std::endl;
 float Delta_Ac1 = 0.0;
  if (sum1 != 0)
  {
    Ac1 = diff1 / sum1;
    float Delta_diff1 = sqrt(sum1);
    float Delta_sum1 = sqrt(sum1);
    float Delta_N_plus1 = sqrt(N_plus1);
    float Delta_N_minus1 = sqrt(N_minus1);
    float Delta_Ac1 = Ac1 * sqrt((Delta_diff1 / diff1) * (Delta_diff1 / diff1) + (Delta_sum1 / sum1) * (Delta_sum1 / sum1));
    int bin = h4->FindBin(NetMass1);
    h4->SetBinContent(bin, Ac1);
    h4->SetBinError(bin, Delta_Ac1);
    std::cout << "Ac1 = " << Ac1 << std::endl;
  }
  else
  {
    //  std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  float Delta_Ac2 =0.0;
  if (sum2 != 0)
  {
    Ac2 = diff2 / sum2;
    float Delta_diff2 = sqrt(sum2);
    float Delta_sum2 = sqrt(sum2);
    float Delta_N_plus2 = sqrt(N_plus2);
    float Delta_N_minus2 = sqrt(N_minus2);
    float Delta_Ac2 = Ac2 * sqrt((Delta_diff2 / diff2) * (Delta_diff2 / diff2) + (Delta_sum2 / sum2) * (Delta_sum2 / sum2));
    int bin = h4->FindBin(NetMass2);
    h4->SetBinContent(bin, Ac2);
    h4->SetBinError(bin, Delta_Ac2);
    std::cout << "Ac2 = " << Ac2 << std::endl;
  }
  else
  {
    //   std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  float Delta_Ac3 = 0.0;
  if (sum3 != 0)
  {
    Ac3 = diff3 / sum3;
    float Delta_diff3 = sqrt(sum3);
    float Delta_sum3 = sqrt(sum3);
    float Delta_N_plus3 = sqrt(N_plus3);
    float Delta_N_minus3 = sqrt(N_minus3);
    float Delta_Ac3 = Ac3 * sqrt((Delta_diff3 / diff3) * (Delta_diff3 / diff3) + (Delta_sum3 / sum3) * (Delta_sum3 / sum3));
    int bin = h4->FindBin(NetMass3);
    h4->SetBinContent(bin, Ac3);
    h4->SetBinError(bin, Delta_Ac3);
    std::cout << "Ac3 = " << Ac3 << std::endl;
  }
  else
  {
    //   std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  float Delta_Ac4 =0.0;
  if (sum4 != 0)
  {
    Ac4 = diff4 / sum4;
    float Delta_diff4 = sqrt(sum4);
    float Delta_sum4 = sqrt(sum4);
    float Delta_N_plus4 = sqrt(N_plus4);
    float Delta_N_minus4 = sqrt(N_minus4);
    float Delta_Ac4 = Ac4 * sqrt((Delta_diff4 / diff4) * (Delta_diff4 / diff4) + (Delta_sum4 / sum4) * (Delta_sum4 / sum4));
    int bin = h4->FindBin(NetMass4);
    h4->SetBinContent(bin, Ac4);
    h4->SetBinError(bin, Delta_Ac4);
    std::cout << "Ac4 = " << Ac4 << std::endl;
  }
  else
  {
    //   std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
float Delta_Ac5 = 0.0;
  if (sum5 != 0)
  {
    Ac5 = diff5 / sum5;
    float Delta_diff5 = sqrt(sum5);
    float Delta_sum5 = sqrt(sum5);
    float Delta_N_plus5 = sqrt(N_plus5);
    float Delta_N_minus5 = sqrt(N_minus5);
    float Delta_Ac5 = Ac5 * sqrt((Delta_diff5 / diff5) * (Delta_diff5 / diff5) + (Delta_sum5 / sum5) * (Delta_sum5 / sum5));
    int bin = h4->FindBin(NetMass5);
    h4->SetBinContent(bin, Ac5);
    h4->SetBinError(bin, Delta_Ac5);
    std::cout << "Ac5 = " << Ac5 << std::endl;
  }
  else
  {
    //   std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
 float Delta_Ac6 = 0.0;
  if (sum6 != 0)
  {
    Ac6 = diff6 / sum6;
    float Delta_diff6 = sqrt(sum6);
    float Delta_sum6 = sqrt(sum6);
    float Delta_N_plus6 = sqrt(N_plus6);
    float Delta_N_minus6 = sqrt(N_minus6);
    float Delta_Ac6 = Ac6 * sqrt((Delta_diff6 / diff6) * (Delta_diff6 / diff6) + (Delta_sum6 / sum6) * (Delta_sum6 / sum6));
    int bin = h4->FindBin(NetMass6);
    h4->SetBinContent(bin, Ac6);
    h4->SetBinError(bin, Delta_Ac6);
    std::cout << "Ac6 = " << Ac6 << std::endl;
  }
  else
  {
    //   std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }
  float Delta_Ac7 = 0.0;
  if (sum7 != 0)
  {
    Ac7 = diff7 / sum7;
    float Delta_diff7 = sqrt(sum7);
    float Delta_sum7 = sqrt(sum7);
    float Delta_N_plus7 = sqrt(N_plus7);
    float Delta_N_minus7 = sqrt(N_minus7);
    float Delta_Ac7 = Ac7 * sqrt((Delta_diff7 / diff7) * (Delta_diff7 / diff7) + (Delta_sum7 / sum7) * (Delta_sum7 / sum7));
    int bin = h4->FindBin(NetMass7);
    h4->SetBinContent(bin, Ac7);
    h4->SetBinError(bin, Delta_Ac7);
    std::cout << "Ac7 = " << Ac7 << std::endl;
  }
  else
  {
    //   std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
  }

  //   std::cout<<"topantitop_mass = "<< NetMass <<std::endl;
  ////////////////filling the histogram////
  h1->Fill(Pt_top);
  h2->Fill(Pt_topantitop);
  h3->Fill(NetMass);
  // float Ac = sum / diff;
  // std::cout << "Ac = " << Ac << std::endl;
  h4->Write();
  Rapidity_diff_Gen->Write();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  h1->Draw();
  c1->Update();
  c1->SaveAs("/eos/user/s/ssnehshu/GenLevel/Top_Pt.png");

  TCanvas *c2 = new TCanvas();
  c2->cd();
  h2->Draw();
  c2->Update();
  c2->SaveAs("/eos/user/s/ssnehshu/GenLevel/TopAntitop_Pt.png");

  TCanvas *c3 = new TCanvas();
  c3->cd();
  h3->Draw();
  c3->Update();
  c3->SaveAs("/eos/user/s/ssnehshu/GenLevel/TopAntitop_mass.png");

  TCanvas *c4 = new TCanvas();
  c4->cd();
  h4->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
  h4->GetXaxis()->SetTitleSize(0.04);
  h4->GetXaxis()->SetLabelSize(0.03);
  h4->GetYaxis()->SetTitle("A_{c_{ele}}");
  h4->GetYaxis()->SetTitleOffset(1.1);
  h4->GetYaxis()->SetTitleSize(0.04);
  h4->GetYaxis()->SetLabelSize(0.03);

  h4->SetMarkerSize(1);
  h4->SetMarkerColor(kRed);
  h4->SetLineWidth(1);
  // gStyle->SetEndErrorSize(1);
  //  gStyle->SetErrorX(.3);

  h4->Draw("E1");
  c4->Update();
  c4->SaveAs("/eos/user/s/ssnehshu/GenLevel/Ac_ele vs ttbarmass.png");

  TCanvas *c5 = new TCanvas();
  c5->cd();
  h5->Draw();
  c5->Update();
  c5->SaveAs("/eos/user/s/ssnehshu/GenLevel/Top_Rapidity.png");

  TCanvas *c6 = new TCanvas();
  c6->cd();
  h6->Draw();
  c6->Update();
  c6->SaveAs("/eos/user/s/ssnehshu/GenLevel/Antitop_Rapidity.png");


  TCanvas *c7 = new TCanvas();
  c7->cd();
  Rapidity_diff_Gen->Draw();
  c7->Update();
  c7->SaveAs("/eos/user/s/ssnehshu/GenLevel/Rapidity_Gen.png");

  fl->Close();
  file1->Close();
}
