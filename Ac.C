

#include "TROOT.h"

void Ac()
{
    TFile *fl1 = new TFile("output_TTGamma.root", "read");
    TH1F *hist1 = (TH1F *)fl1->Get("h14");
    TH1F *hist3 = (TH1F *)fl1->Get("h17");
    TH1F *hist5 = (TH1F *)fl1->Get("h18");
    TH1F *hist7 = (TH1F *)fl1->Get("h19");
    TH1F *hist9 = (TH1F *)fl1->Get("h21");
    TH1F *hist11 = (TH1F *)fl1->Get("h22");
    TTree *tree1 = (TTree *)fl1->Get("tree1"); // tree1 is for TTGamma sample

    TFile *fl2 = new TFile("output_TTbar.root", "read");
    TH1F *hist2 = (TH1F *)fl2->Get("h14");
    TH1F *hist4 = (TH1F *)fl2->Get("h17");
    TH1F *hist6 = (TH1F *)fl2->Get("h18");
    TH1F *hist8 = (TH1F *)fl2->Get("h19");
    TH1F *hist10 = (TH1F *)fl2->Get("h21");
    TH1F *hist12 = (TH1F *)fl2->Get("h22");
    TTree *tree2 = (TTree *)fl2->Get("tree2"); // tree2 is for TTbar sample

    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);
    hist1->SetFillColor(38);
    hist1->SetFillStyle(3144);

    hist3->SetLineColor(kRed);
    hist4->SetLineColor(kBlue);

    hist5->SetLineColor(kRed);
    hist6->SetLineColor(kBlue);

    hist7->SetLineColor(kRed);
    hist8->SetLineColor(kBlue);

    hist9->SetLineColor(kRed);
    hist10->SetLineColor(kBlue);

    hist11->SetLineColor(kRed);
    hist12->SetLineColor(kBlue);

    float N_plus1, N_minus1, N_plus2, N_minus2;
    float N_plus3, N_minus3, N_plus4, N_minus4;

    tree1->SetBranchAddress("N_plus1", &N_plus1); // N+1 na N-1 values are when # of photon is zero
    tree1->SetBranchAddress("N_minus1", &N_minus1);
    tree1->SetBranchAddress("N_plus2", &N_plus2); // N+2 and N-2 values are when # of photon is one.
    tree1->SetBranchAddress("N_minus2", &N_minus2);

    tree2->SetBranchAddress("N_plus3", &N_plus3); // N+3 na N-3 values are when # of photon is zero
    tree2->SetBranchAddress("N_minus3", &N_minus3);
    tree2->SetBranchAddress("N_plus4", &N_plus4); // N+4 and N-4 values are when # of photon is one.
    tree2->SetBranchAddress("N_minus4", &N_minus4);

    int nEntries1 = tree1->GetEntries();
    int nEntries2 = tree2->GetEntries();

    // Loop over the entries and read the value of N_plus
    for (int i = 0; i < nEntries1; i++)
    {
        tree1->GetEntry(i);
    }

    for (int j = 0; j < nEntries2; j++)
    {
        tree2->GetEntry(j);
    }

    ///////////////////////////////////////////////////Calculating charge asymmetry for all 1,2,3,4////////////////////////////////////////////////////////
    float Ac_1, Ac_2, Ac_3, Ac_4;
    float Sum1 = N_plus1 + N_minus1;
    float Diff1 = N_plus1 - N_minus1;
    float Sum2 = N_plus2 + N_minus2;
    float Diff2 = N_plus2 - N_minus2;
    float Sum3 = N_plus3 + N_minus3;
    float Diff3 = N_plus3 - N_minus3;
    float Sum4 = N_plus4 + N_minus4;
    float Diff4 = N_plus4 - N_minus4;

    if (Sum1 != 0)
    {
        Ac_1 = Diff1 / Sum1;
        std::cout << "Ac_1 = " << Ac_1 << std::endl;
    }

    if (Sum2 != 0)
    {
        Ac_2 = Diff2 / Sum2;
        std::cout << "Ac_2 = " << Ac_2 << std::endl;
    }

    if (Sum3 != 0)
    {
        Ac_3 = Diff3 / Sum3;
        std::cout << "Ac_3 = " << Ac_3 << std::endl;
    }

    if (Sum4 != 0)
    {
        Ac_4 = Diff4 / Sum4;
        std::cout << "Ac_4 = " << Ac_4 << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // When number of photon is zero

    Float_t Ac1;

    Float_t N_plus_net1 = N_plus1 + N_plus3;
    Float_t N_minus_net1 = N_minus1 + N_minus3;

    Float_t sum1 = N_plus_net1 + N_minus_net1;
    Float_t diff1 = N_plus_net1 - N_minus_net1;
    Float_t Delta_sum1, Delta_diff1;

    if (sum1 != 0)
    {
        Ac1 = diff1 / sum1;
        Delta_diff1 = sqrt(sum1);
        Delta_sum1 = sqrt(sum1);
        Float_t Delta_N_plus_net1 = sqrt(N_plus1 + N_plus3);
        Float_t Delta_N_minus_net1 = sqrt(N_minus1 + N_minus3);
        Float_t Delta_Ac1 = Ac1 * sqrt((Delta_diff1 / diff1) * (Delta_diff1 / diff1) + (Delta_sum1 / sum1) * (Delta_sum1 / sum1));

        std::cout << "Ac1 = " << Ac1 << " ± " << Delta_Ac1 << std::endl;
    }
    else
    {
        std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
    }
    // when number of photon is one
    Float_t Ac2;

    Float_t N_plus_net2 = N_plus2 + N_plus4;
    Float_t N_minus_net2 = N_minus2 + N_minus4;

    Float_t sum2 = N_plus_net2 + N_minus_net2;
    Float_t diff2 = N_plus_net2 - N_minus_net2;
    Float_t Delta_sum2, Delta_diff2;

    if (sum2 != 0)
    {
        Ac2 = diff2 / sum2;
        Delta_diff2 = sqrt(sum2);
        Delta_sum2 = sqrt(sum2);
        Float_t Delta_N_plus_net2 = sqrt(N_plus2 + N_plus4);
        Float_t Delta_N_minus_net2 = sqrt(N_minus2 + N_minus4);
        Float_t Delta_Ac2 = Ac2 * sqrt((Delta_diff2 / diff2) * (Delta_diff2 / diff2) + (Delta_sum2 / sum2) * (Delta_sum2 / sum2));
        std::cout << "Ac2 = " << Ac2 << " ± " << Delta_Ac2 << std::endl;
    }
    else
    {
        std::cout << "Cannot calculate asymmetry: sum is zero" << std::endl;
    }

    /*std::cout << "N_plus1 for event " << N_plus1 << std::endl;
     std::cout << "N_minus1 for event " << N_minus1 << std::endl;
     std::cout << "N_plus2 for event " << N_plus2 << std::endl;
     std::cout << "N_minus2 for event " << N_minus2 << std::endl;
     std::cout << "N_plus3 for event " << N_plus3 << std::endl;
     std::cout << "N_minus3 for event " << N_minus3 << std::endl;
     std::cout << "N_plus4 for event " << N_plus4 << std::endl;
     std::cout << "N_minus4 for event " << N_minus4 << std::endl;*/

    TCanvas *c1 = new TCanvas("c1", "Charge Asymmetry vs TTbarmass Electron Channel");
    hist1->Draw();

    hist2->Draw(" SAME");
    hist1->SetStats(kFALSE);
    hist2->SetStats(kFALSE);
    hist1->GetXaxis()->SetTitleSize(0.04);
    hist1->GetYaxis()->SetTitleOffset(1.1);
    hist1->GetYaxis()->SetLabelSize(0.03);
    hist1->GetYaxis()->SetTitleSize(0.04);
    hist1->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
    hist1->GetYaxis()->SetTitle("A_{c_{ele}}");
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend1->AddEntry(hist1, "TTGamma", "E1P");
    legend1->AddEntry(hist2, "TTbar", "E1P");
    legend1->Draw("same");
    c1->Update();
    c1->SaveAs("/eos/user/s/ssnehshu/temp/Ac_vs_ttbarmass_compare_Electron_Channel.png");

    TCanvas *c2 = new TCanvas("c2", "Charge Asymmetry vs TTbarmass Muon Channel");
    hist3->Draw();
    hist4->Draw("same");
    hist3->SetStats(kFALSE);
    hist4->SetStats(kFALSE);
    hist3->GetXaxis()->SetTitleSize(0.04);
    hist3->GetYaxis()->SetTitleOffset(1.1);
    hist3->GetYaxis()->SetLabelSize(0.03);
    hist3->GetYaxis()->SetTitleSize(0.04);
    hist3->GetXaxis()->SetTitle("t#bar{t}_{mass} (GeV)");
    hist3->GetYaxis()->SetTitle("A_{c_{mu}}");
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(hist3, "TTGamma", "E1P");
    legend2->AddEntry(hist4, "TTbar", "E1P");
    legend2->Draw("same");
    c2->Update();
    c2->SaveAs("/eos/user/s/ssnehshu/temp/Ac_vs_ttbarmass_compare_Muon_Channel.png");

    TCanvas *c3 = new TCanvas("c3", "Charge Asymmetry vs TTbar Transverse Momentum Electron Channel");
    hist5->Draw();
    hist6->Draw("same");
    hist5->SetStats(kFALSE);
    hist6->SetStats(kFALSE);
    hist5->GetXaxis()->SetTitleSize(0.04);
    hist5->GetYaxis()->SetTitleOffset(1.1);
    hist5->GetYaxis()->SetLabelSize(0.03);
    hist5->GetYaxis()->SetTitleSize(0.04);
    hist5->GetXaxis()->SetTitle("t#bar{t}_{P_{t}}(GeV/c)");
    hist5->GetYaxis()->SetTitle("A_{c_{ele}}");
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend3->AddEntry(hist5, "TTGamma", "E1P");
    legend3->AddEntry(hist6, "TTbar", "E1P");
    legend3->Draw("same");
    c3->Update();
    c3->SaveAs("/eos/user/s/ssnehshu/temp/Ac_vs_ttbar_transversemomentum_compare_Electron_Channel.png");

    TCanvas *c4 = new TCanvas("c4", "Charge Asymmetry vs TTbar Transverse Momentum Muon Channel");
    hist7->Draw();
    hist8->Draw("same");
    hist7->SetStats(kFALSE);
    hist8->SetStats(kFALSE);
    hist7->GetXaxis()->SetTitleSize(0.04);
    hist7->GetYaxis()->SetTitleOffset(1.1);
    hist7->GetYaxis()->SetLabelSize(0.03);
    hist7->GetYaxis()->SetTitleSize(0.04);
    hist7->GetXaxis()->SetTitle("t#bar{t}_{P_{t}}(GeV/c)");
    hist7->GetYaxis()->SetTitle("A_{c_{mu}}");
    TLegend *legend4 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend4->AddEntry(hist7, "TTGamma", "E1P");
    legend4->AddEntry(hist8, "TTbar", "E1P");
    legend4->Draw("same");
    c4->Update();
    c4->SaveAs("/eos/user/s/ssnehshu/temp/Ac_vs_ttbar_transversemomentum_compare_Muon_Channel.png");

    TCanvas *c5 = new TCanvas("c5", "Charge Asymmetry vs Photon transverse momentum Electron Channel");
    hist9->Draw();
    hist10->Draw("same");
    hist9->SetStats(kFALSE);
    hist10->SetStats(kFALSE);
    hist9->GetXaxis()->SetTitleSize(0.04);
    hist9->GetYaxis()->SetTitleOffset(1.1);
    hist9->GetYaxis()->SetLabelSize(0.03);
    hist9->GetYaxis()->SetTitleSize(0.04);
    hist9->GetXaxis()->SetTitle("#gamma_{p_{t}} (GeV/c)");
    hist9->GetYaxis()->SetTitle("A_{c_{ele}}");
    TLegend *legend5 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend5->AddEntry(hist9, "TTGamma", "E1P");
    legend5->AddEntry(hist10, "TTbar", "E1P");
    legend5->Draw("same");
    c5->Update();
    c5->SaveAs("/eos/user/s/ssnehshu/temp/Ac_vs_photon_pt_compare_electron_Channel.png");

    TCanvas *c6 = new TCanvas("c6", "Charge Asymmetry vs  Photon transverse momentum Muon Channel");
    hist11->Draw();
    hist12->Draw("same");
    hist11->SetStats(kFALSE);
    hist12->SetStats(kFALSE);
    hist11->GetXaxis()->SetTitleSize(0.04);
    hist11->GetYaxis()->SetTitleOffset(1.1);
    hist11->GetYaxis()->SetLabelSize(0.03);
    hist11->GetYaxis()->SetTitleSize(0.04);
    hist11->GetXaxis()->SetTitle("#gamma_{p_{t}} (GeV/c)");
    hist11->GetYaxis()->SetTitle("A_{c_{mu}}");
    TLegend *legend6 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend6->AddEntry(hist11, "TTGamma", "E1P");
    legend6->AddEntry(hist12, "TTbar", "E1P");
    legend6->Draw("same");
    c6->Update();
    c6->SaveAs("/eos/user/s/ssnehshu/temp/Ac_vs_photon_pt_compare_Muon_Channel.png");

    // Close the file
    fl1->Close();
    fl2->Close();
}
