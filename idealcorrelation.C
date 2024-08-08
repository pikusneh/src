#include "TROOT.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TPad.h>

void idealcorrelation()
{

    TFile *fl = new TFile("output_GenLevel.root", "read");
    TH1F *hist1 = (TH1F *)fl->Get("Rapidity_diff_Gen");
    TH1F *hist2 = (TH1F *)hist1->Clone();

    TH2F *histresponse = new TH2F("histresponse", "histresponse", 20, -3, 3, 20, -3, 3);
    // Create a new histogram to store the bin-by-bin difference
    TH1F *histDiff = new TH1F("histDiff", "Bin-by-Bin Difference", hist1->GetNbinsX(), hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax());

    // Calculate the bin-by-bin difference
    for (int iBin = 1; iBin <= hist2->GetNbinsX(); ++iBin)
    {
        double diff = hist2->GetBinContent(iBin) - hist1->GetBinContent(iBin);
        histDiff->SetBinContent(iBin, diff);
    }
    Int_t nbins = 20;
    // Loop over Gen level histogram
    for (Int_t k = 1; k <= nbins; k++)
    {
        // Loop over the same Gen level histogram
        for (Int_t j = 1; j <= nbins; j++)
        {
            Float_t genContent1 = hist1->GetBinContent(k);
            Float_t genContent2 = hist2->GetBinContent(j);
            histresponse->SetBinContent(k, j, genContent1 * genContent2);
        }
    }
    
    TCanvas *c1 = new TCanvas("c1", "Rapidity_Diff_of_reco_gen", 800, 800);

    // Create upper pad for stacked histograms
    TPad *upperPad = new TPad("upperPad", "Upper Pad", 0.0, 0.4, 1.0, 1.0);
    upperPad->SetBottomMargin(0.0);
    upperPad->Draw();
    upperPad->cd();

    THStack *hs1 = new THStack("hs1", "Rapidity_Diff_of_reco_gen");
    hs1->Add(hist1);
    hs1->Add(hist2);
    hs1->Draw("nostack");
    hist1->SetLineColor(kBlue);
    hist2->SetLineColor(kRed);
    hist1->Draw();
    hist2->Draw("HIST SAME");

    // Create and customize a legend for the upper pad
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend1->AddEntry(hist1, "Rapidity_diff_reco_ele", "f"); // "f" for filled histograms
    legend1->AddEntry(hist2, "Rapidity_diff_gen", "f");      // "f" for filled histograms
    legend1->SetBorderSize(0);                               // Remove border from legend
    legend1->SetFillColor(0);                                // Set legend background color to white
    legend1->Draw();

    // Create lower pad for the bin-by-bin difference
    c1->cd(); // Go back to the main canvas
    TPad *lowerPad = new TPad("lowerPad", "Lower Pad", 0.0, 0.0, 1.0, 0.3);
    lowerPad->SetTopMargin(0.0);
    // lowerPad->SetGridx();
    // lowerPad->SetGridy();
    lowerPad->Draw();
    lowerPad->cd();

    // Draw the bin-by-bin difference as a separate histogram in the lower pad
    histDiff->SetLineColor(kGreen);
    histDiff->Draw("HIST");

    // Customize axis labels and titles for the lower pad
    histDiff->GetXaxis()->SetTitle("Rapidity_diff");
    histDiff->GetYaxis()->SetTitle("Difference");
    histDiff->GetXaxis()->SetTitleSize(0.04);
    histDiff->GetYaxis()->SetTitleOffset(0.6);
    histDiff->GetYaxis()->SetLabelSize(0.03);
    histDiff->GetYaxis()->SetTitleSize(0.04);

    TFile *file1 = new TFile("output_reco_gen_rapdiff_ideal.root", "RECREATE");
    c1->Write();
    c1->Update();
    c1->SaveAs("/eos/user/s/ssnehshu/reco_gen_rapidity_diff_ideal.png");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas *c2 = new TCanvas("c2", "Correlation_Rapidity_Diff_of_reco_gen", 1000, 1000);
    histresponse->Draw("COLZTEXT");
    histresponse->GetYaxis()->SetTitle("Gen 2 Level");
    histresponse->GetXaxis()->SetTitle("Gen 1 Level");

    histresponse->Write();
    TFile *outputFile = new TFile("output_correlation_raid_diff_gen_reco_ele_ideal.root", "RECREATE");
    c2->Write();
    c2->Update();
    c2->SaveAs("/eos/user/s/ssnehshu/correlation_reco_gen_rapidity_diff_ideal.png");

    fl->Close();
    file1->Close();
    outputFile->Close();
}
