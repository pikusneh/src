#include "TH2D.h"

const Double_t cutdummy = -99999.0;

//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear(Double_t xt)
{
    Double_t xeff = 0.3 + (1.0 - 0.3) / 20 * (xt + 10.0); // efficiency
    Double_t x = gRandom->Rndm();
    if (x > xeff)
        return cutdummy;
    Double_t xsmear = gRandom->Gaus(-2.5, 0.2); // bias and smear
    return xt + xsmear;
}

void unfoldingtest()
{
    TFile *file1 = new TFile("output_GenLevel.root", "READ");
    TH1F *hist1 = (TH1F *)file1->Get("h7");
    TFile *file2 = new TFile("output_TTGamma.root", "READ");
    TH1F *hist2 = (TH1F *)file2->Get("h25");

    // Create the response matrix
    TH2D *responseMatrix = new TH2D("responseMatrix", "Response Matrix", hist2->GetNbinsX(), hist2->GetXaxis()->GetXmin(), hist2->GetXaxis()->GetXmax(), hist1->GetNbinsX(), hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax());

    // Open the subtracted_histogram ROOT file
    TFile *subtracted_file = new TFile("subtracted_histogram.root", "READ");

    // Get the histogram from the subtracted_histogram file
    TH1F *h11 = (TH1F *)subtracted_file->Get("h11");
    // Fill the response matrix
    for (Int_t i = 1; i <= hist2->GetNbinsX(); i++)
    {
        Double_t xt = hist2->GetBinContent(i); // Charge asymmetry value of h14
        Double_t x = smear(xt);                // Assuming you have a smear function
        if (x >= hist1->GetXaxis()->GetXmin() && x <= hist1->GetXaxis()->GetXmax())
        {
            Int_t bin = hist1->FindBin(x);
            responseMatrix->Fill(xt, hist1->GetBinContent(bin)); // Charge asymmetry value of h4
        }
    }

    // Perform unfolding using the response matrix
    // ...

    // Perform unfolding: y = Ax
    TH1F *unfolded_hist = new TH1F("unfolded_hist", "Unfolded Histogram", nbins_gen, h11->GetXaxis()->GetXmin(), h11->GetXaxis()->GetXmax());
    for (int i = 0; i < nbins_gen; i++)
    {
        double x = 0.0;
        for (int j = 0; j < nbins_meas; j++)
        {
            double y = h11->GetBinContent(j + 1);
            double A = (*responseMatrix)(j, i);
            x += A * y;
        }
        unfolded_hist->SetBinContent(i + 1, x);
    }

    // Save the unfolded histogram to a new ROOT file
    TFile *output_file = new TFile("unfolding.root", "RECREATE");
    unfolded_hist->Write();
    output_file->Close();

    // // Save the response matrix to a ROOT file
    // TFile *outputFile = new TFile("response_matrix.root", "RECREATE");
    // responseMatrix->Write("responseMatrix");
    // outputFile->Close();

    // Save the histogram and the response matrix to a ROOT file
    // TFile *outputFile = new TFile("unfolding_output.root", "RECREATE");

    // responseMatrix->Write("responseMatrix");
}
