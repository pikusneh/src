#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TCanvas.h>

void performUnfolding() {
  // Open the response matrix ROOT file
  TFile* response_file = new TFile("response_matrix.root", "READ");

  // Read the response matrix TMatrixD
  TMatrixD* responseMatrix = (TMatrixD*)response_file->Get("responseMatrix");

  // Open the subtracted_histogram ROOT file
  TFile* subtracted_file = new TFile("subtracted_histogram.root", "READ");

  // Get the histogram from the subtracted_histogram file
  TH1F* h11 = (TH1F*)subtracted_file->Get("h11");

  // Check dimensions compatibility
  int nbins_gen = responseMatrix->GetNrows();
  int nbins_meas = responseMatrix->GetNcols();
  if (nbins_meas != h11->GetNbinsX()) {
    cout << "Error: Incompatible dimensions between response matrix and h11 histogram." << endl;
    return;
  }

  // Perform unfolding: y = Ax
  TH1F* unfolded_hist = new TH1F("unfolded_hist", "Unfolded Histogram", nbins_gen, h11->GetXaxis()->GetXmin(), h11->GetXaxis()->GetXmax());
  for (int i = 0; i < nbins_gen; i++) {
    double x = 0.0;
    for (int j = 0; j < nbins_meas; j++) {
      double y = h11->GetBinContent(j+1);
      double A = (*responseMatrix)(j, i);
      x += A * y;
    }
    unfolded_hist->SetBinContent(i+1, x);
  }

  // Save the unfolded histogram to a new ROOT file
  TFile* output_file = new TFile("unfolding.root", "RECREATE");
  unfolded_hist->Write();
  output_file->Close();

  // Clean up
  delete unfolded_hist;
  delete response_file;
  delete subtracted_file;
  delete output_file;
}

void unfoldingExample() {
  performUnfolding();
}

int main() {
  unfoldingExample();
  return 0;
}
