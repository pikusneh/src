#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TCanvas.h>

void createResponseMatrix() {
  // Open the input ROOT files
  TFile* file_TTGamma = new TFile("output_TTGamma.root", "READ");
  TFile* file_GenLevel = new TFile("output_GenLevel.root", "READ");

  // Get the histograms from the input files
  TH1F* h_TTGamma = (TH1F*)file_TTGamma->Get("h25");
  TH1F* h_GenLevel = (TH1F*)file_GenLevel->Get("h7");

  // Create the response matrix histogram
  TH2F* h_response = new TH2F("h_response", "Response Matrix", h_TTGamma->GetNbinsX(), h_TTGamma->GetXaxis()->GetXmin(), h_TTGamma->GetXaxis()->GetXmax(), h_GenLevel->GetNbinsX(), h_GenLevel->GetXaxis()->GetXmin(), h_GenLevel->GetXaxis()->GetXmax());

  // Fill the response matrix
  for (int i = 1; i <= h_GenLevel->GetNbinsX(); i++) {
    for (int j = 1; j <= h_TTGamma->GetNbinsX(); j++) {
      double x = h_GenLevel->GetBinContent(i);
      double y = h_TTGamma->GetBinContent(j);
      h_response->SetBinContent(j, i, y);
    }
  }


  // Save the response matrix histogram to a new ROOT file
  TFile* output_file = new TFile("response_matrix.root", "RECREATE");
  h_response->Write();
  output_file->Close();

  // Clean up
  delete file_TTGamma;
  delete file_GenLevel;
}

void responseMatrix() {
  createResponseMatrix();
}

int main() {
  responseMatrix();
  return 0;
}
