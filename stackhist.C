


void stackhist ()
{

   TFile *fl1 = new TFile ("TTGamma.root","read");
   TFile *fl2 = new TFile ("TTbar.root","read");
   TFile *fl3 = new TFile ("Ele.root","read");
   TFile *fl4 = new TFile ("Mu.root","read");
  
   TH1F *hist1a = (TH1F*)fl1->Get("h9ele");
   TH1F *hist1b = (TH1F*)fl1->Get("h9mu");
   TH1F *hist1c = (TH1F*)fl1->Get("h10ele");
   TH1F *hist1d = (TH1F*)fl1->Get("h10mu");
   TH1F *hist2a = (TH1F*)fl2->Get("h9ele");
   TH1F *hist2b = (TH1F*)fl2->Get("h9mu");
   TH1F *hist2c = (TH1F*)fl2->Get("h10ele");
   TH1F *hist2d = (TH1F*)fl2->Get("h10mu");
   TH1F *hist3a = (TH1F*)fl3->Get("h9ele");
   TH1F *hist3b = (TH1F*)fl3->Get("h9mu");
   TH1F *hist3c = (TH1F*)fl3->Get("h10ele");
   TH1F *hist3d = (TH1F*)fl3->Get("h10mu");
   TH1F *hist4a = (TH1F*)fl4->Get("h9ele");
   TH1F *hist4b = (TH1F*)fl4->Get("h9mu");
   TH1F *hist4c = (TH1F*)fl4->Get("h10ele");
   TH1F *hist4d = (TH1F*)fl4->Get("h10mu");
   
   THStack *hs1 = new THStack("hs1", "Top Electron channel");
   hist1a->SetMarkerStyle(kStar); 
   hist1a->SetMarkerColor(kBlack);
   hs1->Add(hist1a);
   hist2a->SetMarkerColor(kRed);
   hist2a->SetMarkerStyle(kStar); 
   hs1->Add(hist2a);

   THStack *hs2 = new THStack("hs2", "Top Muon channel");
   hist1b->SetMarkerStyle(kStar); 
   hist1b->SetMarkerColor(kBlue);
   hs2->Add(hist1b);
   hist2b->SetMarkerStyle(kStar); 
   hist2b->SetMarkerColor(kGreen);
   hs2->Add(hist2b);
 
   THStack *hs3 = new THStack("hs3", "Antitop Electron channel");
   hist1c->SetMarkerStyle(kStar); 
   hist1c->SetMarkerColor(kMagenta);
   hs3->Add(hist1c);
   hist2c->SetMarkerStyle(kStar); 
   hist2c->SetMarkerColor(kCyan);
   hs3->Add(hist2c);
 
   THStack *hs4 = new THStack("hs4", "Antitop Muon channel");
   hist1d->SetMarkerStyle(kStar); 
   hist1d->SetMarkerColor(kViolet);
   hs4->Add(hist1d);
   hist2d->SetMarkerStyle(kStar); 
   hist2d->SetMarkerColor(kAzure);
   hs4->Add(hist2d);

   TCanvas *c1 = new TCanvas("c1", "Top Electron channel");
   hs1->Draw("hist,nostack");
   hist3a->SetMarkerStyle(kCircle);
   hist3a->Draw("same,p");
//   hist4a->SetMarkerStyle(kCircle);
//   hist4a->Draw("same");
   TLegend *legend1 = new TLegend(0.7,0.7,0.9,0.9);
legend1->AddEntry(hist1a, "TTGamma", "f");
legend1->AddEntry(hist2a, "TTbar", "f");
legend1->Draw();
  
   TCanvas *c2 = new TCanvas("c2", "Top Muon channel");
   hs2->Draw("hist,nostack");
//   hist3b->SetMarkerStyle(kCircle);
//  hist3b->Draw("same");
   hist4b->SetMarkerStyle(kCircle);
   hist4b->Draw("hist,same,p");
   TLegend *legend2 = new TLegend(0.7,0.7,0.9,0.9);
   legend2->AddEntry(hist1b, "TTGamma", "f");
   legend2->AddEntry(hist2b, "TTbar", "f");
   legend2->Draw();

   TCanvas *c3 = new TCanvas("c3", "Antitop Electron channel");
   hs3->Draw("hist,nostack");
   hist3c->SetMarkerStyle(kCircle);
   hist3c->SetLineStyle(1);
   hist3c->SetLineWidth(2);
   hist3c->Draw("hist,same,p");
//   hist4c->SetMarkerStyle(kCircle);
//   hist4c->Draw("same");
    TLegend *legend3 = new TLegend(0.7,0.7,0.9,0.9);
   legend3->AddEntry(hist1c, "TTGamma", "f");
   legend3->AddEntry(hist2c, "TTbar", "f");
   legend3->Draw();
  
   TCanvas *c4 = new TCanvas("c4", "Antitop Muon channel");
   hs4->Draw("hist,nostack");
//   hist3d->SetMarkerStyle(kCircle);
//   hist3d->Draw("same");
   hist4d->SetMarkerStyle(kCircle);
   hist4d->Draw("hist,same,p");
    TLegend *legend4 = new TLegend(0.7,0.7,0.9,0.9);
   legend4->AddEntry(hist1d, "TTGamma", "f");
   legend4->AddEntry(hist2d, "TTbar", "f");
   legend4->Draw();


  
   c1->Update();
   c1->SaveAs("/eos/user/s/ssnehshu/plots/Topele1.png");
   c2->Update();
   c2->SaveAs("/eos/user/s/ssnehshu/plots/Topmu1.png");
   c3->Update();
   c3->SaveAs("/eos/user/s/ssnehshu/plots/AntiTopele1.png");
   c4->Update();
   c4->SaveAs("/eos/user/s/ssnehshu/plots/AntiTopmu1.png");
 


  /* TCanvas *d1 = new TCanvas("d1", "Summed Ele Data");
   hist3a->Draw();
   hist3b->Draw("same");
   hist3c->Draw("same");
   hist3d->Draw("same");
 
   d1->Update();
   d1->SaveAs("/eos/user/s/ssnehshu/plots/EleDataSum.png");
*/
   fl1->Close();
   fl2->Close();
   fl3->Close();
   fl4->Close();
}
