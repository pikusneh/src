#include "VarTree.h"

Event::Event(TChain *trev) : tr(trev)
{
    SetBranchStatuses();
    SetBranchAddresses(); // Link branch addresses to eventData structure
}

Event::~Event()
{
}

void Event::SetBranchStatuses()
{
    tr->SetBranchStatus("*", 0);
    //  chain->SetBranchStatus("elePt", 1);
    tr->SetBranchStatus("elePt", 1);
    tr->SetBranchStatus("muPt", 1);
    tr->SetBranchStatus("jetPt", 1);
    tr->SetBranchStatus("phoEt", 1);
    tr->SetBranchStatus("phoEta", 1);
    tr->SetBranchStatus("phoPhi", 1);
    tr->SetBranchStatus("nEle", 1);
    tr->SetBranchStatus("nMu", 1);
    tr->SetBranchStatus("nGenJet", 1);
    tr->SetBranchStatus("nGenPart", 1);
    tr->SetBranchStatus("genPDGID", 1);
    tr->SetBranchStatus("genMomIdx", 1);
    tr->SetBranchStatus("genMass", 1);
    tr->SetBranchStatus("genPhi", 1);
    tr->SetBranchStatus("genEta", 1);
    tr->SetBranchStatus("genPt", 1);
    tr->SetBranchStatus("evtWeight", 1);
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
}

void Event::SetBranchAddresses()
{
    // Use only branches you need to reduce memory consumption and increase processing speed
    tr->SetBranchAddress("elePt", &elePt);
    tr->SetBranchAddress("muPt", &muPt);
    tr->SetBranchAddress("jetPt", &jetPt);
    tr->SetBranchAddress("phoEt", &phoEt);
    tr->SetBranchAddress("phoEta", &phoEta);
    tr->SetBranchAddress("phoPhi", &phoPhi);
    tr->SetBranchAddress("nEle", &nEle);
    tr->SetBranchAddress("nMu", &nMu);
    tr->SetBranchAddress("nGenJet", &nGenJet);
    tr->SetBranchAddress("nGenPart", &nGenPart);
    tr->SetBranchAddress("genPDGID", &genPDGID);
    tr->SetBranchAddress("genMomIdx", &genMomIdx);
    tr->SetBranchAddress("genMass", &genMass);
    tr->SetBranchAddress("genPhi", &genPhi);
    tr->SetBranchAddress("genEta", &genEta);
    tr->SetBranchAddress("genPt", &genPt);
    tr->SetBranchAddress("evtWeight", &evtWeight);
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
}

void Event::LoopOverEvents();
{
    int numevents = tr->GetEntries();

    for (int j = 0; j < numevents; j++)
    {
        tr->GetEntry(j);
        int topindex = -1;
        int antitopindex = -1;
        int photonindex = -1;
        int electronindex = -1;
        bool passPresel_Ele_gen = false;

        weight_gen = evtWeight;
        weight = evtWeight * PUweight * muEffWeight * eleEffWeight * btagWeight_1a * prefireSF;
        TLorentzVector Top_gen;
        TLorentzVector Anti_Top_gen;
        TLorentzVector Top;
        TLorentzVector AntiTop;
        TLorentzVector Photon_gen;
        TLorentzVector Ele_gen;
        int totalPhotons = 0;
        int totalbjet = 0;

        for (int i = 0; i < nGenPart; i++)
        {
            if (genPDGID->at(i) == 6)
            {
                topindex = i;
            }
            if (genPDGID->at(i) == -6)
            {
                antitopindex = i;
            }
            if (genPDGID->at(i) == 11)
            { 
                electronindex = i;  
            }
            Ele_gen.SetPtEtaPhiM((*genPt)[electronindex], (*genEta)[electronindex], (*genPhi)[electronindex], (*genMass)[electronindex]);
            float EleRapidity_gen = Ele_gen.Rapidity();
            float YEle_gen = TMath::Abs(EleRapidity_gen);
            float Ele_Pt = Ele_gen.Pt();

            if (Ele_Pt >= 35.0 && YEle_gen <= 2.4)
            {
                passPresel_Ele_gen = true;
            }

            if (genPDGID->at(i) == 22)
            {
                photonindex = i;
                totalPhotons++;
            }
            // for (int j = 0; j < nGenJet; j++)
            // {
            if (genPDGID->at(i) == 5 || genPDGID->at(i) == -5) // && genPt->at(i) >= 30.0 && fabs(genEta->at(i)) <= 2.4)
            {
                totalbjet++;
                // }
            }
        }
    }
}
