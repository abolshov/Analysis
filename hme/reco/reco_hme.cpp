#include <iostream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <chrono>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TRandom3.h"

#include "HistManager.hpp"
#include "EstimatorTools.hpp"
#include "CombTools.hpp"

static constexpr int N_RECO_JETS = 12;

int main()
{
    auto start = std::chrono::system_clock::now();

    auto file_pdf = std::make_unique<TFile>("pdf.root", "READ");
    auto pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf")));

    TFile *myFile = TFile::Open("nano_0.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    TRandom3 rg;
    rg.SetSeed(42);

    // reco variables
    Float_t         lep1_pt;
    Float_t         lep1_eta;
    Float_t         lep1_phi;
    Float_t         lep1_mass;

    Float_t         PuppiMET_pt;
    Float_t         PuppiMET_phi;

    Int_t           ncentralJet;
    Float_t         centralJet_pt[N_RECO_JETS];   
    Float_t         centralJet_eta[N_RECO_JETS];   
    Float_t         centralJet_phi[N_RECO_JETS];   
    Float_t         centralJet_mass[N_RECO_JETS]; 
    Float_t         centralJet_btagPNetB[N_RECO_JETS]; 

    myTree->SetBranchAddress("lep1_pt", &lep1_pt);
    myTree->SetBranchAddress("lep1_eta", &lep1_eta);
    myTree->SetBranchAddress("lep1_phi", &lep1_phi);
    myTree->SetBranchAddress("lep1_mass", &lep1_mass);

    myTree->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt);
    myTree->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi);

    myTree->SetBranchAddress("ncentralJet", &ncentralJet);
    myTree->SetBranchAddress("centralJet_pt", centralJet_pt);
    myTree->SetBranchAddress("centralJet_eta", centralJet_eta);
    myTree->SetBranchAddress("centralJet_phi", centralJet_phi);
    myTree->SetBranchAddress("centralJet_mass", centralJet_mass);
    myTree->SetBranchAddress("centralJet_btagPNetB", centralJet_btagPNetB);

    // gen variables
    // V1: leptonic 
    // prod1: lep
    // prod2: neutrino
    Double_t        genV1prod1_pt;
    Double_t        genV1prod1_vis_pt;
    Double_t        genV1prod1_eta;
    Double_t        genV1prod1_vis_eta;
    Double_t        genV1prod1_phi;
    Double_t        genV1prod1_vis_phi;
    Double_t        genV1prod1_mass;
    Double_t        genV1prod1_vis_mass;

    Double_t        genV1prod2_pt;
    Double_t        genV1prod2_vis_pt;
    Double_t        genV1prod2_eta;
    Double_t        genV1prod2_vis_eta;
    Double_t        genV1prod2_phi;
    Double_t        genV1prod2_vis_phi;
    Double_t        genV1prod2_mass;
    Double_t        genV1prod2_vis_mass;

    // V2: hadronic
    Double_t        genV2prod1_pt;
    Double_t        genV2prod1_vis_pt;
    Double_t        genV2prod1_eta;
    Double_t        genV2prod1_vis_eta;
    Double_t        genV2prod1_phi;
    Double_t        genV2prod1_vis_phi;
    Double_t        genV2prod1_mass;
    Double_t        genV2prod1_vis_mass;

    Double_t        genV2prod2_pt;
    Double_t        genV2prod2_vis_pt;
    Double_t        genV2prod2_eta;
    Double_t        genV2prod2_vis_eta;
    Double_t        genV2prod2_phi;
    Double_t        genV2prod2_vis_phi;
    Double_t        genV2prod2_mass;
    Double_t        genV2prod2_vis_mass;

    myTree->SetBranchAddress("genV1prod1_pt", &genV1prod1_pt);
    myTree->SetBranchAddress("genV1prod1_vis_pt", &genV1prod1_vis_pt);
    myTree->SetBranchAddress("genV1prod1_eta", &genV1prod1_eta);
    myTree->SetBranchAddress("genV1prod1_vis_eta", &genV1prod1_vis_eta);
    myTree->SetBranchAddress("genV1prod1_phi", &genV1prod1_phi);
    myTree->SetBranchAddress("genV1prod1_vis_phi", &genV1prod1_vis_phi);
    myTree->SetBranchAddress("genV1prod1_mass", &genV1prod1_mass);
    myTree->SetBranchAddress("genV1prod1_vis_mass", &genV1prod1_vis_mass);

    myTree->SetBranchAddress("genV1prod2_pt", &genV1prod2_pt);
    myTree->SetBranchAddress("genV1prod2_vis_pt", &genV1prod2_vis_pt);
    myTree->SetBranchAddress("genV1prod2_eta", &genV1prod2_eta);
    myTree->SetBranchAddress("genV1prod2_vis_eta", &genV1prod2_vis_eta);
    myTree->SetBranchAddress("genV1prod2_phi", &genV1prod2_phi);
    myTree->SetBranchAddress("genV1prod2_vis_phi", &genV1prod2_vis_phi);
    myTree->SetBranchAddress("genV1prod2_mass", &genV1prod2_mass);
    myTree->SetBranchAddress("genV1prod2_vis_mass", &genV1prod2_vis_mass);

    myTree->SetBranchAddress("genV2prod1_pt", &genV2prod1_pt);
    myTree->SetBranchAddress("genV2prod1_vis_pt", &genV2prod1_vis_pt);
    myTree->SetBranchAddress("genV2prod1_eta", &genV2prod1_eta);
    myTree->SetBranchAddress("genV2prod1_vis_eta", &genV2prod1_vis_eta);
    myTree->SetBranchAddress("genV2prod1_phi", &genV2prod1_phi);
    myTree->SetBranchAddress("genV2prod1_vis_phi", &genV2prod1_vis_phi);
    myTree->SetBranchAddress("genV2prod1_mass", &genV2prod1_mass);
    myTree->SetBranchAddress("genV2prod1_vis_mass", &genV2prod1_vis_mass);

    myTree->SetBranchAddress("genV2prod2_pt", &genV2prod2_pt);
    myTree->SetBranchAddress("genV2prod2_vis_pt", &genV2prod2_vis_pt);
    myTree->SetBranchAddress("genV2prod2_eta", &genV2prod2_eta);
    myTree->SetBranchAddress("genV2prod2_vis_eta", &genV2prod2_vis_eta);
    myTree->SetBranchAddress("genV2prod2_phi", &genV2prod2_phi);
    myTree->SetBranchAddress("genV2prod2_vis_phi", &genV2prod2_vis_phi);
    myTree->SetBranchAddress("genV2prod2_mass", &genV2prod2_mass);
    myTree->SetBranchAddress("genV2prod2_vis_mass", &genV2prod2_vis_mass);

    HistManager hm;

    std::string hme_mass("hme_mass");
    hm.Add(hme_mass, "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    int nEvents = myTree->GetEntries();
    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        if (ncentralJet < 4)
        {
            continue;
        }

        TLorentzVector reco_lep_p4, reco_bj1_p4, reco_bj2_p4, reco_met_p4;
        reco_lep_p4.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_mass);
        reco_met_p4.SetPtEtaPhiM(PuppiMET_pt, 0.0, PuppiMET_phi, 0.0);

        // establish jet with greater pt for computing corrections
        // jets are sorted by btag so first two are b jets
        if (centralJet_pt[0] > centralJet_pt[1])
        {
            reco_bj1_p4.SetPtEtaPhiM(centralJet_pt[0], centralJet_eta[0], centralJet_phi[0], centralJet_mass[0]);
            reco_bj2_p4.SetPtEtaPhiM(centralJet_pt[1], centralJet_eta[1], centralJet_phi[1], centralJet_mass[1]);
        }
        else 
        {
            reco_bj2_p4.SetPtEtaPhiM(centralJet_pt[0], centralJet_eta[0], centralJet_phi[0], centralJet_mass[0]);
            reco_bj1_p4.SetPtEtaPhiM(centralJet_pt[1], centralJet_eta[1], centralJet_phi[1], centralJet_mass[1]);
        }

        std::vector<TLorentzVector> light_jets;
        for (int j = 2; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[i], centralJet_eta[i], centralJet_phi[i], centralJet_mass[i]);
            light_jets.push_back(jet);
        }

        auto best_onshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(80.0 - (p1 + p2).M()); });
        auto best_offshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(40.0 - (p1 + p2).M()); });

        auto [idx1, idx2] = best_onshell_pair;
        std::vector<TLorentzVector> input_onshell = {reco_bj1_p4, reco_bj2_p4, light_jets[idx1], light_jets[idx2], reco_lep_p4, reco_met_p4};

        std::tie(idx1, idx2) = best_offshell_pair;
        std::vector<TLorentzVector> input_offshell = {reco_bj1_p4, reco_bj2_p4, light_jets[idx1], light_jets[idx2], reco_lep_p4, reco_met_p4};

        auto hme = EstimateMass(input_onshell, pdf, rg, i);
        if (hme)
        {
            auto [mass, eff] = hme.value();
            hm.FillWeighted(hme_mass, mass, 0.5);
        }

        hme = EstimateMass(input_offshell, pdf, rg, i);
        if (hme)
        {
            auto [mass, eff] = hme.value();
            hm.FillWeighted(hme_mass, mass, 0.5);
        }
    }

    hm.Draw();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    return 0;
}