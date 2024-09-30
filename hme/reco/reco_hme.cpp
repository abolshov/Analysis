#include <iostream>
#include <memory>
#include <chrono>
#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TCanvas.h"

#include "HistManager.hpp"
#include "EstimatorTools.hpp"
#include "CombTools.hpp"

double InterquantileRange(std::unique_ptr<TH1F> const& h)
{
    int const nq = 50;
    double xq[nq];  // position where to compute the quantiles in [0,1]
    double yq[nq];  // array to contain the quantiles
    for (int i = 0; i < nq; ++i) 
    {
        xq[i] = static_cast<double>(i+1)/nq;
    }
    h->GetQuantiles(nq, yq, xq);
    return yq[41] - yq[7];
}

enum FailureType { Resc, Nu, Both };

static constexpr int N_RECO_JETS = 12;

int main()
{
    auto start = std::chrono::system_clock::now();

    auto file_pdf = std::make_unique<TFile>("pdf.root", "READ");
    auto pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_b1")));

    TFile *myFile = TFile::Open("nano_0.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    TRandom3 rg;
    rg.SetSeed(42);

    std::ofstream output("hme_output.txt");

    ULong64_t event;
    myTree->SetBranchAddress("event", &event);

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
    // H->VV
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

    // H->bb
    Double_t        genb1_pt;
    Double_t        genb1_vis_pt;
    Double_t        genb1_eta;
    Double_t        genb1_vis_eta;
    Double_t        genb1_phi;
    Double_t        genb1_vis_phi;
    Double_t        genb1_mass;
    Double_t        genb1_vis_mass;
    Double_t        genb2_pt;
    Double_t        genb2_vis_pt;
    Double_t        genb2_eta;
    Double_t        genb2_vis_eta;
    Double_t        genb2_phi;
    Double_t        genb2_vis_phi;
    Double_t        genb2_mass;
    Double_t        genb2_vis_mass;

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

    myTree->SetBranchAddress("genb1_pt", &genb1_pt);
    myTree->SetBranchAddress("genb1_vis_pt", &genb1_vis_pt);
    myTree->SetBranchAddress("genb1_eta", &genb1_eta);
    myTree->SetBranchAddress("genb1_vis_eta", &genb1_vis_eta);
    myTree->SetBranchAddress("genb1_phi", &genb1_phi);
    myTree->SetBranchAddress("genb1_vis_phi", &genb1_vis_phi);
    myTree->SetBranchAddress("genb1_mass", &genb1_mass);
    myTree->SetBranchAddress("genb1_vis_mass", &genb1_vis_mass);
    myTree->SetBranchAddress("genb2_pt", &genb2_pt);
    myTree->SetBranchAddress("genb2_vis_pt", &genb2_vis_pt);
    myTree->SetBranchAddress("genb2_eta", &genb2_eta);
    myTree->SetBranchAddress("genb2_vis_eta", &genb2_vis_eta);
    myTree->SetBranchAddress("genb2_phi", &genb2_phi);
    myTree->SetBranchAddress("genb2_vis_phi", &genb2_vis_phi);
    myTree->SetBranchAddress("genb2_mass", &genb2_mass);
    myTree->SetBranchAddress("genb2_vis_mass", &genb2_vis_mass);

    HistManager hm;

    // std::string hme_mass("hme_mass");
    // hm.Add(hme_mass, "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 2000}, 100);

    std::string W_mass_onshell("W_mass_onshell");
    hm.Add(W_mass_onshell, "Mass of picked W for onshell pair", {"W mass, [GeV]", "Count"}, {0, 250}, 100);

    std::string W_mass_offshell("W_mass_offshell");
    hm.Add(W_mass_offshell, "Mass of picked W for offshell pair", {"W mass, [GeV]", "Count"}, {0, 250}, 100);

    std::string err_onshell("err_onshell");
    hm.Add(err_onshell, "Errors of HME with onshell pair", {err_names[Error::RescFact], err_names[Error::NuEqn]}, {0, 100}, {0, 100}, {10, 10});

    std::string err_offshell("err_offshell");
    hm.Add(err_offshell, "Errors of HME with offshell pair", {err_names[Error::RescFact], err_names[Error::NuEqn]}, {0, 100}, {0, 100}, {10, 10});

    std::vector<char const*> labels = {"only resc", "only nu", "resc + nu"};

    std::vector<int> onshell_counts(labels.size(), 0);
    std::vector<int> offshell_counts(labels.size(), 0);

    int hme_events = 0;
    int hme_worked = 0;

    int nEvents = myTree->GetEntries();
    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        if (event % 2 != 1)
        {
            continue;
        }

        if (ncentralJet < 4)
        {
            continue;
        }

        TLorentzVector genb1_p4, genb2_p4, genq1_p4, genq2_p4; // quarks
        genb1_p4.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        genb2_p4.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);
        genq1_p4.SetPtEtaPhiM(genV2prod1_pt, genV2prod1_eta, genV2prod1_phi, genV2prod1_mass);
        genq2_p4.SetPtEtaPhiM(genV2prod2_pt, genV2prod2_eta, genV2prod2_phi, genV2prod2_mass);

        if (genb1_p4.DeltaR(genb2_p4) < 0.4)
        {
            continue;
        }

        if (genq1_p4.DeltaR(genq2_p4) < 0.4)
        {
            continue;
        }

        TLorentzVector gen_bj1_p4, gen_bj2_p4, gen_lj1_p4, gen_lj2_p4; // gen jets matched to quarks
        gen_bj1_p4.SetPtEtaPhiM(genb1_vis_pt, genb1_vis_eta, genb1_vis_phi, genb1_vis_mass);
        gen_bj2_p4.SetPtEtaPhiM(genb2_vis_pt, genb2_vis_eta, genb2_vis_phi, genb2_vis_mass);
        gen_lj1_p4.SetPtEtaPhiM(genV2prod1_vis_pt, genV2prod1_vis_eta, genV2prod1_vis_phi, genV2prod1_vis_mass);
        gen_lj2_p4.SetPtEtaPhiM(genV2prod2_vis_pt, genV2prod2_vis_eta, genV2prod2_vis_phi, genV2prod2_vis_mass);

        TLorentzVector gen_lep_p4;
        gen_lep_p4.SetPtEtaPhiM(genV1prod1_vis_pt, genV1prod1_vis_eta, genV1prod1_vis_phi, genV1prod1_vis_mass);

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
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);
            light_jets.push_back(jet);
        }

        auto best_onshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(80.0 - (p1 + p2).M()); });
        auto best_offshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(40.0 - (p1 + p2).M()); });

        auto [i1, i2] = best_onshell_pair;
        std::vector<TLorentzVector> input_onshell = {reco_bj1_p4, reco_bj2_p4, light_jets[i1], light_jets[i2], reco_lep_p4, reco_met_p4};

        auto [j1, j2] = best_offshell_pair;
        std::vector<TLorentzVector> input_offshell = {reco_bj1_p4, reco_bj2_p4, light_jets[j1], light_jets[j2], reco_lep_p4, reco_met_p4};

        ++hme_events;

        std::vector<int> err_cnt(err_names.size(), 0);
        [[maybe_unused]] auto hme_onshell = EstimateMass(input_onshell, pdf, rg, i, err_cnt);
        if (!hme_onshell)
        {
            hm.Fill(err_onshell, err_cnt[Error::RescFact], err_cnt[Error::NuEqn]);

            bool both = std::all_of(err_cnt.begin(), err_cnt.end(), [](int c){ return c > 0; });
            bool resc = err_cnt[Error::NuEqn] == 0;
            bool nu = err_cnt[Error::RescFact] == 0;

            onshell_counts[FailureType::Resc] += resc; 
            onshell_counts[FailureType::Nu] += nu; 
            onshell_counts[FailureType::Both] += both; 

            hm.Fill(W_mass_onshell, (light_jets[i1] + light_jets[i2]).M());
        }

        err_cnt = std::vector<int>(err_names.size(), 0);
        [[maybe_unused]] auto hme_offshell = EstimateMass(input_offshell, pdf, rg, i, err_cnt);
        if (!hme_offshell)
        {
            hm.Fill(err_offshell, err_cnt[Error::RescFact], err_cnt[Error::NuEqn]);
            
            bool both = std::all_of(err_cnt.begin(), err_cnt.end(), [](int c){ return c > 0; });
            bool resc = err_cnt[Error::NuEqn] == 0;
            bool nu = err_cnt[Error::RescFact] == 0;

            offshell_counts[FailureType::Resc] += resc; 
            offshell_counts[FailureType::Nu] += nu; 
            offshell_counts[FailureType::Both] += both; 

            hm.Fill(W_mass_offshell, (light_jets[j1] + light_jets[j2]).M());
        }
    }

    int nx = labels.size();

    auto onshell_summary = std::make_unique<TH1F>("onshell_summary", "Summary of HME failures for onshell pair", nx, 0, nx);
    onshell_summary->SetStats(0);
    onshell_summary->SetFillStyle(3544);
    onshell_summary->SetLineWidth(2);
    onshell_summary->SetFillColorAlpha(kBlue, 0.75);

    auto offshell_summary = std::make_unique<TH1F>("offshell_summary", "Summary of HME failures for offshell pair", nx, 0, nx);
    offshell_summary->SetStats(0);
    offshell_summary->SetFillStyle(3544);
    offshell_summary->SetLineWidth(2);
    offshell_summary->SetFillColorAlpha(kBlue, 0.75);

    auto canvas = std::make_unique<TCanvas>("canvas", "canvas");
    canvas->SetGrid();

    for (int b = 1; b <= nx; ++b)
    {
        onshell_summary->SetBinContent(b, onshell_counts[b-1]);
        onshell_summary->GetXaxis()->SetBinLabel(b, labels[b-1]);

        offshell_summary->SetBinContent(b, offshell_counts[b-1]);
        offshell_summary->GetXaxis()->SetBinLabel(b, labels[b-1]);
    }

    onshell_summary->Draw();
    canvas->SaveAs("histograms/onshell_summary.png");

    offshell_summary->Draw();
    canvas->SaveAs("histograms/offshell_summary.png");

    hm.Draw();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Finished processing, total events = " << nEvents << "\n";
    std::cout << "Events passed to HME = " << hme_events << "\n"; 
    std::cout << "HME successful = " << hme_worked << "\n"; 
    std::cout << "Processing time = " << elapsed.count() << " s\n";

    return 0;
}