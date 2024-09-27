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

    std::string hme_mass_onshell("hme_mass_onshell");
    hm.Add(hme_mass_onshell, "HME X->HH mass wiht onshell jets for W->qq", {"X->HH mass, [GeV]", "Count"}, {0, 2000}, 100);

    std::string hme_mass_offshell("hme_mass_offshell");
    hm.Add(hme_mass_offshell, "HME X->HH mass wiht offshell jets for W->qq", {"X->HH mass, [GeV]", "Count"}, {0, 2000}, 100);

    std::string hme_succ_onshell("hme_succ_onshell");
    hm.Add(hme_succ_onshell, "HME success rate with onshell jets for W->qq", {"X->HH mass, [GeV]", "Count"}, {0, 1}, 10);

    std::string hme_succ_offshell("hme_succ_offshell");
    hm.Add(hme_succ_offshell, "HME success rate with offshell jets for W->qq", {"X->HH mass, [GeV]", "Count"}, {0, 1}, 10);

    std::string numLightJets_onshell_fail("numLightJets_onshell_fail");
    hm.Add(numLightJets_onshell_fail, "Number of light jets in events where HME with onshell pair fails", {"Number of light jets", "Count"}, {-0.5, 12.5}, 13);

    std::string numLightJets_offshell_fail("numLightJets_offshell_fail");
    hm.Add(numLightJets_offshell_fail, "Number of light jets in events where HME with offshell pair fails", {"Number of light jets", "Count"}, {-0.5, 12.5}, 13);

    std::string numLightJets_both_fail("numLightJets_both_fail");
    hm.Add(numLightJets_both_fail, "Number of light jets in events where HME fails", {"Number of light jets", "Count"}, {-0.5, 12.5}, 13);

    std::string debug_hist("debug_hist");
    hm.Add(debug_hist, "HME failures summary", {"Number of problems", "Count"}, {-0.5, 3.5}, 4);

    std::string nJets_vs_massDiff("nJets_vs_massDiff");
    hm.Add(nJets_vs_massDiff, "HME mass difference vs number of light jets", {"Number of light jets", "Mass difference, [GeV]"}, {-0.5, 12.5}, {-10.5, 120.5}, {13, 13});

    std::string debug_off("debug_off");
    hm.Add(debug_off, "HME failures with offshell pair summary", {"Number of problems", "Count"}, {-0.5, 3.5}, 4);

    std::string debug_on("debug_on");
    hm.Add(debug_on, "HME failures with onshell pair summary", {"Number of problems", "Count"}, {-0.5, 3.5}, 4);

    std::string same_pair("same_pair");
    hm.Add(same_pair, "Mass difference for same pair of light jets", {"Mass difference, [GeV]", "Count"}, {-1, 100}, 30);

    std::string diff_pairs("diff_pairs");
    hm.Add(diff_pairs, "Mass difference for different pairs of light jets", {"Mass difference, [GeV]", "Count"}, {-1, 100}, 30);

    auto h = std::make_unique<TH1F>("h", "h", 200, 0, 2000);

    int hme_events = 0;
    int hme_worked = 0;

    int has_bad_lep = 0;
    int has_bad_bjets = 0;
    int has_bad_light_jets = 0;

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

        auto [idx1, idx2] = best_onshell_pair;
        std::vector<TLorentzVector> input_onshell = {reco_bj1_p4, reco_bj2_p4, light_jets[idx1], light_jets[idx2], reco_lep_p4, reco_met_p4};

        std::tie(idx1, idx2) = best_offshell_pair;
        std::vector<TLorentzVector> input_offshell = {reco_bj1_p4, reco_bj2_p4, light_jets[idx1], light_jets[idx2], reco_lep_p4, reco_met_p4};

        ++hme_events;

        auto hme_onshell = EstimateMass(input_onshell, pdf, rg, i);
        double m1 = -1.0;
        if (hme_onshell)
        {
            auto [mass, sr] = hme_onshell.value();
            hm.Fill(hme_mass_onshell, mass);
            hm.Fill(hme_succ_onshell, sr);
            // hm.FillWeighted(hme_mass, mass, 0.5);
            m1 = mass;
            output << mass << " ";
        }
        else
        {
            output << -1.0 << " "; 
            hm.Fill(numLightJets_onshell_fail, light_jets.size());
        }

        auto hme_offshell = EstimateMass(input_offshell, pdf, rg, i);
        double m2 = -1.0;
        if (hme_offshell)
        {
            auto [mass, sr] = hme_offshell.value();
            hm.Fill(hme_mass_offshell, mass);
            hm.Fill(hme_succ_offshell, sr);
            // hm.FillWeighted(hme_mass, mass, 0.5);
            m2 = mass;
            output << mass << " ";
        }
        else
        {
            output << -1.0 << " "; 
            hm.Fill(numLightJets_offshell_fail, light_jets.size());
        }

        if (!hme_offshell)
        {
            bool bad_lep = gen_lep_p4.DeltaR(reco_lep_p4) > 0.4;

            bool bad_bjet_1 = genb1_p4.DeltaR(reco_bj1_p4) > 0.4 && genb1_p4.DeltaR(reco_bj2_p4) > 0.4;
            bool bad_bjet_2 = genb2_p4.DeltaR(reco_bj1_p4) > 0.4 && genb2_p4.DeltaR(reco_bj2_p4) > 0.4;
            bool bad_bjets = bad_bjet_1 || bad_bjet_2;

            bool bad_ljet_1 = genq1_p4.DeltaR(light_jets[best_offshell_pair.first]) > 0.4 && genq1_p4.DeltaR(light_jets[best_offshell_pair.second]) > 0.4;
            bool bad_ljet_2 = genq2_p4.DeltaR(light_jets[best_offshell_pair.first]) > 0.4 && genq2_p4.DeltaR(light_jets[best_offshell_pair.second]) > 0.4;
            bool bad_light_jets = bad_ljet_1 || bad_ljet_2;

            int n_problems = bad_lep + bad_bjets + bad_light_jets;
            hm.Fill(debug_off, n_problems);
        }

        if (!hme_onshell)
        {
            bool bad_lep = gen_lep_p4.DeltaR(reco_lep_p4) > 0.4;

            bool bad_bjet_1 = genb1_p4.DeltaR(reco_bj1_p4) > 0.4 && genb1_p4.DeltaR(reco_bj2_p4) > 0.4;
            bool bad_bjet_2 = genb2_p4.DeltaR(reco_bj1_p4) > 0.4 && genb2_p4.DeltaR(reco_bj2_p4) > 0.4;
            bool bad_bjets = bad_bjet_1 || bad_bjet_2;

            bool bad_ljet_1 = genq1_p4.DeltaR(light_jets[best_onshell_pair.first]) > 0.4 && genq1_p4.DeltaR(light_jets[best_onshell_pair.second]) > 0.4;
            bool bad_ljet_2 = genq2_p4.DeltaR(light_jets[best_onshell_pair.first]) > 0.4 && genq2_p4.DeltaR(light_jets[best_onshell_pair.second]) > 0.4;
            bool bad_light_jets = bad_ljet_1 || bad_ljet_2;

            int n_problems = bad_lep + bad_bjets + bad_light_jets;
            hm.Fill(debug_on, n_problems);
        }

        if (!hme_offshell && !hme_onshell)
        {
            hm.Fill(numLightJets_both_fail, light_jets.size());

            bool bad_lep = gen_lep_p4.DeltaR(reco_lep_p4) > 0.4;
            if (bad_lep)
            {
                ++has_bad_lep;
            }

            bool bad_bjet_1 = genb1_p4.DeltaR(reco_bj1_p4) > 0.4 && genb1_p4.DeltaR(reco_bj2_p4) > 0.4;
            bool bad_bjet_2 = genb2_p4.DeltaR(reco_bj1_p4) > 0.4 && genb2_p4.DeltaR(reco_bj2_p4) > 0.4;
            bool bad_bjets = bad_bjet_1 || bad_bjet_2;
            if (bad_bjets)
            {
                ++has_bad_bjets;
            }

            std::vector<TLorentzVector> slj = {light_jets[best_onshell_pair.first], 
                                               light_jets[best_onshell_pair.second],
                                               light_jets[best_offshell_pair.first],
                                               light_jets[best_offshell_pair.second]};

            bool bad_light_jet_1 = std::all_of(slj.begin(), slj.end(), [&genq1_p4](TLorentzVector const& v){ return v.DeltaR(genq1_p4) > 0.4; });
            bool bad_light_jet_2 = std::all_of(slj.begin(), slj.end(), [&genq2_p4](TLorentzVector const& v){ return v.DeltaR(genq2_p4) > 0.4; });
            bool bad_light_jets = bad_light_jet_1 || bad_light_jet_2;
            if (bad_light_jets)
            {
                ++has_bad_light_jets;
            }

            int n_problems = bad_lep + bad_bjets + bad_light_jets;
            hm.Fill(debug_hist, n_problems);   
        }

        if (m1 > 0.0 && m2 > 0.0)
        {
            double dm = std::abs(m1 - m2);
            hm.Fill(nJets_vs_massDiff, light_jets.size(), dm);

            if (best_offshell_pair == best_onshell_pair)
            {
                hm.Fill(same_pair, dm);
            }
            else
            {
                hm.Fill(diff_pairs, dm);
            }

            h->Fill(m1, 0.5);
            h->Fill(m2, 0.5);
        }
        else if (m1 > 0.0 && m2 < 0.0)
        {
            h->Fill(m1);
        }
        else if (m1 < 0.0 && m2 > 0.0)
        {
            h->Fill(m2);
        }

        hme_worked += (hme_onshell || hme_offshell);

        output << event << "\n";
    }

    std::vector<int> counts{ has_bad_light_jets,
                             has_bad_bjets,
                             has_bad_lep };

    std::vector<char const*> labels{ "wrong light jets", 
                                     "wrong b jets",
                                     "wrong lepton" };

    int nx = labels.size();
    auto debug_breakdown = std::make_unique<TH1F>("debug_breakdown", "Breakdown of events where HME failed to categories", nx, 0, nx);
    debug_breakdown->SetStats(0);
    debug_breakdown->SetFillStyle(3544);
    debug_breakdown->SetLineWidth(2);
    debug_breakdown->SetFillColorAlpha(kBlue, 0.75);

    auto canvas = std::make_unique<TCanvas>("canvas", "canvas");
    canvas->SetGrid();

    for (int b = 1; b <= nx; ++b)
    {
        debug_breakdown->SetBinContent(b, counts[b-1]);
        debug_breakdown->GetXaxis()->SetBinLabel(b, labels[b-1]);
    }

    debug_breakdown->Draw();
    canvas->SaveAs("histograms/debug_breakdown.png");

    hm.Draw();

    hm.DrawStack({hme_mass_onshell, hme_mass_offshell}, "Impact of onshell/offshell light jet candidates on HME", "hme_offshell_vs_hme_onshell.png");
    hm.DrawStack({same_pair, diff_pairs}, "Impact of light jet selection on HME mass difference", "mass_diff_cmp.png");

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Finished processing, total events = " << nEvents << "\n";
    std::cout << "Events passed to HME = " << hme_events << "\n"; 
    std::cout << "HME successful = " << hme_worked << "\n"; 
    std::cout << "Combned width = " << InterquantileRange(h) << "\n";
    std::cout << "Processing time = " << elapsed.count() << " s\n";

    return 0;
}