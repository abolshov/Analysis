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

    Int_t           lep1_type;
    Int_t           lep1_genLep_kind;

    Float_t         PuppiMET_pt;
    Float_t         PuppiMET_phi;

    Int_t           ncentralJet;
    Float_t         centralJet_pt[N_RECO_JETS];   
    Float_t         centralJet_eta[N_RECO_JETS];   
    Float_t         centralJet_phi[N_RECO_JETS];   
    Float_t         centralJet_mass[N_RECO_JETS]; 
    Float_t         centralJet_btagPNetB[N_RECO_JETS]; 
    Float_t         centralJet_PNetRegPtRawCorr[N_RECO_JETS];    
    Float_t         centralJet_PNetRegPtRawRes[N_RECO_JETS]; 
    // Float_t         centralJet_PNetRegPtRawCorrNeutrino[N_RECO_JETS];  

    myTree->SetBranchAddress("lep1_pt", &lep1_pt);
    myTree->SetBranchAddress("lep1_eta", &lep1_eta);
    myTree->SetBranchAddress("lep1_phi", &lep1_phi);
    myTree->SetBranchAddress("lep1_mass", &lep1_mass);

    myTree->SetBranchAddress("lep1_type", &lep1_type);
    myTree->SetBranchAddress("lep1_genLep_kind", &lep1_genLep_kind);

    myTree->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt);
    myTree->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi);

    myTree->SetBranchAddress("ncentralJet", &ncentralJet);
    myTree->SetBranchAddress("centralJet_pt", centralJet_pt);
    myTree->SetBranchAddress("centralJet_eta", centralJet_eta);
    myTree->SetBranchAddress("centralJet_phi", centralJet_phi);
    myTree->SetBranchAddress("centralJet_mass", centralJet_mass);
    myTree->SetBranchAddress("centralJet_btagPNetB", centralJet_btagPNetB);
    myTree->SetBranchAddress("centralJet_PNetRegPtRawCorr", centralJet_PNetRegPtRawCorr);
    myTree->SetBranchAddress("centralJet_PNetRegPtRawRes", centralJet_PNetRegPtRawRes);
    // myTree->SetBranchAddress("centralJet_PNetRegPtRawCorrNeutrino", centralJet_PNetRegPtRawCorrNeutrino);

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

    std::string hme_mass("hme_mass");
    hm.Add(hme_mass, "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 2000}, 100);

    auto h = std::make_unique<TH1F>("h", "h", 200, 0, 2000);

    int hme_events = 0;
    int hme_worked = 0;

    // int lep_reco = 0;
    // int jet_mult = 0;
    // int res_top = 0;
    // int quark_accept = 0;
    // int matching = 0;
    // int tot = 0;

    int nEvents = myTree->GetEntries();
    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        if (event % 2 != 1)
        {
            continue;
        }
        // ++tot;

        if (ncentralJet < 4)
        {
            continue;
        }
        // ++jet_mult;

        // bool reco_lep1_mu = (lep1_type == 2);
        // bool reco_lep1_ele = (lep1_type == 1);

        // bool gen_lep1_mu = ((lep1_genLep_kind == 2) || (lep1_genLep_kind == 4));
        // bool gen_lep1_ele = ((lep1_genLep_kind == 1) || (lep1_genLep_kind == 3));

        // bool corr_lep_reco = ((reco_lep1_mu && gen_lep1_mu) || (reco_lep1_ele && gen_lep1_ele));
        
        // if (!corr_lep_reco)
        // {
        //     continue;
        // }
        // ++lep_reco;

        TLorentzVector genb1_p4, genb2_p4, genq1_p4, genq2_p4; // quarks
        genb1_p4.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        genb2_p4.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);
        genq1_p4.SetPtEtaPhiM(genV2prod1_pt, genV2prod1_eta, genV2prod1_phi, genV2prod1_mass);
        genq2_p4.SetPtEtaPhiM(genV2prod2_pt, genV2prod2_eta, genV2prod2_phi, genV2prod2_mass);

        if (genb1_p4.DeltaR(genb2_p4) < 0.4 || genq1_p4.DeltaR(genq2_p4) < 0.4)
        {
            continue;
        }
        // ++res_top;

        // bool bq_accept = (genb1_p4.Pt() > 20.0 && std::abs(genb1_p4.Eta()) < 2.5) && (genb2_p4.Pt() > 20.0 && std::abs(genb2_p4.Eta()) < 2.5);
        // bool lq_accept = (genq1_p4.Pt() > 20.0 && std::abs(genq1_p4.Eta()) < 5.0) && (genq2_p4.Pt() > 20.0 && std::abs(genq2_p4.Eta()) < 5.0);
        // bool quarks_accept = bq_accept && lq_accept;
        // if (!quarks_accept)
        // {
        //     continue;
        // }
        // ++quark_accept;

        TLorentzVector gen_bj1_p4, gen_bj2_p4, gen_lj1_p4, gen_lj2_p4; // gen jets matched to quarks
        gen_bj1_p4.SetPtEtaPhiM(genb1_vis_pt, genb1_vis_eta, genb1_vis_phi, genb1_vis_mass);
        gen_bj2_p4.SetPtEtaPhiM(genb2_vis_pt, genb2_vis_eta, genb2_vis_phi, genb2_vis_mass);
        gen_lj1_p4.SetPtEtaPhiM(genV2prod1_vis_pt, genV2prod1_vis_eta, genV2prod1_vis_phi, genV2prod1_vis_mass);
        gen_lj2_p4.SetPtEtaPhiM(genV2prod2_vis_pt, genV2prod2_vis_eta, genV2prod2_vis_phi, genV2prod2_vis_mass);

        // std::vector<TLorentzVector> gen_jets = {gen_bj1_p4, gen_bj2_p4, gen_lj1_p4, gen_lj2_p4};
        // if (!std::all_of(gen_jets.begin(), gen_jets.end(), [](TLorentzVector const& p){ return p != TLorentzVector{}; }))
        // {
        //     continue;
        // }
        // ++matching;

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
        std::vector<double> resolutions;
        for (int j = 2; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);

            // use PNet correction to correct p4 of light jets
            jet *= centralJet_PNetRegPtRawCorr[j];
            light_jets.push_back(jet);

            // save resolutions of pt corrections of light jets
            resolutions.push_back(centralJet_pt[j]*centralJet_PNetRegPtRawRes[j]);
        }

        auto best_onshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(80.0 - (p1 + p2).M()); });
        auto best_offshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(40.0 - (p1 + p2).M()); });

        auto [i1, i2] = best_onshell_pair;
        std::vector<TLorentzVector> input_onshell = {reco_bj1_p4, reco_bj2_p4, light_jets[i1], light_jets[i2], reco_lep_p4, reco_met_p4};
        std::pair<double, double> res_onshell = {resolutions[i1], resolutions[i2]};

        auto [j1, j2] = best_offshell_pair;
        std::vector<TLorentzVector> input_offshell = {reco_bj1_p4, reco_bj2_p4, light_jets[j1], light_jets[j2], reco_lep_p4, reco_met_p4};
        std::pair<double, double> res_offshell = {resolutions[j1], resolutions[j2]};

        ++hme_events;

        [[maybe_unused]] auto hme_onshell = EstimateMass(input_onshell, pdf, rg, i, res_onshell);
        [[maybe_unused]] auto hme_offshell = EstimateMass(input_offshell, pdf, rg, i, res_offshell);

        hme_worked += (hme_onshell || hme_offshell);

        if (hme_offshell && hme_onshell)
        {
            auto hme = rg.Uniform(0, 1) > 0.27 ? hme_onshell : hme_offshell;
            if (hme)
            {
                auto [mass, succ_rate] = hme.value();
                hm.Fill(hme_mass, mass);
                h->Fill(mass);
                output << mass << "\n";
            }
        }
        else if (hme_offshell)
        {
            auto [mass, succ_rate] = hme_offshell.value();
            hm.Fill(hme_mass, mass);
            h->Fill(mass);
            output << mass << "\n";
        }
        else if (hme_onshell)
        {
            auto [mass, succ_rate] = hme_onshell.value();
            hm.Fill(hme_mass, mass);
            h->Fill(mass);
            output << mass << "\n";
        }
        else 
        {
            output << -1.0 << "\n";
        }
    }

    // std::vector<char const*> labels = {"total",
    //                                    ">= 4 jets",
    //                                    "lep reco",
    //                                    "resolved",
    //                                    "acceptance",
    //                                    "matching"};

    // std::vector<double> count = {100.0*tot/tot,
    //                              100.0*jet_mult/tot,
    //                              100.0*lep_reco/tot,
    //                              100.0*res_top/tot,
    //                              100.0*quark_accept/tot,
    //                              100.0*matching/tot};

    // int nx = labels.size();

    // auto cutflow = std::make_unique<TH1F>("cutflow", "HME benchmark cutflow", nx, 0, nx);
    // cutflow->SetStats(0);
    // cutflow->SetFillStyle(3544);
    // cutflow->SetLineWidth(2);
    // cutflow->SetFillColorAlpha(kBlue, 0.75);

    // auto canvas = std::make_unique<TCanvas>("canvas", "canvas");
    // canvas->SetGrid();

    // for (int b = 1; b <= nx; ++b)
    // {
    //     cutflow->SetBinContent(b, count[b-1]);
    //     cutflow->GetXaxis()->SetBinLabel(b, labels[b-1]);
    // }

    // cutflow->Draw();
    // canvas->SaveAs("histograms/cutflow.png");

    hm.Draw();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Finished processing, total events = " << nEvents << "\n";
    std::cout << "Events passed to HME = " << hme_events << "\n"; 
    std::cout << "HME successful = " << 100.0*hme_worked/hme_events << "%\n"; 
    std::cout << "Combined width = " << InterquantileRange(h) << "\n";
    std::cout << "Processing time = " << elapsed.count() << " s\n";

    return 0;
}