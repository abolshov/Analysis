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

double InterquantileRange(std::unique_ptr<TH1F> const& h, unsigned l = 16, unsigned r = 84)
{
    int const nq = 100;
    double xq[nq];  // position where to compute the quantiles in [0,1]
    double yq[nq];  // array to contain the quantiles
    for (int i = 0; i < nq; ++i) 
    {
        xq[i] = static_cast<double>(i+1)/nq;
    }
    h->GetQuantiles(nq, yq, xq);
    return yq[r - 1] - yq[l - 1];
}

static constexpr int N_RECO_JETS = 12;

int main()
{
    std::cout << std::setprecision(4);
    auto start = std::chrono::system_clock::now();

    TH1::AddDirectory(false);

    auto file_pdf = std::make_unique<TFile>("pdf.root", "READ");
    auto pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_b1")));

    auto pdf_b1b2 = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get("pdf_b1b2")));
    auto pdf_hh_dEtadPhi = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get("pdf_hh_dEtadPhi")));
    auto pdf_hh_pt_e = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get("pdf_hh_pt_e")));

    std::vector<std::unique_ptr<TH2F>> pdfs_2d;
    pdfs_2d.push_back(std::move(pdf_b1b2));
    pdfs_2d.push_back(std::move(pdf_hh_dEtadPhi));
    pdfs_2d.push_back(std::move(pdf_hh_pt_e));

    auto pdf_numet_pt = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_numet_pt")));
    auto pdf_numet_dphi = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_numet_dphi")));
    auto pdf_nulep_deta = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_nulep_deta")));
    auto pdf_hh_dphi = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_hh_dphi")));
    auto pdf_mbb = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_mbb")));
    auto pdf_mww = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_mbb"))); // using anything but mbb here makes everyhting worse!
    auto pdf_hh_deta = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf_hh_deta")));

    std::vector<std::unique_ptr<TH1F>> pdfs_1d;
    pdfs_1d.push_back(std::move(pdf_numet_pt));
    pdfs_1d.push_back(std::move(pdf_numet_dphi));
    pdfs_1d.push_back(std::move(pdf_nulep_deta));
    pdfs_1d.push_back(std::move(pdf_hh_dphi));
    pdfs_1d.push_back(std::move(pdf_mbb));
    pdfs_1d.push_back(std::move(pdf_mww));
    pdfs_1d.push_back(std::move(pdf_hh_deta));

    file_pdf->Close();

    TFile *myFile = TFile::Open("nano_0.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    TRandom3 rg;
    rg.SetSeed(42);

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
    hm.Add(hme_mass, "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 2500}, 100);

    auto h = std::make_unique<TH1F>("h", "h", 200, 0, 2500);

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

        bool reco_lep1_mu = (lep1_type == 2);
        bool reco_lep1_ele = (lep1_type == 1);

        bool gen_lep1_mu = ((lep1_genLep_kind == 2) || (lep1_genLep_kind == 4));
        bool gen_lep1_ele = ((lep1_genLep_kind == 1) || (lep1_genLep_kind == 3));

        bool corr_lep_reco = ((reco_lep1_mu && gen_lep1_mu) || (reco_lep1_ele && gen_lep1_ele));
        
        if (!corr_lep_reco)
        {
            continue;
        }

        TLorentzVector genb1_p4, genb2_p4, genq1_p4, genq2_p4; // quarks
        genb1_p4.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
        genb2_p4.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);
        genq1_p4.SetPtEtaPhiM(genV2prod1_pt, genV2prod1_eta, genV2prod1_phi, genV2prod1_mass);
        genq2_p4.SetPtEtaPhiM(genV2prod2_pt, genV2prod2_eta, genV2prod2_phi, genV2prod2_mass);

        if (genb1_p4.DeltaR(genb2_p4) < 0.4 || genq1_p4.DeltaR(genq2_p4) < 0.4)
        {
            continue;
        }

        bool bq_accept = (genb1_p4.Pt() > 20.0 && std::abs(genb1_p4.Eta()) < 2.5) && (genb2_p4.Pt() > 20.0 && std::abs(genb2_p4.Eta()) < 2.5);
        bool lq_accept = (genq1_p4.Pt() > 20.0 && std::abs(genq1_p4.Eta()) < 5.0) && (genq2_p4.Pt() > 20.0 && std::abs(genq2_p4.Eta()) < 5.0);
        bool quarks_accept = bq_accept && lq_accept;
        if (!quarks_accept)
        {
            continue;
        }

        bool gen_b_match_exists = genb1_vis_pt > 0.0 && genb2_vis_pt > 0.0;
        bool gen_light_match_exists = genV2prod1_vis_pt > 0.0 && genV2prod2_vis_pt > 0.0;
        bool gen_match_exists = gen_b_match_exists && gen_light_match_exists;
        if (!gen_match_exists)
        {
            continue;
        }

        ++hme_events;

        TLorentzVector reco_lep_p4, reco_bj1_p4, reco_bj2_p4, reco_met_p4;
        reco_lep_p4.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_mass);
        reco_met_p4.SetPtEtaPhiM(PuppiMET_pt, 0.0, PuppiMET_phi, 0.0);

        std::vector<TLorentzVector> jets;
        std::vector<double> resolutions;
        for (int j = 0; j < ncentralJet; ++j)
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(centralJet_pt[j], centralJet_eta[j], centralJet_phi[j], centralJet_mass[j]);
            jets.push_back(jet);

            jet *= centralJet_PNetRegPtRawCorr[j];
            resolutions.push_back(centralJet_pt[j]*centralJet_PNetRegPtRawRes[j]);
        }

        bool worked = false;
        // int n_jets = jets.size();
        int n_jets = 3;
        std::vector<double> estimations;
        std::vector<double> efficiencies;
        int comb_id = 1;
        for (int bj1_idx = 0; bj1_idx < n_jets; ++bj1_idx)
        {
            for (int bj2_idx = bj1_idx + 1; bj2_idx < n_jets; ++bj2_idx)
            {
                TLorentzVector const& bj1_p4 = jets[bj1_idx].Pt() > jets[bj2_idx].Pt() ? jets[bj1_idx] : jets[bj2_idx];
                TLorentzVector const& bj2_p4 = jets[bj1_idx].Pt() > jets[bj2_idx].Pt() ? jets[bj2_idx] : jets[bj1_idx];

                auto pairs = MakePairs(jets.size(), bj1_idx, bj2_idx);
                for (auto [lj1_idx, lj2_idx]: pairs)
                {
                    TLorentzVector const& lj1_p4 = jets[lj1_idx];
                    TLorentzVector const& lj2_p4 = jets[lj2_idx];

                    double lj1_res = resolutions[lj1_idx];
                    double lj2_res = resolutions[lj2_idx];

                    std::vector<TLorentzVector> input = {bj1_p4, bj2_p4, lj1_p4, lj2_p4, reco_lep_p4, reco_met_p4};
                    std::pair<double, double> lj_resolutions = {lj1_res, lj2_res};

                    // auto hme = EstimateMass(input, pdf, rg, i, lj_resolutions);
                    auto hme = EstimateMass(input, pdfs_2d, pdfs_1d, rg, i, comb_id++, lj_resolutions);
                    if (hme)
                    {
                        worked = true;
                        auto [mass, succ_rate] = hme.value();
                        estimations.push_back(mass);
                        efficiencies.push_back(succ_rate);
                    }
                }

                // auto [lj1_idx, lj2_idx] = FindByAngle(jets, bj1_idx, bj2_idx);
                // TLorentzVector lj1_p4 = jets[lj1_idx];
                // TLorentzVector lj2_p4 = jets[lj2_idx];

                // double lj1_res = resolutions[lj1_idx];
                // double lj2_res = resolutions[lj2_idx];

                // std::vector<TLorentzVector> input = {bj1_p4, bj2_p4, lj1_p4, lj2_p4, reco_lep_p4, reco_met_p4};
                // std::pair<double, double> lj_resolutions = {lj1_res, lj2_res};

                // // auto hme = EstimateMass(input, pdf, rg, i, lj_resolutions);
                // auto hme = EstimateMass(input, pdfs_2d, pdfs_1d, rg, i, comb_id++, lj_resolutions);
                // if (hme)
                // {
                //     worked = true;
                //     auto [mass, succ_rate] = hme.value();
                //     estimations.push_back(mass);
                //     efficiencies.push_back(succ_rate);
                // }
            }
        }
        
        if (worked)
        {
            ++hme_worked;
            
            // auto it = std::max_element(estimations.begin(), estimations.end());
            // hm.Fill(hme_mass, *it);
            // h->Fill(*it);

            auto it = std::max_element(efficiencies.begin(), efficiencies.end());
            int idx = it - efficiencies.begin();
            hm.Fill(hme_mass, estimations[idx]);
            // hm.Fill(hme_integral, *it);
            h->Fill(estimations[idx]);

            // hm.Fill(mass_vs_int_2d, estimations[idx], *it);

            // double mass = std::accumulate(estimations.begin(), estimations.end(), 0.0)/estimations.size(); 
            // hm.Fill(hme_mass, mass);
            // h->Fill(mass);
        }
    }

    hm.Draw();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "Finished processing, total events = " << nEvents << "\n";
    std::cout << "Events passed to HME = " << hme_events << "\n"; 
    std::cout << "HME successful = " << 100.0*hme_worked/hme_events << "%\n"; 
    std::cout << "HME width = " << InterquantileRange(h) << "\n";
    std::cout << "HME value = " << h->GetXaxis()->GetBinCenter(h->GetMaximumBin()) << "\n";
    std::cout << "Processing time = " << elapsed.count() << " s\n";

    myFile->Close();

    return 0;
}