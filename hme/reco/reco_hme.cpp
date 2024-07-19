#include <iostream>
#include <algorithm>
#include <numeric>
#include <memory>
#include <fstream>
#include <optional>
#include <chrono>
#include <cassert>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TString.h"
#include "TRandom3.h"

#include "HistManager.hpp"

using OptPair = std::optional<std::pair<double, double>>;

inline constexpr int N_ATTEMPTS = 10;
inline constexpr int N_ITER = 10000;

inline constexpr double HIGGS_MASS = 125.03;
inline constexpr double HIGGS_WIDTH = 0.004;

inline constexpr double MAX_MASS = 1500;
inline constexpr int N_BINS = 1500;

inline constexpr double TOL = 10e-7;

inline constexpr bool USE_RECO = true;
inline constexpr bool CHEAT = false;

template <typename T>
std::vector<int> sort_indices(T* v, int n) 
{
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&v](int i1, int i2) {return v[i1] > v[i2];});
  return idx;
}

OptPair JetRescFact(TLorentzVector& j1, TLorentzVector& j2, std::unique_ptr<TH1F>& pdf, double mass)
{
    if (j1.Pt() < j2.Pt()) 
    {
        std::swap(j1, j2);
    }

    for (int i = 0; i < N_ATTEMPTS; ++i)
    {
        double c1 = pdf->GetRandom();
    
        double x1 = j2.M2();
        double x2 = 2*c1*(j1*j2);
        double x3 = c1*c1*j1.M2() - mass*mass;

        double discriminant = (x2*x2 - 4.0*x1*x3);

        if (x2 < 0.0)
        {
            continue;
        }
        if (discriminant < 0.0 || std::abs(x1) < TOL)
        {
            continue;
        }

        double c2 = (-1.0*x2 + sqrt(discriminant))/(2*x1);

        if (c2 < 0.0)
        {
            continue;
        }

        return std::make_optional<std::pair<double, double>>(c1, c2);
    }

    return std::nullopt;
}

int main()
{
    auto start = std::chrono::system_clock::now();
    std::cout << std::boolalpha;
    std::cout << "Regime:\n";
    std::cout << "  CHEAT: " << CHEAT << "\n";
    std::cout << "  USE_RECO: " << USE_RECO << "\n";

    std::unique_ptr<TFile> file_pdf;
    std::unique_ptr<TH1F> pdf;

    if constexpr (USE_RECO)
    {
        file_pdf = std::make_unique<TFile>("pdf.root", "READ");
        pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("pdf")));
    }

    if constexpr (!USE_RECO)
    {
        file_pdf = std::make_unique<TFile>("pdfs.root", "READ");
        pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("lead_bjet_pdf")));    
    }

    TFile *myFile = TFile::Open("GluGlutoRadiontoHHto2B2Vto2B2JLNu_M-450/nano_0.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    Double_t        genHVV_pt;
    Double_t        genHVV_eta;
    Double_t        genHVV_phi;
    Double_t        genHVV_mass;

    Int_t           ncentralJet;
    Float_t         centralJet_pt[20];   //[ncentralJet]
    Float_t         centralJet_eta[20];   //[ncentralJet]
    Float_t         centralJet_phi[20];   //[ncentralJet]
    Float_t         centralJet_mass[20];   //[ncentralJet]

    Float_t         centralJet_btagPNetB[20];   //[ncentralJet]
    Float_t         centralJet_btagDeepFlavB[20];   //[ncentralJet]

    Double_t        genb1_pt;
    Double_t        genb1_eta;
    Double_t        genb1_phi;
    Double_t        genb1_mass;

    Double_t        genb2_pt;
    Double_t        genb2_eta;
    Double_t        genb2_phi;
    Double_t        genb2_mass;

    Double_t        genb1_vis_pt;
    Double_t        genb1_vis_eta;
    Double_t        genb1_vis_phi;
    Double_t        genb1_vis_mass;

    Double_t        genb2_vis_pt;
    Double_t        genb2_vis_eta;
    Double_t        genb2_vis_phi;
    Double_t        genb2_vis_mass;

    myTree->SetBranchAddress("genHVV_pt", &genHVV_pt);
    myTree->SetBranchAddress("genHVV_eta", &genHVV_eta);
    myTree->SetBranchAddress("genHVV_phi", &genHVV_phi);
    myTree->SetBranchAddress("genHVV_mass", &genHVV_mass);

    myTree->SetBranchAddress("ncentralJet", &ncentralJet);
    myTree->SetBranchAddress("centralJet_pt", &centralJet_pt);
    myTree->SetBranchAddress("centralJet_eta", &centralJet_eta);
    myTree->SetBranchAddress("centralJet_phi", &centralJet_phi);
    myTree->SetBranchAddress("centralJet_mass", &centralJet_mass);

    myTree->SetBranchAddress("centralJet_btagPNetB", &centralJet_btagPNetB);
    myTree->SetBranchAddress("centralJet_btagDeepFlavB", &centralJet_btagDeepFlavB);

    myTree->SetBranchAddress("genb1_pt", &genb1_pt);
    myTree->SetBranchAddress("genb1_eta", &genb1_eta);
    myTree->SetBranchAddress("genb1_phi", &genb1_phi);
    myTree->SetBranchAddress("genb1_mass", &genb1_mass);

    myTree->SetBranchAddress("genb2_pt", &genb2_pt);
    myTree->SetBranchAddress("genb2_eta", &genb2_eta);
    myTree->SetBranchAddress("genb2_phi", &genb2_phi);
    myTree->SetBranchAddress("genb2_mass", &genb2_mass);

    myTree->SetBranchAddress("genb1_vis_pt", &genb1_vis_pt);
    myTree->SetBranchAddress("genb1_vis_eta", &genb1_vis_eta);
    myTree->SetBranchAddress("genb1_vis_phi", &genb1_vis_phi);
    myTree->SetBranchAddress("genb1_vis_mass", &genb1_vis_mass);

    myTree->SetBranchAddress("genb2_vis_pt", &genb2_vis_pt);
    myTree->SetBranchAddress("genb2_vis_eta", &genb2_vis_eta);
    myTree->SetBranchAddress("genb2_vis_phi", &genb2_vis_phi);
    myTree->SetBranchAddress("genb2_vis_mass", &genb2_vis_mass);

    HistManager hm;

    std::string res_mass_default("res_mass_default");
    hm.Add(res_mass_default, "X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::string hme_mass("hme_mass");
    hm.Add(hme_mass, "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    std::string succ_rate("succ_rate");
    hm.Add(succ_rate, "HME Success rate", {"Success rate", "Count"}, {0, 1}, 10);

    TRandom3 rg;
    rg.SetSeed(0);

    int nEvents = myTree->GetEntries();
    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        TLorentzVector bj1, bj2;
        int j1, j2;
        if constexpr (USE_RECO)
        {
            std::vector<int> sorted_by_PNet = sort_indices<Float_t>(centralJet_btagPNetB, ncentralJet);
            // std::vector<int> sorted_by_DeepFlav = sort_indices<Float_t>(centralJet_btagDeepFlavB, ncentralJet);

            j1 = sorted_by_PNet[0];
            j2 = sorted_by_PNet[1];

            bj1.SetPtEtaPhiM(centralJet_pt[j1], centralJet_eta[j1], centralJet_phi[j1], centralJet_mass[j1]);
            bj2.SetPtEtaPhiM(centralJet_pt[j2], centralJet_eta[j2], centralJet_phi[j2], centralJet_mass[j2]);
        }

        if constexpr (!USE_RECO)
        {
            bj1.SetPtEtaPhiM(genb1_vis_pt, genb1_vis_eta, genb2_vis_phi, genb2_vis_mass);
            bj2.SetPtEtaPhiM(genb2_vis_pt, genb2_vis_eta, genb2_vis_phi, genb2_vis_mass);
        }

        bool first_match = true;
        bool second_match = true;

        if constexpr (CHEAT)
        {
            if constexpr (USE_RECO)
            {
                TLorentzVector bq1, bq2;
                bq1.SetPtEtaPhiM(genb1_pt, genb1_eta, genb1_phi, genb1_mass);
                bq2.SetPtEtaPhiM(genb2_pt, genb2_eta, genb2_phi, genb2_mass);

                first_match = bq1.DeltaR(bj1) < 0.4 || bq1.DeltaR(bj2) < 0.4;
                second_match = bq2.DeltaR(bj1) < 0.4 || bq2.DeltaR(bj2) < 0.4;
            }
            else 
            {
                first_match = genb1_vis_pt > 20 && std::abs(genb1_vis_eta) < 2.5;
                second_match = genb2_vis_pt > 20 && std::abs(genb2_vis_eta) < 2.5;
            }
        }

        if (first_match && second_match)
        {
            TLorentzVector HtoVV;
            HtoVV.SetPtEtaPhiM(genHVV_pt, genHVV_eta, genHVV_phi, genHVV_mass);

            auto evt_mass = std::make_unique<TH1F>("X_mass", Form("X->HH mass: event %d", i), N_BINS, 0.0, MAX_MASS);
            int failed_iter = 0;
            for (int iter = 0; iter < N_ITER; ++iter)
            {
                double mh = rg.Gaus(HIGGS_MASS, HIGGS_WIDTH);
                OptPair b_jet_resc = JetRescFact(bj1, bj2, pdf, mh);
                if (b_jet_resc)
                {
                    auto [c1, c2] = b_jet_resc.value();
                    TLorentzVector bb1 = c1*bj1;
                    TLorentzVector bb2 = c2*bj2;

                    TLorentzVector tmp(HtoVV);
                    HtoVV.SetPtEtaPhiM(genHVV_pt, genHVV_eta, genHVV_phi, genHVV_mass);
                    tmp += bb1;
                    tmp += bb2;

                    double X_mass = tmp.M();
                    evt_mass->Fill(X_mass);
                }
                else 
                {
                    ++failed_iter;
                    continue;
                }
            }

            if (evt_mass->GetEntries())
            {
                int binmax = evt_mass->GetMaximumBin(); 
                double mass = evt_mass->GetXaxis()->GetBinCenter(binmax);
                double success_rate = 1.0 - static_cast<double>(failed_iter)/N_ITER;

                hm.Fill(hme_mass, mass);
                hm.Fill(succ_rate, success_rate);

                hm.Fill(res_mass_default, (HtoVV + bj1 + bj2).M());
            }
        }
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n";

    std::string level = USE_RECO ? "reco" : "gen";

    hm.Draw();
    hm.DrawStack({res_mass_default, hme_mass}, Form("X->HH: HME (gen H->WW & corr %s jet H->bb) vs default calculation", level.data()), "hme_vs_default.png");

    return 0;
}