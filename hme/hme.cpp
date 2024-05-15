#include <iostream>
#include <numeric>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom3.h"

#include "MatchingTools.hpp"
#include "HistManager.hpp"
#include "EstimatorTools.hpp"

static constexpr int MAX_AK4_GENJET = 35;
static constexpr int MAX_GENPART = 350;

int main()
{
    auto start = std::chrono::system_clock::now();
    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);
    TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-600-2.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    auto file_pdf = std::make_unique<TFile>("pdfs.root", "READ");
    auto lead_bjet_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("lead_bjet_pdf")));

    TRandom3 rg;
    rg.SetSeed(0);

    UInt_t          nGenPart;
    Float_t         GenPart_eta[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_mass[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_phi[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_pt[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_pdgId[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_status[MAX_GENPART];   //[nGenPart]

    // ak4 jets
    UInt_t          nGenJetAK4;
    Float_t         GenJetAK4_eta[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_mass[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_phi[MAX_AK4_GENJET];   //[nGenJetAK4]
    Float_t         GenJetAK4_pt[MAX_AK4_GENJET];   //[nGenJetAK4]
    Int_t           GenJetAK4_partonFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]
    UChar_t         GenJetAK4_hadronFlavour[MAX_AK4_GENJET];   //[nGenJetAK4]

    Float_t         GenMET_phi;
    Float_t         GenMET_pt;

    // gen particles 
    myTree->SetBranchAddress("nGenPart", &nGenPart);
    myTree->SetBranchAddress("GenPart_eta", &GenPart_eta);
    myTree->SetBranchAddress("GenPart_mass", &GenPart_mass);
    myTree->SetBranchAddress("GenPart_phi", &GenPart_phi);
    myTree->SetBranchAddress("GenPart_pt", &GenPart_pt);
    myTree->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    myTree->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    myTree->SetBranchAddress("GenPart_status", &GenPart_status);
    myTree->SetBranchAddress("GenPart_phi", &GenPart_phi);
    // myTree->SetBranchAddress("GenPart_statusFlags", &GenPart_statusFlags);

    // ak4 gen jets
    myTree->SetBranchAddress("nGenJet", &nGenJetAK4);
    myTree->SetBranchAddress("GenJet_eta", &GenJetAK4_eta);
    myTree->SetBranchAddress("GenJet_mass", &GenJetAK4_mass);
    myTree->SetBranchAddress("GenJet_phi", &GenJetAK4_phi);
    myTree->SetBranchAddress("GenJet_pt", &GenJetAK4_pt);
    myTree->SetBranchAddress("GenJet_partonFlavour", &GenJetAK4_partonFlavour);
    myTree->SetBranchAddress("GenJet_hadronFlavour", &GenJetAK4_hadronFlavour);

    myTree->SetBranchAddress("GenMET_phi", &GenMET_phi);
    myTree->SetBranchAddress("GenMET_pt", &GenMET_pt);

    int nEvents = myTree->GetEntries();

    // counters
    int has_bbWW_decay = 0;
    int has_tau = 0;
    int has_only_emu = 0;
    int has_at_least_4_primary_jets = 0;
    int has_at_least_4_clean_jets = 0;

    int has_at_least_1_bflav_jet = 0;
    int has_at_least_2_bflav_jet = 0;
    int has_at_least_2_bflav_jet_passing_pt = 0;
    int has_at_least_2_bflav_jet_passing_eta = 0;
    int has_accept_lep = 0;

    // hists
    HistManager hm;

    std::string res_mass_default("res_mass_default");
    std::string succ_rate_default("succ_rate_default");
    hm.Add(res_mass_default, "X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);
    hm.Add(succ_rate_default, "HME Success rate", {"Success rate", "Count"}, {0, 1}, 10);

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        if (!sig.empty())
        {
            ++has_bbWW_decay;
            X_mass = static_cast<int>(GenPart_mass[sig[SIG::X]]);

            if (HasOnlyEleMu(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                ++has_only_emu;

                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};

                TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                // int b_idx = sig[SIG::b];
                // int bbar_idx = sig[SIG::bbar];
                int q1_idx = sig[SIG::q1];
                int q2_idx = sig[SIG::q2];

                // CUTS START
                if (!PassLeptonCut(l_p4))
                {
                    continue;
                }
                ++has_accept_lep;

                std::vector<int> selected_jets = PrimaryJetSelection(genjet_ak4);
                if (selected_jets.size() < 4)
                {
                    continue;
                }
                ++has_at_least_4_primary_jets;   

                // perform jet cleaning: reject jets overlapping with the lepton
                auto JetLepOverlap = [&l_p4, &genjet_ak4](int i)
                { 
                    TLorentzVector p = GetP4(genjet_ak4, i);
                    return p.DeltaR(l_p4) < DR_THRESH;
                };
                selected_jets.erase(std::remove_if(selected_jets.begin(), selected_jets.end(), JetLepOverlap), selected_jets.end());
                if (selected_jets.size() < 4)
                {
                    continue;
                }
                ++has_at_least_4_clean_jets;

                std::vector<TLorentzVector> b_jets, light_jets;
                for (auto j: selected_jets)
                {
                    unsigned flav = static_cast<unsigned>(GenJetAK4_hadronFlavour[j]);
                    if (flav == 5)
                    {
                        b_jets.push_back(GetP4(genjet_ak4, j));
                    }
                    else
                    {
                        light_jets.push_back(GetP4(genjet_ak4, j));
                    }
                }

                int num_bflav_jets = b_jets.size();
                if (num_bflav_jets < 1)
                {
                    continue;
                }
                ++has_at_least_1_bflav_jet;

                if (num_bflav_jets < 2)
                {
                    continue;
                }
                ++has_at_least_2_bflav_jet;

                auto PtCut = [](TLorentzVector const& p) { return p.Pt() > MIN_B_GENJET_PT; };
                int num_bjets_passing_pt = std::count_if(b_jets.begin(), b_jets.end(), PtCut);
                if (num_bjets_passing_pt < 2)
                {
                    continue;
                }
                ++has_at_least_2_bflav_jet_passing_pt;

                auto EtaCut = [](TLorentzVector const& p) { return std::abs(p.Eta()) < MAX_GENJET_ETA; };
                int num_bjets_passing_eta = std::count_if(b_jets.begin(), b_jets.end(), EtaCut);
                if (num_bjets_passing_eta < 2)
                {
                    continue;
                }
                ++has_at_least_2_bflav_jet_passing_eta;

                if (light_jets.size() < 2)
                {
                    continue;
                }

                TLorentzVector lq1_p4 = GetP4(genpart, q1_idx);
                TLorentzVector lq2_p4 = GetP4(genpart, q2_idx);

                TLorentzVector bj1_p4, bj2_p4;
                if (b_jets.size() == 2)
                {
                    bj1_p4 = b_jets[0];
                    bj2_p4 = b_jets[1];
                }
                else
                {
                    continue;
                    // double min_delta_m = GenPart_mass[h_idx];
                    // for (size_t f = 0; f < b_jets.size(); ++f)
                    // {
                    //     for (size_t s = 1; s < b_jets.size(); ++s)
                    //     {
                    //         double m_jj = (b_jets[f] + b_jets[s]).M();
                    //         if (std::abs(m_jj - GenPart_mass[h_idx]) < min_delta_m)
                    //         {
                    //             bj1_p4 = b_jets[f];
                    //             bj2_p4 = b_jets[s];
                    //             min_delta_m = m_jj;
                    //         }
                    //     }
                    // }
                }

                TLorentzVector genmet;
                genmet.SetPtEtaPhiM(GenMET_pt, 0, GenMET_phi, 0);

                std::vector<TLorentzVector> input{ bj1_p4, bj2_p4, lq1_p4, lq2_p4, genmet };
                auto res = EstimateMass(input, lead_bjet_pdf, rg, i);
                if (res)
                {
                    auto [mass, eff] = res.value();
                    hm.Fill(res_mass_default, mass);
                    hm.Fill(succ_rate_default, eff);
                }

                // loop over all possible pairs of light jets
                // int num_light_jets = light_jets.size();
                // for (int j1 = 0; j1 < num_light_jets; ++j1)
                // {
                //     for (int j2 = j1 + 1; j2 < num_light_jets; ++j2)
                //     {
                //         std::vector<TLorentzVector> input{ bj1_p4, bj2_p4, light_jets[j1], light_jets[j2], l_p4, genmet };
                //         auto res = EstimateMass(input, lead_bjet_pdf, rg, i);
                //         if (res)
                //         {
                //             auto [mass, eff] = res.value();
                //             hm.Fill(res_mass_default, mass);
                //             hm.Fill(succ_rate_default, eff);
                //         }
                //     }
                // }
            }
            else
            {
                ++has_tau;
            }
        }
    }

    std::cout << "Finished processing\n";
    hm.Draw();

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << std::setprecision(3);
    std::cout << "nEvents = " << nEvents << ", processing time = " << elapsed.count() << " s\n"
              << "\tHave bbWW decay: " << has_bbWW_decay << "/" << nEvents << " (" << 100.0*has_bbWW_decay/nEvents << "%)\n"
              << "\tHave tau leptons: " << has_tau << "/" << nEvents << " (" << 100.0*has_tau/nEvents << "%)\n" 
              << "\tHave only electrons/muons: " << has_only_emu << "/" << nEvents << " (" << 100.0*has_only_emu/nEvents << "%)\n"
              << "\tHave lepton in acceptance region: " << has_accept_lep << "/" << has_only_emu << " (" << 100.0*has_accept_lep/has_only_emu << "%)\n"
              << "\tHave at least 4 primary jets: " << has_at_least_4_primary_jets << "/" << has_accept_lep << " (" << 100.0*has_at_least_4_primary_jets/has_accept_lep << "%)\n"
              << "\tHave at least 4 clean jets: " << has_at_least_4_clean_jets << "/" << has_at_least_4_primary_jets << " (" << 100.0*has_at_least_4_clean_jets/has_at_least_4_primary_jets << "%)\n"
              << "\tHave at least 1 b-flavored jet: " << has_at_least_1_bflav_jet << "/" << has_at_least_4_clean_jets << " (" << 100.0*has_at_least_1_bflav_jet/has_at_least_4_clean_jets << "%)\n"
              << "\tHave at least 2 b-flavored jets: " << has_at_least_2_bflav_jet << "/" << has_at_least_1_bflav_jet << " (" << 100.0*has_at_least_2_bflav_jet/has_at_least_1_bflav_jet << "%)\n"
              << "\tHave at least 2 b-flavored jets passing pt cut: " << has_at_least_2_bflav_jet_passing_pt << "/" << has_at_least_2_bflav_jet << " (" << 100.0*has_at_least_2_bflav_jet_passing_pt/has_at_least_2_bflav_jet << "%)\n"
              << "\tHave at least 2 b-flavored jets passing eta cut: " << has_at_least_2_bflav_jet_passing_eta << "/" << has_at_least_2_bflav_jet_passing_pt << " (" << 100.0*has_at_least_2_bflav_jet_passing_eta/has_at_least_2_bflav_jet_passing_pt << "%)\n";
    std::cout << "============================================================================\n";

    return 0;
}
