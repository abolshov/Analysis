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
#include "TH2.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "MatchingTools.hpp"
#include "HistManager.hpp"
#include "EstimatorTools.hpp"

static constexpr int MAX_AK4_GENJET = 35;
static constexpr int MAX_GENPART = 350;

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

int main()
{
    auto start = std::chrono::system_clock::now();
    auto gStyle = std::make_unique<TStyle>();
    gStyle->SetPalette(kRainBow);
    TFile *myFile = TFile::Open("GluGluToRadionToHHTo2B2WToLNu2J_M-450.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    auto file_pdf = std::make_unique<TFile>("pdfs.root", "READ");
    // auto lead_bjet_pt_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("lead_bjet_pt_pdf")));
    auto lead_bjet_pt_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_pdf->Get("lead_b_resc_picked_by_j")));
    // auto pdf2d = std::unique_ptr<TH2F>(static_cast<TH2F*>(file_pdf->Get("bjet_2d_pt_pdf")));

    // auto file_weight_pdf = std::make_unique<TFile>("pdfs_test.root", "READ");
    // auto Hbb_corr_weight_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_weight_pdf->Get("Hbb_corr_vs_uncorr_pt")));
    // auto Hbb_mass_weight_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_weight_pdf->Get("Hbb_mass_w_true_corr")));
    // auto met_pdf = std::unique_ptr<TH1F>(static_cast<TH1F*>(file_weight_pdf->Get("nu_frac_met")));

    TRandom3 rg;
    rg.SetSeed(0);

    std::cout << std::setprecision(3);
    std::cout << std::boolalpha;

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

    // hists
    HistManager hm;

    std::string hme_mass("hme_mass");
    hm.Add(hme_mass, "HME X->HH mass", {"X->HH mass, [GeV]", "Count"}, {0, 1000}, 250);

    int selected_events = 0;

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);

        if (!sig.empty())
        {
            if (HasOnlyEleMu(sig, GenPart_genPartIdxMother, GenPart_pdgId)) 
            {
                KinematicData genpart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, static_cast<int>(nGenPart)};
                KinematicData genjet_ak4{GenJetAK4_pt, GenJetAK4_eta, GenJetAK4_phi, GenJetAK4_mass, static_cast<int>(nGenJetAK4)};

                TLorentzVector l_p4 = GetP4(genpart, sig[SIG::l]);
                TLorentzVector nu_p4 = GetP4(genpart, sig[SIG::nu]);

                TLorentzVector bq1_p4 = GetP4(genpart, sig[SIG::b]);
                TLorentzVector bq2_p4 = GetP4(genpart, sig[SIG::bbar]);

                TLorentzVector q1_p4 = GetP4(genpart, sig[SIG::q1]);
                TLorentzVector q2_p4 = GetP4(genpart, sig[SIG::q2]);

                if (q1_p4.DeltaR(q2_p4) < DR_THRESH)
                {
                    continue;
                }

                if (bq1_p4.DeltaR(bq2_p4) < DR_THRESH)
                {
                    continue;
                }

                if (q1_p4.Pt() < MIN_GENJET_PT || q2_p4.Pt() < MIN_GENJET_PT)
                {
                    continue;
                }

                if (std::abs(q1_p4.Eta()) > MAX_GENJET_ETA || std::abs(q2_p4.Eta()) > MAX_GENJET_ETA)
                {
                    continue;
                }

                if (q1_p4.DeltaR(l_p4) < DR_THRESH || q2_p4.DeltaR(l_p4) < DR_THRESH)
                {
                    continue;
                }

                TLorentzVector Hbb_p4 = GetP4(genpart, sig[SIG::H_bb]);
            
                // CUTS START
                if (!PassLeptonCut(l_p4))
                {
                    continue;
                }

                std::vector<int> selected_jets = PrimaryJetSelection(genjet_ak4);
                if (selected_jets.size() < 4)
                {
                    continue;
                }

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

                if (num_bflav_jets < 2)
                {
                    continue;
                }

                auto PtCut = [](TLorentzVector const& p) { return p.Pt() > MIN_GENJET_PT; };
                int num_bjets_passing_pt = std::count_if(b_jets.begin(), b_jets.end(), PtCut);
                if (num_bjets_passing_pt < 2)
                {
                    continue;
                }

                auto EtaCut = [](TLorentzVector const& p) { return std::abs(p.Eta()) < MAX_GENJET_ETA; };
                int num_bjets_passing_eta = std::count_if(b_jets.begin(), b_jets.end(), EtaCut);
                if (num_bjets_passing_eta < 2)
                {
                    continue;
                }

                if (light_jets.size() < 2)
                {
                    continue;
                }

                TLorentzVector bj1_p4, bj2_p4;
                if (b_jets.size() == 2)
                {
                    bj1_p4 = b_jets[0].DeltaR(bq1_p4) < b_jets[1].DeltaR(bq1_p4) ? b_jets[0] : b_jets[1];
                    bj2_p4 = b_jets[0].DeltaR(bq2_p4) < b_jets[1].DeltaR(bq2_p4) ? b_jets[0] : b_jets[1];
                }
                else
                {
                    continue;
                }

                if (bj1_p4 == bj2_p4)
                {
                    continue;
                }

                TLorentzVector genmet;
                genmet.SetPtEtaPhiM(GenMET_pt, 0, GenMET_phi, 0);

                // start hme here
                ++selected_events;

                // choose 2 pairs of light jet candidates
                auto best_onshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(80.0 - (p1 + p2).M()); });
                auto best_offshell_pair = ChooseBestPair(light_jets, [](TLorentzVector const& p1, TLorentzVector const& p2){ return std::abs(40.0 - (p1 + p2).M()); });

                auto [idx1, idx2] = best_onshell_pair;
                std::vector<TLorentzVector> input_onshell = {bj1_p4, bj2_p4, light_jets[idx1], light_jets[idx2], l_p4, genmet};

                std::tie(idx1, idx2) = best_offshell_pair;
                std::vector<TLorentzVector> input_offshell = {bj1_p4, bj2_p4, light_jets[idx1], light_jets[idx2], l_p4, genmet};

                auto hme = EstimateMass(input_onshell, lead_bjet_pt_pdf, rg, i);
                if (hme)
                {
                    auto [mass, eff] = hme.value();
                    hm.FillWeighted(hme_mass, mass, 0.5);
                }

                hme = EstimateMass(input_offshell, lead_bjet_pt_pdf, rg, i);
                if (hme)
                {
                    auto [mass, eff] = hme.value();
                    hm.FillWeighted(hme_mass, mass, 0.5);
                }
            }
        }
    }

    hm.Draw();

    std::cout << "Finished processing\n";
    std::cout << "Selected events = " << selected_events << "\n";

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    return 0;
}
