#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

#include "Utils.hpp"
#include "Matching.hpp"
#include "Clustering.hpp"

constexpr bool debug = false;
constexpr bool matched_debug = false;

constexpr int MAX_GENJET = 16;
constexpr int MAX_GENPART = 300;

int main(int argc, char* argv[])
{
    bool interactive = false;
    if (argc == 2)
    {
        std::string arg(argv[1]);
        if (arg == "-i")
        {
            interactive = true;
            std::cout << "Running in interactive mode\n";
            std::cout << "Possible commands:\n"
                      << "\tc for continue\n"
                      << "\tf for finish\n";
        }
    }
    // TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    TFile *myFile = TFile::Open("NanoAOD_600GeV_1000Events.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    UInt_t          nGenPart;
    Float_t         GenPart_eta[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_mass[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_phi[MAX_GENPART];   //[nGenPart]
    Float_t         GenPart_pt[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_pdgId[MAX_GENPART];   //[nGenPart]
    Int_t           GenPart_status[MAX_GENPART];   //[nGenPart]
    // Int_t           GenPart_statusFlags[MAX_GENPART];   //[nGenPart]

    UInt_t          nGenJet;
    Float_t         GenJet_eta[MAX_GENJET];   //[nGenJet]
    Float_t         GenJet_mass[MAX_GENJET];   //[nGenJet]
    Float_t         GenJet_phi[MAX_GENJET];   //[nGenJet]
    Float_t         GenJet_pt[MAX_GENJET];   //[nGenJet]

    Int_t           GenJet_partonFlavour[MAX_GENJET];   //[nGenJet]
    // UChar_t         GenJet_hadronFlavour[MAX_GENJET];   //[nGenJet] doesn't work for some reason - prints empty spaces

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

    myTree->SetBranchAddress("nGenJet", &nGenJet);
    myTree->SetBranchAddress("GenJet_eta", &GenJet_eta);
    myTree->SetBranchAddress("GenJet_mass", &GenJet_mass);
    myTree->SetBranchAddress("GenJet_phi", &GenJet_phi);
    myTree->SetBranchAddress("GenJet_pt", &GenJet_pt);

    myTree->SetBranchAddress("GenJet_partonFlavour", &GenJet_partonFlavour);
    // myTree->SetBranchAddress("GenJet_hadronFlavour", &GenJet_hadronFlavour);

    auto tot_gen_jets = std::make_unique<TH1I>("tot_gen_jets", "Number of gen jets", 19, 0, 19);
    auto bad_gen_jets = std::make_unique<TH1I>("bad_gen_jets", "Number of bad gen jets", 19, 0, 19);
    auto tot_cand = std::make_unique<TH1I>("tot_cand", "Number of potential jet constituents", 60, 0, 60);
    auto unused_cand = std::make_unique<TH1I>("unused_cand", "Number of unused candidates", 60, 0, 60);
    auto empty_jets = std::make_unique<TH1I>("empty_jets", "Number of empty jets", 19, 0, 19);

    int diHiggsSL_cnt = 0;
    int good_matching = 0;

    int nEvents = myTree->GetEntries();
    std::cout << "nEvents = " << nEvents << "\n";

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        if (interactive)
        {
            std::cout << "Event " << i << "\n";
            char c;
            while(std::cin.get(c))
            {
                if (c == 'c') break;
                if (c == 'f')
                {
                    std::cout << "Finish\n";
                    return 0;
                }
            }
        }

        GenPartIndex idx = IsDiHiggsSL(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, nGenPart);
        if (!idx) continue;

        PtEtaPhiMArray genPart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, nGenPart};
        PtEtaPhiMArray genJet{GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, nGenJet};

        std::vector<Bool_t> candidates = PossibleJetConstituents(GenPart_genPartIdxMother, GenPart_status, nGenPart, idx);
        Int_t n_cand = std::count(candidates.begin(), candidates.end(), true);
        Int_t n_unused = n_cand;
        Int_t n_bad_jets = 0;
        Int_t n_empty_jets = 0;
        // std::cout << "n_cand = " << n_cand << "\n";
        // std::cout << "nGenJet = " << nGenJet << "\n";

        for (Int_t idx = 0; idx < static_cast<Int_t>(nGenPart); ++idx)
        {
            std::cout << std::boolalpha;
            if (candidates[idx]) std::cout << "idx = " << idx
                                           << ", pdgId = " << GenPart_pdgId[idx] 
                                           << ", status = " << GenPart_status[idx] 
                                           << ", mother_idx = " << GenPart_genPartIdxMother[idx] 
                                           << ", mother_pdgId = " << GenPart_pdgId[GenPart_genPartIdxMother[idx]] << "\n";
            if (GenPart_pdgId[idx] == 21 && candidates[idx])
            {
                // auto it = std::find(GenPart_genPartIdxMother, GenPart_genPartIdxMother + nGenPart, idx);
                // std::cout << "\t" << "nGenPart = " << nGenPart << ", it = " << it - GenPart_genPartIdxMother << "\n";
                TLorentzVector g;
                g.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
                std::cout << "\tg = ";
                Print(g);
            }
            if (GenPart_pdgId[idx] == 2212 && candidates[idx])
            {
                TLorentzVector p;
                p.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
                std::cout << "\tp = ";
                Print(p);
            }
            // std::cout << "pdgId = " << GenPart_pdgId[idx] << ", status = " << GenPart_status[idx] << ", mother_idx = " << GenPart_genPartIdxMother[idx] << ": " << candidates[idx] << "\n";
        }   
        // break;

        tot_gen_jets->Fill(nGenJet);
        tot_cand->Fill(n_cand);

        if constexpr (debug) std::cout << "----------------------EVENT #" << i << " START-----------------------------\n";
        std::vector<Overlap> overlaps = FindOverlaps(genJet, genPart, candidates);
        if constexpr (debug)
        {
            std::cout << "Found overlaps:\n";
            for (auto const& ov: overlaps)
            {
                std::cout << ov << "\n";
            }
        }

        std::vector<Int_t> used_parts;
        for (Int_t jetIdx = 0; jetIdx < static_cast<Int_t>(nGenJet); ++jetIdx) // loop over all jets in event
        {
            TLorentzVector jet, jet_cand;
            jet.SetPtEtaPhiM(GenJet_pt[jetIdx], GenJet_eta[jetIdx], GenJet_phi[jetIdx], GenJet_mass[jetIdx]);
            
            if constexpr (debug)
            {
                std::cout << "----------------------JET #" << jetIdx << " START-------------------------------\n";
                std::cout << "GenJet = ";
                Print(jet);
                std::cout << "\n";
            }

            Int_t constituent_count = 0;

            for (Int_t partIdx = 0; partIdx < static_cast<Int_t>(nGenPart); ++partIdx)
            {
                assert(nGenPart == candidates.size());
                if (!candidates[partIdx]) continue; // skip all particles that can't potentially constitue jet

                TLorentzVector part;
                part.SetPtEtaPhiM(GenPart_pt[partIdx], GenPart_eta[partIdx], GenPart_phi[partIdx], GenPart_mass[partIdx]);

                float dR = part.DeltaR(jet);
                if (dR < 0.4)
                {   
                    // if the particle has NOT been used before, add it to jet_cand and mark as used
                    if (std::find(used_parts.begin(), used_parts.end(), partIdx) != used_parts.end()) 
                    {
                        if constexpr (debug)
                        {
                            std::cout << "\n\tWARNING! Particle " << partIdx << " has been used before, skipping in " << jetIdx << " with dR = " << dR << "\n";
                        }  
                        continue;
                    }

                    jet_cand += part;
                    --n_unused;
                    used_parts.push_back(partIdx);
                    ++constituent_count;

                    Int_t motherIdx = GenPart_genPartIdxMother[partIdx];
                    Int_t motherPdgId = (motherIdx == -1) ? 0 : GenPart_pdgId[motherIdx];
                    Int_t partStatus = GenPart_status[partIdx];

                    if constexpr (debug)
                    {
                        std::cout << "\tpartIdx = " << partIdx << "\n"
                                  << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
                                  << "\tmotherPdgId = " << ((motherPdgId == 0) ? "No mother" : std::to_string(motherPdgId)) << "\n"
                                  << "\tpartStatus = " << partStatus << "\n"
                                  << "\tdR = " << dR << "\n";
                        std::cout << "\t";
                        Print(part);
                        std::cout << "\n";
                    }
                }
            }

            if (constituent_count == 0) ++n_empty_jets;
            constituent_count = 0;

            if constexpr (debug)
            {
                std::cout << "jet_cand = ";
                Print(jet_cand);
                std::cout << "dR = " << jet.DeltaR(jet_cand) << "\n";
            }
            
            Double_t pt_ratio = jet_cand.Pt()/jet.Pt();
            if (pt_ratio > 1.1 || pt_ratio < 0.9) ++n_bad_jets;

            if constexpr (debug) std::cout << "----------------------JET #" << jetIdx << " END-------------------------------\n";
        }

        empty_jets->Fill(n_empty_jets);

        used_parts.clear();
        // if (overlaps.size() > 1) break;

        unused_cand->Fill(n_unused);
        bad_gen_jets->Fill(n_bad_jets);

        // std::cout << "n_unused = " << n_unused << "\n";
        // std::cout << "n_bad_jets = " << n_bad_jets << "\n";
        if constexpr (debug) std::cout << "----------------------EVENT #" << i << " END-----------------------------\n";
                
        // break;
        ++diHiggsSL_cnt;

        GenJetIndex match = Match(idx, genPart, genJet, GenPart_pdgId, GenJet_partonFlavour);

        if (!match) continue;
        ++good_matching;

        auto [bj1Idx, bj2Idx, lj1Idx, lj2Idx] = match;
        TLorentzVector bj1, bj2, lj1, lj2;
        bj1.SetPtEtaPhiM(GenJet_pt[bj1Idx], GenJet_eta[bj1Idx], GenJet_phi[bj1Idx], GenJet_mass[bj1Idx]);
        bj2.SetPtEtaPhiM(GenJet_pt[bj2Idx], GenJet_eta[bj2Idx], GenJet_phi[bj2Idx], GenJet_mass[bj2Idx]);
        lj1.SetPtEtaPhiM(GenJet_pt[lj1Idx], GenJet_eta[lj1Idx], GenJet_phi[lj1Idx], GenJet_mass[lj1Idx]);
        lj2.SetPtEtaPhiM(GenJet_pt[lj2Idx], GenJet_eta[lj2Idx], GenJet_phi[lj2Idx], GenJet_mass[lj2Idx]);

        TLorentzVector bq1, bq2, lq1, lq2;
        bq1.SetPtEtaPhiM(GenPart_pt[idx.b1], GenPart_eta[idx.b1], GenPart_phi[idx.b1], GenPart_mass[idx.b1]);
        bq2.SetPtEtaPhiM(GenPart_pt[idx.b2], GenPart_eta[idx.b2], GenPart_phi[idx.b2], GenPart_mass[idx.b2]);
        lq1.SetPtEtaPhiM(GenPart_pt[idx.q1], GenPart_eta[idx.q1], GenPart_phi[idx.q1], GenPart_mass[idx.q1]);
        lq2.SetPtEtaPhiM(GenPart_pt[idx.q2], GenPart_eta[idx.q2], GenPart_phi[idx.q2], GenPart_mass[idx.q2]);

        if constexpr (matched_debug)
        {
            std::cout << "----------------------EVENT #" << i << " MATCHED DEBUG START-----------------------------\n";
            std::cout << "bj1 = ";
            Print(bj1);
            std::cout << "bq1 = ";
            Print(bq1);
            std::cout << "bj2 = ";
            Print(bj2);
            std::cout << "bq2 = ";
            Print(bq2);
            std::cout << "lj1 = ";
            Print(lj1);
            std::cout << "lq1 = ";
            Print(lq1);
            std::cout << "lj2 = ";
            Print(lj2);
            std::cout << "lq2 = ";
            Print(lq2);
            std::cout << "\n";
        }

        TLorentzVector cand_b1, cand_b2, cand_j1, cand_j2;
        for (Int_t partIdx = 0; partIdx < static_cast<Int_t>(nGenPart); ++partIdx)
        {
            if (!candidates[partIdx]) continue;

            TLorentzVector part;
            part.SetPtEtaPhiM(GenPart_pt[partIdx], GenPart_eta[partIdx], GenPart_phi[partIdx], GenPart_mass[partIdx]);

            Int_t motherIdx = GenPart_genPartIdxMother[partIdx];
            Int_t motherPdgId = (motherIdx == -1) ? 0 : GenPart_pdgId[motherIdx];
            Int_t partStatus = GenPart_status[partIdx];

            Float_t dR_b1, dR_b2, dR_j1, dR_j2;
            dR_b1 = bj1.DeltaR(part);
            dR_b2 = bj2.DeltaR(part);
            dR_j1 = lj1.DeltaR(part);
            dR_j2 = lj2.DeltaR(part);

            if (dR_b1 < 0.4) 
            {
                cand_b1 += part;
                if constexpr (matched_debug)
                {
                    std::cout << "\tParticle added to b jet #1\n";
                    std::cout << "\tpartIdx = " << partIdx << "\n"
                            << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
                            << "\tmotherPdgId = " << ((motherPdgId == 0) ? "No mother" : std::to_string(motherPdgId)) << "\n"
                            << "\tpartStatus = " << partStatus << "\n"
                            << "\tdR = " << dR_b1 << "\n";
                    std::cout << "\t";
                    Print(part);
                    std::cout << "\n";
                }
            }
            if (dR_b2 < 0.4) 
            {
                cand_b2 += part;
                if constexpr (matched_debug)
                {
                    std::cout << "\tParticle added to b jet #2\n";
                    std::cout << "\tpartIdx = " << partIdx << "\n"
                            << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
                            << "\tmotherPdgId = " << ((motherPdgId == 0) ? "No mother" : std::to_string(motherPdgId)) << "\n"
                            << "\tpartStatus = " << partStatus << "\n"
                            << "\tdR = " << dR_b2 << "\n";
                    std::cout << "\t";
                    Print(part);
                    std::cout << "\n";
                }
            }
            if (dR_j1 < 0.4) 
            {
                cand_j1 += part;
                if constexpr (matched_debug)
                {
                    std::cout << "\tParticle added to light jet #1\n";
                    std::cout << "\tpartIdx = " << partIdx << "\n"
                            << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
                            << "\tmotherPdgId = " << ((motherPdgId == 0) ? "No mother" : std::to_string(motherPdgId)) << "\n"
                            << "\tpartStatus = " << partStatus << "\n"
                            << "\tdR = " << dR_j1 << "\n";
                    std::cout << "\t";
                    Print(part);
                    std::cout << "\n";
                }
            }
            if (dR_j2 < 0.4) 
            {
                cand_j2 += part;
                if constexpr (matched_debug)
                {
                    std::cout << "\tParticle added to light jet #2\n";
                    std::cout << "\tpartIdx = " << partIdx << "\n"
                            << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
                            << "\tmotherPdgId = " << ((motherPdgId == 0) ? "No mother" : std::to_string(motherPdgId)) << "\n"
                            << "\tpartStatus = " << partStatus << "\n"
                            << "\tdR = " << dR_j2 << "\n";
                    std::cout << "\t";
                    Print(part);
                    std::cout << "\n";
                }
            }
        }

        if constexpr (matched_debug)
        {
            std::cout << "cand_b1 = ";
            Print(cand_b1);
            std::cout << "cand_b2 = ";
            Print(cand_b2);
            std::cout << "cand_j1 = ";
            Print(cand_j1);
            std::cout << "cand_j2 = ";
            Print(cand_j2);
            std::cout << "----------------------EVENT #" << i << " MATCHED DEBUG END-----------------------------\n";
        }
    }

    // auto c1 = std::make_unique<TCanvas>("c1", "c1");
    // c1->SetGrid();
    // c1->SetTickx();
    // c1->SetTicky();

    // auto stack = std::make_unique<THStack>("stack", "Total # of gen jets vs # of bad jets");
    // auto legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);

    // bad_gen_jets->SetLineWidth(3);
    // bad_gen_jets->SetLineColor(2);
    // stack->Add(bad_gen_jets.get());
    // legend->AddEntry(bad_gen_jets.get(), "bad gen jets", "l");

    // tot_gen_jets->SetLineWidth(3);
    // tot_gen_jets->SetLineColor(4);
    // stack->Add(tot_gen_jets.get());
    // legend->AddEntry(tot_gen_jets.get(), "total gen jets", "l");

    // stack->Draw("nostack");
    // stack->GetXaxis()->SetTitle("Number of gen jets");
    // legend->Draw();
    // c1->SaveAs("jet_numbers.png");

    // auto stack = std::make_unique<THStack>("stack", "All potential jet constituents vs # of unused constituents");
    // auto legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);

    // unused_cand->SetLineWidth(3);
    // unused_cand->SetLineColor(2);
    // stack->Add(unused_cand.get());
    // legend->AddEntry(unused_cand.get(), "unused", "l");

    // tot_cand->SetLineWidth(3);
    // tot_cand->SetLineColor(4);
    // stack->Add(tot_cand.get());
    // legend->AddEntry(tot_cand.get(), "all constit", "l");

    // stack->Draw("nostack");
    // stack->GetXaxis()->SetTitle("Number of gen jet constituents");
    // legend->Draw();
    // c1->SaveAs("jet_constituents.png");

    // empty_jets->SetLineWidth(3);
    // empty_jets->SetLineColor(4);
    // empty_jets->GetXaxis()->SetTitle("Number of empty jets");
    // empty_jets->Draw();
    // c1->SaveAs("empty_jets.png");

    std::cout << "diHiggsSL_cnt = " << diHiggsSL_cnt << "\n";
    std::cout << "good_matching = " << good_matching << "\n";
    
    return 0;
}