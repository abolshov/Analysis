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

constexpr bool debug = true;

int main()
{
    TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    TTree *myTree = static_cast<TTree*>(myFile->Get("Events"));

    UInt_t          nGenPart;
    Float_t         GenPart_eta[115];   //[nGenPart]
    Float_t         GenPart_mass[115];   //[nGenPart]
    Float_t         GenPart_phi[115];   //[nGenPart]
    Float_t         GenPart_pt[115];   //[nGenPart]
    Int_t           GenPart_genPartIdxMother[115];   //[nGenPart]
    Int_t           GenPart_pdgId[115];   //[nGenPart]
    Int_t           GenPart_status[115];   //[nGenPart]
    // Int_t           GenPart_statusFlags[115];   //[nGenPart]

    UInt_t          nGenJet;
    Float_t         GenJet_eta[19];   //[nGenJet]
    Float_t         GenJet_mass[19];   //[nGenJet]
    Float_t         GenJet_phi[19];   //[nGenJet]
    Float_t         GenJet_pt[19];   //[nGenJet]

    Int_t           GenJet_partonFlavour[19];   //[nGenJet]
    // UChar_t         GenJet_hadronFlavour[19];   //[nGenJet] doesn't work for some reason - prints empty spaces

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
    auto tot_cand = std::make_unique<TH1I>("tot_cand", "Number of potential jet constituents", 115, 0, 115);
    auto unused_cand = std::make_unique<TH1I>("unused_cand", "Number of unused candidates", 115, 0, 115);

    int diHiggsSL_cnt = 0;
    int good_matching = 0;

    int nEvents = myTree->GetEntries();
    std::cout << "nEvents = " << nEvents << "\n";

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        GenPartIndex idx = IsDiHiggsSL(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, nGenPart);
        if (!idx) continue;

        std::vector<Bool_t> candidates = PossibleJetConstituents(GenPart_genPartIdxMother, GenPart_status, nGenPart, idx);
        Int_t n_cand = std::count(candidates.begin(), candidates.end(), true);
        Int_t n_unused = n_cand;
        Int_t n_bad_jets = 0;
        // std::cout << "n_cand = " << n_cand << "\n";
        // std::cout << "nGenJet = " << nGenJet << "\n";

        tot_gen_jets->Fill(nGenJet);
        tot_cand->Fill(n_cand);

        for (Int_t jetIdx = 0; jetIdx < static_cast<Int_t>(nGenJet); ++jetIdx) // loop over all jets in event
        {
            TLorentzVector jet, jet_cand;
            jet.SetPtEtaPhiM(GenJet_pt[jetIdx], GenJet_eta[jetIdx], GenJet_phi[jetIdx], GenJet_mass[jetIdx]);
            
            // if constexpr (debug)
            // {
            //     std::cout << "Event " << i << ":\n";
            //     std::cout << "Jet = ";
            //     Print(jet);
            // }

            for (Int_t partIdx = 0; partIdx < static_cast<Int_t>(nGenPart); ++partIdx)
            {
                assert(nGenPart == candidates.size());
                if (!candidates[partIdx]) continue; // skip all particles that can't potentially constitue jet

                TLorentzVector part;
                part.SetPtEtaPhiM(GenPart_pt[partIdx], GenPart_eta[partIdx], GenPart_phi[partIdx], GenJet_mass[partIdx]);

                float dR = part.DeltaR(jet);
                if (dR < 0.4)
                {   
                    jet_cand += part;
                    --n_unused;
                    // Int_t motherIdx = GenPart_genPartIdxMother[partIdx];
                    // Int_t motherPdgId = (motherIdx == -1) ? 0 : GenPart_pdgId[motherIdx];
                    // Int_t partStatus = GenPart_status[partIdx];

                    // if constexpr (debug)
                    // {
                    //     std::cout << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
                    //               << "\tmotherPdgId = " << motherPdgId << "\n"
                    //               << "\tpartStatus = " << partStatus << "\n"
                    //               << "\tdR = " << dR << "\n";
                    //     std::cout << "\t";
                    //     Print(part);
                    //     std::cout << "\n";
                    // }
                }
            }

            // std::cout << "jet_cand = ";
            // Print(jet_cand);
            // std::cout << "dR = " << jet.DeltaR(jet_cand) << "\n";
            
            Double_t pt_ratio = jet_cand.Pt()/jet.Pt();
            if (pt_ratio > 1.1 || pt_ratio < 0.9) ++n_bad_jets;

            // if constexpr (debug) std::cout << "----------------------JET #" << jetIdx << " END-------------------------------\n";
        }

        unused_cand->Fill(n_unused);
        bad_gen_jets->Fill(n_bad_jets);

        // std::cout << "n_unused = " << n_unused << "\n";
        // std::cout << "n_bad_jets = " << n_bad_jets << "\n";
        // if constexpr (debug) std::cout << "----------------------EVENT #" << i << " END-----------------------------\n";
                
        // break;
        ++diHiggsSL_cnt;

        PtEtaPhiMArray genPart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, nGenPart};
        PtEtaPhiMArray genJet{GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, nGenJet};
        GenJetIndex match = Match(idx, genPart, genJet, GenPart_pdgId, GenJet_partonFlavour);

        if (!match) continue;
        ++good_matching;
    }

    auto c1 = std::make_unique<TCanvas>("c1", "c1");
    c1->SetGrid();
    c1->SetTickx();
    c1->SetTicky();

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

    auto stack = std::make_unique<THStack>("stack", "All potential jet constituents vs # of unused constituents");
    auto legend = std::make_unique<TLegend>(0.7, 0.7, 0.9, 0.9);

    unused_cand->SetLineWidth(3);
    unused_cand->SetLineColor(2);
    stack->Add(unused_cand.get());
    legend->AddEntry(unused_cand.get(), "unused", "l");

    tot_cand->SetLineWidth(3);
    tot_cand->SetLineColor(4);
    stack->Add(tot_cand.get());
    legend->AddEntry(tot_cand.get(), "all constit", "l");

    stack->Draw("nostack");
    stack->GetXaxis()->SetTitle("Gen jet constituents");
    legend->Draw();
    c1->SaveAs("jet_constituents.png");

    std::cout << "diHiggsSL_cnt = " << diHiggsSL_cnt << "\n";
    std::cout << "good_matching = " << good_matching << "\n";
    
    return 0;
}