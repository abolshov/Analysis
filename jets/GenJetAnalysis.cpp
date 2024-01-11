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

constexpr int MAX_GENJET = 16;
constexpr int MAX_GENPART = 300;
constexpr double MATCH_PT_THRESH = 3.0;

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

        GenPartIndex index = IsDiHiggsSL(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, nGenPart);
        if (!index) continue;

        PtEtaPhiMArray genPart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, nGenPart};
        PtEtaPhiMArray genJet{GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, nGenJet};

        // std::vector<Bool_t> candidates = PossibleJetConstituents(GenPart_genPartIdxMother, GenPart_status, nGenPart, index);

        // std::cout << std::boolalpha;
        // std::cout << "b1 at " << index.b1 << "\n";
        // std::cout << "mother is " << GenPart_pdgId[GenPart_genPartIdxMother[index.b1]] << " at " << GenPart_genPartIdxMother[index.b1] << "\n";
        // TLorentzVector p;
        // for (Int_t idx = 0; idx < static_cast<Int_t>(nGenPart); ++idx)
        // {
            // if (IsFinal(idx, GenPart_genPartIdxMother, nGenPart) && DecayProductOf(idx, index.b1, GenPart_genPartIdxMother))
            // {
            //     std::cout << "\t" << GenPart_pdgId[idx] << ": at " << idx << " with pt = " << GenPart_pt[idx] << "\n";
            //     std::cout << "\t\t";
            //     std::cout << GenPart_pdgId[idx] << "(" << idx << ") ";
            //     Int_t cur_mother = GenPart_genPartIdxMother[idx];
            //     while (cur_mother != -1)
            //     {
            //         std::cout << GenPart_pdgId[cur_mother] << "(" << cur_mother << ") ";
            //         if (cur_mother == index.b1) break;
            //         cur_mother = GenPart_genPartIdxMother[cur_mother];
            //     }
            //     std::cout << "\n";
            //     TLorentzVector tmp;
            //     tmp.SetPtEtaPhiM(GenPart_pt[idx], GenPart_eta[idx], GenPart_phi[idx], GenPart_mass[idx]);
            //     p += tmp;
            // }
            // std::cout << GenPart_pdgId[idx] << ": " 
            //           << candidates[idx] << ", " 
            //           << IsFinal(idx, GenPart_genPartIdxMother, nGenPart) << ", "
            //           << GenPart_pt[idx] << ", "
            //           << DecayProductOf(idx, index.b2, GenPart_genPartIdxMother) << "\n";
        // }
        // std::cout << "b1 pt = " << GenPart_pt[index.b1] << "\n";
        // std::cout << "b2 pt = " << GenPart_pt[index.b2] << "\n";
        // std::cout << "products pt = " << p.Pt() << "\n";

        // TLorentzVector b1;
        // b1.SetPtEtaPhiM(GenPart_pt[index.b1], GenPart_eta[index.b1], GenPart_phi[index.b1], GenPart_mass[index.b1]);

        // for (Int_t jetIdx = 0; jetIdx < static_cast<Int_t>(nGenJet); ++jetIdx)
        // {
        //     TLorentzVector jet;
        //     jet.SetPtEtaPhiM(GenJet_pt[jetIdx], GenJet_eta[jetIdx], GenJet_phi[jetIdx], GenJet_mass[jetIdx]);
        //     std::cout << "dR = " << jet.DeltaR(b1) << ", pt = " << GenJet_pt[jetIdx] << "\n";
        // }
        // break;
        // Double_t pt_sum = 0.0;
        // for (Int_t idx = 0; idx < static_cast<Int_t>(nGenPart); ++idx)
        // {
        //     if (candidates[idx])
        //     {
        //         std::vector<Int_t> decay;
        //         Int_t part_id = GenPart_pdgId[idx];
        //         // std::cout << part_id << " ";
        //         decay.push_back(idx);
        //         Int_t mother_idx = GenPart_genPartIdxMother[idx];
        //         while (mother_idx != -1)
        //         {
        //             // std::cout << GenPart_pdgId[mother_idx] << " ";
        //             decay.push_back(mother_idx);
        //             mother_idx = GenPart_genPartIdxMother[mother_idx];
        //         }

        //         auto IsQuark = [&GenPart_pdgId](Int_t i) { return GenPart_pdgId[i] == -5; };
        //         if (std::find_if(decay.begin(), decay.end(), IsQuark) != decay.end())
        //         {
        //             if (candidates[decay[0]])
        //             {
        //                 pt_sum += GenPart_pt[decay[0]];
        //             }
        //             for (auto const& part: decay)
        //             {
        //                 std::cout << GenPart_pdgId[part] << " ";
        //             }
        //             std::cout << "\n";
        //             for (auto const& part: decay)
        //             {
        //                 std::cout << GenPart_status[part] << " ";
        //             }
        //             std::cout << "\n";
        //             for (auto const& part: decay)
        //             {
        //                 // std::cout << GenPart_pt[part] << " ";
        //                 std::cout << part << " ";
        //             }
        //             std::cout << "\n";
        //             std::cout << "===========================\n";
        //         }
        //     }
        // }
        // if (pt_sum != 0.0)
        // {
        //     std::cout << "Event " << i << ":\n";
        //     std::cout << "pt_sum = " << pt_sum << "\n";
        //     std::cout << "index.b1 = " << index.b1 << "\n";
        //     std::cout << "index.b2 = " << index.b2 << "\n";
        //     std::cout << "b1_pt = " << GenPart_pt[index.b1] << "\n";
        //     std::cout << "b2_pt = " << GenPart_pt[index.b2] << "\n";
        //     std::cout << "===================\n";
        // }

        ++diHiggsSL_cnt;

        GenJetIndex match = Match(index, genPart, genJet, GenPart_pdgId, GenJet_partonFlavour);

        if (!match) continue;
        ++good_matching;

        auto [bj1Idx, bj2Idx, lj1Idx, lj2Idx] = match;
        TLorentzVector bj1, bj2, lj1, lj2;
        bj1.SetPtEtaPhiM(GenJet_pt[bj1Idx], GenJet_eta[bj1Idx], GenJet_phi[bj1Idx], GenJet_mass[bj1Idx]);
        bj2.SetPtEtaPhiM(GenJet_pt[bj2Idx], GenJet_eta[bj2Idx], GenJet_phi[bj2Idx], GenJet_mass[bj2Idx]);
        lj1.SetPtEtaPhiM(GenJet_pt[lj1Idx], GenJet_eta[lj1Idx], GenJet_phi[lj1Idx], GenJet_mass[lj1Idx]);
        lj2.SetPtEtaPhiM(GenJet_pt[lj2Idx], GenJet_eta[lj2Idx], GenJet_phi[lj2Idx], GenJet_mass[lj2Idx]);
        
        TLorentzVector bq1, bq2, lq1, lq2;
        bq1.SetPtEtaPhiM(GenPart_pt[index.b1], GenPart_eta[index.b1], GenPart_phi[index.b1], GenPart_mass[index.b1]);
        bq2.SetPtEtaPhiM(GenPart_pt[index.b2], GenPart_eta[index.b2], GenPart_phi[index.b2], GenPart_mass[index.b2]);
        lq1.SetPtEtaPhiM(GenPart_pt[index.q1], GenPart_eta[index.q1], GenPart_phi[index.q1], GenPart_mass[index.q1]);
        lq2.SetPtEtaPhiM(GenPart_pt[index.q2], GenPart_eta[index.q2], GenPart_phi[index.q2], GenPart_mass[index.q2]);

        // printing problematic event info
        bool problematic_b = std::abs(bj1.Pt()/bq1.Pt()) > MATCH_PT_THRESH || std::abs(bj2.Pt()/bq2.Pt()) > MATCH_PT_THRESH;
        bool problematic_l = std::abs(lj1.Pt()/lq1.Pt()) > MATCH_PT_THRESH || std::abs(lj2.Pt()/lq2.Pt()) > MATCH_PT_THRESH;
        bool problematic_match = problematic_b || problematic_l;
        if (problematic_match)
        {
            std::cout << "------------------------------------------\n";
            std::cout << "Event " << i << "\n";
            std::cout << index << "\n";

            if (problematic_b)
            {
                std::cout << "Bad b quark/jet match\n";
                std::cout << "bj1 = ";
                Print(bj1);
                std::cout << "bj2 = ";
                Print(bj2);
                std::cout << "bq1 = ";
                Print(bq1);
                std::cout << "bq2 = ";
                Print(bq2);

                // for (Int_t idx = 0; idx < static_cast<Int_t>(nGenPart); ++idx)
                // {
                //     if (IsFinal(idx, GenPart_genPartIdxMother, nGenPart) && DecayProductOf(idx, index.b1, GenPart_genPartIdxMother))
                //     {
                //         // std::cout << "\t" << GenPart_pdgId[idx] << ": at " << idx << " with pt = " << GenPart_pt[idx] << "\n";
                //         std::cout << "\t\t";
                //         std::cout << GenPart_pdgId[idx] << "(" << idx << ") ";
                //         Int_t cur_mother = GenPart_genPartIdxMother[idx];
                //         while (cur_mother != -1)
                //         {
                //             std::cout << GenPart_pdgId[cur_mother] << "(" << cur_mother << ") ";
                //             if (cur_mother == index.b1) break;
                //             cur_mother = GenPart_genPartIdxMother[cur_mother];
                //         }
                //         std::cout << "\n";
                //     }
                // }
            }

            if (problematic_l)
            {
                std::cout << "Bad light quark/jet match\n";
                std::cout << "lj1 = ";
                Print(lj1);
                std::cout << "lj2 = ";
                Print(lj2);
                std::cout << "lq1 = ";
                Print(lq1);
                std::cout << "lq2 = ";
                Print(lq2);

                std::cout << "q1: pdgId = " << GenPart_pdgId[index.q1] << " at " << index.q1 << "\n";
                std::cout << "mother is " << GenPart_pdgId[GenPart_genPartIdxMother[index.q1]] << " at " << GenPart_genPartIdxMother[index.q1] << "\n";
                for (Int_t idx = 0; idx < static_cast<Int_t>(nGenPart); ++idx)
                {
                    if (IsFinal(idx, GenPart_genPartIdxMother, nGenPart) && DecayProductOf(idx, index.q1, GenPart_genPartIdxMother))
                    {
                        std::cout << "\tpdgId = " << GenPart_pdgId[idx] << ": at " << idx << " with pt = " << GenPart_pt[idx] << "\n";
                        std::cout << "\t\t";
                        std::cout << GenPart_pdgId[idx] << "(" << idx << ") ";
                        Int_t cur_mother = GenPart_genPartIdxMother[idx];
                        while (cur_mother != -1)
                        {
                            std::cout << GenPart_pdgId[cur_mother] << "(" << cur_mother << ") ";
                            if (cur_mother == index.q1) break;
                            cur_mother = GenPart_genPartIdxMother[cur_mother];
                        }
                        std::cout << "\n";
                    }
                }

                std::cout << "q2: pdgId = " << GenPart_pdgId[index.q2] << " at " << index.q2 << "\n";
                std::cout << "mother is " << GenPart_pdgId[GenPart_genPartIdxMother[index.q2]] << " at " << GenPart_genPartIdxMother[index.q2] << "\n";
                for (Int_t idx = 0; idx < static_cast<Int_t>(nGenPart); ++idx)
                {
                    if (IsFinal(idx, GenPart_genPartIdxMother, nGenPart) && DecayProductOf(idx, index.q2, GenPart_genPartIdxMother))
                    {
                        std::cout << "\tpdgId = " << GenPart_pdgId[idx] << ": at " << idx << " with pt = " << GenPart_pt[idx] << "\n";
                        std::cout << "\t\t";
                        std::cout << GenPart_pdgId[idx] << "(" << idx << ") ";
                        Int_t cur_mother = GenPart_genPartIdxMother[idx];
                        while (cur_mother != -1)
                        {
                            std::cout << GenPart_pdgId[cur_mother] << "(" << cur_mother << ") ";
                            if (cur_mother == index.q2) break;
                            cur_mother = GenPart_genPartIdxMother[cur_mother];
                        }
                        std::cout << "\n";
                    }
                }
            }
        }
    }

    std::cout << "diHiggsSL_cnt = " << diHiggsSL_cnt << "\n";
    std::cout << "good_matching = " << good_matching << "\n";
    
    return 0;
}