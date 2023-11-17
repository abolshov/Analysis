#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "Utils.hpp"
#include "Matching.hpp"

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

    int diHiggsSL_cnt = 0;
    int good_matching = 0;

    int nEvents = myTree->GetEntries();
    std::cout << "nEvents = " << nEvents << "\n";

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);

        GenPartIndex idx = IsDiHiggsSL(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, nGenPart);
        if (!idx) continue;

        // std::cout << "Event " << i << "\n";
        // std::cout << "parton flavors: ";
        // std::copy(GenJet_partonFlavour, GenJet_partonFlavour + nGenJet, std::ostream_iterator<Int_t>(std::cout, " "));
        // std::cout << "\n";
        ++diHiggsSL_cnt;

        PtEtaPhiMArray genPart{GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, nGenPart};
        PtEtaPhiMArray genJet{GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, nGenJet};
        GenJetIndex match = Match(idx, genPart, genJet, GenPart_pdgId, GenJet_partonFlavour);

        if (!match) continue;
        ++good_matching;

        // for (Int_t jetIdx = 0; jetIdx < static_cast<Int_t>(nGenJet); ++jetIdx)
        // {
        //     TLorentzVector jet;
        //     jet.SetPtEtaPhiM(GenJet_pt[jetIdx], GenJet_eta[jetIdx], GenJet_phi[jetIdx], GenJet_mass[jetIdx]);
            
        //     if constexpr (debug)
        //     {
        //         std::cout << "Event " << i << ":\n";
        //         std::cout << "Jet = ";
        //         Print(jet);
        //     }

        //     for (Int_t partIdx = 0; partIdx < static_cast<Int_t>(nGenPart); ++partIdx)
        //     {
        //         TLorentzVector part;
        //         part.SetPtEtaPhiM(GenPart_pt[partIdx], GenPart_eta[partIdx], GenPart_phi[partIdx], GenJet_mass[partIdx]);

        //         float dR = part.DeltaR(jet);
        //         if (dR < 0.4)
        //         {   
        //             Int_t motherIdx = GenPart_genPartIdxMother[partIdx];
        //             Int_t motherPdgId = GenPart_pdgId[motherIdx];
        //             Int_t partStatus = GenPart_status[partIdx];

        //             if constexpr (debug)
        //             {
        //                 std::cout << "\tpdgId = " << GenPart_pdgId[partIdx] << "\n"
        //                           << "\tmotherPdgId = " << motherPdgId << "\n"
        //                           << "\tpartStatus = " << partStatus << "\n"
        //                           << "\tdR = " << dR << "\n";
        //                         //   << "\tbrute_force_dR = " << sqrt((jet.Phi() - part.Phi())*(jet.Phi() - part.Phi()) + (jet.Eta() - part.Eta())*(jet.Eta() - part.Eta())) << "\n";
        //                 std::cout << "\t";
        //                 Print(part);
        //                 std::cout << "\n";
        //             }
        //         }
        //     }
        //     if constexpr (debug) std::cout << "----------------------JET #" << jetIdx << " END-----------------------------\n";
        // }
        // if constexpr (debug) std::cout << "----------------------EVENT #" << i << " END-----------------------------\n";
    }

    std::cout << "diHiggsSL_cnt = " << diHiggsSL_cnt << "\n";
    std::cout << "good_matching = " << good_matching << "\n";
    
    return 0;
}