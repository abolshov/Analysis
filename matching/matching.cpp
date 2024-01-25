#include <iostream>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <memory>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"

#include "MatchingTools.hpp"

static constexpr int MAX_GENJET = 16;
static constexpr int MAX_GENPART = 270;

// static constexpr double MATCH_PT_THRESH = 3.0;

int main()
{
    // TFile *myFile = TFile::Open("NanoAODproduction_2017_cfg_NANO_1_RadionM400_HHbbWWToLNu2J.root");
    // TFile *myFile = TFile::Open("Run3_450GeV.root");
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

    int nEvents = myTree->GetEntries();
    std::cout << "nEvents = " << nEvents << "\n";

    for (int i = 0; i < nEvents; ++i)
    {
        myTree->GetEntry(i);
        
        int rad_idx = FindLast(RADION_ID, GenPart_pdgId, nGenPart);
        std::cout << "Last radion occurs at " << rad_idx << "\n";
        auto descendants = GetDescendants(rad_idx, GenPart_genPartIdxMother, nGenPart);
        // std::cout << descendants.size() << " " << descendants.capacity() << " " << std::boolalpha << descendants.empty() << "\n";
        PrintDecay(descendants, GenPart_pdgId, GenPart_genPartIdxMother, true);
        // auto test = FindSpecificDescendants({25}, -1, GenPart_genPartIdxMother, GenPart_pdgId, nGenPart);
        // std::cout << test.size() << " " << test.capacity() << " " << std::boolalpha << test.empty() << "\n";
        auto sig  = GetSignal(GenPart_pdgId, GenPart_genPartIdxMother, nGenPart);
        for (auto const& s: sig)
        {
            std::cout << s << " ";
        }
        std::cout << "\n";
        break;
    }    
    
    return 0;
}