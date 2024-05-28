#ifndef EVENTDATA_HPP
#define EVENTDATA_HPP

#include "TTree.h"
#include "TLorentzVector.h"

#include "Constants.hpp"

class EventData
{
    private:
        TTree&          m_tree;
    
    public:
        UInt_t          nGenPart;
        Float_t         GenPart_eta[MAX_GENPART];   
        Float_t         GenPart_mass[MAX_GENPART];   
        Float_t         GenPart_phi[MAX_GENPART];   
        Float_t         GenPart_pt[MAX_GENPART];   
        Int_t           GenPart_genPartIdxMother[MAX_GENPART];  
        Int_t           GenPart_pdgId[MAX_GENPART];   
        Int_t           GenPart_status[MAX_GENPART];   

        UInt_t          nGenJetAK4;
        Float_t         GenJetAK4_eta[MAX_AK4_GENJET];   
        Float_t         GenJetAK4_mass[MAX_AK4_GENJET];   
        Float_t         GenJetAK4_phi[MAX_AK4_GENJET];  
        Float_t         GenJetAK4_pt[MAX_AK4_GENJET];  
        Int_t           GenJetAK4_partonFlavour[MAX_AK4_GENJET];   
        UChar_t         GenJetAK4_hadronFlavour[MAX_AK4_GENJET];  

        Float_t         GenMET_phi;
        Float_t         GenMET_pt;

        EventData(TTree& tree);
        inline TLorentzVector ComputeP4(int i, Float_t const* pt, Float_t const* eta, Float_t const* phi, Float_t const* mass)
        {
            TLorentzVector p;
            p.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
            return p;
        }
};

#endif