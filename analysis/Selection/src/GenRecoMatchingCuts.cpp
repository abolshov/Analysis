#include "GenRecoMatchingCuts.hpp"
#include "Utils.hpp"

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4D.h"
using LorentzVectorF_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Float_t>>;

#include "Math/GenVector/VectorUtil.h" 
using ROOT::Math::VectorUtil::DeltaR;

AK4BJetMatchingCut::AK4BJetMatchingCut(Float_t dr_thresh)
:   Specification<Event>("ak4_matching", "each gen b quark has a ak4 reco jet match")
,   m_match_thresh(dr_thresh)
{}

bool AK4BJetMatchingCut::IsSatisfied(Event const& event)
{
    using Quark = ReadConst::Quark;
    LorentzVectorF_t b1{event.gen_quark_pt[static_cast<size_t>(Quark::b1)],
                        event.gen_quark_eta[static_cast<size_t>(Quark::b1)],
                        event.gen_quark_phi[static_cast<size_t>(Quark::b1)],
                        event.gen_quark_mass[static_cast<size_t>(Quark::b1)]};
    LorentzVectorF_t b2{event.gen_quark_pt[static_cast<size_t>(Quark::b2)],
                        event.gen_quark_eta[static_cast<size_t>(Quark::b2)],
                        event.gen_quark_phi[static_cast<size_t>(Quark::b2)],
                        event.gen_quark_mass[static_cast<size_t>(Quark::b2)]};

    Int_t nj = event.n_reco_jet;
    
    Int_t b1_match = -1;
    Int_t b2_match = -1;

    Float_t min_dr1 = m_match_thresh;
    Float_t min_dr2 = m_match_thresh;

    for (Int_t i = 0; i < nj; ++i)
    {
        LorentzVectorF_t jet{event.reco_jet_pt[i],
                             event.reco_jet_eta[i],
                             event.reco_jet_phi[i],
                             event.reco_jet_mass[i]};
        
        Float_t dr1 = DeltaR(jet, b1);
        if (dr1 < m_match_thresh && dr1 < min_dr1)
        {   
            min_dr1 = dr1;
            b1_match = i;
        }

        Float_t dr2 = DeltaR(jet, b2);
        if (dr2 < m_match_thresh && dr2 < min_dr2)
        {   
            min_dr2 = dr2;
            b1_match = i;
        }
    }

    return b1_match != -1 && b2_match != -1 && b1_match != b2_match;
}

AK8BJetMatchingCut::AK8BJetMatchingCut(Float_t dr_thresh)
:   Specification<Event>("ak8_matching", "gen H->bb has a reco ak8 fatjet match")
,   m_match_thresh(dr_thresh)
{}

bool AK8BJetMatchingCut::IsSatisfied(Event const& event)
{
    LorentzVectorF_t hbb{event.hbb_pt, event.hbb_eta, event.hbb_phi, event.hbb_mass};

    Int_t nfj = event.n_reco_fatjet;
    Int_t match = -1;
    Float_t min_dr = m_match_thresh;

    for (Int_t i = 0; i < nfj; ++i)
    {
        LorentzVectorF_t fatjet{event.reco_fatjet_pt[i],
                                event.reco_fatjet_eta[i],
                                event.reco_fatjet_phi[i],
                                event.reco_fatjet_mass[i]};
        
        Float_t dr = DeltaR(fatjet, hbb);
        if (dr < m_match_thresh && dr < min_dr)
        {   
            min_dr = dr;
            match = i;
        }
    }

    return match != -1;
}