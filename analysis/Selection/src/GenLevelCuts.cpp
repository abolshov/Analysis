#include "GenLevelCuts.hpp"
#include "ReadoutConstants.hpp"

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4D.h"
using LorentzVectorF_t = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Float_t>>;

#include "Math/GenVector/VectorUtil.h" 
using ROOT::Math::VectorUtil::DeltaR;

BQuarkAcceptCut::BQuarkAcceptCut(Float_t pt, Float_t eta)
    :   Specification<Event>("bquark_accept_cut", StrCat("pt > ", std::to_string(pt), " && abs(eta) < ", std::to_string(eta)))
    ,   m_pt(pt)
    ,   m_eta(eta)
    {}

bool BQuarkAcceptCut::IsSatisfied(Event const& event)
{
    using Quark = ReadConst::Quark;
    bool b1_cut = event.gen_quark_pt[static_cast<size_t>(Quark::b1)] > m_pt && std::abs(event.gen_quark_eta[static_cast<size_t>(Quark::b1)]) < m_eta;
    bool b2_cut = event.gen_quark_pt[static_cast<size_t>(Quark::b2)] > m_pt && std::abs(event.gen_quark_eta[static_cast<size_t>(Quark::b2)]) < m_eta;
    return b1_cut && b2_cut;
}

Resolved2bCut::Resolved2bCut(Float_t dr_thresh)
    :   Specification<Event>("res2b_cut ", StrCat("dR(b1, b2) > ", std::to_string(dr_thresh)))
    ,   m_threshold(dr_thresh)
    {}

bool Resolved2bCut::IsSatisfied(Event const& event)
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
    return DeltaR(b1, b2) >= m_threshold;
}

LeptonAcceptCutDL::LeptonAcceptCutDL(Float_t pt)
    :   Specification<Event>("lepton_accept_DL", StrCat("pt(l1) > ", std::to_string(pt), ", pt(l2) > ", std::to_string(pt)))
    ,   m_pt(pt)
    {}

bool LeptonAcceptCutDL::IsSatisfied(Event const& event)
{
    using Lep = ReadConst::Lep;
    return event.gen_lep_pt[static_cast<size_t>(Lep::lep1)] > m_pt && event.gen_lep_pt[static_cast<size_t>(Lep::lep2)];
}