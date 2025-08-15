#include "Cuts.hpp"
#include "ReadoutConstants.hpp"
#include "Utils.hpp"

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