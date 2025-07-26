#include "Cuts.hpp"
#include "ReadoutConstants.hpp"

bool BQuarkAcceptCut::IsSatisfied(Event const& event)
{
    using Quark = ReadConst::Quark;
    bool b1_cut = event.gen_quark_pt[static_cast<size_t>(Quark::b1)] > m_pt && std::abs(event.gen_quark_eta[static_cast<size_t>(Quark::b1)]) < m_eta;
    bool b2_cut = event.gen_quark_pt[static_cast<size_t>(Quark::b2)] > m_pt && std::abs(event.gen_quark_eta[static_cast<size_t>(Quark::b2)]) < m_eta;
    return b1_cut && b2_cut;
}