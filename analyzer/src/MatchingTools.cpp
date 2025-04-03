#include "MatchingTools.hpp"
#include "Constants.hpp"

#include <unordered_set>

#include "Math/GenVector/VectorUtil.h" 
using ROOT::Math::VectorUtil::DeltaR;

int MatchIdx(LorentzVectorF_t const& parton, VecLVF_t const& jets)
{
    auto Cmp = [&parton](LorentzVectorF_t const& j1, LorentzVectorF_t const& j2)
    {
        return DeltaR(parton, j1) < DeltaR(parton, j2); 
    };
    auto it = std::min_element(jets.begin(), jets.end(), Cmp);
    return DeltaR(parton, *it) < 0.4 ? it - jets.begin() : -1;
}

TString MakeTrueLabel(VecLVF_t const& gen, VecLVF_t const& reco)
{
    std::unordered_set<int> matches;
    int b1_match = MatchIdx(gen[static_cast<size_t>(Quark::b1)], reco);
    int b2_match = MatchIdx(gen[static_cast<size_t>(Quark::b2)], reco);
    matches.insert(b1_match);
    matches.insert(b2_match);

    if (gen.size() == NUM_BQ)
    {
        if (matches.size() != NUM_BQ || matches.count(-1))
        {
            return {};
        }
        return Form("b%db%d", std::min(b1_match, b2_match), std::max(b1_match, b2_match));
    }

    int q1_match = MatchIdx(gen[static_cast<size_t>(Quark::q1)], reco);
    int q2_match = MatchIdx(gen[static_cast<size_t>(Quark::q2)], reco);
    matches.insert(q1_match);
    matches.insert(q2_match);

    if (gen.size() == static_cast<size_t>(Quark::count))
    {
        if (matches.size() != static_cast<size_t>(Quark::count) || matches.count(-1))
        {
            return {};
        }
        return Form("b%db%dq%dq%d", std::min(b1_match, b2_match), std::max(b1_match, b2_match), std::min(q1_match, q2_match), std::max(q1_match, q2_match));
    }

    return {};
}

JetComb FindMatch(VecLVF_t const& quarks, VecLVF_t const& jets,  Channel ch)
{
    JetComb res{};
    res.b1 = MatchIdx(quarks[static_cast<size_t>(Quark::b1)], jets);
    res.b2 = MatchIdx(quarks[static_cast<size_t>(Quark::b2)], jets);
    if (ch == Channel::SL)
    {
        res.q1 = MatchIdx(quarks[static_cast<size_t>(Quark::q1)], jets);
        res.q2 = MatchIdx(quarks[static_cast<size_t>(Quark::q2)], jets);
    }
    return res;
}