#include "MatchingTools.hpp"
#include "Constants.hpp"

#include <unordered_set>

#include "Math/GenVector/VectorUtil.h" 
using ROOT::Math::VectorUtil::DeltaR;

int MatchIdx(LorentzVectorF_t const& parton, VecLVF_t const& jets)
{
    int res = -1;
    int sz = jets.size();
    Float_t min_dr = 10.0;
    for (int i = 0; i < sz; ++i)
    {
        Float_t dr = DeltaR(parton, jets[i]);
        if (dr < min_dr)
        {
            min_dr = dr;
            res = i;
        }
    }
    return res;
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