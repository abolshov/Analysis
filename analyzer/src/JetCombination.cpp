#include "JetCombination.hpp"

JetComb::JetComb(int b1_idx, int b2_idx, std::array<Float_t, MAX_RECO_JET> const& btag_score)
:   b1(b1_idx)
,   b2(b2_idx)
,   q1(-1)
,   q2(-1)
,   btag_proba(btag_score[b1_idx]*btag_score[b2_idx])
,   qvg_proba(1.0)
{}

JetComb::JetComb(int b1_idx, int b2_idx, int q1_idx, int q2_idx,
                 std::array<Float_t, MAX_RECO_JET> const& btag_score,
                 std::array<Float_t, MAX_RECO_JET> const& qvg_score)
:   b1(b1_idx)
,   b2(b2_idx)
,   q1(q1_idx)
,   q2(q2_idx)
,   btag_proba(btag_score[b1_idx]*btag_score[b2_idx])
,   qvg_proba(qvg_score[q1_idx]*qvg_score[q2_idx])
{}

#ifdef DEV
    #include <unordered_set>

    bool JetComb::HasUniqueJets(Channel ch) const
    {
        std::unordered_set<int> indices;
        indices.insert(b1);
        indices.insert(b2);
        if (ch == Channel::SL)
        {
            indices.insert(q1);
            indices.insert(q2);
        }
        // return indices.size() == QuarkCount(ch) && !indices.contains(-1); // c++20
        return indices.size() == QuarkCount(ch) && !indices.count(-1);
    }
#endif