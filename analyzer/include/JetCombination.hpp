#ifndef JET_COMB_HPP
#define JET_COMB_HPP

#include "TString.h"

#include "Constants.hpp"

// redesign this struct to distinguish between:
// - DL channel, slim jets
// - DL channel, fat jet
// - SL channel, slim H->bb, slim W->qq
// - SL channel, slim H->bb, fat W->qq
// - SL channel, fat H->bb, slim W->qq
// - SL channel, fat H->bb, fat W->qq
struct JetComb
{
    JetComb() = default;
    JetComb(int b1_idx, int b2_idx, std::array<Float_t, MAX_RECO_JET> const& btag_score); // constructor for DL
    JetComb(int b1_idx, int b2_idx, int q1_idx, int q2_idx,
            std::array<Float_t, MAX_RECO_JET> const& btag_score,
            std::array<Float_t, MAX_RECO_JET> const& qvg_score); // constructor for SL

    int b1 = -1;
    int b2 = -1;
    int q1 = -1;
    int q2 = -1;
    Float_t btag_proba = 1.0f;
    Float_t qvg_proba = 1.0f;

    inline Float_t GetProbability() const { return btag_proba*qvg_proba; } 
    inline TString ToString(Channel ch) const
    {
        return ch == Channel::DL ? Form("b%db%d", b1, b2) : Form("b%db%dq%dq%d", b1, b2, q1, q2);
    }

    #ifdef DEV
        bool HasUniqueJets(Channel ch) const;
    #endif
};

#endif