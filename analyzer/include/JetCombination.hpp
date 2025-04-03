#ifndef JET_COMB_HPP
#define JET_COMB_HPP

#include "TString.h"

#include "Constants.hpp"

struct JetComb
{
    int b1 = -1;
    int b2 = -1;
    int q1 = -1;
    int q2 = -1;

    inline TString ToString(Channel ch) const
    {
        return ch == Channel::DL ? Form("b%db%d", b1, b2) : Form("b%db%dq%dq%d", b1, b2, q1, q2);
    }

    #ifdef DEV
        bool HasUniqueJets(Channel ch) const;
    #endif
};

#endif